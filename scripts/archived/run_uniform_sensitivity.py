"""Run the image's uniform refractor example and one-parameter sweeps.

This thin experiment driver reuses the repository QMC sampler, Sinkhorn solver,
notebook artifact builder, and NPZ writer. The default output is
``results/uniform_sensitivity``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from functools import partial
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from refracter.cost import get_cost_function, set_kappa, validate_transport_possible
from refracter.qmc import gen_spherical_patch
from refracter.sinkhorn import entropic_sinkhorn_divergence, run_sinkhorn_divergence
from scripts.archived.generate_sinkhorn_refracter_examples import (
    patch_indicator,
    save_surface_plot,
)
from scripts.generate_results import build_refraction_artifacts, save_refraction_npz


DEG = math.pi / 180.0
NK = 1600
CHUNK_SIZE = 512
GRID_RESOLUTION = 256
GRID_LIMIT = 0.6
PUSHFORWARD_EPSILON = 5e-3
PUSHFORWARD_MAX_ITER = 2000
PUSHFORWARD_TOLERANCE = 1e-9
L2_COST = get_cost_function("l2")
SOURCE = {"theta_deg": [15.0, 60.0], "phi_deg": [15.0, 45.0]}
TARGET = {"theta_deg": [18.0, 36.0], "phi_deg": [18.0, 36.0]}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "results" / "uniform_sensitivity",
        help="folder for NPZ runs, surface PNGs, and L2 cost summaries",
    )
    return parser.parse_args()


def sample_patch(bounds: dict[str, list[float]]) -> np.ndarray:
    theta_min, theta_max = np.asarray(bounds["theta_deg"]) * DEG
    phi_min, phi_max = np.asarray(bounds["phi_deg"]) * DEG
    return gen_spherical_patch(NK, theta_min, theta_max, phi_min, phi_max, skip=0)


def pushforward_sinkhorn_3d_l2(
    pushed_points: np.ndarray,
    target_points: np.ndarray,
    p: np.ndarray,
    q: np.ndarray,
) -> tuple[float, dict[str, object]]:
    """Measure the 3D empirical pushforward discrepancy with shared EOT code."""
    return entropic_sinkhorn_divergence(
        pushed_points,
        target_points,
        p,
        q,
        epsilon=PUSHFORWARD_EPSILON,
        chunk_size=CHUNK_SIZE,
        max_iter=PUSHFORWARD_MAX_ITER,
        tolerance=PUSHFORWARD_TOLERANCE,
        cost_fn=L2_COST,
    )


def cases() -> list[dict[str, object]]:
    output = [{"label": "original_uniform", "category": "original", "kappa": 0.6,
               "source": SOURCE.copy(), "parameter": "baseline"}]
    for kappa in (0.25, 0.50, 0.90):
        output.append({
            "label": f"kappa_{kappa:.2f}".replace(".", "p"),
            "category": "kappa_sweep",
            "kappa": kappa,
            "source": SOURCE.copy(),
            "parameter": f"kappa={kappa:.2f}",
        })
    for shift in (5, 15, 25):
        output.append({
            "label": f"source_phi_plus_{shift:02d}deg",
            "category": "source_phi_translation",
            "kappa": 0.6,
            "source": {
                "theta_deg": SOURCE["theta_deg"].copy(),
                "phi_deg": [v + shift for v in SOURCE["phi_deg"]],
            },
            "parameter": f"source phi bounds +{shift} degrees",
        })
    for shift in (10, 25, 40):
        output.append({
            "label": f"source_theta_plus_{shift:02d}deg",
            "category": "source_theta_translation",
            "kappa": 0.6,
            "source": {
                "theta_deg": [v + shift for v in SOURCE["theta_deg"]],
                "phi_deg": SOURCE["phi_deg"].copy(),
            },
            "parameter": f"source theta bounds +{shift} degrees",
        })
    for case in output:
        case["target"] = TARGET.copy()
    return output


def run_case(case: dict[str, object], output_dir: Path) -> dict[str, object]:
    started = time.perf_counter()
    label = str(case["label"])
    kappa = float(case["kappa"])
    source = case["source"]
    target = case["target"]
    x, y = sample_patch(source), sample_patch(target)
    p = np.full(NK, 1.0 / NK, dtype=np.float64)
    q = np.full(NK, 1.0 / NK, dtype=np.float64)

    set_kappa(kappa)
    validate_transport_possible(x, y, chunk_size=CHUNK_SIZE)
    max_dot = float(np.max(x @ y.T))
    print(f"\n{label}: kappa={kappa:g}, source={source}")

    solution = run_sinkhorn_divergence(
        x, y, p, q, chunk_size=CHUNK_SIZE, verbose=False,
    )
    artifacts = build_refraction_artifacts(
        x,
        y,
        p,
        q,
        solution["f_raw"],
        solution["g_raw"],
        solution["f"],
        solution["g"],
        partial(patch_indicator, bounds=source),
        partial(patch_indicator, bounds=target),
        chunk=CHUNK_SIZE,
        grid_res=GRID_RESOLUTION,
        grid_limit=GRID_LIMIT,
        kappa=kappa,
    )
    sinkhorn_divergence, sinkhorn_diagnostics = pushforward_sinkhorn_3d_l2(
        artifacts["y_push_3d"], y, p, q,
    )
    sinkhorn_w2 = math.sqrt(max(sinkhorn_divergence, 0.0))
    save_surface_plot(
        output_dir,
        label,
        kappa,
        artifacts["Ref"],
        artifacts["R"],
        "sqrt(Sε) 3D",
        sinkhorn_w2,
        filename=f"{label}_surface.png",
        subtitle=f"source θ={source['theta_deg']}, φ={source['phi_deg']}",
    )
    mass_l2_squared = float(np.sum((artifacts["pushed_mass"] - q) ** 2))
    mass_l2 = math.sqrt(mass_l2_squared)

    config = {
        "label": label,
        "category": case["category"],
        "parameter": case["parameter"],
        "kappa": kappa,
        "source_bounds_deg": source,
        "target_bounds_deg": target,
        "nk": NK,
        "k_final": 8 * int(np.floor(np.sqrt(NK))),
        "pushforward_metric": "3D squared Euclidean Sinkhorn divergence",
        "pushforward_epsilon": PUSHFORWARD_EPSILON,
    }
    metadata = {
        "src_density": "uniform",
        "tgt_density": "uniform",
        "source_bounds_json": json.dumps(source, sort_keys=True),
        "target_bounds_json": json.dumps(target, sort_keys=True),
        "pushforward_target_sinkhorn_divergence_3d_l2": np.float64(sinkhorn_divergence),
        "pushforward_target_sinkhorn_w2_3d_l2": np.float64(sinkhorn_w2),
        "pushforward_target_sinkhorn_epsilon": np.float64(PUSHFORWARD_EPSILON),
        "pushforward_target_sinkhorn_max_iter": np.int64(PUSHFORWARD_MAX_ITER),
        "pushforward_target_sinkhorn_tolerance": np.float64(PUSHFORWARD_TOLERANCE),
        "pushforward_target_sinkhorn_converged": np.bool_(sinkhorn_diagnostics["all_converged"]),
        "pushforward_target_sinkhorn_diagnostics_json": json.dumps(
            sinkhorn_diagnostics, sort_keys=True,
        ),
        "pushforward_target_mass_l2_squared": np.float64(mass_l2_squared),
        "pushforward_target_mass_l2": np.float64(mass_l2),
        "max_source_target_dot": np.float64(max_dot),
        "hard_map_unique_targets": np.int64(np.count_nonzero(artifacts["pushed_mass"])),
        "run_config_json": json.dumps(config, sort_keys=True),
    }
    npz_path = output_dir / f"{label}.npz"
    save_refraction_npz(
        npz_path,
        x=x,
        y=y,
        p=p,
        q=q,
        f_raw=solution["f_raw"],
        g_raw=solution["g_raw"],
        f_id=solution["f_id"],
        g_id=solution["g_id"],
        f=solution["f"],
        g=solution["g"],
        artifacts=artifacts,
        kappa=kappa,
        metadata=metadata,
    )
    result = {
        **config,
        "npz": npz_path.name,
        "max_source_target_dot": max_dot,
        "pushforward_target_sinkhorn_divergence_3d_l2": sinkhorn_divergence,
        "pushforward_target_sinkhorn_w2_3d_l2": sinkhorn_w2,
        "pushforward_target_sinkhorn_converged": sinkhorn_diagnostics["all_converged"],
        "pushforward_target_sinkhorn_diagnostics": sinkhorn_diagnostics,
        "pushforward_target_mass_l2_squared": mass_l2_squared,
        "pushforward_target_mass_l2": mass_l2,
        "unique_targets_hit": int(metadata["hard_map_unique_targets"]),
        "runtime_seconds": time.perf_counter() - started,
    }
    print(
        f"  saved {npz_path.name}; S_eps(3D L2)={sinkhorn_divergence:.10g}; "
        f"sqrt(S_eps)={sinkhorn_w2:.10g}; "
        f"converged={sinkhorn_diagnostics['all_converged']}; "
        f"unique targets={result['unique_targets_hit']}/{NK}; "
        f"time={result['runtime_seconds']:.1f}s"
    )
    return result


def write_summary(output_dir: Path, results: list[dict[str, object]]) -> None:
    csv_fields = [
        "case", "category", "parameter", "kappa", "NK", "k_final",
        "source_theta_min_deg", "source_theta_max_deg",
        "source_phi_min_deg", "source_phi_max_deg",
        "target_theta_min_deg", "target_theta_max_deg",
        "target_phi_min_deg", "target_phi_max_deg", "max_kappa_dot",
        "epsilon", "ground_cost", "sinkhorn_divergence_3d_l2",
        "sqrt_sinkhorn_divergence_3d_l2", "sinkhorn_converged",
        "cross_iterations", "source_self_iterations", "target_self_iterations",
        "mass_L2_squared", "mass_L2", "unique_targets", "total_targets", "npz",
    ]
    csv_rows = []
    lines = [
        "Uniform refraction example and one-parameter sensitivity sweeps",
        "",
        f"NK: {NK}",
        "Original setup: kappa=0.6; source theta=[15, 60] deg, phi=[15, 45] deg;",
        "target theta=[18, 36] deg, phi=[18, 36] deg.",
        "Surface solver/artifacts: existing refracter Sinkhorn and generate_results workflow.",
        "Uniform discrete source and target weights; original Halton sampling convention.",
        "",
        "Pushforward metric: debiased entropic Sinkhorn divergence S_epsilon between",
        "the pushed 3D sphere points and target 3D sphere points, using squared Euclidean",
        f"cost and epsilon={PUSHFORWARD_EPSILON:g}. sqrt(S_epsilon) is a regularized",
        "W2-like value, not exact unregularized W2. Convergence is reported per case.",
        "Mass-L2 separately compares pushed target-atom masses with target masses.",
        "All cases pass validate_transport_possible.",
        "",
        "case | category | kappa | source theta deg | source phi deg | max kappa dot |",
        "S_epsilon (3D L2) | sqrt(S_epsilon) | converged | EOT iterations cross/self/self |",
        "mass L2^2 | mass L2 | unique targets | NPZ",
        "-" * 190,
    ]
    for result in results:
        source = result["source_bounds_deg"]
        target = result["target_bounds_deg"]
        diag = result["pushforward_target_sinkhorn_diagnostics"]
        cross_iter = diag["cross"]["iterations"]
        source_iter = diag["source_self"]["iterations"]
        target_iter = diag["target_self"]["iterations"]
        lines.append(
            f"{result['label']} | {result['category']} | {result['kappa']:.6f} | "
            f"{source['theta_deg']} | {source['phi_deg']} | "
            f"{result['kappa'] * result['max_source_target_dot']:.10f} | "
            f"{result['pushforward_target_sinkhorn_divergence_3d_l2']:.12g} | "
            f"{result['pushforward_target_sinkhorn_w2_3d_l2']:.12g} | "
            f"{result['pushforward_target_sinkhorn_converged']} | "
            f"{cross_iter}/{source_iter}/{target_iter} | "
            f"{result['pushforward_target_mass_l2_squared']:.12g} | "
            f"{result['pushforward_target_mass_l2']:.12g} | "
            f"{result['unique_targets_hit']}/{NK} | {result['npz']}"
        )
        csv_rows.append({
            "case": result["label"],
            "category": result["category"],
            "parameter": result["parameter"],
            "kappa": result["kappa"],
            "NK": NK,
            "k_final": result["k_final"],
            "source_theta_min_deg": source["theta_deg"][0],
            "source_theta_max_deg": source["theta_deg"][1],
            "source_phi_min_deg": source["phi_deg"][0],
            "source_phi_max_deg": source["phi_deg"][1],
            "target_theta_min_deg": target["theta_deg"][0],
            "target_theta_max_deg": target["theta_deg"][1],
            "target_phi_min_deg": target["phi_deg"][0],
            "target_phi_max_deg": target["phi_deg"][1],
            "max_kappa_dot": result["kappa"] * result["max_source_target_dot"],
            "epsilon": PUSHFORWARD_EPSILON,
            "ground_cost": "squared_euclidean_3d",
            "sinkhorn_divergence_3d_l2": result["pushforward_target_sinkhorn_divergence_3d_l2"],
            "sqrt_sinkhorn_divergence_3d_l2": result["pushforward_target_sinkhorn_w2_3d_l2"],
            "sinkhorn_converged": result["pushforward_target_sinkhorn_converged"],
            "cross_iterations": cross_iter,
            "source_self_iterations": source_iter,
            "target_self_iterations": target_iter,
            "mass_L2_squared": result["pushforward_target_mass_l2_squared"],
            "mass_L2": result["pushforward_target_mass_l2"],
            "unique_targets": result["unique_targets_hit"],
            "total_targets": NK,
            "npz": result["npz"],
        })
    (output_dir / "L2_costs.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    with (output_dir / "L2_costs_numeric.csv").open(
        "w", newline="", encoding="utf-8",
    ) as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=csv_fields)
        writer.writeheader()
        writer.writerows(csv_rows)
    metadata = {
        "nk": NK,
        "k_final": 8 * int(np.floor(np.sqrt(NK))),
        "solver": "refracter.sinkhorn.run_sinkhorn_divergence",
        "artifact_builder": "scripts.generate_results.build_refraction_artifacts",
        "npz_writer": "scripts.generate_results.save_refraction_npz",
        "l2_metric": "sqrt of debiased entropic Sinkhorn divergence using squared 3D Euclidean ground cost",
        "pushforward_metric": "debiased entropic Sinkhorn divergence with squared 3D Euclidean cost",
        "pushforward_epsilon": PUSHFORWARD_EPSILON,
        "pushforward_max_iter": PUSHFORWARD_MAX_ITER,
        "pushforward_tolerance": PUSHFORWARD_TOLERANCE,
        "pushforward_solver": "refracter.sinkhorn.entropic_sinkhorn_divergence",
        "results": results,
    }
    (output_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = [run_case(case, args.output_dir) for case in cases()]
    write_summary(args.output_dir, results)
    print(f"\nWrote cost summaries and surface images to {args.output_dir}")


if __name__ == "__main__":
    main()
