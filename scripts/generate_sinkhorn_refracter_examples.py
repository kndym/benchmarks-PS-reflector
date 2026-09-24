"""
Generate radius-coloured Sinkhorn refracter plots and paper-style pushforward
error diagnostics.

The output contains one default refraction case and ten deterministic random
cases.  Every case uses equal-sized clouds sampled uniformly with respect to
spherical area on source and target patches, uniform discrete source/target
weights, and the repository Sinkhorn-divergence solver.

For each case this script saves:

* ``*_surface_3d.png``: the refracter surface ``2 R x`` in 3-D, coloured by R;
* ``*_source_xy.png``: source unit vectors projected to their x/y coordinates,
  coloured by the same radius values;
* ``*.npz``: notebook-compatible clouds, potentials, refractor/c-transform
  arrays, projected distributions, and hard-map diagnostics.

The hard push-forward map is

    x_i -> y_argmin_j [ c(x_i, y_j) - g_raw[j] ].

The push-forward diagnostics use independent log-domain Sinkhorn-divergence
solves, matching the paper's small-regularisation error estimator.  The primary
error uses squared Euclidean cost after north-pole stereographic projection;
angular and ambient 3-D costs are retained as supplementary diagnostics.

Run from the repository root:

    python scripts/generate_sinkhorn_refracter_examples.py

Useful options:

    --nk 256 --random-cases 10 --seed 20260914
    --output-dir results/sinkhorn_refracter_examples
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import logsumexp

from refracter.build import c_transform_fc, c_transform_gc
from refracter.cost import (
    cost_matrix_chunk,
    get_kappa,
    set_kappa,
    validate_transport_possible,
)
from refracter.distributions import stereo_north
from refracter.sinkhorn import run_sinkhorn_divergence


DEG = math.pi / 180.0
DEFAULT_BOUNDS = {
    "source": {
        "theta_deg": [15.0, 60.0],
        "phi_deg": [15.0, 45.0],
    },
    "target": {
        "theta_deg": [18.0, 36.0],
        "phi_deg": [18.0, 36.0],
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nk", type=int, default=400,
                        help="points in each cloud (default: 400)")
    parser.add_argument("--random-cases", type=int, default=10,
                        help="number of random cases in addition to default (default: 10)")
    parser.add_argument("--seed", type=int, default=20260914,
                        help="seed for random kappa and patch bounds")
    parser.add_argument("--chunk-size", type=int, default=256,
                        help="chunk size for cost/Sinkhorn operations")
    parser.add_argument("--pushforward-eps", type=float, default=1e-6,
                        help="regularisation epsilon for paper-style pushforward Sinkhorn")
    parser.add_argument("--pushforward-max-iter", type=int, default=2000,
                        help="iteration cap for each pushforward Sinkhorn solve")
    parser.add_argument("--pushforward-tol", type=float, default=1e-10,
                        help="stopping tolerance for each pushforward Sinkhorn solve")
    parser.add_argument("--output-dir", default=None,
                        help="output directory (default: results/sinkhorn_refracter_examples)")
    return parser.parse_args()


def _halton(index: int, base: int) -> float:
    result = 0.0
    factor = 1.0
    current = index
    while current > 0:
        factor /= base
        result += factor * (current % base)
        current //= base
    return result


def gen_spherical_patch(n: int, bounds: dict[str, list[float]], skip: int = 0) -> np.ndarray:
    """Sample a spherical patch uniformly in solid angle using a Halton cloud."""
    theta_min, theta_max = np.asarray(bounds["theta_deg"], dtype=float) * DEG
    phi_min, phi_max = np.asarray(bounds["phi_deg"], dtype=float) * DEG

    # phi is the polar angle; equal-area sampling is uniform in cos(phi).
    cos_phi_max = math.cos(phi_min)
    cos_phi_min = math.cos(phi_max)
    points = np.empty((n, 3), dtype=np.float64)
    for row in range(n):
        index = skip + row
        u_phi = _halton(index, 2)
        u_theta = _halton(index, 3)
        cos_phi = cos_phi_min + u_phi * (cos_phi_max - cos_phi_min)
        sin_phi = math.sqrt(max(0.0, 1.0 - cos_phi * cos_phi))
        theta = theta_min + u_theta * (theta_max - theta_min)
        points[row] = [
            sin_phi * math.cos(theta),
            sin_phi * math.sin(theta),
            cos_phi,
        ]
    return points


def _angles_deg(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    phi = np.arccos(np.clip(points[:, 2], -1.0, 1.0)) / DEG
    theta = np.arctan2(points[:, 1], points[:, 0]) / DEG
    theta = np.mod(theta, 360.0)
    return theta, phi


def patch_coverage(points: np.ndarray, bounds: dict[str, list[float]]) -> dict[str, object]:
    """Check that every generated point lies within its requested patch."""
    theta, phi = _angles_deg(points)
    theta_min, theta_max = bounds["theta_deg"]
    phi_min, phi_max = bounds["phi_deg"]
    tolerance = 2e-10
    inside = (
        (theta >= theta_min - tolerance)
        & (theta <= theta_max + tolerance)
        & (phi >= phi_min - tolerance)
        & (phi <= phi_max + tolerance)
    )
    return {
        "count": int(len(points)),
        "inside_count": int(inside.sum()),
        "fully_covered": bool(inside.all()),
        "theta_observed_deg": [float(theta.min()), float(theta.max())],
        "phi_observed_deg": [float(phi.min()), float(phi.max())],
    }


def random_patch(rng: np.random.Generator) -> dict[str, list[float]]:
    """Draw a compact upper-hemisphere patch with non-wrapping azimuth."""
    theta_min = float(rng.uniform(0.0, 280.0))
    theta_width = float(rng.uniform(12.0, 72.0))
    phi_min = float(rng.uniform(5.0, 55.0))
    phi_width = float(rng.uniform(12.0, 32.0))
    return {
        "theta_deg": [theta_min, theta_min + theta_width],
        "phi_deg": [phi_min, min(phi_min + phi_width, 88.0)],
    }


def random_case(rng: np.random.Generator) -> dict[str, object]:
    return {
        "kappa": float(rng.uniform(0.25, 0.90)),
        "source": random_patch(rng),
        "target": random_patch(rng),
    }


def patch_indicator(points: np.ndarray, bounds: dict[str, list[float]]) -> np.ndarray:
    """Evaluate a uniform spherical-patch density on arbitrary points."""
    theta, phi = _angles_deg(np.asarray(points, dtype=np.float64))
    theta_min, theta_max = bounds["theta_deg"]
    phi_min, phi_max = bounds["phi_deg"]
    tolerance = 2e-10
    inside = (
        (theta >= theta_min - tolerance)
        & (theta <= theta_max + tolerance)
        & (phi >= phi_min - tolerance)
        & (phi <= phi_max + tolerance)
    )
    return inside.astype(np.float64)


def _sinkhorn_ot_cost(
    cost: np.ndarray,
    source_weights: np.ndarray,
    target_weights: np.ndarray,
    epsilon: float,
    max_iter: int,
    tolerance: float,
) -> tuple[float, int, bool]:

    if epsilon <= 0.0:
        raise ValueError(f"epsilon must be positive, got {epsilon}")
    if max_iter <= 0 or tolerance <= 0.0:
        raise ValueError("max_iter and tolerance must be positive")

    cost = np.asarray(cost, dtype=np.float64)
    source_weights = np.asarray(source_weights, dtype=np.float64)
    target_weights = np.asarray(target_weights, dtype=np.float64)
    if cost.ndim != 2 or cost.shape != (len(source_weights), len(target_weights)):
        raise ValueError("cost shape must match the source and target weights")
    if not np.isfinite(cost).all():
        raise ValueError("cost matrix must contain only finite values")
    if np.any(source_weights <= 0.0) or np.any(target_weights <= 0.0):
        raise ValueError("Sinkhorn error weights must be strictly positive")

    source_weights = source_weights / source_weights.sum()
    target_weights = target_weights / target_weights.sum()
    log_source = np.log(source_weights)
    log_target = np.log(target_weights)
    inv_epsilon = 1.0 / epsilon
    f = np.zeros(len(source_weights), dtype=np.float64)
    g = np.zeros(len(target_weights), dtype=np.float64)
    converged = False

    for iteration in range(1, max_iter + 1):
        f_new = -epsilon * logsumexp(
            (g[None, :] - cost) * inv_epsilon + log_target[None, :],
            axis=1,
        )
        g_new = -epsilon * logsumexp(
            (f_new[:, None] - cost) * inv_epsilon + log_source[:, None],
            axis=0,
        )
        change = max(
            float(np.max(np.abs(f_new - f))),
            float(np.max(np.abs(g_new - g))),
        )
        f, g = f_new, g_new
        if change <= tolerance:
            converged = True
            break

    value = float(np.dot(source_weights, f) + np.dot(target_weights, g))
    return value, iteration, converged


def _sinkhorn_divergence(
    source_points: np.ndarray,
    target_points: np.ndarray,
    source_weights: np.ndarray,
    target_weights: np.ndarray,
    cost_fn,
    epsilon: float,
    max_iter: int,
    tolerance: float,
) -> tuple[float, dict[str, object]]:
    """Return paper-style Sinkhorn divergence for one ground cost."""
    cross = np.asarray(cost_fn(source_points, target_points), dtype=np.float64)
    source_self = np.asarray(cost_fn(source_points, source_points), dtype=np.float64)
    target_self = np.asarray(cost_fn(target_points, target_points), dtype=np.float64)

    ot_cross, cross_iter, cross_converged = _sinkhorn_ot_cost(
        cross, source_weights, target_weights, epsilon, max_iter, tolerance
    )
    ot_source, source_iter, source_converged = _sinkhorn_ot_cost(
        source_self, source_weights, source_weights, epsilon, max_iter, tolerance
    )
    ot_target, target_iter, target_converged = _sinkhorn_ot_cost(
        target_self, target_weights, target_weights, epsilon, max_iter, tolerance
    )
    raw_divergence = ot_cross - 0.5 * (ot_source + ot_target)
    return max(0.0, float(raw_divergence)), {
        "cross_iterations": cross_iter,
        "source_self_iterations": source_iter,
        "target_self_iterations": target_iter,
        "all_converged": bool(cross_converged and source_converged and target_converged),
        "raw_divergence": float(raw_divergence),
    }


def sinkhorn_pushforward_errors(
    pushed_points: np.ndarray,
    target_points: np.ndarray,
    source_weights: np.ndarray,
    target_weights: np.ndarray,
    epsilon: float,
    max_iter: int,
    tolerance: float,
) -> dict[str, object]:
    """Estimate pushforward errors with the paper's independent Sinkhorn solves.

    The paper compares the ray-traced output and target in a common planar
    stereographic chart with squared Euclidean ground cost.  The angle and
    ambient 3-D metrics are retained as supplementary diagnostics, but all
    three are now Sinkhorn divergences rather than linear assignments.
    """
    pushed_u, pushed_v = stereo_north(pushed_points)
    target_u, target_v = stereo_north(target_points)
    pushed_plane = np.column_stack([pushed_u, pushed_v])
    target_plane = np.column_stack([target_u, target_v])

    def angle_cost(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return np.arccos(np.clip(a @ b.T, -1.0, 1.0))

    def xy_squared_cost(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        delta = a[:, None, :] - b[None, :, :]
        return np.sum(delta * delta, axis=2)

    def xyz_squared_cost(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        delta = a[:, None, :] - b[None, :, :]
        return np.sum(delta * delta, axis=2)

    angle_divergence, angle_diag = _sinkhorn_divergence(
        pushed_points, target_points, source_weights, target_weights,
        angle_cost, epsilon, max_iter, tolerance,
    )
    xy_divergence, xy_diag = _sinkhorn_divergence(
        pushed_plane, target_plane, source_weights, target_weights,
        xy_squared_cost, epsilon, max_iter, tolerance,
    )
    xyz_divergence, xyz_diag = _sinkhorn_divergence(
        pushed_points, target_points, source_weights, target_weights,
        xyz_squared_cost, epsilon, max_iter, tolerance,
    )
    return {
        "sinkhorn_angle_cost": angle_divergence,
        "sinkhorn_w2_xy_l2": float(np.sqrt(xy_divergence)),
        "sinkhorn_w2_3d_l2": float(np.sqrt(xyz_divergence)),
        "sinkhorn_angle_diagnostics": angle_diag,
        "sinkhorn_xy_diagnostics": xy_diag,
        "sinkhorn_3d_diagnostics": xyz_diag,
        "epsilon": float(epsilon),
        "max_iter": int(max_iter),
        "tolerance": float(tolerance),
    }


def hard_c_transform_pushforward(
    x: np.ndarray,
    y: np.ndarray,
    g_raw: np.ndarray,
    source_weights: np.ndarray,
    chunk_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return hard target indices, pushed atoms, and accumulated target mass."""
    finite_target = np.isfinite(g_raw)
    if not finite_target.any():
        raise ValueError("hard c-transform requires at least one finite target potential")

    target_indices = np.flatnonzero(finite_target)
    y_supported = y[finite_target]
    g_supported = g_raw[finite_target]
    map_indices = np.empty(len(x), dtype=np.int64)

    for start in range(0, len(x), chunk_size):
        stop = min(start + chunk_size, len(x))
        costs = cost_matrix_chunk(x[start:stop], y_supported)
        map_indices[start:stop] = target_indices[np.argmin(
            costs - g_supported[None, :], axis=1
        )]

    pushed_mass = np.bincount(
        map_indices,
        weights=source_weights,
        minlength=len(y),
    ).astype(np.float64)
    return map_indices, y[map_indices], pushed_mass


def _set_equal_axes_3d(ax, points: np.ndarray) -> None:
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    centre = 0.5 * (mins + maxs)
    half_range = max(0.5 * float(np.max(maxs - mins)), 1e-12)
    ax.set_xlim(centre[0] - half_range, centre[0] + half_range)
    ax.set_ylim(centre[1] - half_range, centre[1] + half_range)
    ax.set_zlim(centre[2] - half_range, centre[2] + half_range)


def save_plots(
    case_dir: Path,
    case_label: str,
    case: dict[str, object],
    x: np.ndarray,
    ref: np.ndarray,
    radii: np.ndarray,
    pushforward_errors: dict[str, object],
) -> tuple[str, str]:
    w2_xy = float(pushforward_errors["sinkhorn_w2_xy_l2"])
    w2_3d = float(pushforward_errors["sinkhorn_w2_3d_l2"])
    kappa = float(case["kappa"])
    cmap = "viridis"

    surface_path = case_dir / f"{case_label}_surface_3d.png"
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")
    scatter = ax.scatter(
        ref[:, 0], ref[:, 1], ref[:, 2],
        c=radii, cmap=cmap, s=18, alpha=0.9, depthshade=True,
    )
    fig.colorbar(scatter, ax=ax, shrink=0.68, pad=0.10, label="Refracter radius R")
    ax.set_title(
        f"{case_label}: Sinkhorn refracter surface\n"
        f"κ={kappa:.5f}; Sinkhorn W₂={w2_xy:.5f}"
    )
    ax.set_xlabel("surface X")
    ax.set_ylabel("surface Y")
    ax.set_zlabel("surface Z")
    _set_equal_axes_3d(ax, ref)
    fig.tight_layout()
    fig.savefig(surface_path, dpi=160, bbox_inches="tight")
    plt.close(fig)

    source_xy_path = case_dir / f"{case_label}_source_xy.png"
    fig, ax = plt.subplots(figsize=(8, 7))
    scatter = ax.scatter(
        x[:, 0], x[:, 1], c=radii, cmap=cmap, s=18, alpha=0.9,
    )
    fig.colorbar(scatter, ax=ax, label="Refracter radius R")
    ax.set_title(
        f"{case_label}: source unit vectors projected to xy\n"
        f"κ={kappa:.5f}; Sinkhorn W₂={w2_xy:.5f}"
    )
    ax.set_xlabel("source unit-vector x")
    ax.set_ylabel("source unit-vector y")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(source_xy_path, dpi=160, bbox_inches="tight")
    plt.close(fig)

    return surface_path.name, source_xy_path.name


def _json_ready(value):
    if isinstance(value, dict):
        return {str(k): _json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(v) for v in value]
    if isinstance(value, np.ndarray):
        return [_json_ready(v) for v in value.tolist()]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def _fmt_range(values: list[float], digits: int = 5) -> str:
    return f"[{values[0]:.{digits}f}, {values[1]:.{digits}f}]"


def write_report(output_dir: Path, metadata: dict[str, object]) -> None:
    report_path = output_dir / "README.md"
    lines = [
        "# Sinkhorn refracter examples",
        "",
        "Generated by `scripts/generate_sinkhorn_refracter_examples.py`.",
        "",
        "## Method",
        "",
        f"- Seed: `{metadata['seed']}`; cases: `{metadata['case_count']}` "
        f"(1 default + {metadata['random_cases']} random); NK: `{metadata['nk']}`.",
        "- Every cloud is sampled uniformly in spherical area: azimuth is uniform and "
        "the cosine of the polar angle is uniform.",
        "- Source and target weights are both discrete uniform measures, with one atom "
        "of mass `1/NK` per point.",
        "- The surface is `Ref_i = 2 R_i x_i`, with `R_i = exp(f_i)` from the "
        "Sinkhorn-divergence-corrected source potential.",
        "- The hard c-transform map uses the raw OT target potential: "
        "`j*(i) = argmin_j [c(x_i,y_j) - g_raw_j]`.",
        f"- Push-forward errors use an independent paper-style Sinkhorn divergence "
        f"with ε=`{metadata['pushforward_epsilon']:.3g}`. The primary planar "
        "metric compares north-pole stereographic coordinates with squared "
        "Euclidean cost and reports `sqrt(S_ε)`; angle and 3-D metrics are "
        "supplementary Sinkhorn diagnostics.",
        "- The NPZ files also include the legacy notebook fields `R`, `Ref`, `gc`, "
        "`fc`, `Refc`, meshgrids, and projected distributions so each case can be "
        "loaded by `notebooks/Benchmark_refraction.ipynb`.",
        "",
        "The new `validate_transport_possible` check was run for every case before "
        "Sinkhorn. `target_fully_covered` means every generated target point lies "
        "inside the requested target theta/phi bounds; `finite_cost_reachable` "
        "means every source-target pair passes the strict finite-cost condition.",
        "",
        "## Sinkhorn push-forward summary",
        "",
        "| Case | κ | source θ° | source φ° | target θ° | target φ° | target covered | finite/reachable | Sε angle | sqrt(Sε) xy | sqrt(Sε) 3-D | hard-map target hit rate |",
        "|---|---:|---|---|---|---|:---:|:---:|---:|---:|---:|---:|",
    ]

    for result in metadata["cases"]:
        bounds = result["bounds"]
        lines.append(
            f"| {result['label']} | {result['kappa']:.6f} | "
            f"{_fmt_range(bounds['source']['theta_deg'], 2)} | "
            f"{_fmt_range(bounds['source']['phi_deg'], 2)} | "
            f"{_fmt_range(bounds['target']['theta_deg'], 2)} | "
            f"{_fmt_range(bounds['target']['phi_deg'], 2)} | "
            f"{str(result['target_coverage']['fully_covered']).lower()} | "
            f"{str(result['finite_cost_reachable']).lower()} | "
            f"{result['sinkhorn_angle_cost']:.8g} | "
            f"{result['sinkhorn_w2_xy_l2']:.8g} | "
            f"{result['sinkhorn_w2_3d_l2']:.8g} | "
            f"{result['hard_map_target_hit_rate']:.4f} |"
        )

    lines.extend(["", "## Per-case plots and diagnostics", ""])
    for result in metadata["cases"]:
        lines.extend([
            f"### {result['label']}",
            "",
            f"Sinkhorn angle cost = **{result['sinkhorn_angle_cost']:.8g}**; "
            f"sqrt(Sε) planar L2 = **{result['sinkhorn_w2_xy_l2']:.8g}**; "
            f"sqrt(Sε) 3-D L2 = **{result['sinkhorn_w2_3d_l2']:.8g}**.",
            "",
            f"- [3D surface coloured by radius]({result['surface_plot']})",
            f"- [Source unit vectors projected to xy, coloured by radius]({result['source_xy_plot']})",
            f"- `target_coverage`: `{result['target_coverage']['inside_count']}/{result['target_coverage']['count']}` points inside bounds",
            f"- `finite_cost_reachable`: `{result['finite_cost_reachable']}`; maximum κ·(x·y) = `{result['max_kappa_dot']:.8g}`",
            f"- radius range/mean: `{result['radius_min']:.8g}` / `{result['radius_max']:.8g}` / `{result['radius_mean']:.8g}`",
            f"- raw c-transform residual max `|f_raw - c-transform(g_raw)|`: `{result['c_transform_residual_max']:.8g}`",
            f"- hard-map target mass hit rate: `{result['hard_map_target_hit_rate']:.4f}`; unique target points hit: `{result['hard_map_unique_targets']}/{result['nk']}`",
            f"- Sinkhorn divergence total cost: `{result['total_cost']:.8g}`; "
            f"pushforward solves converged: `{result['pushforward_sinkhorn_converged']}`; "
            f"runtime: `{result['runtime_seconds']:.3f}` s",
            "",
        ])

    report_path.write_text("\n".join(lines), encoding="utf-8")
    (output_dir / "metadata.json").write_text(
        json.dumps(_json_ready(metadata), indent=2), encoding="utf-8"
    )


def run_case(case_index: int, case: dict[str, object], args: argparse.Namespace,
             output_dir: Path) -> dict[str, object]:
    label = "case_%02d_%s" % (
        case_index,
        "default" if case_index == 0 else "random",
    )
    case_dir = output_dir / label
    case_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()

    set_kappa(float(case["kappa"]))
    x = gen_spherical_patch(args.nk, case["source"], skip=0)
    y = gen_spherical_patch(args.nk, case["target"], skip=17 + case_index * 101)

    source_coverage = patch_coverage(x, case["source"])
    target_coverage = patch_coverage(y, case["target"])
    if not source_coverage["fully_covered"] or not target_coverage["fully_covered"]:
        raise RuntimeError(f"{label}: generated point fell outside its patch bounds")

    # This is the repository's new strict finite-cost reachability check.
    validate_transport_possible(x, y, chunk_size=args.chunk_size)
    max_kappa_dot = float(float(case["kappa"]) * np.max(x @ y.T))

    p = np.full(args.nk, 1.0 / args.nk, dtype=np.float64)
    q = np.full(args.nk, 1.0 / args.nk, dtype=np.float64)
    sinkhorn = run_sinkhorn_divergence(
        x, y, p, q, chunk_size=args.chunk_size, verbose=False,
    )

    f = sinkhorn["f"]
    g_raw = sinkhorn["g_raw"]
    radii = np.exp(f)
    ref = 2.0 * x * radii[:, None]

    # Check the raw pair used to define the hard map.
    f_from_g_raw = c_transform_gc(x, y, g_raw, chunk_size=args.chunk_size)
    supported = np.isfinite(f_from_g_raw) & np.isfinite(sinkhorn["f_raw"])
    residual = np.abs(sinkhorn["f_raw"][supported] - f_from_g_raw[supported])

    map_indices, pushed_points, pushed_mass = hard_c_transform_pushforward(
        x, y, g_raw, p, args.chunk_size,
    )
    pushforward_errors = sinkhorn_pushforward_errors(
        pushed_points,
        y,
        p,
        q,
        epsilon=args.pushforward_eps,
        max_iter=args.pushforward_max_iter,
        tolerance=args.pushforward_tol,
    )
    pushforward_converged = bool(
        pushforward_errors["sinkhorn_angle_diagnostics"]["all_converged"]
        and pushforward_errors["sinkhorn_xy_diagnostics"]["all_converged"]
        and pushforward_errors["sinkhorn_3d_diagnostics"]["all_converged"]
    )
    unique_targets = int(np.count_nonzero(pushed_mass > 0.0))
    plots = save_plots(
        case_dir, label, case, x, ref, radii, pushforward_errors,
    )

    # The notebook's historical NPZ schema stores all derived arrays needed by
    # its figures.  Keep the richer example diagnostics alongside those fields.
    grid_res = 256
    grid_side = np.linspace(-0.6, 0.6, grid_res)
    UU, VV = np.meshgrid(grid_side, grid_side, indexing="ij")
    N2 = UU * UU + VV * VV
    denom = 1.0 + N2
    x_grid = np.stack([
        2.0 * UU / denom,
        2.0 * VV / denom,
        (1.0 - N2) / denom,
    ], axis=-1).reshape(-1, 3)
    X_MeshGrid = patch_indicator(x_grid, case["source"]).reshape(grid_res, grid_res)
    Y_MeshGrid = patch_indicator(x_grid, case["target"]).reshape(grid_res, grid_res)
    gc = c_transform_gc(x, y, g_raw, chunk_size=args.chunk_size)
    fc = c_transform_fc(x, y, sinkhorn["f_raw"], chunk_size=args.chunk_size)
    Refc = 2.0 * x * np.exp(gc)[:, None]
    u_y, v_y = stereo_north(y)
    Y_projected = np.column_stack([u_y, v_y, q])
    u_push, v_push = stereo_north(pushed_points)
    Y_Pushed_projected = np.column_stack([
        u_push,
        v_push,
        patch_indicator(pushed_points, case["target"]),
    ])

    np.savez_compressed(
        case_dir / f"{label}.npz",
        x=x,
        y=y,
        x_s=np.zeros((0, 3), dtype=np.float64),
        y_s=np.zeros((0, 3), dtype=np.float64),
        p=p,
        q=q,
        f_raw=sinkhorn["f_raw"],
        g_raw=g_raw,
        f=f,
        g=sinkhorn["g"],
        f_id=sinkhorn["f_id"],
        g_id=sinkhorn["g_id"],
        radii=radii,
        R=radii,
        Ref=ref,
        gc=gc,
        fc=fc,
        Refc=Refc,
        X_MeshGrid=X_MeshGrid,
        Y_MeshGrid=Y_MeshGrid,
        grid_side=grid_side,
        Y_projected=Y_projected,
        Y_Pushed_projected=Y_Pushed_projected,
        y_push_3d=pushed_points,
        hard_map_indices=map_indices,
        pushed_points=pushed_points,
        pushed_mass=pushed_mass,
        kappa=np.float64(case["kappa"]),
        src_density="uniform",
        tgt_density="uniform",
        pushforward_sinkhorn_epsilon=np.float64(args.pushforward_eps),
        sinkhorn_angle_cost=np.float64(pushforward_errors["sinkhorn_angle_cost"]),
        sinkhorn_w2_xy_l2=np.float64(pushforward_errors["sinkhorn_w2_xy_l2"]),
        sinkhorn_w2_3d_l2=np.float64(pushforward_errors["sinkhorn_w2_3d_l2"]),
        pushforward_sinkhorn_converged=np.bool_(pushforward_converged),
        pushforward_sinkhorn_max_iter=np.int64(args.pushforward_max_iter),
        pushforward_sinkhorn_tolerance=np.float64(args.pushforward_tol),
    )

    result = {
        "label": label,
        "index": case_index,
        "kappa": float(case["kappa"]),
        "bounds": case,
        "nk": args.nk,
        "source_coverage": source_coverage,
        "target_coverage": target_coverage,
        "finite_cost_reachable": True,
        "max_kappa_dot": max_kappa_dot,
        "radius_min": float(radii.min()),
        "radius_max": float(radii.max()),
        "radius_mean": float(radii.mean()),
        "c_transform_residual_max": float(residual.max()) if len(residual) else float("nan"),
        "hard_map_unique_targets": unique_targets,
        "hard_map_target_hit_rate": float(unique_targets / args.nk),
        "pushed_mass_min": float(pushed_mass.min()),
        "pushed_mass_max": float(pushed_mass.max()),
        "total_cost": float(sinkhorn["total_cost"]),
        "sinkhorn_angle_cost": pushforward_errors["sinkhorn_angle_cost"],
        "sinkhorn_w2_xy_l2": pushforward_errors["sinkhorn_w2_xy_l2"],
        "sinkhorn_w2_3d_l2": pushforward_errors["sinkhorn_w2_3d_l2"],
        "pushforward_sinkhorn_converged": pushforward_converged,
        "pushforward_sinkhorn_diagnostics": {
            "angle": pushforward_errors["sinkhorn_angle_diagnostics"],
            "xy": pushforward_errors["sinkhorn_xy_diagnostics"],
            "3d": pushforward_errors["sinkhorn_3d_diagnostics"],
        },
        "surface_plot": f"{label}/{plots[0]}",
        "source_xy_plot": f"{label}/{plots[1]}",
        "npz": f"{label}/{label}.npz",
        "runtime_seconds": time.perf_counter() - started,
    }
    print(
        f"{label}: kappa={result['kappa']:.6f}, "
        f"S_angle={result['sinkhorn_angle_cost']:.8g}, "
        f"sqrt(S)_xy={result['sinkhorn_w2_xy_l2']:.8g}, "
        f"sqrt(S)_3d={result['sinkhorn_w2_3d_l2']:.8g}, "
        f"target={target_coverage['inside_count']}/{target_coverage['count']}, "
        f"unique hard targets={unique_targets}/{args.nk}, "
        f"time={result['runtime_seconds']:.2f}s"
    )
    return result


def main() -> None:
    args = parse_args()
    if (
        args.nk <= 0
        or args.random_cases < 0
        or args.chunk_size <= 0
        or args.pushforward_eps <= 0.0
        or args.pushforward_max_iter <= 0
        or args.pushforward_tol <= 0.0
    ):
        raise ValueError(
            "--nk, --chunk-size, --pushforward-eps, --pushforward-max-iter, "
            "and --pushforward-tol must be positive; --random-cases cannot be negative"
        )

    output_dir = Path(args.output_dir) if args.output_dir else (
        REPO_ROOT / "results" / "sinkhorn_refracter_examples"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    cases = [{
        "kappa": 0.6,
        "source": DEFAULT_BOUNDS["source"],
        "target": DEFAULT_BOUNDS["target"],
    }]
    cases.extend(random_case(rng) for _ in range(args.random_cases))

    print(
        f"Generating {len(cases)} Sinkhorn refracter cases with NK={args.nk}, "
        f"chunk_size={args.chunk_size}, seed={args.seed}"
    )
    results = [run_case(i, case, args, output_dir) for i, case in enumerate(cases)]
    metadata = {
        "script": str(Path(__file__).relative_to(REPO_ROOT)),
        "seed": args.seed,
        "nk": args.nk,
        "random_cases": args.random_cases,
        "case_count": len(cases),
        "chunk_size": args.chunk_size,
        "pushforward_epsilon": args.pushforward_eps,
        "pushforward_max_iter": args.pushforward_max_iter,
        "pushforward_tolerance": args.pushforward_tol,
        "solver": "refracter.sinkhorn.run_sinkhorn_divergence",
        "pushforward_solver": "paper-style log-domain Sinkhorn divergence",
        "cost": "-log(1 - kappa * (x dot y))",
        "reachability_check": "refracter.cost.validate_transport_possible",
        "ground_metric": "arccos(y dot y')",
        "results": results,
        "cases": results,
        "final_kappa_in_process": get_kappa(),
    }
    write_report(output_dir, metadata)
    print(f"\nWrote report: {output_dir / 'README.md'}")
    print(f"Wrote metadata: {output_dir / 'metadata.json'}")


if __name__ == "__main__":
    main()
