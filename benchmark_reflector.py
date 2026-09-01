"""
benchmark_reflector.py

Reflector benchmark pipeline (κ = 1.0, SquareToCircle / SquareToTwoGaussSide).

This module provides run_benchmark(), which replicates the full C++ reflector pipeline:

  1. Load QMC point clouds.
  2. Evaluate source/target density functions P and Q.
  3. Run the Sinkhorn divergence algorithm.
  4. Build the reflector surface (R, Ref).
  5. Compute c-transforms (gc, fc, Refc).
  6. Build the regular grid and compute the reflector on it.
  7. Ray-trace the push-forward.
  8. Save all outputs (same filenames as C++).
  9. Print total cost and timing.
 10. Return a dict of all results.

Usage
-----
    from benchmark_reflector import run_benchmark
    results = run_benchmark(benchmark='SquareToCircle', output_dir='Output')

or from the command line:

    python benchmark_reflector.py --benchmark SquareToCircle --output_dir Output
"""

import os
import time
import numpy as np

from refracter.qmc import load_main_cloud, load_small_cloud, load_push_cloud
from refracter.distributions import BENCHMARKS, stereo_south
from refracter.cost import set_kappa
from refracter.sinkhorn import _run_sinkhorn_divergence_inner
from refracter.build import (
    build_reflector,
    c_transform_gc,
    c_transform_fc,
    build_regular_grid,
    reflector_on_regular_grid,
)
from refracter.pushforward import ray_trace


# ---------------------------------------------------------------------------
# Save helpers (matching C++ output format)
# ---------------------------------------------------------------------------

def _save_vector(path: str, arr: np.ndarray):
    """Save a 1-D array in C++ format: first line N, then one value per line."""
    arr = np.asarray(arr, dtype=np.float64)
    with open(path, "w") as fh:
        fh.write(f"{len(arr)}  \n")
        for v in arr:
            fh.write(f"{v:.18e} \n")
        fh.write("\n\n")


def _save_matrix(path: str, arr: np.ndarray):
    """Save a 2-D array in C++ format: first line N, then space-separated rows."""
    arr = np.asarray(arr, dtype=np.float64)
    n_rows, n_cols = arr.shape
    with open(path, "w") as fh:
        fh.write(f"{n_rows}  \n")
        for row in arr:
            fh.write(" ".join(f"{v:.18e}" for v in row) + "\n")
        fh.write("\n\n")


def _save_projected(path: str, rows: np.ndarray):
    """Save projected push-forward: u v density, one row per valid ray."""
    with open(path, "w") as fh:
        for row in rows:
            fh.write(f"{row[0]:.18e} {row[1]:.18e} {row[2]:.18e} \n")


def _save_grid(path: str, grid: np.ndarray):
    """Save a C++ density mesh: one space-separated row per grid row."""
    with open(path, "w") as fh:
        for row in np.asarray(grid, dtype=np.float64):
            fh.write(" ".join(f"{v:.18e}" for v in row) + "\n")


def _save_discrete_mesh(path: str, side: np.ndarray, density_fn, upper: bool):
    """Save supported mesh points as X Y density*Jacobian, like C++."""
    with open(path, "w") as fh:
        for X in side:
            for Y in side:
                N2 = X * X + Y * Y
                denom = 1.0 + N2
                z = (1.0 - N2) / denom
                if not upper:
                    z = -z
                point = np.array([[2.0 * X / denom, 2.0 * Y / denom, z]])
                density = float(density_fn(point)[0])
                if density == 0.0:
                    continue
                jac = 4.0 / (denom * denom)
                fh.write(f"{X:.18e} {Y:.18e} {density * jac:.18e} \n")


# ---------------------------------------------------------------------------
# Main benchmark function
# ---------------------------------------------------------------------------

def run_benchmark(benchmark: str = 'SquareToCircle',
                  output_dir: str = None,
                  chunk_size: int = 512,
                  verbose: bool = True) -> dict:
    """Run the full reflector benchmark pipeline.

    Parameters
    ----------
    benchmark : str
        Name of the benchmark.  One of 'SquareToCircle' or
        'SquareToTwoGaussSide'.
    output_dir : str or None
        Directory to save output files.  If None a timestamped directory is
        created in the current working directory (same behaviour as C++).
    chunk_size : int
        Block size for chunked matrix operations.
    verbose : bool
        Print progress information.

    Returns
    -------
    dict with keys:
        x, y           : main QMC point clouds
        p, q           : normalised source/target weights
        f, g           : Kantorovich potentials (corrected)
        f_id, g_id     : identity potentials
        R, Ref         : reflector scale and surface points
        gc, fc         : c-transforms
        Refc           : reflector from c-transform of g
        x_regular      : regular grid points
        Regular_side   : regular grid axis
        f_regular      : reflector potential on regular grid
        Ref_regular    : reflector surface on regular grid
        push_result    : (K, 3) array of (u, v, density) for push-forward
        total_cost     : float
        elapsed        : float, total wall-clock time in seconds
        output_dir     : str, path to output directory
    """
    t0 = time.perf_counter()

    # --- Benchmark configuration ---
    if benchmark not in BENCHMARKS:
        raise ValueError(
            f"Unknown benchmark '{benchmark}'.  "
            f"Available: {list(BENCHMARKS.keys())}"
        )
    cfg = BENCHMARKS[benchmark]
    if benchmark == "Refraction":
        raise ValueError(
            "The main QMC cloud and ray-tracing pipeline are reflector-only; "
            "use scripts/generate_results.py for the refraction patch benchmark."
        )
    set_kappa(cfg.get("kappa", 1.0))
    P_func = cfg["P"]
    Q_func = cfg["Q"]
    testname = cfg["testname"]

    if verbose:
        print(f"=== Benchmark: {benchmark} ({testname}) ===\n")

    # --- Output directory ---
    if output_dir is None:
        ts = time.strftime("_%Y_%m_%d_%H-%M-%S")
        output_dir = f"Output{ts}_{testname}"
    os.makedirs(output_dir, exist_ok=True)

    if verbose:
        print(f"Output directory: {output_dir}\n")

    # --- 1. Load QMC clouds ---
    t1 = time.perf_counter()
    if verbose:
        print("Loading QMC point clouds...")
    x, y = load_main_cloud()
    x_s, y_s = load_small_cloud()
    push_cloud = load_push_cloud()
    if verbose:
        print(f"  Main cloud: x={x.shape}, y={y.shape}")
        print(f"  Small cloud: x_s={x_s.shape}, y_s={y_s.shape}")
        print(f"  Push cloud: {push_cloud.shape}")
        print(f"  Load time: {time.perf_counter()-t1:.2f}s\n")

    # --- 2. Evaluate densities ---
    t2 = time.perf_counter()
    if verbose:
        print("Evaluating source/target densities...")

    p_raw = P_func(x)
    q_raw = Q_func(y)

    p = p_raw / p_raw.sum()
    q = q_raw / q_raw.sum()

    # Small-grid densities
    p_s_raw = P_func(x_s)
    q_s_raw = Q_func(y_s)

    # Guard: if all zeros, fall back to uniform
    if p_s_raw.sum() == 0:
        p_s_raw = np.ones(len(x_s))
    if q_s_raw.sum() == 0:
        q_s_raw = np.ones(len(y_s))

    if verbose:
        print(f"  Source support: {int((p_raw > 0).sum())} / {len(p_raw)}")
        print(f"  Target support: {int((q_raw > 0).sum())} / {len(q_raw)}")
        print(f"  Density eval time: {time.perf_counter()-t2:.2f}s\n")

    # --- 3. Sinkhorn divergence ---
    t3 = time.perf_counter()
    if verbose:
        print("Running Sinkhorn divergence algorithm...")

    sink_result = _run_sinkhorn_divergence_inner(
        x, y, p_raw, q_raw,
        x_s, y_s, p_s_raw, q_s_raw,
        chunk_size=chunk_size,
        verbose=verbose,
    )
    f     = sink_result["f"]
    g     = sink_result["g"]
    f_id  = sink_result["f_id"]
    g_id  = sink_result["g_id"]

    if verbose:
        print(f"  Sinkhorn time: {time.perf_counter()-t3:.2f}s\n")

    # --- 4. Build reflector ---
    t4 = time.perf_counter()
    if verbose:
        print("Building reflector surface...")

    R, Ref = build_reflector(x, f, f_id)
    # Note: build_reflector treats f as the corrected potential (already has
    # f_id subtracted by _run_sinkhorn_divergence_inner), so R = exp(f).

    if verbose:
        print(f"  R range: [{R.min():.4e}, {R.max():.4e}]")
        print(f"  Reflector build time: {time.perf_counter()-t4:.2f}s\n")

    # --- 5. C-transforms ---
    t5 = time.perf_counter()
    if verbose:
        print("Computing c-transforms...")

    gc = c_transform_gc(x, y, g, chunk_size=chunk_size)
    fc = c_transform_fc(x, y, f, chunk_size=chunk_size)
    Refc = 2.0 * x * np.exp(gc)[:, None]   # reflector from gc

    dif_f = f - gc
    dif_g = g - fc

    if verbose:
        print(f"  |f - gc|_inf = {np.max(np.abs(dif_f[p>0])):.4e}")
        print(f"  |g - fc|_inf = {np.max(np.abs(dif_g[q>0])):.4e}")
        print(f"  C-transform time: {time.perf_counter()-t5:.2f}s\n")

    # --- 6. Regular grid and reflector ---
    t6 = time.perf_counter()
    if verbose:
        print("Building regular grid and computing reflector on it...")

    x_regular, Regular_side = build_regular_grid(final_grid_res=1025)
    f_regular, Ref_regular = reflector_on_regular_grid(
        x_regular, y, g, chunk_size=chunk_size
    )

    if verbose:
        print(f"  Regular grid: {x_regular.shape}")
        print(f"  Regular grid time: {time.perf_counter()-t6:.2f}s\n")

    # --- 7. Push-forward / ray tracing ---
    t7 = time.perf_counter()
    if verbose:
        print("Ray tracing push-forward...")

    push_result = ray_trace(push_cloud, Ref_regular, Regular_side, P_func)
    own_push_result = ray_trace(x, Ref_regular, Regular_side, P_func)

    if verbose:
        print(f"  Valid rays: {len(push_result)} / {len(push_cloud)}")
        print(f"  Ray trace time: {time.perf_counter()-t7:.2f}s\n")

    # --- 8. Total cost ---
    # C++ Compute_TotalCost() evaluates the same corrected potentials:
    # f_sinkhorn - f_id and g_sinkhorn - g_id.
    total_cost = float(
        np.sum(p[p > 0] * f[p > 0]) +
        np.sum(q[q > 0] * g[q > 0])
    )

    elapsed = time.perf_counter() - t0
    if verbose:
        print(f"Total cost: {total_cost:.6e}")
        print(f"Total elapsed: {elapsed:.2f}s\n")

    # --- 9. Save outputs ---
    if verbose:
        print("Saving outputs...")

    _save_matrix(os.path.join(output_dir, "x_MY.txt"),    x)
    _save_matrix(os.path.join(output_dir, "y_MY.txt"),    y)
    _save_vector(os.path.join(output_dir, "p_MY.txt"),    p)
    _save_vector(os.path.join(output_dir, "q_MY.txt"),    q)
    _save_vector(os.path.join(output_dir, "f_MY.txt"),    f)
    _save_vector(os.path.join(output_dir, "g_MY.txt"),    g)
    _save_vector(os.path.join(output_dir, "F_MY.txt"),    np.ones(len(x)))
    _save_vector(os.path.join(output_dir, "G_MY.txt"),    np.ones(len(y)))
    _save_vector(os.path.join(output_dir, "f_id_MY.txt"), f_id)
    _save_vector(os.path.join(output_dir, "g_id_MY.txt"), g_id)
    _save_vector(os.path.join(output_dir, "F_id_MY.txt"), np.ones(len(x)))
    _save_vector(os.path.join(output_dir, "G_id_MY.txt"), np.ones(len(y)))
    _save_vector(os.path.join(output_dir, "R_MY.txt"),    R)
    _save_matrix(os.path.join(output_dir, "Ref_MY.txt"),  Ref)
    _save_matrix(os.path.join(output_dir, "Refc_MY.txt"), Refc)
    _save_vector(os.path.join(output_dir, "gc_MY.txt"),   gc)
    _save_vector(os.path.join(output_dir, "fc_MY.txt"),   fc)
    _save_vector(os.path.join(output_dir, "dif_fc_MY.txt"), dif_f)
    _save_vector(os.path.join(output_dir, "dif_gc_MY.txt"), dif_g)

    # Target distribution support (Y_projected)
    with open(os.path.join(output_dir, "Y_projected.txt"), "w") as fh:
        for j in range(len(y)):
            if q[j] == 0.0:
                continue
            u, v = y[j, 0] / (1.0 - y[j, 2]), y[j, 1] / (1.0 - y[j, 2])
            fh.write(f"{u:.18e} {v:.18e} {float(Q_func(y[j:j+1])):.18e}\n")

    # Push-forward (Y_Pushed_projected)
    _save_projected(os.path.join(output_dir, "Y_Pushed_projected.txt"), push_result)
    _save_projected(
        os.path.join(output_dir, "Y_Pushed_projected_owndisc.txt"),
        own_push_result,
    )

    # Source support in the north-pole projection, matching C++ X_projected.
    x_supported = x[p_raw > 0]
    u_x = x_supported[:, 0] / (1.0 + x_supported[:, 2])
    v_x = x_supported[:, 1] / (1.0 + x_supported[:, 2])
    _save_projected(
        os.path.join(output_dir, "X_projected.txt"),
        np.column_stack([u_x, v_x, P_func(x_supported)]),
    )

    # C++ main.cpp emits 513x513 stereographic meshes with the sphere-area
    # Jacobian included, plus supported-point discrete versions.
    mesh_res = 512 + 1
    mesh_side = np.linspace(-0.6, 0.6, mesh_res)
    MM_U, MM_V = np.meshgrid(mesh_side, mesh_side, indexing="ij")
    MM_N2 = MM_U * MM_U + MM_V * MM_V
    MM_D = 1.0 + MM_N2
    x_mesh = np.stack([2.0 * MM_U / MM_D, 2.0 * MM_V / MM_D,
                       (1.0 - MM_N2) / MM_D], axis=-1).reshape(-1, 3)
    y_mesh = x_mesh.copy()
    y_mesh[:, 2] *= -1.0
    jac = 4.0 / (MM_D * MM_D)
    X_MeshGrid = P_func(x_mesh).reshape(mesh_res, mesh_res) * jac
    Y_MeshGrid = Q_func(y_mesh).reshape(mesh_res, mesh_res) * jac
    _save_grid(os.path.join(output_dir, "X_MeshGrid.txt"), X_MeshGrid)
    _save_grid(os.path.join(output_dir, "Y_MeshGrid.txt"), Y_MeshGrid)
    _save_discrete_mesh(
        os.path.join(output_dir, "X_DiscreteMesh.txt"), mesh_side, P_func, True
    )
    _save_discrete_mesh(
        os.path.join(output_dir, "Y_DiscreteMesh.txt"), mesh_side, Q_func, False
    )

    if verbose:
        print(f"  Saved {len(os.listdir(output_dir))} files to {output_dir}")
        print(f"\n=== Done. Total time: {elapsed:.2f}s ===")

    return {
        "x": x,
        "y": y,
        "p": p,
        "q": q,
        "f": f,
        "g": g,
        "f_id": f_id,
        "g_id": g_id,
        "R": R,
        "Ref": Ref,
        "gc": gc,
        "fc": fc,
        "Refc": Refc,
        "dif_f": dif_f,
        "dif_g": dif_g,
        "x_regular": x_regular,
        "Regular_side": Regular_side,
        "f_regular": f_regular,
        "Ref_regular": Ref_regular,
        "push_result": push_result,
        "own_push_result": own_push_result,
        "X_MeshGrid": X_MeshGrid,
        "Y_MeshGrid": Y_MeshGrid,
        "total_cost": total_cost,
        "elapsed": elapsed,
        "output_dir": output_dir,
        "benchmark": benchmark,
        "testname": testname,
    }


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Run the reflector benchmark (Python rewrite of C++ code)."
    )
    parser.add_argument(
        "--benchmark", "-b",
        default="SquareToCircle",
        choices=list(BENCHMARKS.keys()),
        help="Which benchmark to run (default: SquareToCircle).",
    )
    parser.add_argument(
        "--output_dir", "-o",
        default=None,
        help="Output directory (default: auto-generated timestamped directory).",
    )
    parser.add_argument(
        "--chunk_size", "-c",
        type=int,
        default=512,
        help="Chunk size for matrix operations (default: 512).",
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress progress output.",
    )

    args = parser.parse_args()

    run_benchmark(
        benchmark=args.benchmark,
        output_dir=args.output_dir,
        chunk_size=args.chunk_size,
        verbose=not args.quiet,
    )
