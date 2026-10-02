"""
generate_results.py — Refractor Sinkhorn benchmark, NK=1600

Generates a 1600-point Halton quasi-Monte Carlo cloud on spherical patches,
runs the repository's multi-scale Sinkhorn-divergence pipeline,
and saves a comprehensive results bundle to results/results_refraction_NK{NK}.npz.

Setup:
  * κ = 0.6  (refraction cost:  c(x,y) = -log(1 - 0.6·(x·y)))
  * Both source and target are on the upper hemisphere (spherical patches)
  * Source  Ω  : θ ∈ [π/12, π/3],  φ ∈ [π/12, π/4]
  * Target  Ω* : θ ∈ [π/10, π/5],  φ ∈ [π/10, π/5]
  * North-pole stereographic projection is used for both source and target

Run:  python scripts/generate_results.py [NK]
      (default NK=1600)
"""

import os, sys, time
import numpy as np

# Ensure we can import the refracter package from the repo root
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, REPO_ROOT)

# The cost parameter is set before running the benchmark.
from refracter.cost import set_kappa, cost_matrix_chunk
from refracter.distributions import P_refraction_patch, Q_refraction_patch
from refracter.sinkhorn import run_sinkhorn_divergence
from refracter.build import c_transform_gc, c_transform_fc
from refracter.distributions import stereo_north, stereo_north_inverse
from refracter.qmc import gen_spherical_patch

SRC_THETA = (np.pi / 12, np.pi / 3)
SRC_PHI = (np.pi / 12, np.pi / 4)
TGT_THETA = (np.pi / 10, np.pi / 5)
TGT_PHI = (np.pi / 10, np.pi / 5)
DEFAULT_CHUNK_SIZE = 512
REGULARIZATION_MULTIPLIER = 8


def build_refraction_artifacts(x, y, p, q, f_raw, g_raw, f, g, P, Q,
                               chunk=512, grid_res=256, grid_limit=0.6,
                               kappa=0.6):
    """Build the arrays shared by the result generator and its NPZ helpers."""
    set_kappa(kappa)
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    mask_p = p > 0
    mask_q = q > 0

    R = np.exp(f)
    Ref = 2.0 * x * R[:, None]

    g_raw_masked = np.where(mask_q, g_raw, -1e300)
    f_raw_masked = np.where(mask_p, f_raw, -1e300)
    gc = c_transform_gc(x, y, g_raw_masked, chunk)
    fc = c_transform_fc(x, y, f_raw_masked, chunk)
    Refc = 2.0 * x * np.exp(gc)[:, None]

    grid_side = np.linspace(-grid_limit, grid_limit, grid_res)
    UU, VV = np.meshgrid(grid_side, grid_side, indexing="ij")
    x_grid = stereo_north_inverse(UU, VV).reshape(-1, 3)
    X_MeshGrid = P(x_grid).reshape(grid_res, grid_res)
    Y_MeshGrid = Q(x_grid).reshape(grid_res, grid_res)

    u_y, v_y = stereo_north(y[mask_q])
    Y_projected = np.column_stack([u_y, v_y, q[mask_q]])

    src_idx = np.flatnonzero(mask_p)
    tgt_idx = np.flatnonzero(mask_q)
    if not len(src_idx) or not len(tgt_idx):
        raise ValueError("source and target distributions must have positive support")
    y_tgt = y[tgt_idx]
    g_tgt = np.asarray(g_raw)[tgt_idx]
    hard_map_indices = np.full(len(x), -1, dtype=np.int64)
    for i0 in range(0, len(src_idx), chunk):
        rows = src_idx[i0:i0 + chunk]
        C_block = cost_matrix_chunk(x[rows], y_tgt)
        j_star = np.argmin(C_block - g_tgt[None, :], axis=1)
        hard_map_indices[rows] = tgt_idx[j_star]
    pushed_y = y[hard_map_indices[src_idx]]
    pushed_mass = np.bincount(
        hard_map_indices[src_idx], weights=p[src_idx], minlength=len(y)
    ).astype(np.float64)

    u_push, v_push = stereo_north(pushed_y)
    q_push = Q(pushed_y)
    Y_Pushed_projected = np.column_stack([u_push, v_push, q_push])
    return {
        "R": R,
        "Ref": Ref,
        "gc": gc,
        "fc": fc,
        "Refc": Refc,
        "X_MeshGrid": X_MeshGrid,
        "Y_MeshGrid": Y_MeshGrid,
        "grid_side": grid_side,
        "Y_projected": Y_projected,
        "Y_Pushed_projected": Y_Pushed_projected,
        "y_push_3d": pushed_y,
        "hard_map_indices": hard_map_indices,
        "pushed_mass": pushed_mass,
    }


def save_refraction_npz(path, *, x, y, p, q, f_raw, g_raw, f_id, g_id,
                        f, g, artifacts, kappa=0.6, metadata=None):
    """Save the stable NPZ schema consumed by Benchmark_refraction.ipynb."""
    payload = {
        "x": x,
        "y": y,
        "x_s": np.zeros((0, 3), dtype=np.float64),
        "y_s": np.zeros((0, 3), dtype=np.float64),
        "p": p,
        "q": q,
        "f_raw": f_raw,
        "g_raw": g_raw,
        "f_id": f_id,
        "g_id": g_id,
        "f": f,
        "g": g,
        **artifacts,
        "kappa": np.float64(kappa),
    }
    if metadata:
        overlap = payload.keys() & metadata.keys()
        if overlap:
            raise ValueError(f"metadata cannot replace standard NPZ fields: {sorted(overlap)}")
        payload.update(metadata)
    np.savez(path, **payload)


def main(nk=None):
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    # ---------------------------------------------------------------------------
    # Problem setup
    # ---------------------------------------------------------------------------

    NK = int(nk if nk is not None else (sys.argv[1] if len(sys.argv) > 1 else 1600))
    chunk = DEFAULT_CHUNK_SIZE
    set_kappa(0.6)

    print(f"=== Refraction benchmark (κ=0.6, NK={NK}) ===")
    t_start = time.time()

    print("Generating Halton QMC clouds on spherical patches...")
    x = gen_spherical_patch(NK, *SRC_THETA, *SRC_PHI, skip=0)
    y = gen_spherical_patch(NK, *TGT_THETA, *TGT_PHI, skip=0)

    print(f"  Source: {len(x)} pts, z ∈ [{x[:,2].min():.4f}, {x[:,2].max():.4f}]")
    print(f"  Target: {len(y)} pts, z ∈ [{y[:,2].min():.4f}, {y[:,2].max():.4f}]")

    # Densities — all points are inside the patch by construction → uniform
    # you can import a different distribution
    p_raw = P_refraction_patch(x)
    q_raw = Q_refraction_patch(y)
    print(f"  Source support: {int(p_raw.sum()+0.5)} / {NK}")
    print(f"  Target support: {int(q_raw.sum()+0.5)} / {NK}")

    p = p_raw / p_raw.sum()
    q = q_raw / q_raw.sum()
    k_final = REGULARIZATION_MULTIPLIER * int(np.floor(np.sqrt(NK)))
    print(f"  k_final={k_final}")
    solution = run_sinkhorn_divergence(
        x, y, p, q, chunk_size=chunk, verbose=True,
    )
    f_raw = solution["f_raw"]
    g_raw = solution["g_raw"]
    f_id = solution["f_id"]
    g_id = solution["g_id"]
    f = solution["f"]
    g = solution["g"]

    # ---------------------------------------------------------------------------
    # Steps 5–9 — Build the notebook-compatible refractor result fields
    # ---------------------------------------------------------------------------
    print("\nBuilding refractor, c-transforms, density grids, and push-forward...")
    t0 = time.time()
    artifacts = build_refraction_artifacts(
        x, y, p, q, f_raw, g_raw, f, g,
        P_refraction_patch, Q_refraction_patch, chunk=chunk, kappa=0.6,
    )
    R = artifacts["R"]
    Ref = artifacts["Ref"]
    mask_p = p > 0
    print(f"  done ({time.time()-t0:.2f}s)")
    print(f"R (supported): min={R[mask_p].min():.4f}, max={R[mask_p].max():.4f}, mean={R[mask_p].mean():.4f}")

    # ---------------------------------------------------------------------------
    # Save
    # ---------------------------------------------------------------------------
    results_dir = os.path.join(REPO_ROOT, "results")
    os.makedirs(results_dir, exist_ok=True)
    out = os.path.join(results_dir, f"results_refraction_NK{NK}.npz")
    save_refraction_npz(
        out, x=x, y=y, p=p, q=q,
        f_raw=f_raw, g_raw=g_raw, f_id=f_id, g_id=g_id, f=f, g=g,
        artifacts=artifacts, kappa=0.6,
    )
    sz = os.path.getsize(out) / 1024
    print(f"\nSaved {out}  ({sz:.0f} KB)   total time {time.time()-t_start:.1f}s")

    # ---------------------------------------------------------------------------
    # Figure — Refractor surface (coloured by R)
    # ---------------------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    fig = plt.figure(figsize=(7, 6))
    ax  = fig.add_subplot(111, projection='3d')
    sc  = ax.scatter(Ref[:, 0], Ref[:, 1], Ref[:, 2],
                     c=R, cmap='viridis', s=15,
                     vmin=R.min(), vmax=R.max(), depthshade=True)
    plt.colorbar(sc, ax=ax, label='R (refractor radius)', shrink=0.65)
    ax.set_title(f'Refractor Surface — Python  (κ=0.6, NK={NK}, k_final={k_final})')
    ax.set_xlabel('X');  ax.set_ylabel('Y');  ax.set_zlabel('Z')
    fig.tight_layout()
    fig_path = os.path.join(REPO_ROOT, "figures", f"fig_refractor_3d_NK{NK}.png")
    os.makedirs(os.path.dirname(fig_path), exist_ok=True)
    fig.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f"Saved {fig_path}")


if __name__ == "__main__":
    main()
