

import os, sys, time
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, REPO_ROOT)

# This helper always runs the refraction cost with kappa=0.6.
from refracter.cost import set_kappa
set_kappa(0.6)

from refracter.distributions import (
    make_patch_uniform,
    make_patch_gaussian,
    make_patch_donut,
    make_patch_cross,
)
from refracter.sinkhorn import solve_cold_start_sinkhorn
from generate_results import (
    build_refraction_artifacts,
    gen_spherical_patch,
    save_refraction_npz,
    SRC_THETA,
    SRC_PHI,
    TGT_THETA,
    TGT_PHI,
    DEFAULT_CHUNK_SIZE,
    REGULARIZATION_MULTIPLIER,
)

# ---------------------------------------------------------------------------
# Problem geometry (unchanged from generate_results.py)
# ---------------------------------------------------------------------------

NK    = int(sys.argv[1]) if len(sys.argv) > 1 else 1600
chunk = DEFAULT_CHUNK_SIZE

# Generate clouds ONCE — same QMC points for all pairs
print(f"Generating Halton QMC clouds (NK={NK})...")
x = gen_spherical_patch(NK, *SRC_THETA, *SRC_PHI, skip=0)
y = gen_spherical_patch(NK, *TGT_THETA, *TGT_PHI, skip=0)
print(f"  Source: {len(x)} pts,  Target: {len(y)} pts")

# Regularisation schedule
k_final = REGULARIZATION_MULTIPLIER * int(np.floor(np.sqrt(NK)))
print(f"  k_final={k_final}\n")

# Density factories (keyed by name)
MAKER_MAP = {
    'uniform':  make_patch_uniform,
    'gaussian': make_patch_gaussian,
    'donut':    make_patch_donut,
    'cross':    make_patch_cross,
}
DENSITY_NAMES = ['uniform', 'gaussian', 'donut', 'cross']

results_dir = os.path.join(REPO_ROOT, 'results')
os.makedirs(results_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

total_pairs = len(DENSITY_NAMES) ** 2
pair_idx    = 0

for src_name in DENSITY_NAMES:
    for tgt_name in DENSITY_NAMES:
        pair_idx += 1
        out_path = os.path.join(
            results_dir, f'results_refraction_{src_name}_{tgt_name}_NK{NK}.npz')
        print(f"[{pair_idx:2d}/{total_pairs}] {src_name} -> {tgt_name}  "
              f"-> {os.path.basename(out_path)}")
        t_pair = time.time()

        # ── Densities ──────────────────────────────────────────────────────
        P = MAKER_MAP[src_name](*SRC_THETA, *SRC_PHI)
        Q = MAKER_MAP[tgt_name](*TGT_THETA, *TGT_PHI)

        p_raw = P(x)
        q_raw = Q(y)

        p_sum = p_raw.sum()
        q_sum = q_raw.sum()
        if p_sum == 0 or q_sum == 0:
            print(f"  WARNING: zero-mass density for {src_name}/{tgt_name}, skipping.")
            continue

        p    = p_raw / p_sum
        q    = q_raw / q_sum
        print(f"  src support={int((p>0).sum())}, tgt support={int((q>0).sum())}")

        # ── Shared cold-start Sinkhorn solve ──────────────────────────────
        solution = solve_cold_start_sinkhorn(
            x, y, p, q, k_final, chunk_size=chunk, using_identity=True,
        )
        f_raw = solution["f_raw"]
        g_raw = solution["g_raw"]
        f_id = solution["f_id"]
        g_id = solution["g_id"]
        f = solution["f"]
        g = solution["g"]

        # ── Shared notebook-compatible fields and NPZ schema ──────────────
        artifacts = build_refraction_artifacts(
            x, y, p, q, f_raw, g_raw, f, g, P, Q,
            chunk=chunk, kappa=0.6,
        )
        save_refraction_npz(
            out_path, x=x, y=y, p=p, q=q,
            f_raw=f_raw, g_raw=g_raw, f_id=f_id, g_id=g_id, f=f, g=g,
            artifacts=artifacts, kappa=0.6,
            metadata={"src_density": src_name, "tgt_density": tgt_name},
        )
        sz = os.path.getsize(out_path) / 1024
        print(f"  Saved ({sz:.0f} KB)  pair time {time.time()-t_pair:.1f}s\n")

print(f"Done. {pair_idx} NPZ files written to {results_dir}/")
