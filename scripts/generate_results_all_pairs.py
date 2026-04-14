"""
generate_results_all_pairs.py — All 16 source×target density-pair NPZs

Runs the refraction Sinkhorn pipeline (κ=0.6) for every ordered pair of the
four density shapes:
    uniform  — flat indicator (all weight equal)
    gaussian — isotropic Gaussian centred at patch centre
    donut    — hard annulus (inner=0.25·d_max, outer=0.75·d_max)
    cross    — sum of 4 Gaussians at N/S/E/W of patch centre

The QMC point cloud and patch geometry are identical to generate_results_refraction.py
(NK=1600, same SRC_THETA/PHI and TGT_THETA/PHI). Only the density weights change.

Outputs:
    results/results_refraction_{src}_{tgt}_NK1600.npz   (16 files)

Run:
    python scripts/generate_results_all_pairs.py [NK]
"""

import os, sys, time
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, REPO_ROOT)

# Set kappa BEFORE importing anything that caches cost computations
from reflector.cost import set_kappa, cost_matrix_chunk
set_kappa(0.6)

from reflector.distributions import (
    make_patch_uniform,
    make_patch_gaussian,
    make_patch_donut,
    make_patch_cross,
    stereo_north,
)
from reflector.sinkhorn import (
    sinkhorn_step,
    sinkhorn_identity_f_step,
    sinkhorn_identity_g_step,
)
from reflector.build import c_transform_gc, c_transform_fc

# ---------------------------------------------------------------------------
# Halton quasi-random sequence (identical to generate_results_refraction.py)
# ---------------------------------------------------------------------------

def _halton(index, base):
    result, f, i = 0.0, 1.0, index
    while i > 0:
        f /= base
        result += f * (i % base)
        i //= base
    return result


def gen_spherical_patch(n, theta_min, theta_max, phi_min, phi_max,
                        base2=2, base3=3, skip=0):
    """Sample n points on a spherical patch using Halton QMC.

    Convention (matching existing script):
        z       = cos(phi_gen),  phi_gen ∈ [phi_min, phi_max]
        azimuth = theta_gen,     theta_gen ∈ [theta_min, theta_max]
    """
    pts = []
    idx = skip
    while len(pts) < n:
        u1 = _halton(idx, base2)
        u2 = _halton(idx, base3)
        phi   = phi_min   + u2 * (phi_max   - phi_min)
        theta = theta_min + u1 * (theta_max - theta_min)
        pts.append([np.sin(phi) * np.cos(theta),
                    np.sin(phi) * np.sin(theta),
                    np.cos(phi)])
        idx += 1
    return np.array(pts, dtype=np.float64)


# ---------------------------------------------------------------------------
# Problem geometry (unchanged from generate_results_refraction.py)
# ---------------------------------------------------------------------------

NK    = int(sys.argv[1]) if len(sys.argv) > 1 else 1600
chunk = 512

SRC_THETA = (np.pi / 12, np.pi / 3)
SRC_PHI   = (np.pi / 12, np.pi / 4)
TGT_THETA = (np.pi / 3,  5 * np.pi / 12)
TGT_PHI   = (np.pi / 10, np.pi / 5)

# Generate clouds ONCE — same QMC points for all pairs
print(f"Generating Halton QMC clouds (NK={NK})...")
x = gen_spherical_patch(NK, *SRC_THETA, *SRC_PHI, skip=0)
y = gen_spherical_patch(NK, *TGT_THETA, *TGT_PHI, skip=0)
print(f"  Source: {len(x)} pts,  Target: {len(y)} pts")

# Regularisation schedule
multiplier  = 8
k_final     = multiplier * int(np.floor(np.sqrt(NK)))
id_step     = int(np.floor(np.sqrt(k_final)))
cap_iter    = 16
cap_iter_id = 16
cap_thr     = 1e-5
print(f"  k_final={k_final}, id_step={id_step}\n")

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
# 256×256 stereographic grid (reused for meshgrid evaluation)
# ---------------------------------------------------------------------------
grid_res = 256
gside    = np.linspace(-0.6, 0.6, grid_res)
UU, VV   = np.meshgrid(gside, gside, indexing='ij')
N2       = UU ** 2 + VV ** 2
denom_g  = 1.0 + N2
x_grid   = np.stack([2 * UU / denom_g, 2 * VV / denom_g, (1 - N2) / denom_g],
                    axis=-1).reshape(-1, 3)  # upper-hemisphere grid points

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
        logp = np.where(p > 0, np.log(p), -np.inf)
        logq = np.where(q > 0, np.log(q), -np.inf)

        print(f"  src support={int((p>0).sum())}, tgt support={int((q>0).sum())}")

        # ── Step 1: Cold-start Sinkhorn ────────────────────────────────────
        f = np.zeros(NK, dtype=np.float64)
        g = np.zeros(NK, dtype=np.float64)
        maxdif = cap_thr + 1.0
        i = 0
        t0 = time.time()
        while maxdif > cap_thr:
            f, g, maxdif = sinkhorn_step(x, y, logp, logq, f, g, k_final, chunk)
            i += 1
            if i >= cap_iter:
                break
        print(f"  Step1: {i} iters, maxdif={maxdif:.3e}  ({time.time()-t0:.1f}s)")

        # ── Step 2: Identity F loop ────────────────────────────────────────
        f_id   = np.zeros(NK, dtype=np.float64)
        regvar = 1;  it = 0
        t0 = time.time()
        while regvar < k_final:
            f_id_new, _ = sinkhorn_identity_f_step(x, y, logp, f_id, regvar, chunk)
            f_id = np.where(p > 0, f_id_new, f_id)
            it += 1
            regvar += id_step
        maxdif = cap_thr + 1.0;  i = 0
        while maxdif > cap_thr:
            f_id_new, maxdif = sinkhorn_identity_f_step(x, y, logp, f_id, k_final, chunk)
            f_id = np.where(p > 0, f_id_new, f_id)
            i += 1
            if i >= cap_iter_id:
                break
        print(f"  Step2: {i} final iters, maxdif={maxdif:.3e}  ({time.time()-t0:.1f}s)")

        # ── Step 3: Identity G loop ────────────────────────────────────────
        g_id   = np.zeros(NK, dtype=np.float64)
        regvar = 1;  it = 0
        t0 = time.time()
        while regvar < k_final:
            g_id_new, _ = sinkhorn_identity_g_step(x, y, logq, g_id, regvar, chunk)
            g_id = np.where(q > 0, g_id_new, g_id)
            it += 1
            regvar += id_step
        maxdif = cap_thr + 1.0;  i = 0
        while maxdif > cap_thr:
            g_id_new, maxdif = sinkhorn_identity_g_step(x, y, logq, g_id, k_final, chunk)
            g_id = np.where(q > 0, g_id_new, g_id)
            i += 1
            if i >= cap_iter_id:
                break
        print(f"  Step3: {i} final iters, maxdif={maxdif:.3e}  ({time.time()-t0:.1f}s)")

        # ── Step 4: Normalise and subtract identity ────────────────────────
        f = np.where(p > 0, f, 0.0)
        g = np.where(q > 0, g, 0.0)

        mx_fid = float(np.max(f_id));  f_id -= mx_fid;  g_id += mx_fid
        mx_f   = float(np.max(f));     f    -= mx_f;     g    += mx_f

        f_raw = f.copy()
        g_raw = g.copy()

        f -= f_id;  g -= g_id

        total_cost = float(np.sum(p[p > 0] * f[p > 0]) + np.sum(q[q > 0] * g[q > 0]))
        print(f"  Total cost: {total_cost:.6e}")

        # ── Step 5: Refractor surface ──────────────────────────────────────
        R   = np.exp(f)
        Ref = 2.0 * x * R[:, np.newaxis]

        mask_p = p > 0
        mask_q = q > 0

        # ── Step 6: C-transforms ──────────────────────────────────────────
        g_raw_masked = np.where(mask_q, g_raw, -1e300)
        f_raw_masked = np.where(mask_p, f_raw, -1e300)
        gc   = c_transform_gc(x, y, g_raw_masked, chunk)
        fc   = c_transform_fc(x, y, f_raw_masked, chunk)
        Refc = 2.0 * x * np.exp(gc)[:, np.newaxis]

        # ── Step 7: Density meshgrids (256×256 stereo, upper hemisphere) ──
        X_MeshGrid = P(x_grid).reshape(grid_res, grid_res)
        Y_MeshGrid = Q(x_grid).reshape(grid_res, grid_res)

        # ── Step 8: Y_projected ───────────────────────────────────────────
        u_y, v_y    = stereo_north(y[mask_q])
        Y_projected = np.column_stack([u_y, v_y, q[mask_q]])

        # ── Step 9: Push-forward (argmin OT map) ──────────────────────────
        src_idx  = np.where(mask_p)[0]
        tgt_idx  = np.where(mask_q)[0]
        y_tgt    = y[tgt_idx]
        g_tgt    = g_raw[tgt_idx]
        pushed_y = []
        for i0 in range(0, len(src_idx), chunk):
            i1     = min(i0 + chunk, len(src_idx))
            rows   = src_idx[i0:i1]
            C_blk  = cost_matrix_chunk(x[rows], y_tgt)
            j_star = np.argmin(C_blk - g_tgt[np.newaxis, :], axis=1)
            pushed_y.append(y_tgt[j_star])
        pushed_y = np.vstack(pushed_y)

        u_push, v_push     = stereo_north(pushed_y)
        q_push             = Q(pushed_y)
        Y_Pushed_projected = np.column_stack([u_push, v_push, q_push])
        y_push_3d          = pushed_y

        # ── Save ──────────────────────────────────────────────────────────
        np.savez(
            out_path,
            x=x, y=y,
            x_s=np.zeros((0, 3), dtype=np.float64),
            y_s=np.zeros((0, 3), dtype=np.float64),
            p=p, q=q,
            f_raw=f_raw, g_raw=g_raw,
            f_id=f_id,   g_id=g_id,
            f=f,         g=g,
            R=R,         Ref=Ref,
            gc=gc,       fc=fc,   Refc=Refc,
            X_MeshGrid=X_MeshGrid, Y_MeshGrid=Y_MeshGrid, grid_side=gside,
            Y_projected=Y_projected,
            Y_Pushed_projected=Y_Pushed_projected,
            y_push_3d=y_push_3d,
            kappa=np.float64(0.6),
            src_density=src_name,
            tgt_density=tgt_name,
        )
        sz = os.path.getsize(out_path) / 1024
        print(f"  Saved ({sz:.0f} KB)  pair time {time.time()-t_pair:.1f}s\n")

print(f"Done. {pair_idx} NPZ files written to {results_dir}/")
