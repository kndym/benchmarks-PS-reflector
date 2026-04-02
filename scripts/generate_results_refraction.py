"""
generate_results_refraction.py — Refraction Sinkhorn benchmark, NK=1600

Generates a 1600-point Halton quasi-Monte Carlo cloud on spherical patches,
runs the full Sinkhorn-divergence pipeline (cold start, f=g=0, cap_iter=16),
and saves a comprehensive results bundle to results/results_refraction_NK1600.npz.

Key differences from generate_results.py (SquareToCircle reflector):
  * κ = 0.6  (refraction cost:  c(x,y) = -log(1 - 0.6·(x·y)))
  * Both source and target are on the upper hemisphere (spherical patches)
  * Source  Ω  : θ ∈ [π/12, π/3],  φ ∈ [π/12, π/4]
  * Target  Ω* : θ ∈ [π/10, π/5],  φ ∈ [π/10, π/5]
  * North-pole stereographic projection is used for both source and target

Run:  python scripts/generate_results_refraction.py [NK]
      (default NK=1600)
"""

import os, sys, math, time
import numpy as np

# Ensure we can import the reflector package from the repo root
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, REPO_ROOT)

# Set kappa BEFORE importing anything that caches cost computations
from reflector.cost import set_kappa, cost_matrix_chunk
set_kappa(0.6)

from reflector.distributions import P_refraction_patch, Q_refraction_patch
from reflector.sinkhorn import (
    sinkhorn_step,
    sinkhorn_identity_f_step,
    sinkhorn_identity_g_step,
)
from reflector.build import c_transform_gc, c_transform_fc
from reflector.distributions import stereo_north

# ---------------------------------------------------------------------------
# Halton quasi-random sequence
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
    """Sample n points uniformly on a spherical patch using Halton QMC.

    Uses the equal solid-angle transform:
        cos(θ) = cos(θ_max) + u1 · (cos(θ_min) - cos(θ_max))
        φ      = φ_min + u2 · (φ_max - φ_min)

    Parameters
    ----------
    theta_min, theta_max : polar angle range in radians
    phi_min, phi_max     : azimuthal angle range in radians
    base2, base3 : Halton bases for u1, u2
    skip : number of Halton indices to skip at the start
    """


    #cos_max = np.cos(theta_min)   # larger cos (smaller θ)
    #cos_min = np.cos(theta_max)   # smaller cos (larger θ)

    pts = []
    idx = skip
    while len(pts) < n:
        u1 = _halton(idx, base2)
        u2 = _halton(idx, base3)
        #sin_theta = np.sqrt(max(0.0, 1.0 - cos_theta ** 2))
        phi = phi_min + u2 * (phi_max - phi_min)
        theta= theta_min + u1 * (theta_max -  theta_min)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        pts.append([sin_phi * np.cos(theta),
                    sin_phi * np.sin(theta),
                    cos_phi])
        idx += 1
    return np.array(pts, dtype=np.float64)


# ---------------------------------------------------------------------------
# Problem setup
# ---------------------------------------------------------------------------

NK    = int(sys.argv[1]) if len(sys.argv) > 1 else 1600
chunk = 512

print(f"=== Refraction benchmark (κ=0.6, NK={NK}) ===")
t_start = time.time()

# Patch bounds (matching C++ test header)
SRC_THETA = (np.pi / 12, np.pi / 3)
SRC_PHI   = (np.pi / 12, np.pi / 4)
TGT_THETA = ((np.pi / 10), (np.pi / 5 ) )
TGT_PHI   = ((np.pi / 10)  - 0.15 , (np.pi / 5) - 0.15 )

print("Generating Halton QMC clouds on spherical patches...")
x = gen_spherical_patch(NK, *SRC_THETA, *SRC_PHI, skip=0)
y = gen_spherical_patch(NK, *TGT_THETA, *TGT_PHI, skip=0)

print(f"  Source: {len(x)} pts, z ∈ [{x[:,2].min():.4f}, {x[:,2].max():.4f}]")
print(f"  Target: {len(y)} pts, z ∈ [{y[:,2].min():.4f}, {y[:,2].max():.4f}]")

# Densities — all points are inside the patch by construction → uniform
p_raw = P_refraction_patch(x)
q_raw = Q_refraction_patch(y)
print(f"  Source support: {int(p_raw.sum()+0.5)} / {NK}")
print(f"  Target support: {int(q_raw.sum()+0.5)} / {NK}")

p = p_raw / p_raw.sum()
q = q_raw / q_raw.sum()
logp = np.where(p > 0, np.log(p), -np.inf)
logq = np.where(q > 0, np.log(q), -np.inf)

# Regularisation schedule
multiplier  = 8
k_final     = multiplier * int(np.floor(np.sqrt(NK)))
id_step     = int(np.floor(np.sqrt(k_final)))
cap_iter    = 16
cap_iter_id = 16
cap_thr     = 1e-5

print(f"  k_final={k_final}, id_step={id_step}")

# ---------------------------------------------------------------------------
# Step 1 — Cold-start Sinkhorn at k_final  (f=g=0, cap_iter=16)
# ---------------------------------------------------------------------------
f = np.zeros(NK, dtype=np.float64)
g = np.zeros(NK, dtype=np.float64)

print(f"\nStep 1 — Sinkhorn at k={k_final} (cold start, cap_iter={cap_iter}):")
maxdif = cap_thr + 1.0
i = 0
t0 = time.time()
while maxdif > cap_thr:
    f, g, maxdif = sinkhorn_step(x, y, logp, logq, f, g, k_final, chunk)
    i += 1
    print(f"  iter {i:3d}, maxdif={maxdif:.4e}  ({time.time()-t0:.2f}s)")
    if i >= cap_iter:
        break
print(f"  Final: {i} iters, last maxdif={maxdif:.4e}")

# ---------------------------------------------------------------------------
# Step 2 — Identity F loop  (source marginal)
# ---------------------------------------------------------------------------
print(f"\nStep 2 — Identity Sinkhorn for source (f_id), id_step={id_step}:")
f_id   = np.zeros(NK, dtype=np.float64)
regvar = 1;  it = 0
t0 = time.time()
while regvar < k_final:
    f_id_new, _ = sinkhorn_identity_f_step(x, y, logp, f_id, regvar, chunk)
    f_id = np.where(p > 0, f_id_new, f_id)
    it += 1
    print(f"  iter {it:4d}, k={regvar:4d}")
    regvar += id_step

maxdif = cap_thr + 1.0;  i = 0
while maxdif > cap_thr:
    f_id_new, maxdif = sinkhorn_identity_f_step(x, y, logp, f_id, k_final, chunk)
    f_id = np.where(p > 0, f_id_new, f_id)
    i += 1
    if i >= cap_iter_id: break
print(f"  Identity F: {i} final iters, last maxdif={maxdif:.4e}  ({time.time()-t0:.1f}s)")

# ---------------------------------------------------------------------------
# Step 3 — Identity G loop  (target marginal)
# ---------------------------------------------------------------------------
print(f"\nStep 3 — Identity Sinkhorn for target (g_id), id_step={id_step}:")
g_id   = np.zeros(NK, dtype=np.float64)
regvar = 1;  it = 0
t0 = time.time()
while regvar < k_final:
    g_id_new, _ = sinkhorn_identity_g_step(x, y, logq, g_id, regvar, chunk)
    g_id = np.where(q > 0, g_id_new, g_id)
    it += 1
    print(f"  iter {it:4d}, k={regvar:4d}")
    regvar += id_step

maxdif = cap_thr + 1.0;  i = 0
while maxdif > cap_thr:
    g_id_new, maxdif = sinkhorn_identity_g_step(x, y, logq, g_id, k_final, chunk)
    g_id = np.where(q > 0, g_id_new, g_id)
    i += 1
    if i >= cap_iter_id: break
print(f"  Identity G: {i} final iters, last maxdif={maxdif:.4e}  ({time.time()-t0:.1f}s)")

# ---------------------------------------------------------------------------
# Step 4 — Normalise and subtract identity terms
# ---------------------------------------------------------------------------
f = np.where(p > 0, f, 0.0)
g = np.where(q > 0, g, 0.0)

mx_fid = float(np.max(f_id));  f_id -= mx_fid;  g_id += mx_fid
mx_f   = float(np.max(f));     f    -= mx_f;     g    += mx_f

f_raw = f.copy()
g_raw = g.copy()

f -= f_id;  g -= g_id

total_cost = float(np.sum(p[p>0]*f[p>0]) + np.sum(q[q>0]*g[q>0]))
print(f"\nTotal cost: {total_cost:.6e}")

# ---------------------------------------------------------------------------
# Step 5 — Refractor surface
# ---------------------------------------------------------------------------
R   = np.exp(f)
Ref = 2.0 * x * R[:, np.newaxis]

mask_p = p > 0;  mask_q = q > 0
print(f"R (supported): min={R[mask_p].min():.4f}, max={R[mask_p].max():.4f}, mean={R[mask_p].mean():.4f}")

# ---------------------------------------------------------------------------
# Step 6 — C-transforms
# ---------------------------------------------------------------------------
print("\nStep 6 — C-transforms...")
t0 = time.time()
g_raw_masked = np.where(mask_q, g_raw, -1e300)
f_raw_masked = np.where(mask_p, f_raw, -1e300)
gc   = c_transform_gc(x, y, g_raw_masked, chunk)
fc   = c_transform_fc(x, y, f_raw_masked, chunk)
Refc = 2.0 * x * np.exp(gc)[:, np.newaxis]
print(f"  done ({time.time()-t0:.2f}s)")
print(f"  max|f_raw-gc(g_raw)| (supported) = {np.abs(f_raw[mask_p]-gc[mask_p]).max():.4e}")

# ---------------------------------------------------------------------------
# Step 7 — Density meshgrids (256×256 over stereographic [-0.6, 0.6]²)
#           Both source and target use NORTH-pole projection (upper hemisphere)
# ---------------------------------------------------------------------------
print("\nStep 7 — Density meshgrids...")
grid_res = 256
gside = np.linspace(-0.6, 0.6, grid_res)
UU, VV = np.meshgrid(gside, gside, indexing='ij')
N2 = UU**2 + VV**2;  denom = 1.0 + N2

# Source: upper hemisphere via inverse north-pole stereo
x_grid = np.stack([2*UU/denom, 2*VV/denom, (1-N2)/denom], axis=-1).reshape(-1, 3)
X_MeshGrid = P_refraction_patch(x_grid).reshape(grid_res, grid_res)

# Target: ALSO upper hemisphere via inverse north-pole stereo (both on upper hem.)
y_grid = np.stack([2*UU/denom, 2*VV/denom, (1-N2)/denom], axis=-1).reshape(-1, 3)
Y_MeshGrid = Q_refraction_patch(y_grid).reshape(grid_res, grid_res)

# ---------------------------------------------------------------------------
# Step 8 — Y_projected: target support in north-pole 2D projection
# ---------------------------------------------------------------------------
u_y, v_y    = stereo_north(y[mask_q])
Y_projected = np.column_stack([u_y, v_y, q[mask_q]])

# ---------------------------------------------------------------------------
# Step 9 — Push-forward: optimal transport map  x_i → y_{j*(i)}
# ---------------------------------------------------------------------------
print("\nStep 9 — Push-forward (argmin OT map)...")
t0 = time.time()
src_idx  = np.where(mask_p)[0]
tgt_idx  = np.where(mask_q)[0]
y_tgt    = y[tgt_idx]
g_tgt    = g_raw[tgt_idx]
pushed_y = []
for i0 in range(0, len(src_idx), chunk):
    i1    = min(i0 + chunk, len(src_idx))
    rows  = src_idx[i0:i1]
    C_blk = cost_matrix_chunk(x[rows], y_tgt)
    j_star = np.argmin(C_blk - g_tgt[np.newaxis, :], axis=1)
    pushed_y.append(y_tgt[j_star])
pushed_y = np.vstack(pushed_y)

# Project push-forward via north-pole stereo (both hemispheres are upper)
u_push, v_push     = stereo_north(pushed_y)
q_push             = Q_refraction_patch(pushed_y)
Y_Pushed_projected = np.column_stack([u_push, v_push, q_push])
y_push_3d          = pushed_y
print(f"  {len(y_push_3d)} pushed rays  ({time.time()-t0:.2f}s)")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
results_dir = os.path.join(REPO_ROOT, "results")
os.makedirs(results_dir, exist_ok=True)
out = os.path.join(results_dir, f"results_refraction_NK{NK}.npz")
np.savez(
    out,
    # Clouds
    x=x, y=y,
    x_s=np.zeros((0, 3), dtype=np.float64),
    y_s=np.zeros((0, 3), dtype=np.float64),
    # Densities
    p=p, q=q,
    # Potentials
    f_raw=f_raw, g_raw=g_raw,
    f_id=f_id, g_id=g_id,
    f=f, g=g,
    # Refractor
    R=R, Ref=Ref,
    # C-transforms
    gc=gc, fc=fc, Refc=Refc,
    # Meshgrids
    X_MeshGrid=X_MeshGrid, Y_MeshGrid=Y_MeshGrid, grid_side=gside,
    # Projected distributions
    Y_projected=Y_projected,
    Y_Pushed_projected=Y_Pushed_projected,
    y_push_3d=y_push_3d,
    # Metadata
    kappa=np.float64(0.6),
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
