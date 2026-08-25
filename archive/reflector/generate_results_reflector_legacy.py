# Archived SquareToCircle reflector result generator.
# Key formulas: R = exp(f), Ref = 2*x*R, and corrected f = f_raw - f_id.
# C-transforms: gc[j] = min_i(C(x_i,y_j)-f_i),
# fc[i] = min_j(C(x_i,y_j)-g_j).
# Run: python scripts/generate_results_reflector.py

import os, sys, math, time
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
FIGURES_DIR = os.path.join(REPO_ROOT, "figures")
sys.path.insert(0, REPO_ROOT)

from refracter.distributions import P_square, Q_circle
from refracter.sinkhorn import (
    sinkhorn_step,
    sinkhorn_identity_f_step,
    sinkhorn_identity_g_step,
)
from refracter.build import c_transform_gc, c_transform_fc
from refracter.cost import cost_matrix_chunk

# Generate the base-2/base-3 Halton sequence used by the archived notebook.

def _halton(index, base):
    result, f, i = 0.0, 1.0, index
    while i > 0:
        f /= base
        result += f * (i % base)
        i //= base
    return result

def _sphere_pt(X, Y, upper=True):
    # Inverse stereographic projection: (2X,2Y,1-X^2-Y^2)/(1+X^2+Y^2).
    N2 = X * X + Y * Y
    d  = 1.0 + N2
    z  = (1.0 - N2) / d
    return [2*X/d, 2*Y/d, z if upper else -z]

def gen_cloud(n, upper, skip=0, half=0.6):
    # Generate n points from the Halton sequence.
    pts, idx = [], skip
    while len(pts) < n:
        X = (_halton(idx, 2) - 0.5) * 2.0 * half
        Y = (_halton(idx, 3) - 0.5) * 2.0 * half
        pts.append(_sphere_pt(X, Y, upper))
        idx += 1
    return np.array(pts, dtype=np.float64)

# Set NK, generate both clouds, and evaluate normalized densities.

NK    = int(sys.argv[1]) if len(sys.argv) > 1 else 1600
chunk = 512

print(f"Generating Halton clouds (NK={NK})...")
t_start = time.time()

x   = gen_cloud(NK, upper=True,  skip=0)    # source: upper hemisphere
y   = gen_cloud(NK, upper=False, skip=0)    # target: lower hemisphere
# Keep x_s/y_s as empty arrays so the npz keys still exist (backward compat)
x_s = np.zeros((0, 3), dtype=np.float64)
y_s = np.zeros((0, 3), dtype=np.float64)

# Densities (normalised)
p_raw = P_square(x);  p = p_raw / p_raw.sum()
q_raw = Q_circle(y);  q = q_raw / q_raw.sum()
logp  = np.where(p > 0, np.log(p), -np.inf)
logq  = np.where(q > 0, np.log(q), -np.inf)

print(f"  Source support: {int(p_raw.sum()+0.5)} / {NK}")
print(f"  Target support: {int(q_raw.sum()+0.5)} / {NK}")

# Regularisation schedule (same as run_compare.py / main_compare.cpp)
multiplier  = 8
k_final     = multiplier * int(np.floor(np.sqrt(NK)))   # 8·40 = 320
id_step     = int(np.floor(np.sqrt(k_final)))            # 17
cap_iter    = 16    # match C++ cap_iteration=16
cap_iter_id = 16    # identity loops
cap_thr     = 1e-5

print(f"  k_final={k_final}, id_step={id_step}")

# ---------------------------------------------------------------------------
# Step 1 — Cold-start Sinkhorn at k_final  (f=g=0, cap_iter=16)
#   Matches run_compare.py / main_compare.cpp exactly — no warmstart,
#   no multi-scale ramp.  cap_iter=16 gives R ∈ [0.958, 1.404].
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
print(f"  (main Sinkhorn converged: {maxdif <= cap_thr})")

# ---------------------------------------------------------------------------
# Step 2 — Identity F loop  (source marginal, with unsupported-point masking)
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
# C++ never updates f/g for unsupported points (p[i]==0 / q[j]==0); zero them out here.
f = np.where(p > 0, f, 0.0)
g = np.where(q > 0, g, 0.0)
# Shift f_id so max = 0, compensate in g_id
mx_fid = float(np.max(f_id));  f_id -= mx_fid;  g_id += mx_fid
# Shift f so max = 0, compensate in g
mx_f   = float(np.max(f));     f    -= mx_f;     g    += mx_f

f_raw = f.copy()   # save raw Sinkhorn potentials before identity subtraction
g_raw = g.copy()

f -= f_id;  g -= g_id

total_cost = float(np.sum(p[p>0]*f[p>0]) + np.sum(q[q>0]*g[q>0]))
print(f"\nTotal cost: {total_cost:.6e}")

# ---------------------------------------------------------------------------
# Step 5 — Reflector surface
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
# Use g_raw / f_raw (raw OT potentials) but mask unsupported points to -1e300
# so the min over j only considers supported target / source points.
# gc[i] = min_{j: q[j]>0}(C(xi,yj) - g_raw[j]) should ≈ f_raw[i] at convergence.
g_raw_masked = np.where(mask_q, g_raw, -1e300)
f_raw_masked = np.where(mask_p, f_raw, -1e300)
gc   = c_transform_gc(x, y, g_raw_masked, chunk)
fc   = c_transform_fc(x, y, f_raw_masked, chunk)
Refc = 2.0 * x * np.exp(gc)[:, np.newaxis]
print(f"  done ({time.time()-t0:.2f}s)")
print(f"  max|f_raw-gc(g_raw)| (supported) = {np.abs(f_raw[mask_p]-gc[mask_p]).max():.4e}")

# ---------------------------------------------------------------------------
# Step 7 — Density meshgrids (256×256 over stereographic [-0.6, 0.6]^2)
# ---------------------------------------------------------------------------
print("\nStep 7 — Density meshgrids...")
grid_res = 256
gside = np.linspace(-0.6, 0.6, grid_res)
UU, VV = np.meshgrid(gside, gside, indexing='ij')   # both (grid_res, grid_res)
N2 = UU**2 + VV**2;  denom = 1.0 + N2

# Source: upper hemisphere via inverse north-pole stereo
x_grid = np.stack([2*UU/denom, 2*VV/denom,  (1-N2)/denom], axis=-1).reshape(-1, 3)
X_MeshGrid = P_square(x_grid).reshape(grid_res, grid_res)

# Target: lower hemisphere via inverse south-pole stereo
y_grid = np.stack([2*UU/denom, 2*VV/denom, -(1-N2)/denom], axis=-1).reshape(-1, 3)
Y_MeshGrid = Q_circle(y_grid).reshape(grid_res, grid_res)

# ---------------------------------------------------------------------------
# Step 8 — Y_projected: target support in south-pole 2D projection
# ---------------------------------------------------------------------------
denom_y       = 1.0 - y[mask_q, 2]
u_y, v_y      = y[mask_q, 0] / denom_y, y[mask_q, 1] / denom_y
Y_projected   = np.column_stack([u_y, v_y, q[mask_q]])

# ---------------------------------------------------------------------------
# Step 9 — Push-forward: optimal transport map  x_i → y_{j*(i)}
#   j*(i) = argmin_{j: q[j]>0}  C(x_i, y_j) - g_j   (supported target only)
# ---------------------------------------------------------------------------
print("\nStep 9 — Push-forward (argmin OT map)...")
t0 = time.time()
src_idx  = np.where(mask_p)[0]
tgt_idx  = np.where(mask_q)[0]
y_tgt    = y[tgt_idx]          # only supported target points
g_tgt    = g_raw[tgt_idx]
pushed_y = []
for i0 in range(0, len(src_idx), chunk):
    i1    = min(i0 + chunk, len(src_idx))
    rows  = src_idx[i0:i1]
    C_blk = cost_matrix_chunk(x[rows], y_tgt)     # (blk, n_tgt)
    j_star = np.argmin(C_blk - g_tgt[np.newaxis, :], axis=1)
    pushed_y.append(y_tgt[j_star])
pushed_y = np.vstack(pushed_y)                      # (n_src_support, 3)

denom_push        = 1.0 - pushed_y[:, 2]
u_push, v_push    = pushed_y[:, 0]/denom_push, pushed_y[:, 1]/denom_push
q_push            = Q_circle(pushed_y)
Y_Pushed_projected = np.column_stack([u_push, v_push, q_push])
y_push_3d         = pushed_y                        # 3D for Plotly scene
print(f"  {len(y_push_3d)} pushed rays  ({time.time()-t0:.2f}s)")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out = f"results_NK{NK}.npz"
np.savez(
    out,
    # Clouds
    x=x, y=y, x_s=x_s, y_s=y_s,
    # Densities
    p=p, q=q,
    # Potentials
    f_raw=f_raw, g_raw=g_raw,
    f_id=f_id, g_id=g_id,
    f=f, g=g,
    # Reflector
    R=R, Ref=Ref,
    # C-transforms
    gc=gc, fc=fc, Refc=Refc,
    # Meshgrids
    X_MeshGrid=X_MeshGrid, Y_MeshGrid=Y_MeshGrid, grid_side=gside,
    # Projected distributions
    Y_projected=Y_projected,
    Y_Pushed_projected=Y_Pushed_projected,
    y_push_3d=y_push_3d,
)
sz = os.path.getsize(out) / 1024
print(f"\nSaved {out}  ({sz:.0f} KB)   total time {time.time()-t_start:.1f}s")

# ---------------------------------------------------------------------------
# Figure — Reflector surface  (all NK points, coloured by R)
# Matches the C++ comparison plots: plot ALL NK points so that the
# unsupported ring (low R) and the supported dome (high R) both appear,
# giving the characteristic non-smooth bumpy shape.
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
plt.colorbar(sc, ax=ax, label='R (reflector radius)', shrink=0.65)
ax.set_title(f'Reflector Surface — Python  (NK={NK}, k_final={k_final})')
ax.set_xlabel('X');  ax.set_ylabel('Y');  ax.set_zlabel('Z')
fig.tight_layout()
fig_path = os.path.join(FIGURES_DIR, f"fig_reflector_3d_NK{NK}.png")
os.makedirs(FIGURES_DIR, exist_ok=True)
fig.savefig(fig_path, dpi=150, bbox_inches='tight')
print(f"Saved {fig_path}")
