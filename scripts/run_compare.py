"""
run_compare.py — Python Sinkhorn reflector for arbitrary NK (Halton QMC).

Matches main_compare.cpp exactly:
  - Same Halton cloud (base-2/base-3, half=0.6, skip=0)
  - Same P/Q (SquareToCircle)
  - No warm-start (f=g=0 init)
  - Same cap_iteration=16, cap_thr=1e-5
  - identity maxdif computed over supported points only

Usage:
    python run_compare.py 300
    python run_compare.py 1600
"""

import os, sys, time
import numpy as np
from scipy.special import logsumexp

REPO_ROOT   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(REPO_ROOT, 'results')
sys.path.insert(0, REPO_ROOT)
from reflector.distributions import P_square, Q_circle
from reflector.sinkhorn import (
    sinkhorn_step,
    sinkhorn_identity_f_step,
    sinkhorn_identity_g_step,
)

# ---------------------------------------------------------------------------
# Halton QMC — identical to generate_results.py
# ---------------------------------------------------------------------------

def _halton(index, base):
    result, f, i = 0.0, 1.0, index
    while i > 0:
        f /= base; result += f * (i % base); i //= base
    return result

def _sphere_pt(X, Y, upper=True):
    N2 = X*X + Y*Y; d = 1.0 + N2
    return [2*X/d, 2*Y/d, (1.0-N2)/d if upper else -(1.0-N2)/d]

def gen_cloud(n, upper, skip=0, half=0.6):
    pts, idx = [], skip
    while len(pts) < n:
        X = (_halton(idx, 2) - 0.5) * 2.0 * half
        Y = (_halton(idx, 3) - 0.5) * 2.0 * half
        pts.append(_sphere_pt(X, Y, upper))
        idx += 1
    return np.array(pts, dtype=np.float64)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

NK = int(sys.argv[1]) if len(sys.argv) > 1 else 300
chunk       = 512
multiplier  = 8
k_final     = multiplier * int(np.floor(np.sqrt(NK)))
id_step     = int(np.floor(np.sqrt(k_final)))
cap_iter    = 16      # match C++ cap_iteration
cap_thr     = 1e-5
out_dir     = os.path.join(RESULTS_DIR, f"output_py_NK{NK}")

os.makedirs(out_dir, exist_ok=True)
print(f"NK={NK}, k_final={k_final}, id_step={id_step}")

# ---------------------------------------------------------------------------
# Cloud generation
# ---------------------------------------------------------------------------

print("Generating Halton clouds...")
x = gen_cloud(NK, upper=True,  skip=0)
y = gen_cloud(NK, upper=False, skip=0)

p_raw = P_square(x);  p = p_raw / p_raw.sum()
q_raw = Q_circle(y);  q = q_raw / q_raw.sum()
logp  = np.where(p > 0, np.log(p), -np.inf)
logq  = np.where(q > 0, np.log(q), -np.inf)

print(f"  Source support: {int(p_raw.sum()+0.5)} / {NK}")
print(f"  Target support: {int(q_raw.sum()+0.5)} / {NK}")

# ---------------------------------------------------------------------------
# Step 1 — Final Sinkhorn at k_final  (no warm-start, f=g=0)
# ---------------------------------------------------------------------------

f = np.zeros(NK, dtype=np.float64)
g = np.zeros(NK, dtype=np.float64)

print(f"\nFinal Sinkhorn at k={k_final}:")
maxdif = cap_thr + 1.0;  i = 0
t0 = time.time()
while maxdif > cap_thr:
    f, g, maxdif = sinkhorn_step(x, y, logp, logq, f, g, k_final, chunk)
    i += 1
    print(f"  iter {i:3d}, maxdif={maxdif:.4e}  ({time.time()-t0:.2f}s)")
    if i >= cap_iter: break
print(f"Final: {i} iters, last maxdif={maxdif:.4e}")

# ---------------------------------------------------------------------------
# Step 2 — Identity F
# ---------------------------------------------------------------------------

print(f"\nIdentity F (id_step={id_step}):")
f_id   = np.zeros(NK, dtype=np.float64)
regvar = 1;  it = 0
t0 = time.time()
while regvar < k_final:
    f_id_new, md = sinkhorn_identity_f_step(x, y, logp, f_id, regvar, chunk)
    f_id = np.where(p > 0, f_id_new, f_id)
    it += 1
    print(f"  iter {it:3d}, k={regvar:4d}, maxdif={md:.4e}")
    regvar += id_step

maxdif = cap_thr + 1.0;  i = 0
while maxdif > cap_thr:
    f_id_new, maxdif = sinkhorn_identity_f_step(x, y, logp, f_id, k_final, chunk)
    f_id = np.where(p > 0, f_id_new, f_id)
    i += 1
    if i >= cap_iter: break
print(f"Identity F: {i} final iters, last maxdif={maxdif:.4e}  ({time.time()-t0:.1f}s)")

# ---------------------------------------------------------------------------
# Step 3 — Identity G
# ---------------------------------------------------------------------------

print(f"\nIdentity G (id_step={id_step}):")
g_id   = np.zeros(NK, dtype=np.float64)
regvar = 1;  it = 0
t0 = time.time()
while regvar < k_final:
    g_id_new, md = sinkhorn_identity_g_step(x, y, logq, g_id, regvar, chunk)
    g_id = np.where(q > 0, g_id_new, g_id)
    it += 1
    print(f"  iter {it:3d}, k={regvar:4d}, maxdif={md:.4e}")
    regvar += id_step

maxdif = cap_thr + 1.0;  i = 0
while maxdif > cap_thr:
    g_id_new, maxdif = sinkhorn_identity_g_step(x, y, logq, g_id, k_final, chunk)
    g_id = np.where(q > 0, g_id_new, g_id)
    i += 1
    if i >= cap_iter: break
print(f"Identity G: {i} final iters, last maxdif={maxdif:.4e}  ({time.time()-t0:.1f}s)")

# ---------------------------------------------------------------------------
# Normalise and build reflector
# ---------------------------------------------------------------------------

# C++ never updates f/g for unsupported points (p[i]==0 / q[j]==0), keeping them at 0.
# Python sinkhorn_step updates ALL points; zero out unsupported values to match C++.
f = np.where(p > 0, f, 0.0)
g = np.where(q > 0, g, 0.0)

mx_fid = float(np.max(f_id));  f_id -= mx_fid;  g_id += mx_fid
mx_f   = float(np.max(f));     f    -= mx_f;     g    += mx_f
f -= f_id;  g -= g_id

total_cost = float(np.sum(p[p>0]*f[p>0]) + np.sum(q[q>0]*g[q>0]))
print(f"\nTotal cost: {total_cost:.6e}")

R_arr  = np.exp(f)
Ref_arr = 2.0 * x * R_arr[:, np.newaxis]

mask_p = p > 0
print(f"R (supported): min={R_arr[mask_p].min():.4f}, max={R_arr[mask_p].max():.4f}, "
      f"mean={R_arr[mask_p].mean():.4f}")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

np.save(f"{out_dir}/f.npy",    f)
np.save(f"{out_dir}/g.npy",    g)
np.save(f"{out_dir}/f_id.npy", f_id)
np.save(f"{out_dir}/g_id.npy", g_id)
np.save(f"{out_dir}/R.npy",    R_arr)
np.save(f"{out_dir}/Ref.npy",  Ref_arr)
print(f"\nSaved to {out_dir}/")
