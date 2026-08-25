# Run the 381-point Python reflector case used by main_fast.cpp.
# Outputs are written to output_python/; run from scripts/reflector/.

import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("TkAgg")          # use a GUI backend; fall back silently
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from refracter.qmc import load_small_cloud
from refracter.sinkhorn import (
    sinkhorn_step,
    sinkhorn_identity_f_step,
    sinkhorn_identity_g_step,
)

# Match the C++ SquareToCircle source and target densities.
def P(v):
    # Source: uniform on the north-pole square [-0.5, 0.5]^2.
    p1 = v[0] / (1.0 + v[2])
    p2 = v[1] / (1.0 + v[2])
    return 1.0 if (-0.5 < p1 < 0.5 and -0.5 < p2 < 0.5) else 0.0


def Q(v):
    # Target: uniform on the south-pole disk u^2 + v^2 <= 0.25.
    p1 = v[0] / (1.0 - v[2])
    p2 = v[1] / (1.0 - v[2])
    return 1.0 if (p1 * p1 + p2 * p2 <= 0.25) else 0.0


# Load the 381-point cloud and evaluate both densities.
print("Loading QMC cloud (NK=381)...")
x, y = load_small_cloud()   # y has z negated by qmc.py (lower hemisphere)
NK = len(x)
assert NK == 381, f"Expected 381, got {NK}"

# Evaluate densities
p_w = np.array([P(x[i]) for i in range(NK)], dtype=np.float64)
q_w = np.array([Q(y[i]) for i in range(NK)], dtype=np.float64)

p_sum, q_sum = p_w.sum(), q_w.sum()
p_w /= p_sum
q_w /= q_sum
logp = np.where(p_w > 0, np.log(p_w), -np.inf)
logq = np.where(q_w > 0, np.log(q_w), -np.inf)

print(f"Source support: {int(p_sum + 0.5)} / {NK} points")
print(f"Target support: {int(q_sum + 0.5)} / {NK} points")

# Match the C++ schedule: k_final = 8 * floor(sqrt(NK)).

multiplier  = 8
k_final     = multiplier * int(np.floor(np.sqrt(NK)))   # 152
id_step     = int(np.floor(np.sqrt(k_final)))            # 12
cap_iter    = 16
cap_thr     = 1e-5
chunk_size  = 512

print(f"\nk_final={k_final}, id_step={id_step}, cap_iter={cap_iter}")

# Start all potentials at zero, as in the C++ fast case.

f = np.zeros(NK, dtype=np.float64)
g = np.zeros(NK, dtype=np.float64)

# Run the multi-scale schedule before the final regularization value.

regvar = multiplier * int(np.floor(np.sqrt(381)))   # k_small = 152
it = 0
while regvar < k_final:
    f, g, maxdif = sinkhorn_step(x, y, logp, logq, f, g, regvar, chunk_size)
    it += 1
    print(f"  multi-scale iter {it:3d}, k={regvar}, maxdif={maxdif:.4e}")
    regvar += int(round(k_final ** (1.0 / 3.0)))

# Iterate the main Sinkhorn solve at k_final.

print(f"\nFinal Sinkhorn at k={k_final}:")
maxdif = cap_thr + 1.0
i = 0
t0 = time.time()
while maxdif > cap_thr:
    f, g, maxdif = sinkhorn_step(x, y, logp, logq, f, g, k_final, chunk_size)
    i += 1
    print(f"  iter {i:3d}, maxdif={maxdif:.4e}  ({time.time()-t0:.2f}s)")
    if i >= cap_iter:
        break
print(f"Final: {i} iterations, last change={maxdif:.4e}")

# Compute the source self-transport correction.

print(f"\nIdentity Sinkhorn for source (f_id), id_step={id_step}:")
f_id = np.zeros(NK, dtype=np.float64)
regvar = 1
it = 0
t0 = time.time()
while regvar < k_final:
    f_id_new, maxdif = sinkhorn_identity_f_step(x, y, logp, f_id, regvar, chunk_size)
    f_id = np.where(p_w > 0, f_id_new, f_id)
    it += 1
    print(f"  iter {it:4d}, k={regvar:4d}, maxdif={maxdif:.4e}")
    regvar += id_step

maxdif = cap_thr + 1.0
i = 0
while maxdif > cap_thr:
    f_id_new, maxdif = sinkhorn_identity_f_step(x, y, logp, f_id, k_final, chunk_size)
    f_id = np.where(p_w > 0, f_id_new, f_id)
    i += 1
    if i >= cap_iter:
        break
print(f"Identity F: {i} final iterations, last change={maxdif:.4e}  ({time.time()-t0:.1f}s)")

# Compute the target self-transport correction.

print(f"\nIdentity Sinkhorn for target (g_id), id_step={id_step}:")
g_id = np.zeros(NK, dtype=np.float64)
regvar = 1
it = 0
t0 = time.time()
while regvar < k_final:
    g_id_new, maxdif = sinkhorn_identity_g_step(x, y, logq, g_id, regvar, chunk_size)
    g_id = np.where(q_w > 0, g_id_new, g_id)
    it += 1
    print(f"  iter {it:4d}, k={regvar:4d}, maxdif={maxdif:.4e}")
    regvar += id_step

maxdif = cap_thr + 1.0
i = 0
while maxdif > cap_thr:
    g_id_new, maxdif = sinkhorn_identity_g_step(x, y, logq, g_id, k_final, chunk_size)
    g_id = np.where(q_w > 0, g_id_new, g_id)
    i += 1
    if i >= cap_iter:
        break
print(f"Identity G: {i} final iterations, last change={maxdif:.4e}  ({time.time()-t0:.1f}s)")

# Normalize potentials and subtract the identity terms.

# Shift f_id so max=0, compensate in g_id
max_f_id = float(np.max(f_id))
f_id -= max_f_id
g_id += max_f_id

# Shift f so max=0, compensate in g
max_f = float(np.max(f))
f -= max_f
g += max_f

# Subtract identity
f -= f_id
g -= g_id

# Total cost (informational)
total_cost = float(
    np.sum(p_w[p_w > 0] * f[p_w > 0]) +
    np.sum(q_w[q_w > 0] * g[q_w > 0])
)
print(f"\nTotal cost: {total_cost:.6e}")

# Build reflector: R = exp(f),  Ref = 2 * x * R
R   = np.exp(f)
Ref = 2.0 * x * R[:, np.newaxis]

# Save the reflector arrays used by the comparison notebook.

os.makedirs("output_python", exist_ok=True)
np.save("output_python/f.npy",    f)
np.save("output_python/g.npy",    g)
np.save("output_python/f_id.npy", f_id)
np.save("output_python/g_id.npy", g_id)
np.save("output_python/R.npy",    R)
np.save("output_python/Ref.npy",  Ref)
print("\nSaved to output_python/")
print(f"R: min={R.min():.4f}, max={R.max():.4f}, mean={R.mean():.4f}")

# Show the reflector surface colored by its radial scale.

try:
    from mpl_toolkits.mplot3d import Axes3D   # noqa: F401

    fig = plt.figure(figsize=(9, 7))
    ax  = fig.add_subplot(111, projection="3d")
    sc  = ax.scatter(Ref[:, 0], Ref[:, 1], Ref[:, 2],
                     c=R, cmap="viridis", s=20)
    plt.colorbar(sc, label="R (reflector radius)")
    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
    ax.set_title(f"Reflector Surface — Python  (NK={NK}, k_final={k_final})")
    plt.tight_layout()
    plt.show()
except Exception as e:
    print(f"[plot skipped: {e}]")
