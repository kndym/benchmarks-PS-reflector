"""
refractor_sinkhorn.py — Sinkhorn divergence for the far-field refraction problem.

Cost: c(x, y) = -log(1 - κ · (x·y))  with κ = 0.6
Points sampled uniformly (Halton QMC) on spherical patches:
  Source: θ ∈ [0°,  60°], φ ∈ [0°, 360°]
  Target: θ ∈ [120°,180°], φ ∈ [0°, 360°]

Usage:
    python refractor_sinkhorn.py
"""

import os
import time
import numpy as np
from scipy.special import logsumexp

# ---------------------------------------------------------------------------
# Refraction cost
# ---------------------------------------------------------------------------

KAPPA = 0.6
EPS_CLIP = 1e-15


def cost_matrix_chunk(x_chunk, y):
    """Refraction cost c(x,y) = -log(1 - κ·(x·y)) for x_chunk (M,3), y (N,3)."""
    dots = x_chunk @ y.T                          # (M, N)
    return -np.log(np.clip(1.0 - KAPPA * dots, EPS_CLIP, None))


# ---------------------------------------------------------------------------
# Halton sequence
# ---------------------------------------------------------------------------

def _halton(index, base):
    result, f, i = 0.0, 1.0, index
    while i > 0:
        f /= base
        result += f * (i % base)
        i //= base
    return result


def gen_spherical_patch(n, theta_min_deg, theta_max_deg,
                        phi_min_deg=0.0, phi_max_deg=360.0,
                        base2=2, base3=3, skip=0):
    """Sample n points uniformly on a spherical patch using Halton QMC.

    Uses the equal solid-angle transform:
        cos(θ) = cos(θ_max) + u1 * (cos(θ_min) - cos(θ_max))
        φ      = φ_min + u2 * (φ_max - φ_min)

    Parameters
    ----------
    theta_min_deg, theta_max_deg : polar angle range in degrees
    phi_min_deg, phi_max_deg     : azimuthal angle range in degrees
    base2, base3 : Halton bases for u1, u2
    skip : number of Halton indices to skip at the start
    """
    cos_max = np.cos(np.radians(theta_min_deg))  # larger cos (smaller θ)
    cos_min = np.cos(np.radians(theta_max_deg))  # smaller cos (larger θ)
    phi_min = np.radians(phi_min_deg)
    phi_max = np.radians(phi_max_deg)

    pts = []
    idx = skip
    while len(pts) < n:
        u1 = _halton(idx, base2)
        u2 = _halton(idx, base3)
        cos_theta = cos_min + u1 * (cos_max - cos_min)
        sin_theta = np.sqrt(max(0.0, 1.0 - cos_theta ** 2))
        phi = phi_min + u2 * (phi_max - phi_min)
        pts.append([sin_theta * np.cos(phi),
                    sin_theta * np.sin(phi),
                    cos_theta])
        idx += 1
    return np.array(pts, dtype=np.float64)


# ---------------------------------------------------------------------------
# Log-domain logsumexp helpers (identical structure to reflector/sinkhorn.py
# but calling the local refraction cost_matrix_chunk)
# ---------------------------------------------------------------------------

def _logsumexp_g_update(x, y, logp, f, g, k, chunk_size):
    N = len(y)
    lse_g = np.full(N, -np.inf)
    a = k * f + logp
    for i_start in range(0, len(x), chunk_size):
        i_end = min(i_start + chunk_size, len(x))
        x_chunk = x[i_start:i_end]
        a_chunk = a[i_start:i_end]
        C_chunk = cost_matrix_chunk(x_chunk, y)
        vals = a_chunk[:, None] - k * C_chunk
        chunk_lse = logsumexp(vals, axis=0)
        lse_g = np.logaddexp(lse_g, chunk_lse)
    return lse_g


def _logsumexp_f_update(x, y, logq, f, g, k, chunk_size):
    N = len(x)
    lse_f = np.full(N, -np.inf)
    b = k * g + logq
    for i_start in range(0, len(x), chunk_size):
        i_end = min(i_start + chunk_size, len(x))
        x_chunk = x[i_start:i_end]
        C_chunk = cost_matrix_chunk(x_chunk, y)
        vals = b[None, :] - k * C_chunk
        chunk_lse = logsumexp(vals, axis=1)
        lse_f[i_start:i_end] = chunk_lse
    return lse_f


# ---------------------------------------------------------------------------
# Sinkhorn step functions (adapted from reflector/sinkhorn.py)
# ---------------------------------------------------------------------------

def sinkhorn_step(x, y, logp, logq, f, g, k, chunk_size=512):
    """Standard log-domain Sinkhorn step (no damping)."""
    lse_g = _logsumexp_g_update(x, y, logp, f, g, k, chunk_size)
    log_G = -k * g - lse_g
    g_new = g + log_G / (2.0 * k)

    g_eff = -lse_g / k
    lse_f = _logsumexp_f_update(x, y, logq, f, g_eff, k, chunk_size)
    log_F = -k * f - lse_f
    f_new = f + log_F / (2.0 * k)

    maxdif = float(np.max(np.abs(log_F / k)))
    return f_new, g_new, maxdif


def sinkhorn_identity_f_step(x, y, logp, f_id, k, chunk_size=512):
    """Damped identity Sinkhorn step for the source marginal."""
    N = len(x)
    b = k * f_id + logp
    lse = np.full(N, -np.inf, dtype=np.float64)

    for i_start in range(0, N, chunk_size):
        i_end = min(i_start + chunk_size, N)
        x_chunk = x[i_start:i_end]
        C_chunk = cost_matrix_chunk(x_chunk, y)
        vals = b[None, :] - k * C_chunk
        chunk_lse = logsumexp(vals, axis=1)
        lse[i_start:i_end] = chunk_lse

    std = -lse / k
    f_id_new = 0.5 * f_id + 0.5 * std

    log_F = -k * f_id - lse
    supported = np.isfinite(logp)
    maxdif = float(np.max(np.abs(log_F[supported] / k))) if supported.any() else 0.0
    return f_id_new, maxdif


def sinkhorn_identity_g_step(x, y, logq, g_id, k, chunk_size=512):
    """Damped identity Sinkhorn step for the target marginal."""
    N = len(y)
    a = k * g_id + logq
    lse = np.full(N, -np.inf, dtype=np.float64)

    for i_start in range(0, N, chunk_size):
        i_end = min(i_start + chunk_size, N)
        x_chunk = x[i_start:i_end]
        a_chunk = a[i_start:i_end]
        C_chunk = cost_matrix_chunk(x_chunk, y)
        vals = a_chunk[:, None] - k * C_chunk
        chunk_lse = logsumexp(vals, axis=0)
        lse = np.logaddexp(lse, chunk_lse)

    std = -lse / k
    g_id_new = 0.5 * g_id + 0.5 * std

    log_G = -k * g_id - lse
    supported = np.isfinite(logq)
    maxdif = float(np.max(np.abs(log_G[supported] / k))) if supported.any() else 0.0
    return g_id_new, maxdif


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

NK         = 400
chunk      = 512
multiplier = 8
k_final    = multiplier * int(np.floor(np.sqrt(NK)))   # 8*20 = 160
id_step    = int(np.floor(np.sqrt(k_final)))            # floor(sqrt(160)) = 12
cap_iter   = 16
cap_thr    = 1e-5
out_dir    = "output_refraction"

os.makedirs(out_dir, exist_ok=True)
print(f"NK={NK}, k_final={k_final}, id_step={id_step}")

# ---------------------------------------------------------------------------
# Cloud generation — uniform on spherical patches
# ---------------------------------------------------------------------------

print("Generating Halton QMC clouds on spherical patches...")
x = gen_spherical_patch(NK, theta_min_deg=0.0,   theta_max_deg=60.0)
y = gen_spherical_patch(NK, theta_min_deg=120.0, theta_max_deg=180.0)

print(f"  Source: {len(x)} points, z range [{x[:,2].min():.4f}, {x[:,2].max():.4f}]")
print(f"  Target: {len(y)} points, z range [{y[:,2].min():.4f}, {y[:,2].max():.4f}]")

# Uniform weights
p = np.full(NK, 1.0 / NK)
q = np.full(NK, 1.0 / NK)
logp = np.log(p)
logq = np.log(q)

# ---------------------------------------------------------------------------
# Step 1 — Final Sinkhorn at k_final  (no warm-start, f=g=0)
# ---------------------------------------------------------------------------

f = np.zeros(NK, dtype=np.float64)
g = np.zeros(NK, dtype=np.float64)

print(f"\nFinal Sinkhorn at k={k_final}:")
maxdif = cap_thr + 1.0
i = 0
t0 = time.time()
while maxdif > cap_thr:
    f, g, maxdif = sinkhorn_step(x, y, logp, logq, f, g, k_final, chunk)
    i += 1
    print(f"  iter {i:3d}, maxdif={maxdif:.4e}  ({time.time()-t0:.2f}s)")
    if i >= cap_iter:
        break
print(f"Final: {i} iters, last maxdif={maxdif:.4e}")

# ---------------------------------------------------------------------------
# Step 2 — Identity F
# ---------------------------------------------------------------------------

print(f"\nIdentity F (id_step={id_step}):")
f_id   = np.zeros(NK, dtype=np.float64)
regvar = 1
it     = 0
t0     = time.time()
while regvar < k_final:
    f_id_new, md = sinkhorn_identity_f_step(x, y, logp, f_id, regvar, chunk)
    f_id = f_id_new          # all points are supported (uniform weights)
    it += 1
    print(f"  iter {it:3d}, k={regvar:4d}, maxdif={md:.4e}")
    regvar += id_step

maxdif = cap_thr + 1.0
i = 0
while maxdif > cap_thr:
    f_id_new, maxdif = sinkhorn_identity_f_step(x, y, logp, f_id, k_final, chunk)
    f_id = f_id_new
    i += 1
    if i >= cap_iter:
        break
print(f"Identity F: {i} final iters, last maxdif={maxdif:.4e}  ({time.time()-t0:.1f}s)")

# ---------------------------------------------------------------------------
# Step 3 — Identity G
# ---------------------------------------------------------------------------

print(f"\nIdentity G (id_step={id_step}):")
g_id   = np.zeros(NK, dtype=np.float64)
regvar = 1
it     = 0
t0     = time.time()
while regvar < k_final:
    g_id_new, md = sinkhorn_identity_g_step(x, y, logq, g_id, regvar, chunk)
    g_id = g_id_new
    it += 1
    print(f"  iter {it:3d}, k={regvar:4d}, maxdif={md:.4e}")
    regvar += id_step

maxdif = cap_thr + 1.0
i = 0
while maxdif > cap_thr:
    g_id_new, maxdif = sinkhorn_identity_g_step(x, y, logq, g_id, k_final, chunk)
    g_id = g_id_new
    i += 1
    if i >= cap_iter:
        break
print(f"Identity G: {i} final iters, last maxdif={maxdif:.4e}  ({time.time()-t0:.1f}s)")

# ---------------------------------------------------------------------------
# Normalise and compute total cost
# ---------------------------------------------------------------------------

mx_fid = float(np.max(f_id))
f_id  -= mx_fid
g_id  += mx_fid

mx_f = float(np.max(f))
f   -= mx_f
g   += mx_f

f -= f_id
g -= g_id

total_cost = float(np.sum(p * f) + np.sum(q * g))
print(f"\nTotal cost: {total_cost:.6e}")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

np.save(f"{out_dir}/f.npy",    f)
np.save(f"{out_dir}/g.npy",    g)
np.save(f"{out_dir}/f_id.npy", f_id)
np.save(f"{out_dir}/g_id.npy", g_id)
np.save(f"{out_dir}/x.npy",    x)
np.save(f"{out_dir}/y.npy",    y)
print(f"Saved to {out_dir}/")
