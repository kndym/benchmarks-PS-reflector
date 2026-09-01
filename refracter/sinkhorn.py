"""
refracter/sinkhorn.py

Log-domain Sinkhorn divergence matching C++ BenchmarkCode/main.cpp.

SD_ε(μ,ν) = OT_ε(μ,ν) - ½(OT_ε(μ,μ') + OT_ε(ν',ν))

All iterations use logsumexp. The cost matrix is never fully materialised;
chunk_size controls the block size. Identity iterations are damped (half-step).
"""

import numpy as np
from scipy.special import logsumexp

from .cost import cost_matrix_chunk

# ---------------------------------------------------------------------------
# Helper: chunked logsumexp over one axis of the cost matrix
# ---------------------------------------------------------------------------

def _logsumexp_g_update(x, y, logp, f, g, k, chunk_size, threshold=None):
    """Chunked logsumexp for g update: lse[j] = logsumexp_i(k*f[i] + logp[i] - k*C[i,j])."""
    N = len(y)
    lse_g = np.full(N, -np.inf)

    # a[i] = k*f[i] + logp[i]
    a = k * f + logp  # shape (M,)

    for i_start in range(0, len(x), chunk_size):
        i_end = min(i_start + chunk_size, len(x))
        x_chunk = x[i_start:i_end]
        a_chunk = a[i_start:i_end]  # (chunk,)

        C_chunk = cost_matrix_chunk(x_chunk, y)  # (chunk, N)

        # For each j: logsumexp_i( a[i] - k*C[i,j] )
        # = logsumexp over rows of  (a_chunk[:,None] - k*C_chunk)
        vals = a_chunk[:, None] - k * C_chunk  # (chunk, N)
        if threshold is not None:
            keep = logp[i_start:i_end, None] + k * (
                f[i_start:i_end, None] + g[None, :] - C_chunk
            ) >= threshold
            vals = np.where(keep, vals, -np.inf)
        # logsumexp over the chunk dimension (axis=0), accumulate
        chunk_lse = logsumexp(vals, axis=0)  # (N,)

        # Accumulate across chunks using log-sum-exp identity
        lse_g = np.logaddexp(lse_g, chunk_lse)

    return lse_g


def _logsumexp_f_update(x, y, logq, f, g, k, chunk_size, threshold=None,
                        threshold_g=None):
    """Chunked logsumexp for f update: lse[i] = logsumexp_j(k*g[j] + logq[j] - k*C[i,j])."""
    N = len(x)
    lse_f = np.full(N, -np.inf)

    b = k * g + logq  # shape (M,)  — over y

    for i_start in range(0, len(x), chunk_size):
        i_end = min(i_start + chunk_size, len(x))
        x_chunk = x[i_start:i_end]

        C_chunk = cost_matrix_chunk(x_chunk, y)  # (chunk, N_y)

        # For each i in chunk: logsumexp_j( b[j] - k*C[i,j] )
        vals = b[None, :] - k * C_chunk  # (chunk, N_y)
        if threshold is not None:
            if threshold_g is None:
                threshold_g = g
            keep = logq[None, :] + k * (
                f[i_start:i_end, None] + threshold_g[None, :] - C_chunk
            ) >= threshold
            vals = np.where(keep, vals, -np.inf)
        chunk_lse = logsumexp(vals, axis=1)  # (chunk,)

        lse_f[i_start:i_end] = chunk_lse

    return lse_f


# ---------------------------------------------------------------------------
# Main Sinkhorn step
# ---------------------------------------------------------------------------

def sinkhorn_step(x, y, logp, logq, f, g, k, chunk_size=512):
    """One log-domain Sinkhorn step matching C++ ``Sinkhorn_axb``.

    The C++ implementation updates ``G`` first, then ``F``, absorbs half of
    each log scaling factor into the potentials, and reports
    ``max(abs(log(F) / k))``. Unsupported target/source points retain their
    current potentials, as they do in the C++ loops.
    """
    source_supported = np.isfinite(logp)
    target_supported = np.isfinite(logq)
    threshold = np.log(1e-6 / len(x))

    # C++ G update: use the old f and g, then interpret G as an updated
    # target potential before computing the F update.
    lse_g = _logsumexp_g_update(
        x, y, logp, f, g, k, chunk_size, threshold=threshold
    )
    log_G = np.where(
        target_supported & np.isfinite(lse_g), -k * g - lse_g, 0.0
    )
    g_star = g + log_G / k

    # C++ F update: use the updated G and the old f.
    lse_f = _logsumexp_f_update(
        x, y, logq, f, g_star, k, chunk_size,
        threshold=threshold, threshold_g=g
    )
    log_F = np.where(
        source_supported & np.isfinite(lse_f), -k * f - lse_f, 0.0
    )

    # C++ absorbtion(): f += log(F)/(2*k), g += log(G)/(2*k).
    f_new = f + log_F / (2.0 * k)
    g_new = g + log_G / (2.0 * k)

    # C++ computes maxdif from log(F)/k. There is no additional /k here.
    maxdif = (
        float(np.max(np.abs(log_F[source_supported] / k)))
        if source_supported.any()
        else 0.0
    )

    return f_new, g_new, maxdif


# ---------------------------------------------------------------------------
# Identity Sinkhorn steps (damped)
# ---------------------------------------------------------------------------

def sinkhorn_identity_f_step(x, y, logp, f_id, k, chunk_size=512):
    """Damped identity Sinkhorn step for the source marginal (matches C++ Sinkhorn_identity_F_axb).

    Uses source weights on both sides; update is f_id_new = 0.5*f_id + 0.5*std.
    Returns (f_id_new, maxdif).
    """
    N = len(x)
    threshold = np.log(1e-6 / N)

    b = k * f_id + logp  # weights on y side (same as x for identity problem)
    lse = np.full(N, -np.inf, dtype=np.float64)

    for i_start in range(0, N, chunk_size):
        i_end = min(i_start + chunk_size, N)
        x_chunk = x[i_start:i_end]

        C_chunk = cost_matrix_chunk(x_chunk, y)  # (chunk, N)
        vals = b[None, :] - k * C_chunk          # (chunk, N)
        keep = logp[None, :] + k * (
            f_id[i_start:i_end, None] + f_id[None, :] - C_chunk
        ) >= threshold
        vals = np.where(keep, vals, -np.inf)
        chunk_lse = logsumexp(vals, axis=1)       # (chunk,)
        lse[i_start:i_end] = chunk_lse

    std = np.where(np.isfinite(lse), -lse / k, f_id)
    f_id_new = 0.5 * f_id + 0.5 * std

    # maxdif over supported points only — unsupported are masked by the caller
    log_F = np.where(np.isfinite(lse), -k * f_id - lse, 0.0)
    supported = np.isfinite(logp)
    maxdif = float(np.max(np.abs(log_F[supported] / k))) if supported.any() else 0.0

    return f_id_new, maxdif


def sinkhorn_identity_g_step(x, y, logq, g_id, k, chunk_size=512):
    """Damped identity Sinkhorn step for the target marginal (mirrors f_step with q).

    Returns (g_id_new, maxdif).
    """
    N = len(y)
    threshold = np.log(1e-6 / N)

    a = k * g_id + logq  # weights on x side (same as y for identity problem)
    lse = np.full(N, -np.inf, dtype=np.float64)

    for i_start in range(0, N, chunk_size):
        i_end = min(i_start + chunk_size, N)
        x_chunk = x[i_start:i_end]
        a_chunk = a[i_start:i_end]

        C_chunk = cost_matrix_chunk(x_chunk, y)   # (chunk, N)
        vals = a_chunk[:, None] - k * C_chunk      # (chunk, N)
        keep = logq[i_start:i_end, None] + k * (
            g_id[i_start:i_end, None] + g_id[None, :] - C_chunk
        ) >= threshold
        vals = np.where(keep, vals, -np.inf)
        chunk_lse = logsumexp(vals, axis=0)        # (N,)
        lse = np.logaddexp(lse, chunk_lse)

    std = np.where(np.isfinite(lse), -lse / k, g_id)
    g_id_new = 0.5 * g_id + 0.5 * std

    log_G = np.where(np.isfinite(lse), -k * g_id - lse, 0.0)
    supported = np.isfinite(logq)
    maxdif = float(np.max(np.abs(log_G[supported] / k))) if supported.any() else 0.0

    return g_id_new, maxdif


# ---------------------------------------------------------------------------
# Small-grid Sinkhorn (with kernel pre-conditioning)
# ---------------------------------------------------------------------------

def run_small_sinkhorn(x_s, y_s, p_s, q_s, k_reg):
    """Small-grid (381-pt) Sinkhorn with kernel pre-conditioning (matches C++ smallsinkhorn).

    Pre-conditions marginals via K = exp(-k_reg*C), then runs multi-scale loop.
    Returns (f_small, g_small).
    """
    NK_s = len(x_s)

    C = cost_matrix_chunk(x_s, y_s)  # full (381×381) cost matrix

    # Kernel pre-conditioning
    K = np.exp(-k_reg * C)
    p_s = (K @ p_s); p_s /= p_s.sum()
    q_s = (K @ q_s); q_s /= q_s.sum()

    F_s = np.ones(NK_s)
    G_s = np.ones(NK_s)
    f_s = np.zeros(NK_s)
    g_s = np.zeros(NK_s)

    # C++ increments an int by pow(...), so the fractional increment is
    # truncated on assignment rather than rounded.
    step = max(1, int(np.floor(k_reg ** (1.0 / 3.0))))

    k = 1
    while k < k_reg:
        # The C++ small-grid routine updates F first and then G.
        Gq = G_s * q_s
        exp_mat = np.exp(-k * (C - f_s[:, None] - g_s[None, :]))
        sums_F = exp_mat @ Gq
        F_s = np.where(sums_F > 0, 1.0 / sums_F, 1.0)
        F_s = np.where(p_s > 0, F_s, 1.0)

        Fp = F_s * p_s
        exp_mat = np.exp(-k * (C - f_s[:, None] - g_s[None, :]))
        sums_G = exp_mat.T @ Fp
        G_s = np.where(sums_G > 0, 1.0 / sums_G, 1.0)
        G_s = np.where(q_s > 0, G_s, 1.0)

        f_s += np.log(np.maximum(F_s, 1e-300)) / k
        g_s += np.log(np.maximum(G_s, 1e-300)) / k
        F_s[:] = 1.0
        G_s[:] = 1.0

        k += step

    return f_s, g_s


# ---------------------------------------------------------------------------
# Warm-start via c-transform
# ---------------------------------------------------------------------------

def warmstart_from_small(x_s, y_s, x, y, f_s, g_s, chunk_size=512):
    """Warm-start main grid via c-transforms of the small-grid potentials (C++ UseSmall).

    Returns (f_init, g_init).
    """
    NK = len(x)

    # g_init[j] = min_i(C(x_s[i], y[j]) - f_s[i])
    g_init = np.full(NK, np.inf)
    for j_start in range(0, NK, chunk_size):
        j_end = min(j_start + chunk_size, NK)
        C_block = cost_matrix_chunk(x_s, y[j_start:j_end])
        g_init[j_start:j_end] = np.min(C_block - f_s[:, None], axis=0)

    # f_init[i] = min_j(C(x[i], y_s[j]) - g_s[j])
    f_init = np.full(NK, np.inf)
    for i_start in range(0, NK, chunk_size):
        i_end = min(i_start + chunk_size, NK)
        C_block = cost_matrix_chunk(x[i_start:i_end], y_s)
        f_init[i_start:i_end] = np.min(C_block - g_s[None, :], axis=1)

    return f_init, g_init


# ---------------------------------------------------------------------------
# Full Sinkhorn divergence
# ---------------------------------------------------------------------------

def run_sinkhorn_divergence(x, y, p, q, chunk_size=512, verbose=True):
    """Run the full divergence when only main-grid weights are available.

    The C++ executable receives separate small-grid density evaluations from
    its benchmark function.  The generic Python API only receives ``p`` and
    ``q``, so it uses the main grid as its own warm-start grid instead of
    silently loading incompatible benchmark data.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    if len(x) != len(y) or len(p) != len(x) or len(q) != len(y):
        raise ValueError("run_sinkhorn_divergence requires equal-sized x/y grids and weights")
    return _run_sinkhorn_divergence_inner(
        x, y, p, q, x, y, p, q,
        chunk_size=chunk_size,
        verbose=verbose,
    )


def _run_sinkhorn_divergence_inner(x, y, p, q, x_s, y_s, p_s, q_s,
                                   chunk_size=512, verbose=True):
    """Full Sinkhorn divergence (C++ do_sinkhorn_subtracted_axb).

    Returns dict with keys: f, g, f_id, g_id, total_cost.
    """
    NK = len(x)

    # --- Normalise weights ---
    p = np.array(p, dtype=np.float64)
    q = np.array(q, dtype=np.float64)
    p = p / p.sum()
    q = q / q.sum()

    logp = np.full_like(p, -np.inf)
    logq = np.full_like(q, -np.inf)
    np.log(p, out=logp, where=p > 0)
    np.log(q, out=logq, where=q > 0)

    p_s = np.array(p_s, dtype=np.float64)
    q_s = np.array(q_s, dtype=np.float64)
    p_s = p_s / p_s.sum()
    q_s = q_s / q_s.sum()

    multiplier = 8
    k_small = multiplier * int(np.floor(np.sqrt(len(x_s))))  # 152
    k_final  = multiplier * int(np.floor(np.sqrt(NK)))        # 1024
    # C++ stores each fractional increment back into an int (truncation).
    step = int(np.floor(k_final ** (1.0 / 3.0)))               # 10
    cap_iter = 16
    cap_thr  = 1e-5

    if verbose:
        print(f"Regularisation: k_small={k_small}, k_final={k_final}, step={step}")
        print(f"NK={NK}, source support={(p>0).sum()}, target support={(q>0).sum()}")

    # --- 1. Small Sinkhorn ---
    if verbose:
        print("\nRunning small Sinkhorn (warm-start)...")
    f_s, g_s = run_small_sinkhorn(x_s, y_s, p_s.copy(), q_s.copy(), k_small)

    # --- 2. Warm-start main grid ---
    if verbose:
        print("Computing warm-start c-transforms...")
    f, g = warmstart_from_small(x_s, y_s, x, y, f_s, g_s, chunk_size=chunk_size)

    # --- 3. Multi-scale Sinkhorn loop ---
    if verbose:
        print("\nMulti-scale Sinkhorn loop:")
    it = 0
    regvar = k_small
    maxdif = -1.0
    while regvar < k_final:
        f, g, maxdif = sinkhorn_step(x, y, logp, logq, f, g, regvar, chunk_size)
        it += 1
        if verbose:
            print(f"  iter {it:4d}, k={regvar:5d}, maxdif={maxdif:.4e}")
        regvar += step

    # --- 4. Final loop at k_final ---
    if verbose:
        print(f"\nFinal Sinkhorn at k={k_final}:")
    i = 0
    maxdif = cap_thr + 1.0  # ensure at least one iteration
    while maxdif > cap_thr:
        f, g, maxdif = sinkhorn_step(x, y, logp, logq, f, g, k_final, chunk_size)
        i += 1
        if verbose:
            print(f"  iter {i:3d}, maxdif={maxdif:.4e}")
        # C++ checks i > cap_iteration after the update, so cap_iter=16
        # permits 17 final iterations in this main.cpp-compatible loop.
        if i > cap_iter:
            break
    if verbose:
        print(f"Final loop: {i} iterations, last change={maxdif:.4e}")

    # --- 5. Identity F loop ---
    if verbose:
        print("\nIdentity Sinkhorn for source (f_id):")
    f_id = np.zeros(NK)
    # C++ uses regvariable += sqrt(k) with an integer regvariable, so this
    # is truncation/floor rather than round().
    id_step = int(np.floor(np.sqrt(k_final)))  # 32
    regvar = 1
    i = 0
    while regvar < k_final:
        f_id_new, maxdif = sinkhorn_identity_f_step(
            x, y, logp, f_id, regvar, chunk_size
        )
        f_id = np.where(np.isfinite(logp), f_id_new, f_id)
        i += 1
        if verbose:
            print(f"  iter {i:4d}, k={regvar:5d}, maxdif={maxdif:.4e}")
        regvar += id_step

    if verbose:
        print(f"  Final {cap_iter} iterations at k={k_final}:")
    i = 0
    maxdif = cap_thr + 1.0
    while maxdif > cap_thr:
        f_id_new, maxdif = sinkhorn_identity_f_step(
            x, y, logp, f_id, k_final, chunk_size
        )
        f_id = np.where(np.isfinite(logp), f_id_new, f_id)
        i += 1
        if verbose:
            print(f"  iter {i:3d}, maxdif={maxdif:.4e}")
        if i > cap_iter:
            break
    if verbose:
        print(f"Identity F: {i} final iterations, last change={maxdif:.4e}")

    # --- 6. Identity G loop ---
    if verbose:
        print("\nIdentity Sinkhorn for target (g_id):")
    g_id = np.zeros(NK)
    regvar = 1
    i = 0
    while regvar < k_final:
        g_id_new, maxdif = sinkhorn_identity_g_step(
            x, y, logq, g_id, regvar, chunk_size
        )
        g_id = np.where(np.isfinite(logq), g_id_new, g_id)
        i += 1
        if verbose:
            print(f"  iter {i:4d}, k={regvar:5d}, maxdif={maxdif:.4e}")
        regvar += id_step

    if verbose:
        print(f"  Final {cap_iter} iterations at k={k_final}:")
    i = 0
    maxdif = cap_thr + 1.0
    while maxdif > cap_thr:
        g_id_new, maxdif = sinkhorn_identity_g_step(
            x, y, logq, g_id, k_final, chunk_size
        )
        g_id = np.where(np.isfinite(logq), g_id_new, g_id)
        i += 1
        if verbose:
            print(f"  iter {i:3d}, maxdif={maxdif:.4e}")
        if i > cap_iter:
            break
    if verbose:
        print(f"Identity G: {i} final iterations, last change={maxdif:.4e}")

    # Normalise: shift max to 0, then subtract identity terms
    max_f_id = np.max(f_id)
    f_id = f_id - max_f_id
    g_id = g_id + max_f_id

    max_f = np.max(f)
    f = f - max_f
    g = g + max_f

    f = f - f_id
    g = g - g_id

    total_cost = float(
        np.sum(p[p > 0] * f[p > 0]) +
        np.sum(q[q > 0] * g[q > 0])
    )

    if verbose:
        print(f"\nTotal cost (approx): {total_cost:.6e}")

    return {
        "f": f,
        "g": g,
        "f_id": f_id,
        "g_id": g_id,
        "total_cost": total_cost,
    }
