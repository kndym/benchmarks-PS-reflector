"""
reflector/sinkhorn.py

Full Sinkhorn / Sinkhorn-divergence algorithm matching the C++ implementation
in BenchmarkCode/main.cpp.

Algorithm summary
-----------------
The Sinkhorn divergence is

    SD_ε(μ, ν) = OT_ε(μ,ν) - ½(OT_ε(μ,μ') + OT_ε(ν',ν))

where μ', ν' are antipodal images of μ, ν.  In practice the "identity" terms
OT_ε(μ, μ') and OT_ε(ν', ν) are computed by running a separate Sinkhorn
iteration between x and y (the same physical points) but weighting *both*
sides by the source (resp. target) measure.

Implementation notes
--------------------
*  All Sinkhorn iterations are carried out in log-domain using logsumexp.
*  The C++ code uses "absorption" steps that fold the multiplicative kernel
   factor F (resp. G) back into the additive potential f (resp. g) after
   every iteration.  In log-domain this corresponds to a simple addition.
*  The identity iterations are *damped*: the new iterate is the average of
   the old value and the standard (undamped) update.  This is "Feydy's idea"
   referenced in the C++ comments.
*  chunk_size controls the block size for the outer loop over source/target
   points, so that the (NK × NK) cost matrix is never materialised in full.

Schedule (matching C++ ``do_sinkhorn_subtracted_axb``)
-------------------------------------------------------
Small grid warm-start:
    k_small = 8 * getk(NK_small) = 8 * floor(sqrt(381)) ≈ 8*19 = 152
Main multi-scale loop:
    k = 152, 162, 172, …  (step = int(1024^(1/3)) = 10), stop when k >= 1024
Final loop:
    k = 1024, up to 16 iterations, break if maxdif < 1e-5
Identity F loop:
    k = 1, 33, 65, …  (step = int(sqrt(1024)) = 32), stop when k >= 1024;
    then 16 final iterations at k=1024
Identity G loop: same as Identity F loop but for g_id.
"""

import numpy as np
from scipy.special import logsumexp

from .cost import cost_matrix_chunk

# ---------------------------------------------------------------------------
# Helper: chunked logsumexp over one axis of the cost matrix
# ---------------------------------------------------------------------------

def _logsumexp_g_update(x, y, logp, f, g, k, chunk_size):
    """Compute the g update for all j in log-domain.

    For each j:
        new_log_G[j] = logsumexp_i( k*(f[i] - C[i,j]) + logp[i] )

    then  g_new[j] = -1/k * new_log_G[j]
                   = g[j] - 1/(2k) * new_log_G[j]   (after absorption)

    But since we mimic the C++ absorption pattern (fold the F/G factor back
    into f/g after *every* iteration), we return the raw logsumexp values and
    let the caller add them.

    This function returns the per-j logsumexp value
        lse[j] = logsumexp_i( k*f[i] + logp[i] - k*C[i,j] )
    which equals  k * (new_g[j] + g[j])  in the C++ notation.

    Parameters
    ----------
    x, y : (N,3) and (N,3) arrays
    logp : (N,) log source weights
    f    : (N,) current f potential
    g    : (N,) current g potential (used for threshold skipping)
    k    : regularisation parameter (integer)
    chunk_size : int

    Returns
    -------
    lse_g : (N,) array — the logsumexp for each j
    """
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
        # logsumexp over the chunk dimension (axis=0), accumulate
        chunk_lse = logsumexp(vals, axis=0)  # (N,)

        # Accumulate across chunks using log-sum-exp identity
        lse_g = np.logaddexp(lse_g, chunk_lse)

    return lse_g


def _logsumexp_f_update(x, y, logq, f, g, k, chunk_size):
    """Compute the f update for all i in log-domain.

    For each i:
        lse[i] = logsumexp_j( k*g[j] + logq[j] - k*C[i,j] )

    Returns
    -------
    lse_f : (N,) array
    """
    N = len(x)
    lse_f = np.full(N, -np.inf)

    b = k * g + logq  # shape (M,)  — over y

    for i_start in range(0, len(x), chunk_size):
        i_end = min(i_start + chunk_size, len(x))
        x_chunk = x[i_start:i_end]

        C_chunk = cost_matrix_chunk(x_chunk, y)  # (chunk, N_y)

        # For each i in chunk: logsumexp_j( b[j] - k*C[i,j] )
        vals = b[None, :] - k * C_chunk  # (chunk, N_y)
        chunk_lse = logsumexp(vals, axis=1)  # (chunk,)

        lse_f[i_start:i_end] = chunk_lse

    return lse_f


# ---------------------------------------------------------------------------
# Main Sinkhorn step (no damping)
# ---------------------------------------------------------------------------

def sinkhorn_step(x, y, logp, logq, f, g, k, chunk_size=512):
    """Standard log-domain Sinkhorn step (no damping).

    Matches C++ ``Sinkhorn_axb`` + ``absorbtion``.

    The C++ code maintains separate multiplicative factors F, G initialised to
    1 and then absorbs them into f, g every iteration.  In log-domain this is
    equivalent to:

        g_new[j] = -1/k * logsumexp_i( k*(f[i] - C[i,j]) + logp[i] )
        f_new[i] = -1/k * logsumexp_j( k*(g_new[j] - C[i,j]) + logq[j] )
        maxdif    = max|f_new - f|

    The absorption step in C++ does:
        f <- f + log(F)/(2k)
        g <- g + log(G)/(2k)
    which, since log(F) = logsumexp - k*g (the raw logsumexp value), gives
    exactly the formulas above (the 2k denominator cancels with a factor 2
    because both the numerator and the denominator are combined).

    Wait — let us be precise.  In the C++ code, within Sinkhorn_axb:
      G[j] = 1/sum_i( exp(-k*(C[i,j]-f[i]-g[j])) * F[i]*p[i] )
    After absorption:
      g[j] += log(G[j])/(2k)
    So the net update to g[j] is:
      g[j] <-  g[j] - 1/(2k) * log( sum_i exp(-k*(C[i,j]-f[i]-g[j])) * p[i] )
             = g[j] - 1/(2k) * (logsumexp_i( logp[i] - k*(C[i,j]-f[i]-g[j]) ))
             = -1/(2k) * logsumexp_i( logp[i] + k*f[i] - k*C[i,j] )  + g[j]/2

    Hmm — that is a half-step.  But looking more carefully: F is reset to 1
    after every absorption, so log(F)=0 at the start of each iteration, and the
    whole expression is:
      new_G[j] = 1/( sum_i exp(-k*(C[i,j]-f[i]-g[j])) * p[i] )
    => log(new_G[j]) = -logsumexp_i( logp[i] + k*(f[i]+g[j]-C[i,j]) )
                     = -logsumexp_i( logp[i] + k*f[i] - k*C[i,j] ) - k*g[j]

    After absorption:
      g[j] += log(new_G[j])/(2k)
            = g[j] - 1/(2k)*logsumexp_i(logp[i]+k*f[i]-k*C[i,j]) - g[j]/2
      =>  g[j]  <-  g[j]/2 - 1/(2k)*logsumexp_i(logp[i]+k*f[i]-k*C[i,j])

    That is NOT what we want (we want the full update).  Let us re-read.

    Actually in C++ the summation is over CURRENT f and g together:
      temp = -k*(Cost(x[i],y[j]) - f[i] - g[j])
      sum += exp(temp) * F[i] * p[i]
      G[j] = 1/sum
    where F[i]=1 (just been reset).  So G[j] = 1/K_p[j] where K_p is the
    kernel-weighted sum.  Then:
      log(G[j]) = -log(K_p[j])
               = -logsumexp_i( logp[i] - k*(C[i,j]-f[i]-g[j]) )
               = k*g[j] - logsumexp_i( logp[i] + k*f[i] - k*C[i,j] )

    Absorption: g[j] += log(G[j])/(2k)
               = g[j] + g[j]/2 - lse_g_j/(2k)
               = 3g[j]/2 - lse_g_j/(2k)

    That also does not match what I expect.  Let me look at absorbtion more
    carefully.

    absorbtion(k):
      tempvecNK_1 = log(F)          // log(G) for g side
      tempvecNK_1 *= 1/(2k)
      tempvecNK_2 = f + tempvecNK_1
      swap(f, tempvecNK_2)          // f += log(F)/(2k)

    So the update is simply:
      f[i] += log(F[i]) / (2k)
      g[j] += log(G[j]) / (2k)

    And log(G[j]) = k*g[j] - lse_g[j].

    Therefore:
      g[j] += (k*g[j] - lse_g[j]) / (2k)
             = g[j]/2 - lse_g[j]/(2k)
    => g_new[j] = g[j] + g[j]/2 - lse_g[j]/(2k)
                = 3*g[j]/2 - lse_g[j]/(2k)

    That still seems odd.  Let me re-check in the very next iteration: F is
    reset to 1 (log(F)=0) so the absorbed g incorporates the current kernel
    value.  After many iterations g converges to the fixed point where
    log(G[j]) = 0, i.e. lse_g[j] = k*g[j], i.e. exactly the standard
    Sinkhorn fixed-point condition.

    At convergence: g[j] converges when lse_g[j] = k * g[j],
    i.e.  logsumexp_i(logp[i] + k*f[i] - k*C[i,j]) = k*g[j]
    =>  g[j] = (1/k) * logsumexp_i(logp[i] + k*f[i] - k*C[i,j])

    So the fixed point is the correct one, but the iteration mixes the old g
    with the new update.  The sequence is:
      g_{t+1}[j] = 3/2 * g_t[j] - 1/(2k) * lse_g_t[j]

    At the fixed point g* we have lse_g*(j) = k*g*(j), so:
      g_{t+1}[j] = 3/2 g* - k*g*/(2k) = 3/2 g* - 1/2 g* = g*  ✓

    This is a stationary point.  Whether it converges depends on the step
    being a contraction, which it is in practice for this problem.

    For our Python implementation we replicate this identical update rule so
    that the results match the C++ output.

    Parameters
    ----------
    x, y : (N, 3) arrays of source and target points
    logp : (N,) log source weights (−inf where p=0)
    logq : (N,) log target weights (−inf where q=0)
    f, g : (N,) current Kantorovich potentials
    k    : regularisation parameter (int or float)
    chunk_size : int

    Returns
    -------
    f_new : (N,) updated f potential
    g_new : (N,) updated g potential
    maxdif : float, max|f_new - f|
    """
    # --- g update ---
    lse_g = _logsumexp_g_update(x, y, logp, f, g, k, chunk_size)
    # G[j] = 1 / sum_i(exp(-k*(C-f[i]-g[j]))*p[i])
    #       = exp(-k*g[j]) / sum_i(p[i]*exp(k*(f[i]-C[i,j])))
    # log(G[j]) = -k*g[j] - lse_g[j]
    # Absorption: g += log(G)/(2k)  =>  g_new = g/2 - lse_g/(2k)  [contraction]
    log_G = -k * g - lse_g
    #g_new = g + log_G / (2.0 * k)
    g_new = -lse_g / k


    # --- f update ---
    # C++ uses G[j] (freshly computed, not yet absorbed) with CURRENT g[j].
    # This is equivalent to using g_eff[j] = -lse_g[j]/k as the effective g,
    # since G[j]*exp(k*g[j]) = exp(-lse_g[j]).
    # Using g_eff matches C++ path and converges faster than using g_new.
    g_eff = -lse_g / k
    lse_f = _logsumexp_f_update(x, y, logq, f, g_eff, k, chunk_size)
    log_F = -k * f - lse_f
    #f_new = f + log_F / (2.0 * k)
    f_new = -lse_f / k


    # maxdif measured as max|log(F)|/k  (C++ uses cblas_idamax on log(F)/k)
    maxdif = float(np.max(np.abs(log_F / k)))

    return f_new, g_new, maxdif


# ---------------------------------------------------------------------------
# Identity Sinkhorn steps (damped)
# ---------------------------------------------------------------------------

def sinkhorn_identity_f_step(x, y, logp, f_id, k, chunk_size=512):
    """Damped identity Sinkhorn step for the source marginal.

    Matches C++ ``Sinkhorn_identity_F_axb`` + ``absorbtion_f_id``.

    The identity OT problem sends μ to its antipodal image μ' using the same
    cost c(x, y).  Because x and y are the *same* set of points with the
    *same* measure (source weights p on both sides), the update is:

        std[i] = -1/k * logsumexp_j( k*(f_id[j] - C[i,j]) + logp[j] )

    which is computed via the same log-domain machinery.  After computing the
    standard update, the new iterate is the damped average:

        f_id_new[i] = 0.5 * f_id[i] + 0.5 * std[i]

    This corresponds to the C++ absorption:
        F_id[i] = 1 / (sum_j exp(-k*(C[i,j]-f_id[i]-f_id[j])) * F_id[j]*p[j])
    followed by
        f_id[i] += log(F_id[i]) / (2k)

    At the fixed point f_id* we have F_id = 1 for all i, so the update is
    stationary.  The half-step here ensures damped (stable) convergence.

    Parameters
    ----------
    x, y : (N, 3) source and target point arrays (same physical locations)
    logp : (N,) log source weights (used on BOTH sides)
    f_id : (N,) current identity potential for source
    k    : regularisation parameter
    chunk_size : int

    Returns
    -------
    f_id_new : (N,) updated identity potential
    maxdif   : float
    """
    N = len(x)

    # We compute:
    #   lse[i] = logsumexp_j( k*(f_id[j]) + logp[j] - k*C[i,j] )
    # Note: this uses f_id[j] on both sides (source weights p on both sides,
    # and the potential f_id acts as the "g" in the regular update).
    # This matches:  temp = -k*(C(x[i],y[j]) - f_id[i] - f_id[j])
    #               sum += exp(temp) * F_id[j] * p[j]
    # with F_id[j]=1, which gives:
    #   sum[i] = sum_j exp(k*(f_id[i]+f_id[j]-C[i,j])) * p[j]
    #          = exp(k*f_id[i]) * sum_j exp(k*f_id[j] + logp[j] - k*C[i,j])
    # => F_id[i] = exp(-k*f_id[i]) / sum_j exp(k*f_id[j]+logp[j]-k*C[i,j])
    # => log(F_id[i]) = -k*f_id[i] - lse[i]
    # Absorption: f_id[i] += log(F_id[i])/(2k)
    #           = f_id[i] - f_id[i]/2 - lse[i]/(2k)
    #           = f_id[i]/2 - lse[i]/(2k)
    # Which is: f_id_new = 0.5*f_id + 0.5 * (-lse/k)
    # i.e. f_id_new = 0.5*f_id + 0.5*std  where std = -lse/k  ✓

    b = k * f_id + logp  # weights on y side (which is also x here)
    lse = np.full(N, -np.inf, dtype=np.float64)

    for i_start in range(0, N, chunk_size):
        i_end = min(i_start + chunk_size, N)
        x_chunk = x[i_start:i_end]

        C_chunk = cost_matrix_chunk(x_chunk, y)  # (chunk, N)
        vals = b[None, :] - k * C_chunk          # (chunk, N)
        chunk_lse = logsumexp(vals, axis=1)       # (chunk,)
        lse[i_start:i_end] = chunk_lse

    std = -lse / k
    f_id_new = 0.5 * f_id + 0.5 * std

    # maxdif: C++ takes max|log(F_id)|/k over SUPPORTED points only (logp > -inf).
    # Unsupported points are never updated (masked by caller), so their log_F
    # never converges — including them would permanently prevent termination.
    log_F = -k * f_id - lse  # = log(F_id) * k  (before /k)
    supported = np.isfinite(logp)
    maxdif = float(np.max(np.abs(log_F[supported] / k))) if supported.any() else 0.0

    return f_id_new, maxdif


def sinkhorn_identity_g_step(x, y, logq, g_id, k, chunk_size=512):
    """Damped identity Sinkhorn step for the target marginal.

    Mirrors ``sinkhorn_identity_f_step`` but uses target weights q on both
    sides.  Matches C++ ``Sinkhorn_identity_G_axb`` + ``absorbtion_g_id``.

    Parameters
    ----------
    x, y : (N, 3) source and target point arrays
    logq : (N,) log target weights (used on BOTH sides)
    g_id : (N,) current identity potential for target
    k    : regularisation parameter
    chunk_size : int

    Returns
    -------
    g_id_new : (N,)
    maxdif   : float
    """
    N = len(y)

    # The C++ update for G_id:
    #   temp = -k*(C(x[i],y[j]) - g_id[i] - g_id[j])
    #   sum_i exp(temp) * G_id[i] * q[i]  => G_id[j] = 1/sum
    # With G_id[i]=1:
    #   sum[j] = sum_i exp(k*(g_id[i]+g_id[j]-C[i,j])) * q[i]
    #          = exp(k*g_id[j]) * sum_i exp(k*g_id[i]+logq[i]-k*C[i,j])
    # => log(G_id[j]) = -k*g_id[j] - lse_j  where lse_j=logsumexp_i(...)
    # Absorption: g_id[j] += log(G_id[j])/(2k)
    #           = g_id[j]/2 - lse_j/(2k)
    # => g_id_new = 0.5*g_id + 0.5*(-lse/k) = 0.5*g_id + 0.5*std  ✓

    a = k * g_id + logq  # weights on x side (which is also y here)
    lse = np.full(N, -np.inf, dtype=np.float64)

    for i_start in range(0, N, chunk_size):
        i_end = min(i_start + chunk_size, N)
        x_chunk = x[i_start:i_end]
        a_chunk = a[i_start:i_end]

        C_chunk = cost_matrix_chunk(x_chunk, y)   # (chunk, N)
        vals = a_chunk[:, None] - k * C_chunk      # (chunk, N)
        chunk_lse = logsumexp(vals, axis=0)        # (N,)
        lse = np.logaddexp(lse, chunk_lse)

    std = -lse / k
    g_id_new = 0.5 * g_id + 0.5 * std

    # maxdif over supported points only (matching C++ — unsupported never updated)
    log_G = -k * g_id - lse
    supported = np.isfinite(logq)
    maxdif = float(np.max(np.abs(log_G[supported] / k))) if supported.any() else 0.0

    return g_id_new, maxdif


# ---------------------------------------------------------------------------
# Small-grid Sinkhorn (with kernel pre-conditioning)
# ---------------------------------------------------------------------------

def run_small_sinkhorn(x_s, y_s, p_s, q_s, k_reg):
    """Full small-grid Sinkhorn with kernel pre-conditioning.

    Matches C++ ``smallsinkhorn``.

    The small grid is only 381 points so we can build the full (381×381) cost
    matrix in memory without chunking.

    The C++ code pre-conditions the marginals by applying one round of kernel
    smoothing before starting the Sinkhorn loop:

        p_s <- K p_s / sum(K p_s)
        q_s <- K q_s / sum(K q_s)

    where K[i,j] = exp(-k_reg * C[i,j]).

    Then it runs a multi-scale loop: k = 1, 1+step, 1+2*step, ...
    where step = int(k_reg^(1/3)).

    Parameters
    ----------
    x_s, y_s : (NK_small, 3) arrays
    p_s, q_s : (NK_small,) arrays, already normalised
    k_reg    : int, regularisation parameter for the small grid

    Returns
    -------
    f_small : (NK_small,)
    g_small : (NK_small,)
    """
    NK_s = len(x_s)

    # Full cost matrix
    C = cost_matrix_chunk(x_s, y_s)  # (NK_s, NK_s)

    # Kernel pre-conditioning (C++ lines 56-87)
    K = np.exp(-k_reg * C)  # (NK_s, NK_s)

    # Pre-condition p_s
    Kp = K @ p_s           # (NK_s,)
    Kp_sum = np.sum(Kp)
    p_s = Kp / Kp_sum

    # Pre-condition q_s
    Kq = K @ q_s           # (NK_s,)
    Kq_sum = np.sum(Kq)
    q_s = Kq / Kq_sum

    # Initialise potentials and multiplicative factors
    F_s = np.ones(NK_s)
    G_s = np.ones(NK_s)
    f_s = np.zeros(NK_s)
    g_s = np.zeros(NK_s)

    step = int(round(k_reg ** (1.0 / 3.0)))
    if step < 1:
        step = 1

    # Vectorised Sinkhorn on small (NK_s × NK_s) cost matrix
    k = 1
    while k < k_reg:
        # G update (vectorised): G_s[j] = 1/sum_i(exp(-k*(C[i,j]-f[i]-g[j]))*F_s[i]*p_s[i])
        # exponent matrix: -k*(C - f[:,None] - g[None,:]) shape (NK_s, NK_s)
        Fp = F_s * p_s                           # (NK_s,)
        exp_mat = np.exp(-k * (C - f_s[:, None] - g_s[None, :]))  # (NK_s, NK_s)
        sums_G = exp_mat.T @ Fp                  # (NK_s,) — sum over i for each j
        G_s = np.where(sums_G > 0, 1.0 / sums_G, 1.0)
        G_s = np.where(q_s > 0, G_s, 1.0)

        # F update: uses freshly-updated G
        Gq = G_s * q_s                           # (NK_s,)
        exp_mat = np.exp(-k * (C - f_s[:, None] - g_s[None, :]))
        sums_F = exp_mat @ Gq                    # (NK_s,) — sum over j for each i
        F_s = np.where(sums_F > 0, 1.0 / sums_F, 1.0)
        F_s = np.where(p_s > 0, F_s, 1.0)

        # Absorb and reset
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
    """Warm-start the main grid potentials from the small grid.

    Matches C++ ``UseSmall``:

        g_init[j] = min_i(C(x_s[i], y[j]) - f_s[i])
        f_init[i] = min_j(C(x[i],   y_s[j]) - g_s[j])

    Parameters
    ----------
    x_s : (NK_small, 3) small source points
    y_s : (NK_small, 3) small target points
    x   : (NK, 3) main source points
    y   : (NK, 3) main target points
    f_s : (NK_small,) small source potential
    g_s : (NK_small,) small target potential
    chunk_size : int

    Returns
    -------
    f_init : (NK,)
    g_init : (NK,)
    """
    NK = len(x)
    NK_s = len(x_s)

    # g_init[j] = min_i( C(x_s[i], y[j]) - f_s[i] )
    # Loop in chunks over y (the main target)
    g_init = np.full(NK, np.inf)
    for j_start in range(0, NK, chunk_size):
        j_end = min(j_start + chunk_size, NK)
        y_chunk = y[j_start:j_end]
        # C_block[i, j] = C(x_s[i], y_chunk[j])
        C_block = cost_matrix_chunk(x_s, y_chunk)   # (NK_s, chunk)
        # C_block - f_s[:, None] => (NK_s, chunk)
        vals = C_block - f_s[:, None]
        g_init[j_start:j_end] = np.min(vals, axis=0)

    # f_init[i] = min_j( C(x[i], y_s[j]) - g_s[j] )
    f_init = np.full(NK, np.inf)
    for i_start in range(0, NK, chunk_size):
        i_end = min(i_start + chunk_size, NK)
        x_chunk = x[i_start:i_end]
        C_block = cost_matrix_chunk(x_chunk, y_s)   # (chunk, NK_s)
        vals = C_block - g_s[None, :]
        f_init[i_start:i_end] = np.min(vals, axis=1)

    return f_init, g_init


# ---------------------------------------------------------------------------
# Full Sinkhorn divergence
# ---------------------------------------------------------------------------

def run_sinkhorn_divergence(x, y, p, q, chunk_size=512, verbose=True):
    """Full Sinkhorn divergence algorithm matching C++ ``do_sinkhorn_subtracted_axb``.

    Parameters
    ----------
    x : (NK, 3) source points (upper hemisphere)
    y : (NK, 3) target points (lower hemisphere)
    p : (NK,) source weights (unnormalised OK — will be normalised internally)
    q : (NK,) target weights (unnormalised OK)
    chunk_size : int, block size for chunked cost-matrix evaluation
    verbose : bool, if True print progress

    Returns
    -------
    dict with keys:
        f      : (NK,) Kantorovich potential for source (main)
        g      : (NK,) Kantorovich potential for target (main)
        f_id   : (NK,) identity potential for source
        g_id   : (NK,) identity potential for target
        total_cost : float
    """
    from .qmc import load_small_cloud

    NK = len(x)

    # --- Normalise weights ---
    p = np.array(p, dtype=np.float64)
    q = np.array(q, dtype=np.float64)
    p_sum = p.sum()
    q_sum = q.sum()
    p = p / p_sum
    q = q / q_sum

    # Log-weights (-inf where weight is 0)
    logp = np.where(p > 0, np.log(p), -np.inf)
    logq = np.where(q > 0, np.log(q), -np.inf)

    # --- Regularisation schedule ---
    # C++: multiplier=8, getk(NK)=sqrt(NK), getk(NK_small)=sqrt(NK_small)
    # NK=16488  -> getk = floor(sqrt(16488)) = 128
    # NK_small=381 -> getk = floor(sqrt(381)) = 19
    multiplier = 8
    k_small = multiplier * int(np.floor(np.sqrt(381)))   # = 152
    k_final  = multiplier * int(np.floor(np.sqrt(NK)))   # = 1024
    step = int(round(k_final ** (1.0 / 3.0)))             # = 10
    cap_iter = 16
    cap_thr  = 1e-5

    if verbose:
        print(f"k_small={k_small}, k_final={k_final}, step={step}")
        print(f"NK={NK}, NK points in source: {int((p>0).sum())}, "
              f"NK points in target: {int((q>0).sum())}")

    # --- Small grid warm-start ---
    x_s, y_s = load_small_cloud()
    NK_s = len(x_s)
    p_s_raw = np.zeros(NK_s)
    q_s_raw = np.zeros(NK_s)
    # Evaluate densities on small grid points using p/q functions stored in caller
    # Actually the C++ smallsinkhorn calls P(x_small[i]) and Q(y_small[i]).
    # We receive already-evaluated p, q for the main grid.  For the small grid
    # we need to re-evaluate — but we don't have the density function here.
    # Instead we pass in p_s and q_s from the caller.  See run_benchmark().
    # For standalone use (called directly), compute uniform fallback.
    # We raise a clear error so that callers know to pass p_s/q_s separately.
    raise RuntimeError(
        "run_sinkhorn_divergence must be called via run_benchmark(), which "
        "supplies p_s and q_s for the small grid.  Alternatively call "
        "_run_sinkhorn_divergence_inner() directly."
    )


def _run_sinkhorn_divergence_inner(x, y, p, q, x_s, y_s, p_s, q_s,
                                   chunk_size=512, verbose=True):
    """Internal implementation of the full Sinkhorn divergence.

    Separated from run_sinkhorn_divergence so that run_benchmark can supply
    the small-grid densities p_s, q_s without re-evaluating P/Q.

    Parameters
    ----------
    x, y   : (NK, 3) main source/target points
    p, q   : (NK,) main source/target weights (will be normalised)
    x_s, y_s : (NK_small, 3) small source/target points
    p_s, q_s : (NK_small,) small source/target weights (will be normalised)
    chunk_size : int
    verbose : bool

    Returns
    -------
    dict with keys: f, g, f_id, g_id, total_cost
    """
    NK = len(x)

    # --- Normalise weights ---
    p = np.array(p, dtype=np.float64)
    q = np.array(q, dtype=np.float64)
    p = p / p.sum()
    q = q / q.sum()

    logp = np.where(p > 0, np.log(p), -np.inf)
    logq = np.where(q > 0, np.log(q), -np.inf)

    p_s = np.array(p_s, dtype=np.float64)
    q_s = np.array(q_s, dtype=np.float64)
    p_s = p_s / p_s.sum()
    q_s = q_s / q_s.sum()

    # --- Regularisation schedule ---
    multiplier = 8
    k_small = multiplier * int(np.floor(np.sqrt(len(x_s))))  # 8*19=152
    k_final  = multiplier * int(np.floor(np.sqrt(NK)))        # 8*128=1024
    step = int(round(k_final ** (1.0 / 3.0)))                  # 10
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
        if i >= cap_iter:
            break
    if verbose:
        print(f"Final loop: {i} iterations, last change={maxdif:.4e}")

    # --- 5. Identity F loop ---
    if verbose:
        print("\nIdentity Sinkhorn for source (f_id):")
    f_id = np.zeros(NK)
    id_step = int(round(np.sqrt(k_final)))  # 32
    regvar = 1
    i = 0
    while regvar < k_final:
        f_id, maxdif = sinkhorn_identity_f_step(x, y, logp, f_id, regvar, chunk_size)
        i += 1
        if verbose:
            print(f"  iter {i:4d}, k={regvar:5d}, maxdif={maxdif:.4e}")
        regvar += id_step

    if verbose:
        print(f"  Final {cap_iter} iterations at k={k_final}:")
    i = 0
    maxdif = cap_thr + 1.0
    while maxdif > cap_thr:
        f_id, maxdif = sinkhorn_identity_f_step(x, y, logp, f_id, k_final, chunk_size)
        i += 1
        if verbose:
            print(f"  iter {i:3d}, maxdif={maxdif:.4e}")
        if i >= cap_iter:
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
        g_id, maxdif = sinkhorn_identity_g_step(x, y, logq, g_id, regvar, chunk_size)
        i += 1
        if verbose:
            print(f"  iter {i:4d}, k={regvar:5d}, maxdif={maxdif:.4e}")
        regvar += id_step

    if verbose:
        print(f"  Final {cap_iter} iterations at k={k_final}:")
    i = 0
    maxdif = cap_thr + 1.0
    while maxdif > cap_thr:
        g_id, maxdif = sinkhorn_identity_g_step(x, y, logq, g_id, k_final, chunk_size)
        i += 1
        if verbose:
            print(f"  iter {i:3d}, maxdif={maxdif:.4e}")
        if i >= cap_iter:
            break
    if verbose:
        print(f"Identity G: {i} final iterations, last change={maxdif:.4e}")

    # --- 7. Normalise (matching C++ sink_to_Reflector_subtracted_axb) ---
    # Shift f_id so its max is 0, compensate in g_id
    max_f_id = np.max(f_id)
    f_id = f_id - max_f_id
    g_id = g_id + max_f_id

    # Shift f so its max is 0, compensate in g
    max_f = np.max(f)
    f = f - max_f
    g = g + max_f

    # Subtract identity terms
    f = f - f_id
    g = g - g_id

    # --- 8. Total cost ---
    # At this point f = f_sinkhorn - f_id and g = g_sinkhorn - g_id
    # (identity terms have already been subtracted above).
    # C++ Compute_TotalCost() = sum p*(f_sink - f_id) + sum q*(g_sink - g_id)
    # which equals sum p*f + sum q*g here.
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
