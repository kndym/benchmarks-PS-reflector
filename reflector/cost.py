"""
reflector/cost.py

Cost function for the Point Source Far-Field Reflector / Refractor Problem.

The generalised cost is:

    c(x, y) = -log(1 - κ · (x·y))

where x and y are unit vectors on the sphere and κ ∈ (0, 1] is the
refraction parameter.

* κ = 1.0  →  standard Wang (2004) reflector cost:  c = -log(1 - x·y)
* κ = 0.6  →  refraction cost from Figure 4 of the paper

When x·y approaches 1/κ, the cost diverges to +∞.  To avoid numerical
issues we clip the argument of the logarithm away from 0.

Public API
----------
set_kappa(k)                         -> None  (set the module-level κ)
get_kappa()                          -> float (query the current κ)
cost_vec(x_vec, y_vec)               -> scalar
cost_matrix_chunk(x_chunk, y)        -> ndarray (len(x_chunk), len(y))
"""

import numpy as np

# ---------------------------------------------------------------------------
# Module-level refraction parameter
# ---------------------------------------------------------------------------
_KAPPA = 1.0          # default: standard reflector cost

# Clipping value: keep (1 - κ·dot) >= _EPS so that -log(...) is finite.
_EPS_CLIP = 1e-15


def set_kappa(k: float) -> None:
    """Set the refraction parameter κ used by the cost functions.

    Parameters
    ----------
    k : float
        Must be in (0, 1].  Use 1.0 for the standard reflector cost and
        0.6 for the refraction benchmark.
    """
    global _KAPPA
    if not (0.0 < k <= 1.0):
        raise ValueError(f"kappa must be in (0, 1], got {k}")
    _KAPPA = float(k)


def get_kappa() -> float:
    """Return the current refraction parameter κ."""
    return _KAPPA


def cost_vec(x_vec: np.ndarray, y_vec: np.ndarray) -> float:
    """Compute the cost for a single pair of unit vectors.

    Parameters
    ----------
    x_vec : array_like, shape (3,)
    y_vec : array_like, shape (3,)

    Returns
    -------
    float
        c(x, y) = -log(1 - κ·(x·y)), or +inf if the argument is ≤ 0.
    """
    x_vec = np.asarray(x_vec, dtype=np.float64)
    y_vec = np.asarray(y_vec, dtype=np.float64)
    dot = np.dot(x_vec, y_vec)
    arg = 1.0 - _KAPPA * dot
    arg = max(arg, _EPS_CLIP)
    return float(-np.log(arg))


def cost_matrix_chunk(x_chunk: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Compute a block of the cost matrix.

    Computes c(x_chunk[i], y[j]) = -log(1 - κ·(x_chunk[i]·y[j])) for all i, j.

    Parameters
    ----------
    x_chunk : ndarray, shape (M, 3)
        A chunk of source points.
    y : ndarray, shape (N, 3)
        All target points.

    Returns
    -------
    C : ndarray, shape (M, N)
        Block of the cost matrix.
    """
    x_chunk = np.asarray(x_chunk, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    # Dot products: (M, N)
    dots = x_chunk @ y.T
    # Clip to ensure argument of log is positive
    arg = np.clip(1.0 - _KAPPA * dots, _EPS_CLIP, None)
    return -np.log(arg)
