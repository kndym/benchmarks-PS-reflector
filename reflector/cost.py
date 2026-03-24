"""
reflector/cost.py

Cost function for the Point Source Far-Field Reflector Problem.

The cost is the Wang (2004) reflector cost:

    c(x, y) = -log(1 - x·y)

where x and y are unit vectors on the sphere.  When x·y approaches 1 (i.e.
when x = y), the cost diverges to +∞.  To avoid numerical issues we clip the
dot product away from 1 before taking the logarithm.

Public API
----------
cost_vec(x_vec, y_vec)                -> scalar
cost_matrix_chunk(x_chunk, y)         -> ndarray (len(x_chunk), len(y))
"""

import numpy as np

# Clipping value: keep (1 - dot) >= _EPS so that -log(1-dot) is finite.
# 1 - _EPS_CLIP is the maximum dot product we accept.
_EPS_CLIP = 1e-15


def cost_vec(x_vec: np.ndarray, y_vec: np.ndarray) -> float:
    """Compute the reflector cost for a single pair of unit vectors.

    Parameters
    ----------
    x_vec : array_like, shape (3,)
    y_vec : array_like, shape (3,)

    Returns
    -------
    float
        c(x, y) = -log(1 - x·y), or +inf if the vectors are (nearly) equal.
    """
    x_vec = np.asarray(x_vec, dtype=np.float64)
    y_vec = np.asarray(y_vec, dtype=np.float64)
    dot = np.dot(x_vec, y_vec)
    arg = 1.0 - dot
    arg = max(arg, _EPS_CLIP)
    return float(-np.log(arg))


def cost_matrix_chunk(x_chunk: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Compute a block of the cost matrix.

    Computes c(x_chunk[i], y[j]) = -log(1 - x_chunk[i]·y[j]) for all i, j.

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
    arg = np.clip(1.0 - dots, _EPS_CLIP, None)
    return -np.log(arg)
