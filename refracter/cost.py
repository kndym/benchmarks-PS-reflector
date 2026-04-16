"""
refracter/cost.py

Cost function c(x,y) = -log(1 - κ·(x·y)) for unit vectors x, y on the sphere.
κ=1.0 → reflector,  κ=0.6 → refractor.

Call set_kappa() before importing anything that caches cost computations.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Module-level refraction parameter
# ---------------------------------------------------------------------------
_KAPPA = 1.0          # default: standard reflector cost

# Clipping value: keep (1 - κ·dot) >= _EPS so that -log(...) is finite.
_EPS_CLIP = 1e-15


def set_kappa(k: float) -> None:
    """Set the module-level κ; must be in (0, 1]."""
    global _KAPPA
    if not (0.0 < k <= 1.0):
        raise ValueError(f"kappa must be in (0, 1], got {k}")
    _KAPPA = float(k)


def get_kappa() -> float:
    """Return the current refraction parameter κ."""
    return _KAPPA


def cost_vec(x_vec: np.ndarray, y_vec: np.ndarray) -> float:
    """Scalar cost c(x,y) = -log(1 - κ·(x·y)) for a single pair of unit vectors."""
    x_vec = np.asarray(x_vec, dtype=np.float64)
    y_vec = np.asarray(y_vec, dtype=np.float64)
    dot = np.dot(x_vec, y_vec)
    arg = 1.0 - _KAPPA * dot
    arg = max(arg, _EPS_CLIP)
    return float(-np.log(arg))


def cost_matrix_chunk(x_chunk: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return cost matrix block of shape (M, N) for x_chunk (M,3) and y (N,3)."""
    x_chunk = np.asarray(x_chunk, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    # Dot products: (M, N)
    dots = x_chunk @ y.T
    # Clip to ensure argument of log is positive
    arg = np.clip(1.0 - _KAPPA * dots, _EPS_CLIP, None)
    return -np.log(arg)
