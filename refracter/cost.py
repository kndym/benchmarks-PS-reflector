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

def set_kappa(k: float) -> None:
    """Set the module-level κ; must be in (0, 1]."""
    global _KAPPA
    if not (0.0 < k <= 1.0):
        raise ValueError(f"kappa must be in (0, 1], got {k}")
    _KAPPA = float(k)


def get_kappa() -> float:
    """Return the current refraction parameter κ."""
    return _KAPPA


def validate_transport_possible(x: np.ndarray, y: np.ndarray,
                                chunk_size: int = 512) -> None:
    """Raise if any pair violates the finite-cost condition kappa*(x dot y) < 1."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.ndim != 2 or y.ndim != 2 or x.shape[1] != y.shape[1]:
        raise ValueError("transport point clouds must have shape (N, dim) with matching dimensions")
    if len(x) == 0 or len(y) == 0:
        raise ValueError("transport point clouds must be non-empty")
    if chunk_size <= 0:
        raise ValueError(f"chunk_size must be positive, got {chunk_size}")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("transport point clouds must contain only finite values")

    max_dot = -np.inf
    for i_start in range(0, len(x), chunk_size):
        max_dot = max(max_dot, np.max(x[i_start:i_start + chunk_size] @ y.T))
        if _KAPPA * max_dot >= 1.0:
            raise ValueError(
                "transport cost is undefined: require kappa*(x dot y) < 1 "
                f"for every pair (kappa={_KAPPA:g}, max dot={max_dot:.17g})"
            )


def cost_vec(x_vec: np.ndarray, y_vec: np.ndarray) -> float:
    """Scalar cost c(x,y) = -log(1 - κ·(x·y)) for a single pair of unit vectors."""
    x_vec = np.asarray(x_vec, dtype=np.float64)
    y_vec = np.asarray(y_vec, dtype=np.float64)
    dot = np.dot(x_vec, y_vec)
    arg = 1.0 - _KAPPA * dot
    if arg <= 0.0:
        raise ValueError(
            "transport cost is undefined: require kappa*(x dot y) < 1 "
            f"(kappa={_KAPPA:g}, dot={dot:.17g})"
        )
    return float(-np.log(arg))


def cost_matrix_chunk(x_chunk: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return cost matrix block of shape (M, N) for x_chunk (M,3) and y (N,3)."""
    x_chunk = np.asarray(x_chunk, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    # Dot products: (M, N)
    dots = x_chunk @ y.T
    arg = 1.0 - _KAPPA * dots
    if np.any(arg <= 0.0):
        raise ValueError(
            "transport cost is undefined: require kappa*(x dot y) < 1 "
            f"(kappa={_KAPPA:g}, max dot={np.max(dots):.17g})"
        )
    return -np.log(arg)
