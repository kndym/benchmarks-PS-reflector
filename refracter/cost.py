# Cost function for unit vectors: c(x,y) = -log(1 - κ * x·y).
# Use κ = 0.6 for refraction and κ = 1.0 for reflection.

import numpy as np

# The cost parameter is shared by all cost-matrix calls.
_KAPPA = 1.0          # default: standard reflector cost

# Keep the logarithm finite when roundoff makes the cost argument too small.
_EPS_CLIP = 1e-15


def set_kappa(k: float) -> None:
    # Set the cost parameter before running a benchmark.
    global _KAPPA
    if not (0.0 < k <= 1.0):
        raise ValueError(f"kappa must be in (0, 1], got {k}")
    _KAPPA = float(k)


def get_kappa() -> float:
    # Return the current cost parameter.
    return _KAPPA


def cost_vec(x_vec: np.ndarray, y_vec: np.ndarray) -> float:
    # Evaluate the cost for one pair of directions.
    x_vec = np.asarray(x_vec, dtype=np.float64)
    y_vec = np.asarray(y_vec, dtype=np.float64)
    dot = np.dot(x_vec, y_vec)
    arg = 1.0 - _KAPPA * dot
    arg = max(arg, _EPS_CLIP)
    return float(-np.log(arg))


def cost_matrix_chunk(x_chunk: np.ndarray, y: np.ndarray) -> np.ndarray:
    # Evaluate a cost-matrix block without materializing the full matrix.
    x_chunk = np.asarray(x_chunk, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    # Dot products: (M, N)
    dots = x_chunk @ y.T
    # Clip to ensure argument of log is positive
    arg = np.clip(1.0 - _KAPPA * dots, _EPS_CLIP, None)
    return -np.log(arg)
