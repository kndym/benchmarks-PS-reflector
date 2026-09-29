"""Cost functions used by the reflector and refraction OT routines."""

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


def validate_point_clouds(x: np.ndarray, y: np.ndarray,
                          chunk_size: int = 512) -> None:
    """Validate two non-empty, finite point clouds with matching dimensions."""
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


def validate_transport_possible(x: np.ndarray, y: np.ndarray,
                                chunk_size: int = 512) -> None:
    """Raise if any pair violates the finite-cost condition kappa*(x dot y) < 1."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    validate_point_clouds(x, y, chunk_size)

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
    return cost_matrix_refraction_chunk(x_chunk, y)


def cost_matrix_refraction_chunk(x_chunk: np.ndarray, y: np.ndarray,
                                 kappa: float | None = None) -> np.ndarray:
    """Return refraction costs ``-log(1-kappa*(x dot y))`` for a point block."""
    x_chunk = np.asarray(x_chunk, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    kappa = _KAPPA if kappa is None else float(kappa)
    if not (0.0 < kappa <= 1.0):
        raise ValueError(f"kappa must be in (0, 1], got {kappa}")
    # Dot products: (M, N)
    dots = x_chunk @ y.T
    arg = 1.0 - kappa * dots
    if np.any(arg <= 0.0):
        raise ValueError(
            "transport cost is undefined: require kappa*(x dot y) < 1 "
            f"(kappa={kappa:g}, max dot={np.max(dots):.17g})"
        )
    return -np.log(arg)


def cost_matrix_l2_chunk(x_chunk: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return squared Euclidean (L2) costs ``||x-y||_2^2`` for a point block."""
    x_chunk = np.asarray(x_chunk, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x_chunk.ndim != 2 or y.ndim != 2 or x_chunk.shape[1] != y.shape[1]:
        raise ValueError("x_chunk and y must have shape (N, dim) with matching dimensions")
    x_norm2 = np.einsum("ij,ij->i", x_chunk, x_chunk)
    y_norm2 = np.einsum("ij,ij->i", y, y)
    costs = x_norm2[:, None] + y_norm2[None, :] - 2.0 * (x_chunk @ y.T)
    # Roundoff can make an exact zero distance slightly negative.
    np.maximum(costs, 0.0, out=costs)
    return costs


def get_cost_function(name: str = "l2", *, kappa: float | None = None):
    """Return a chunked OT cost function by name.

    ``l2`` (the default) is squared Euclidean distance and does not use kappa.
    ``refraction`` selects ``-log(1-kappa*(x dot y))``; an explicit kappa can
    be supplied without changing the module-level refraction setting.
    """
    if not isinstance(name, str):
        raise TypeError("cost function name must be a string")
    key = name.strip().lower().replace("-", "_")
    if key in {"l2", "squared_l2", "squared_euclidean"}:
        return cost_matrix_l2_chunk
    if key in {"refraction", "refraction_cost"}:
        if kappa is None:
            return cost_matrix_chunk
        if not (0.0 < kappa <= 1.0):
            raise ValueError(f"kappa must be in (0, 1], got {kappa}")
        return lambda x_chunk, y: cost_matrix_refraction_chunk(
            x_chunk, y, kappa=kappa
        )
    raise ValueError("unknown cost function; choose 'l2' or 'refraction'")
