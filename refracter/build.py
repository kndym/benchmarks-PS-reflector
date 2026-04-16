"""
refracter/build.py

Build the reflector surface from Sinkhorn potentials and compute c-transforms.
"""

import numpy as np
from .cost import cost_matrix_chunk


# ---------------------------------------------------------------------------
# Reflector from potentials
# ---------------------------------------------------------------------------

def build_reflector(x: np.ndarray, f: np.ndarray, f_id: np.ndarray):
    """Return (R, Ref) where R = exp(f) and Ref = 2*x*R.

    Expects f already corrected (f_id subtracted) by _run_sinkhorn_divergence_inner.
    """
    x = np.asarray(x, dtype=np.float64)
    f = np.asarray(f, dtype=np.float64)

    R = np.exp(f)                         # shape (NK,)
    Ref = 2.0 * x * R[:, None]            # shape (NK, 3)

    return R, Ref


# ---------------------------------------------------------------------------
# C-transforms
# ---------------------------------------------------------------------------

def c_transform_gc(x: np.ndarray, y: np.ndarray, g: np.ndarray,
                   chunk_size: int = 512) -> np.ndarray:
    """gc[i] = min_j(C(x[i], y[j]) - g[j])  (c-conjugate of g, C++ Get_gc)."""
    NK = len(x)
    gc = np.full(NK, np.inf, dtype=np.float64)

    finite_mask = np.isfinite(g)
    y_finite = y[finite_mask]
    g_finite = g[finite_mask]

    if len(y_finite) == 0:
        return gc

    for i_start in range(0, NK, chunk_size):
        i_end = min(i_start + chunk_size, NK)
        C_block = cost_matrix_chunk(x[i_start:i_end], y_finite)
        gc[i_start:i_end] = np.min(C_block - g_finite[None, :], axis=1)

    return gc


def c_transform_fc(x: np.ndarray, y: np.ndarray, f: np.ndarray,
                   chunk_size: int = 512) -> np.ndarray:
    """fc[j] = min_i(C(x[i], y[j]) - f[i])  (c-conjugate of f, C++ Get_fc)."""
    NK = len(y)
    fc = np.full(NK, np.inf, dtype=np.float64)

    for j_start in range(0, NK, chunk_size):
        j_end = min(j_start + chunk_size, NK)
        C_block = cost_matrix_chunk(x, y[j_start:j_end])
        fc[j_start:j_end] = np.min(C_block - f[:, None], axis=0)

    return fc


# ---------------------------------------------------------------------------
# Regular grid
# ---------------------------------------------------------------------------

def build_regular_grid(final_grid_res: int = 1025):
    """Regular stereographic grid on [-0.6,0.6]² mapped to upper hemisphere.

    Returns (x_regular, Regular_side) with shapes (res², 3) and (res,).
    """
    res = final_grid_res
    side = np.linspace(-0.6, 0.6, res)
    Regular_side = side

    X, Y = np.meshgrid(side, side, indexing='ij')
    N2    = X * X + Y * Y
    denom = 1.0 + N2

    x_regular = np.stack([
        (2.0 * X / denom).ravel(),
        (2.0 * Y / denom).ravel(),
        ((1.0 - N2) / denom).ravel(),
    ], axis=1)

    return x_regular, Regular_side


def reflector_on_regular_grid(x_regular: np.ndarray, y: np.ndarray,
                               g: np.ndarray, chunk_size: int = 512):
    """Compute c-transform of g on the regular grid and build Ref_regular.

    Returns (f_regular, Ref_regular).
    """
    FinalGrid = len(x_regular)

    # Finite-g mask
    finite_mask = np.isfinite(g)
    y_finite = y[finite_mask]
    g_finite = g[finite_mask]

    f_regular = np.full(FinalGrid, np.inf, dtype=np.float64)

    for i_start in range(0, FinalGrid, chunk_size):
        i_end = min(i_start + chunk_size, FinalGrid)
        xr_chunk = x_regular[i_start:i_end]

        C_block = cost_matrix_chunk(xr_chunk, y_finite)   # (chunk, N_finite)
        vals = C_block - g_finite[None, :]
        f_regular[i_start:i_end] = np.min(vals, axis=1)

    Ref_regular = 2.0 * x_regular * np.exp(f_regular)[:, None]   # (FinalGrid, 3)

    return f_regular, Ref_regular
