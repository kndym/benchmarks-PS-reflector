"""
reflector/build.py

Build the reflector surface from Sinkhorn potentials and compute the
c-transforms needed for the Sinkhorn divergence correction.

Functions
---------
build_reflector(x, f, f_id)
    Compute the reflector scale R and reflector points Ref from the potentials.

c_transform_gc(x, y, g, chunk_size)
    Compute gc[i] = min_j(C(x[i], y[j]) - g[j])   (c-transform of g)

c_transform_fc(x, y, f, chunk_size)
    Compute fc[j] = min_i(C(x[i], y[j]) - f[i])   (c-conjugate of f)

build_regular_grid(final_grid_res)
    Build the regular stereographic grid on the upper hemisphere.

reflector_on_regular_grid(x_regular, y, g, chunk_size)
    Evaluate the reflector on the regular grid via the c-transform of g.
"""

import numpy as np
from .cost import cost_matrix_chunk


# ---------------------------------------------------------------------------
# Reflector from potentials
# ---------------------------------------------------------------------------

def build_reflector(x: np.ndarray, f: np.ndarray, f_id: np.ndarray):
    """Compute the reflector scale R and reflector points Ref.

    Matches C++ ``sink_to_Reflector_subtracted_axb`` after the normalisation
    and identity subtraction:

        R[i]   = exp(f[i] - f_id[i])      (but f has already had f_id subtracted
                                             in _run_sinkhorn_divergence_inner)
        Ref[i] = 2 * x[i] * R[i]

    In practice ``_run_sinkhorn_divergence_inner`` already computes
    ``f <- f - f_id`` before returning, so we simply compute:

        R[i]   = exp(f[i])
        Ref[i] = 2 * x[i] * R[i]

    where the passed ``f`` is the *corrected* (Sinkhorn divergence) potential.

    Parameters
    ----------
    x    : (NK, 3) source points
    f    : (NK,) corrected f potential  (= f_sinkhorn - f_id after normalisation)
    f_id : (NK,) identity potential for source (used only for documentation here;
           in the C++ the subtraction happens before this function is called)

    Returns
    -------
    R   : (NK,) reflector radial scale
    Ref : (NK, 3) reflector point directions (unnormalised)
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
    """C-transform of g from the target side to the source side.

    For each source point i:

        gc[i] = min_j ( C(x[i], y[j]) - g[j] )

    This is the c-conjugate of g, mapping target potential to source side.
    Matches C++ ``Get_gc``.

    Parameters
    ----------
    x  : (NK, 3) source points
    y  : (NK, 3) target points
    g  : (NK,) target Kantorovich potential
    chunk_size : int

    Returns
    -------
    gc : (NK,) c-transform evaluated at source points
    """
    NK = len(x)
    gc = np.full(NK, np.inf, dtype=np.float64)

    # Finite-g mask (C++ skips j where g[j] == inf)
    finite_mask = np.isfinite(g)
    y_finite = y[finite_mask]
    g_finite = g[finite_mask]

    if len(y_finite) == 0:
        return gc

    for i_start in range(0, NK, chunk_size):
        i_end = min(i_start + chunk_size, NK)
        x_chunk = x[i_start:i_end]

        C_block = cost_matrix_chunk(x_chunk, y_finite)   # (chunk, N_finite)
        # C_block - g_finite[None, :]  =>  (chunk, N_finite)
        vals = C_block - g_finite[None, :]
        gc[i_start:i_end] = np.min(vals, axis=1)

    return gc


def c_transform_fc(x: np.ndarray, y: np.ndarray, f: np.ndarray,
                   chunk_size: int = 512) -> np.ndarray:
    """C-transform of f from the source side to the target side.

    For each target point j:

        fc[j] = min_i ( C(x[i], y[j]) - f[i] )

    Matches C++ ``Get_fc``.

    Parameters
    ----------
    x  : (NK, 3) source points
    y  : (NK, 3) target points
    f  : (NK,) source Kantorovich potential
    chunk_size : int

    Returns
    -------
    fc : (NK,) c-transform evaluated at target points
    """
    NK = len(y)
    fc = np.full(NK, np.inf, dtype=np.float64)

    for j_start in range(0, NK, chunk_size):
        j_end = min(j_start + chunk_size, NK)
        y_chunk = y[j_start:j_end]

        C_block = cost_matrix_chunk(x, y_chunk)   # (NK_x, chunk)
        # C_block - f[:, None]  =>  (NK_x, chunk)
        vals = C_block - f[:, None]
        fc[j_start:j_end] = np.min(vals, axis=0)

    return fc


# ---------------------------------------------------------------------------
# Regular grid
# ---------------------------------------------------------------------------

def build_regular_grid(final_grid_res: int = 1025):
    """Build the regular stereographic grid on the upper hemisphere.

    The grid covers the square [-0.6, 0.6]² in the north-pole stereographic
    plane and is mapped back to 3-D unit sphere coordinates.

    Matches C++ ``discretization()`` (the x_regular and Regular_side parts).

    Grid indexing: point (i, j) maps to flat index i*final_grid_res + j.

    Parameters
    ----------
    final_grid_res : int
        Number of grid lines per axis.  Default 1025 = 8*128+1.

    Returns
    -------
    x_regular : (final_grid_res², 3) array of 3-D unit sphere points
    Regular_side : (final_grid_res,) array of grid line positions in [-0.6, 0.6]
    """
    res = final_grid_res

    # 1-D axis in [-0.6, 0.6]
    side = np.linspace(-0.6, 0.6, res)          # (res,)
    Regular_side = side                          # Regular_side[j] = Y value at column j

    # Meshgrid: X varies over rows (i), Y over columns (j)
    X, Y = np.meshgrid(side, side, indexing='ij')   # both (res, res)
    N2    = X * X + Y * Y                            # (res, res)
    denom = 1.0 + N2                                 # (res, res)

    x0 = (2.0 * X / denom).ravel()
    x1 = (2.0 * Y / denom).ravel()
    x2 = ((1.0 - N2) / denom).ravel()

    x_regular = np.stack([x0, x1, x2], axis=1)      # (FinalGrid, 3)

    return x_regular, Regular_side


def reflector_on_regular_grid(x_regular: np.ndarray, y: np.ndarray,
                               g: np.ndarray, chunk_size: int = 512):
    """Evaluate the reflector surface on the regular grid.

    For each regular grid point i compute the c-transform of g:

        f_regular[i] = min_j ( C(x_regular[i], y[j]) - g[j] )

    and then build the reflector point:

        Ref_regular[i] = 2 * x_regular[i] * exp(f_regular[i])

    Matches C++ ``sink_to_Reflector_subtracted_axb`` (the f_regular and
    Ref_regular computation block).

    Parameters
    ----------
    x_regular : (FinalGrid, 3) regular grid points on upper hemisphere
    y         : (NK, 3) target points
    g         : (NK,) target Kantorovich potential (corrected)
    chunk_size : int

    Returns
    -------
    f_regular   : (FinalGrid,) c-transform on regular grid
    Ref_regular : (FinalGrid, 3) reflector surface on regular grid
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
