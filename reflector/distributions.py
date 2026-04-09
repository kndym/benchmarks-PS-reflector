"""
reflector/distributions.py

Source and target density functions for the reflector benchmark problems,
together with stereographic projection helpers.

Stereographic projections
-------------------------
North pole (upper hemisphere → plane):
    (x0, x1, x2)  ->  (x0/(1+x2),  x1/(1+x2))

South pole (lower hemisphere → plane):
    (y0, y1, y2)  ->  (y0/(1-y2),  y1/(1-y2))

Inverse north pole (plane → upper hemisphere):
    (u, v)  ->  (2u/(1+N²),  2v/(1+N²),  (1-N²)/(1+N²))   where N²=u²+v²

Benchmark test cases (matching C++ headers)
-------------------------------------------
SquareToCircle
    P(x) = 1  if  stereo_north(x) in [-0.5, 0.5]²
    Q(y) = 1  if  |stereo_south(y)| <= 0.5

SquareToTwoGaussSide
    P(x) = same as above
    Q(y) = exp(-16*((u-0.25)²+(v-0.25)²)) + exp(-16*((u-0.25)²+(v+0.25)²))
           where (u,v) = stereo_south(y)
"""

import numpy as np

# ---------------------------------------------------------------------------
# Stereographic projections
# ---------------------------------------------------------------------------

def stereo_north(pts: np.ndarray) -> tuple:
    """Stereographic projection from the upper hemisphere to the plane.

    Parameters
    ----------
    pts : ndarray, shape (..., 3)
        Points on the unit sphere with x2 > -1.

    Returns
    -------
    u, v : ndarrays matching the leading dimensions of pts.
    """
    pts = np.asarray(pts, dtype=np.float64)
    denom = 1.0 + pts[..., 2]
    u = pts[..., 0] / denom
    v = pts[..., 1] / denom
    return u, v


def stereo_south(pts: np.ndarray) -> tuple:
    """Stereographic projection from the lower hemisphere to the plane.

    Parameters
    ----------
    pts : ndarray, shape (..., 3)
        Points on the unit sphere with x2 < 1.

    Returns
    -------
    u, v : ndarrays matching the leading dimensions of pts.
    """
    pts = np.asarray(pts, dtype=np.float64)
    denom = 1.0 - pts[..., 2]
    u = pts[..., 0] / denom
    v = pts[..., 1] / denom
    return u, v


def stereo_north_inverse(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Inverse of the north-pole stereographic projection.

    Maps a point (u, v) in the plane back to the upper hemisphere.

    Parameters
    ----------
    u, v : array_like
        Plane coordinates.

    Returns
    -------
    pts : ndarray, shape (..., 3)
    """
    u = np.asarray(u, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    N2 = u * u + v * v
    denom = 1.0 + N2
    x0 = 2.0 * u / denom
    x1 = 2.0 * v / denom
    x2 = (1.0 - N2) / denom
    return np.stack([x0, x1, x2], axis=-1)


# ---------------------------------------------------------------------------
# Density functions
# ---------------------------------------------------------------------------

def P_square(x: np.ndarray) -> np.ndarray:
    """Uniform density on the square [-0.5, 0.5]² in stereographic coordinates.

    Equivalent to P(x) in the C++ test headers: returns 1 if the north-pole
    stereographic projection of x falls inside the open square, 0 otherwise.

    Parameters
    ----------
    x : ndarray, shape (N, 3)

    Returns
    -------
    density : ndarray, shape (N,), values in {0, 1}.
    """
    u, v = stereo_north(x)
    inside = (u > -0.5) & (u < 0.5) & (v > -0.5) & (v < 0.5)
    return inside.astype(np.float64)


def Q_circle(y: np.ndarray) -> np.ndarray:
    """Uniform density on the disc of radius 0.5 in south-pole coordinates.

    Equivalent to Q(y) in test_3D_SquareToCircle: returns 1 if
    |stereo_south(y)|² <= 0.25, else 0.

    Parameters
    ----------
    y : ndarray, shape (N, 3)

    Returns
    -------
    density : ndarray, shape (N,), values in {0, 1}.
    """
    u, v = stereo_south(y)
    inside = (u * u + v * v) <= 0.25
    return inside.astype(np.float64)


def Q_two_gaussians(y: np.ndarray) -> np.ndarray:
    """Two Gaussians side by side in south-pole stereographic coordinates.

    Equivalent to Q(y) in test_3D_SquareToTwoGaussSide:

        Q = exp(-16*((u-0.25)²+(v-0.25)²)) + exp(-16*((u-0.25)²+(v+0.25)²))

    where (u, v) = stereo_south(y).

    Parameters
    ----------
    y : ndarray, shape (N, 3)

    Returns
    -------
    density : ndarray, shape (N,), non-negative floats.
    """
    u, v = stereo_south(y)
    g1 = np.exp(-16.0 * ((u - 0.25) ** 2 + (v - 0.25) ** 2))
    g2 = np.exp(-16.0 * ((u - 0.25) ** 2 + (v + 0.25) ** 2))
    return g1 + g2


# ---------------------------------------------------------------------------
# Refraction patch densities (flat indicator on spherical patches)
#
# Both source and target live on the *upper* hemisphere.
# Bounds are given in spherical coordinates (θ = polar, φ = azimuthal):
#
#   Source  Ω  : θ ∈ [π/12, π/3],  φ ∈ [π/12, π/4]
#   Target  Ω* : θ ∈ [π/10, π/5],  φ ∈ [π/10, π/5]
#
# Matches the C++ header test_3D_FlatPatchRefraction_logcost_MonteCarlo.h.
# ---------------------------------------------------------------------------

# Default patch bounds (multiples of π)
_SRC_THETA = (np.pi / 12, np.pi / 3)
_SRC_PHI   = (np.pi / 12, np.pi / 4)
_TGT_THETA = (np.pi / 10, np.pi / 5)         # polar angle range
_TGT_PHI   = (np.pi / 3, 5 * np.pi / 12)   # azimuthal range — non-overlapping with source


def P_refraction_patch(x: np.ndarray) -> np.ndarray:
    """Flat indicator on the source refraction patch.

    Returns 1 if (θ, φ) of x falls inside
        θ ∈ [π/12, π/3],  φ ∈ [π/12, π/4]
    and 0 otherwise.  φ is taken in [0, 2π).

    Parameters
    ----------
    x : ndarray, shape (N, 3)

    Returns
    -------
    density : ndarray, shape (N,), values in {0, 1}.
    """
    x = np.asarray(x, dtype=np.float64)
    theta = np.arccos(np.clip(x[..., 2], -1.0, 1.0))
    phi   = np.arctan2(x[..., 1], x[..., 0])
    phi   = np.where(phi < 0, phi + 2 * np.pi, phi)

    inside = (
        (theta >= _SRC_THETA[0]) & (theta <= _SRC_THETA[1]) &
        (phi   >= _SRC_PHI[0])   & (phi   <= _SRC_PHI[1])
    )
    return inside.astype(np.float64)


def Q_refraction_patch(y: np.ndarray) -> np.ndarray:
    """Flat indicator on the target refraction patch.

    Returns 1 if (θ, φ) of y falls inside
        θ ∈ [π/10, π/5],  φ ∈ [π/10, π/5]
    and 0 otherwise.  φ is taken in [0, 2π).

    Parameters
    ----------
    y : ndarray, shape (N, 3)

    Returns
    -------
    density : ndarray, shape (N,), values in {0, 1}.
    """
    y = np.asarray(y, dtype=np.float64)
    theta = np.arccos(np.clip(y[..., 2], -1.0, 1.0))
    phi   = np.arctan2(y[..., 1], y[..., 0])
    phi   = np.where(phi < 0, phi + 2 * np.pi, phi)

    inside = (
        (theta >= _TGT_THETA[0]) & (theta <= _TGT_THETA[1]) &
        (phi   >= _TGT_PHI[0])   & (phi   <= _TGT_PHI[1])
    )
    return inside.astype(np.float64)


# ---------------------------------------------------------------------------
# Benchmark registry
# ---------------------------------------------------------------------------

BENCHMARKS = {
    "SquareToCircle": {
        "testname": "3D_SquareToCircle_logcost_MonteCarlo",
        "P": P_square,
        "Q": Q_circle,
    },
    "SquareToTwoGaussSide": {
        "testname": "3D_SquareToTwoGaussSide_logcost_MonteCarlo",
        "P": P_square,
        "Q": Q_two_gaussians,
    },
    "Refraction": {
        "testname": "3D_FlatPatchRefraction_logcost_MonteCarlo",
        "P": P_refraction_patch,
        "Q": Q_refraction_patch,
        "kappa": 0.6,
    },
}
