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
# Density factories for spherical-patch benchmarks
#
# Each factory accepts the same (theta_min, theta_max, phi_min, phi_max) that
# would be passed to gen_spherical_patch — where phi acts as the polar angle
# (z = cos(phi)) and theta acts as the azimuthal angle.  It returns a callable
#   density(pts: ndarray (N, 3)) -> ndarray (N,)
# that evaluates the density at each unit-sphere point.
#
# All densities are defined using angular distance from the patch centre:
#   d_i = arccos(clip(pts[i] @ centre, -1, 1))
# and are scaled relative to d_max, the maximum angular distance from the
# centre to any corner of the patch.
# ---------------------------------------------------------------------------

def _patch_centre_and_dmax(theta_min, theta_max, phi_min, phi_max):
    """Compute the Cartesian centre and max corner angular distance for a patch.

    The patch is parameterised in gen_spherical_patch convention:
        z       = cos(phi),  phi ∈ [phi_min, phi_max]
        azimuth = theta,     theta ∈ [theta_min, theta_max]

    Returns
    -------
    centre : ndarray, shape (3,), unit vector at patch centre.
    d_max  : float, max angular distance (radians) from centre to any corner.
    """
    phi_c   = 0.5 * (phi_min   + phi_max)
    theta_c = 0.5 * (theta_min + theta_max)
    centre  = np.array([
        np.sin(phi_c) * np.cos(theta_c),
        np.sin(phi_c) * np.sin(theta_c),
        np.cos(phi_c),
    ])
    d_max = 0.0
    for t in (theta_min, theta_max):
        for p in (phi_min, phi_max):
            corner = np.array([np.sin(p) * np.cos(t),
                               np.sin(p) * np.sin(t),
                               np.cos(p)])
            d = np.arccos(np.clip(centre @ corner, -1.0, 1.0))
            if d > d_max:
                d_max = d
    return centre, d_max


def make_patch_uniform(theta_min, theta_max, phi_min, phi_max):
    """Factory: flat (uniform) density — returns 1.0 for every point."""
    def density(pts):
        pts = np.asarray(pts, dtype=np.float64)
        return np.ones(len(pts), dtype=np.float64)
    return density


def make_patch_gaussian(theta_min, theta_max, phi_min, phi_max):
    """Factory: isotropic Gaussian centred at the patch centre.

    sigma = d_max / 3, so the density is ~0.01 at the patch corners.
    """
    centre, d_max = _patch_centre_and_dmax(theta_min, theta_max, phi_min, phi_max)
    sigma = d_max / 3.0

    def density(pts):
        pts = np.asarray(pts, dtype=np.float64)
        d   = np.arccos(np.clip(pts @ centre, -1.0, 1.0))
        return np.exp(-0.5 * (d / sigma) ** 2)
    return density


def make_patch_donut(theta_min, theta_max, phi_min, phi_max):
    """Factory: soft Gaussian annulus centred at the patch centre.

    Peaks at angular distance r0 = 0.5 * d_max from the centre,
    with ring width sigma = d_max / 8.
    """
    centre, d_max = _patch_centre_and_dmax(theta_min, theta_max, phi_min, phi_max)
    r0    = 0.5  * d_max
    sigma = d_max / 8.0

    def density(pts):
        pts = np.asarray(pts, dtype=np.float64)
        d   = np.arccos(np.clip(pts @ centre, -1.0, 1.0))
        return np.exp(-0.5 * ((d - r0) / sigma) ** 2)
    return density


def make_patch_cross(theta_min, theta_max, phi_min, phi_max):
    """Factory: cross of 4 Gaussians placed at N/S/E/W from the patch centre.

    The four arm centres are at angular offset delta = d_max / 3 from the
    patch centre along two orthogonal tangent directions.
    Each Gaussian has sigma = d_max / 6.
    """
    centre, d_max = _patch_centre_and_dmax(theta_min, theta_max, phi_min, phi_max)
    delta = d_max / 3.0
    sigma = d_max / 6.0

    # Two orthonormal tangent vectors to 'centre'
    ref = np.array([0.0, 0.0, 1.0])
    if abs(centre @ ref) > 0.9:          # centre near a pole
        ref = np.array([1.0, 0.0, 0.0])
    e1 = np.cross(centre, ref);  e1 /= np.linalg.norm(e1)
    e2 = np.cross(centre, e1);   e2 /= np.linalg.norm(e2)

    arm_dirs = [centre + delta * e1,
                centre - delta * e1,
                centre + delta * e2,
                centre - delta * e2]
    arm_centres = [v / np.linalg.norm(v) for v in arm_dirs]

    def density(pts):
        pts  = np.asarray(pts, dtype=np.float64)
        vals = np.zeros(len(pts), dtype=np.float64)
        for ac in arm_centres:
            d = np.arccos(np.clip(pts @ ac, -1.0, 1.0))
            vals += np.exp(-0.5 * (d / sigma) ** 2)
        return vals
    return density


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
