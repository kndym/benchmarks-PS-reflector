"""
refracter/distributions.py

Stereographic projections and density functions for the reflector/refractor benchmarks.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Stereographic projections
# ---------------------------------------------------------------------------

def stereo_north(pts: np.ndarray) -> tuple:
    """North-pole stereographic projection: (x0,x1,x2) → (x0/(1+x2), x1/(1+x2))."""
    pts = np.asarray(pts, dtype=np.float64)
    denom = 1.0 + pts[..., 2]
    u = pts[..., 0] / denom
    v = pts[..., 1] / denom
    return u, v


def stereo_south(pts: np.ndarray) -> tuple:
    """South-pole stereographic projection: (y0,y1,y2) → (y0/(1-y2), y1/(1-y2))."""
    pts = np.asarray(pts, dtype=np.float64)
    denom = 1.0 - pts[..., 2]
    u = pts[..., 0] / denom
    v = pts[..., 1] / denom
    return u, v


def stereo_north_inverse(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Inverse north-pole projection: (u,v) → upper hemisphere point (3,)."""
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
    """Uniform density on the square [-0.5,0.5]² in north-pole stereographic coords."""
    u, v = stereo_north(x)
    inside = (u > -0.5) & (u < 0.5) & (v > -0.5) & (v < 0.5)
    return inside.astype(np.float64)


def Q_circle(y: np.ndarray) -> np.ndarray:
    """Uniform density on the disc |stereo_south(y)| ≤ 0.5."""
    u, v = stereo_south(y)
    inside = (u * u + v * v) <= 0.25
    return inside.astype(np.float64)


def Q_two_gaussians(y: np.ndarray) -> np.ndarray:
    """Two side-by-side Gaussians in south-pole coordinates (SquareToTwoGaussSide)."""
    u, v = stereo_south(y)
    g1 = np.exp(-16.0 * ((u - 0.25) ** 2 + (v - 0.25) ** 2))
    g2 = np.exp(-16.0 * ((u - 0.25) ** 2 + (v + 0.25) ** 2))
    return g1 + g2


# Refraction patch densities — flat indicators on spherical patches (upper hemisphere)

# Default patch bounds (multiples of π)
_SRC_THETA = (np.pi / 12, np.pi / 3)
_SRC_PHI   = (np.pi / 12, np.pi / 4)
_TGT_THETA = (np.pi / 10, np.pi / 5)         # polar angle range
_TGT_PHI   = (np.pi / 3, 5 * np.pi / 12)   # azimuthal range — non-overlapping with source


def P_refraction_patch(x: np.ndarray) -> np.ndarray:
    """Indicator for source patch: θ ∈ [π/12, π/3], φ ∈ [π/12, π/4]."""
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
    """Indicator for target patch: θ ∈ [π/10, π/5], φ ∈ [π/3, 5π/12]."""
    y = np.asarray(y, dtype=np.float64)
    theta = np.arccos(np.clip(y[..., 2], -1.0, 1.0))
    phi   = np.arctan2(y[..., 1], y[..., 0])
    phi   = np.where(phi < 0, phi + 2 * np.pi, phi)

    inside = (
        (theta >= _TGT_THETA[0]) & (theta <= _TGT_THETA[1]) &
        (phi   >= _TGT_PHI[0])   & (phi   <= _TGT_PHI[1])
    )
    return inside.astype(np.float64)


# Density factories: each returns density(pts (N,3)) -> (N,) using angular distance

def _patch_centre_and_dmax(theta_min, theta_max, phi_min, phi_max):
    """Return (centre, d_max): unit-vector patch centre and max corner angular distance."""
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
    """Isotropic Gaussian centred at the patch centre (sigma = d_max/3)."""
    centre, d_max = _patch_centre_and_dmax(theta_min, theta_max, phi_min, phi_max)
    sigma = d_max / 3.0

    def density(pts):
        pts = np.asarray(pts, dtype=np.float64)
        d   = np.arccos(np.clip(pts @ centre, -1.0, 1.0))
        return np.exp(-0.5 * (d / sigma) ** 2)
    return density


def make_patch_donut(theta_min, theta_max, phi_min, phi_max):
    """Gaussian annulus peaking at r0 = 0.5*d_max with sigma = d_max/8."""
    centre, d_max = _patch_centre_and_dmax(theta_min, theta_max, phi_min, phi_max)
    r0    = 0.5  * d_max
    sigma = d_max / 8.0

    def density(pts):
        pts = np.asarray(pts, dtype=np.float64)
        d   = np.arccos(np.clip(pts @ centre, -1.0, 1.0))
        return np.exp(-0.5 * ((d - r0) / sigma) ** 2)
    return density


def make_patch_cross(theta_min, theta_max, phi_min, phi_max):
    """Four Gaussians at N/S/E/W offsets (delta = d_max/3, sigma = d_max/6)."""
    centre, d_max = _patch_centre_and_dmax(theta_min, theta_max, phi_min, phi_max)
    delta = d_max / 3.0
    sigma = d_max / 6.0

    ref = np.array([0.0, 0.0, 1.0])
    if abs(centre @ ref) > 0.9:  # avoid parallel vectors near pole
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
