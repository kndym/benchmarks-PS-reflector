"""
refracter/pushforward.py

Ray-tracing pushforward on the regular reflector grid.

For each source direction in the push cloud:
  1. Find which patch in the regular grid contains the projected source point.
  2. Compute bilinear interpolation coefficients for the reflector surface
     in the local patch.
  3. Solve the quadratic equation to find the reflector scale h along the ray.
  4. Compute the outward surface normal from the bilinear coefficients.
  5. Apply the Snell's law (mirror) reflection formula.
  6. Project the reflected direction onto the south-pole stereographic plane.
  7. Return (u_reflected, v_reflected, source_density) for valid points.

Matching C++ ``Pushforward_Ref_regular`` / ``Do_RegularPush`` /
``IntCoef`` / ``rightrootofquadraticequation`` / ``binsearch``.
"""

import numpy as np


# ---------------------------------------------------------------------------
# Binary search
# ---------------------------------------------------------------------------

def _binsearch(val: float, arr: np.ndarray) -> int:
    """Return index i such that arr[i] <= val < arr[i+1].

    Matches C++ ``binsearch(0, FinalGridResolution, val, Regular_side)``.

    If val is outside [arr[0], arr[-1]], returns the closest valid index
    (0 or len(arr)-2).

    Parameters
    ----------
    val : float
    arr : 1-D array, monotonically increasing

    Returns
    -------
    i : int, 0 <= i <= len(arr)-2
    """
    a = 0
    b = len(arr)
    while True:
        mid = (a + b) // 2
        if mid == a:
            return a
        if val < arr[mid]:
            b = mid
        else:
            a = mid


# ---------------------------------------------------------------------------
# Bilinear interpolation coefficients
# ---------------------------------------------------------------------------

def _bilinear_coefs(Ref_regular: np.ndarray, i: int, j: int, k_res: int) -> np.ndarray:
    """Compute bilinear interpolation coefficients for the reflector z-component.

    We fit a bilinear surface through the four corners of the grid patch:

        z = a0 + a1*x + a2*y + a3*x*y

    where (x, y, z) are the 3-D coordinates of the reflector points at the
    four corners (i,j), (i+1,j), (i,j+1), (i+1,j+1).

    The system is solved by Gauss elimination following C++ ``IntCoef``.

    Parameters
    ----------
    Ref_regular : (FinalGrid, 3) reflector surface on regular grid
    i, j : int, patch indices (row, col)
    k_res : int, FinalGridResolution (number of grid points per side)

    Returns
    -------
    coef : (4,) array [a0, a1, a2, a3]
    """
    k = k_res

    # Corner flat indices: C++ indexation is:
    #  1=ij, 2=(i+1)j, 3=i(j+1), 4=(i+1)(j+1)
    idx_ij   = i * k + j
    idx_i1j  = (i + 1) * k + j
    idx_ij1  = i * k + (j + 1)
    idx_i1j1 = (i + 1) * k + (j + 1)

    r_ij   = Ref_regular[idx_ij]
    r_i1j  = Ref_regular[idx_i1j]
    r_ij1  = Ref_regular[idx_ij1]
    r_i1j1 = Ref_regular[idx_i1j1]

    # Differences from corner (i,j)
    xdif2 = r_i1j[0]  - r_ij[0]
    xdif3 = r_ij1[0]  - r_ij[0]
    xdif4 = r_i1j1[0] - r_ij[0]

    ydif2 = r_i1j[1]  - r_ij[1]
    ydif3 = r_ij1[1]  - r_ij[1]
    ydif4 = r_i1j1[1] - r_ij[1]

    fdif2 = r_i1j[2]  - r_ij[2]
    fdif3 = r_ij1[2]  - r_ij[2]
    fdif4 = r_i1j1[2] - r_ij[2]

    # xy products at each corner (for the bilinear a3 term)
    xy_ij   = r_ij[0]   * r_ij[1]
    xy_i1j  = r_i1j[0]  * r_i1j[1]
    xy_ij1  = r_ij1[0]  * r_ij1[1]
    xy_i1j1 = r_i1j1[0] * r_i1j1[1]

    pdif2 = xy_i1j  - xy_ij
    pdif3 = xy_ij1  - xy_ij
    pdif4 = xy_i1j1 - xy_ij

    # Gauss elimination (matching C++ exactly)
    if abs(xdif2) < 1e-300:
        # degenerate patch — return zeros to signal failure
        return np.zeros(4)

    by3 = ydif3 - ydif2 * xdif3 / xdif2
    by4 = ydif4 - ydif2 * xdif4 / xdif2

    bp3 = pdif3 - pdif2 * xdif3 / xdif2
    bp4 = pdif4 - pdif2 * xdif4 / xdif2

    bf3 = fdif3 - fdif2 * xdif3 / xdif2
    bf4 = fdif4 - fdif2 * xdif4 / xdif2

    denom_a3 = bp4 - bp3 * by4 / by3 if abs(by3) > 1e-300 else 0.0
    if abs(denom_a3) < 1e-300:
        return np.zeros(4)

    a3 = (bf4 - bf3 * by4 / by3) / denom_a3
    a2 = (bf3 - a3 * bp3) / by3  if abs(by3) > 1e-300 else 0.0
    a1 = (fdif2 - a3 * pdif2 - a2 * ydif2) / xdif2
    a0 = r_ij[2] - a3 * xy_ij - a2 * r_ij[1] - a1 * r_ij[0]

    return np.array([a0, a1, a2, a3])


# ---------------------------------------------------------------------------
# Quadratic root
# ---------------------------------------------------------------------------

def _quadratic_root(a: float, b: float, c: float, avg: float) -> float:
    """Stable quadratic formula; return root closest to avg.

    Solves a*h² + b*h + c = 0 using numerically stable Citardauq formula
    (matches C++ ``rightrootofquadraticequation``).

    Parameters
    ----------
    a, b, c : coefficients
    avg     : reference value for selecting the right root

    Returns
    -------
    h : float, the root closest to avg, or nan on failure.
    """
    diskr = b * b - 4.0 * a * c
    if diskr < 0.0:
        return float('nan')

    diskr = np.sqrt(diskr)

    if b < 0.0:
        root1 = (2.0 * c) / (-b + diskr)
        root2 = (-b + diskr) / (2.0 * a) if abs(a) > 1e-300 else float('nan')
    else:
        root1 = (-b - diskr) / (2.0 * a) if abs(a) > 1e-300 else float('nan')
        root2 = (2.0 * c) / (-b - diskr) if abs(-b - diskr) > 1e-300 else float('nan')

    # Choose root closest to avg
    d1 = abs(root1 - avg) if not np.isnan(root1) else np.inf
    d2 = abs(root2 - avg) if not np.isnan(root2) else np.inf

    if d1 <= d2:
        return root1
    return root2


# ---------------------------------------------------------------------------
# Single-ray push
# ---------------------------------------------------------------------------

def _do_regular_push(sx, sy, sz, Ref_regular, Regular_side, k_res):
    """Trace a single ray and return the reflected direction + validity flag.

    Matches C++ ``Do_RegularPush``.

    Parameters
    ----------
    sx, sy, sz : float, source direction (unit vector on upper hemisphere)
    Ref_regular : (FinalGrid, 3) reflector surface
    Regular_side : (k_res,) grid axis coordinates
    k_res : int, FinalGridResolution

    Returns
    -------
    (rx, ry, rz, ok) where (rx,ry,rz) is the reflected direction and
    ok is True if the ray was successfully traced.
    """
    # Project source onto stereographic plane (north pole)
    denom_proj = 1.0 + sz
    if abs(denom_proj) < 1e-300:
        return 0.0, 0.0, 0.0, False

    X = sx / denom_proj
    Y = sy / denom_proj

    # Find patch
    i = _binsearch(X, Regular_side)
    j = _binsearch(Y, Regular_side)

    # Guard against boundary issues
    if i >= k_res - 1 or j >= k_res - 1:
        return 0.0, 0.0, 0.0, False

    # Bilinear coefficients
    coef = _bilinear_coefs(Ref_regular, i, j, k_res)
    if np.all(coef == 0.0):
        return 0.0, 0.0, 0.0, False

    a0, a1, a2, a3 = coef

    # Average reflector length for initial guess
    idx_ij   = i * k_res + j
    idx_i1j  = (i + 1) * k_res + j
    idx_ij1  = i * k_res + (j + 1)
    idx_i1j1 = (i + 1) * k_res + (j + 1)

    avg_len = 0.25 * (
        np.linalg.norm(Ref_regular[idx_ij]) +
        np.linalg.norm(Ref_regular[idx_i1j]) +
        np.linalg.norm(Ref_regular[idx_ij1]) +
        np.linalg.norm(Ref_regular[idx_i1j1])
    )

    # Solve h*z = a0 + a1*h*x + a2*h*y + a3*h²*x*y
    # => a3*x*y * h² + (a1*x + a2*y - z) * h + a0 = 0
    qa = a3 * sx * sy
    qb = a1 * sx + a2 * sy - sz
    qc = a0

    h = _quadratic_root(qa, qb, qc, avg_len)
    if np.isnan(h) or h <= 0.0:
        return 0.0, 0.0, 0.0, False

    # Surface normal components (unnormalised)
    # n1 = -(a1 + a3*h*y),  n2 = -(a2 + a3*h*x),  n3 = 1  (pointing up)
    # C++ uses:
    #   tempvar1 = a1 + a3*h*y   (= -n1)
    #   tempvar2 = a2 + a3*h*x   (= -n2)
    #   tempvar3 = 2*(a0/h - a3*h*x*y) / (tempvar1²+tempvar2²+1)
    tempvar1 = a1 + a3 * h * sy
    tempvar2 = a2 + a3 * h * sx
    norm_sq = tempvar1 * tempvar1 + tempvar2 * tempvar2 + 1.0
    tempvar3 = 2.0 * (a0 / h - a3 * h * sx * sy) / norm_sq

    # Reflected direction: ray - 2*(ray·n)*n  but in C++ notation:
    #   vec[0] = x + tempvar1*tempvar3
    #   vec[1] = y + tempvar2*tempvar3
    #   vec[2] = z - tempvar3
    rx = sx + tempvar1 * tempvar3
    ry = sy + tempvar2 * tempvar3
    rz = sz - tempvar3

    return rx, ry, rz, True


# ---------------------------------------------------------------------------
# Batch push-forward
# ---------------------------------------------------------------------------

def ray_trace(push_cloud: np.ndarray, Ref_regular: np.ndarray,
              Regular_side: np.ndarray, P_func) -> np.ndarray:
    """Compute the push-forward of the source measure via ray tracing.

    For each source point in push_cloud:
      1. Find the grid patch containing the stereographic projection.
      2. Compute bilinear interpolation of the reflector surface.
      3. Solve for the reflector scale h along the ray.
      4. Compute the surface normal and apply Snell's law reflection.
      5. Project the reflected ray via south-pole stereographic projection.

    Matches C++ ``Pushforward_Ref_regular`` (the Y_Pushed_projected output).

    Parameters
    ----------
    push_cloud : (M, 3) source ray directions (unit vectors, upper hemisphere)
    Ref_regular : (FinalGrid, 3) reflector surface points on regular grid
    Regular_side : (k_res,) grid axis positions
    P_func : callable, P_func(pts) -> density array; evaluated on push_cloud

    Returns
    -------
    result : (K, 3) array of (u_south, v_south, density) for the K valid rays.
    """
    k_res = len(Regular_side)
    M = len(push_cloud)

    # Evaluate source density
    densities = P_func(push_cloud)   # (M,)

    rows = []
    for ind in range(M):
        sx, sy, sz = push_cloud[ind]
        rho = float(densities[ind])

        rx, ry, rz, ok = _do_regular_push(sx, sy, sz, Ref_regular, Regular_side, k_res)
        if not ok:
            continue

        # South-pole stereographic projection of reflected direction
        denom_s = 1.0 - rz
        if abs(denom_s) < 1e-300:
            continue

        u_s = rx / denom_s
        v_s = ry / denom_s

        rows.append((u_s, v_s, rho))

    if len(rows) == 0:
        return np.zeros((0, 3), dtype=np.float64)

    return np.array(rows, dtype=np.float64)
