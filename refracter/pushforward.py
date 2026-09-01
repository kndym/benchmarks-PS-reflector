"""
refracter/pushforward.py

Ray-tracing push-forward on the regular reflector grid.
Matches C++ Pushforward_Ref_regular / Do_RegularPush / IntCoef.
"""

import numpy as np


# ---------------------------------------------------------------------------
# Binary search
# ---------------------------------------------------------------------------

def _binsearch(val: float, arr: np.ndarray) -> int:
    """Return i such that arr[i] <= val < arr[i+1]; clamped to valid range."""
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
    """Bilinear coefficients [a0,a1,a2,a3] for z = a0+a1*x+a2*y+a3*x*y over patch (i,j)."""
    k = k_res
    idx_ij   = i * k + j
    idx_i1j  = (i + 1) * k + j
    idx_ij1  = i * k + (j + 1)
    idx_i1j1 = (i + 1) * k + (j + 1)

    r_ij   = Ref_regular[idx_ij]
    r_i1j  = Ref_regular[idx_i1j]
    r_ij1  = Ref_regular[idx_ij1]
    r_i1j1 = Ref_regular[idx_i1j1]

    xdif2 = r_i1j[0]  - r_ij[0]
    xdif3 = r_ij1[0]  - r_ij[0]
    xdif4 = r_i1j1[0] - r_ij[0]

    ydif2 = r_i1j[1]  - r_ij[1]
    ydif3 = r_ij1[1]  - r_ij[1]
    ydif4 = r_i1j1[1] - r_ij[1]

    fdif2 = r_i1j[2]  - r_ij[2]
    fdif3 = r_ij1[2]  - r_ij[2]
    fdif4 = r_i1j1[2] - r_ij[2]

    xy_ij   = r_ij[0]   * r_ij[1]
    xy_i1j  = r_i1j[0]  * r_i1j[1]
    xy_ij1  = r_ij1[0]  * r_ij1[1]
    xy_i1j1 = r_i1j1[0] * r_i1j1[1]

    pdif2 = xy_i1j  - xy_ij
    pdif3 = xy_ij1  - xy_ij
    pdif4 = xy_i1j1 - xy_ij

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
    """Numerically stable quadratic root closest to avg (Citardauq formula)."""
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
    """Trace a single ray; return (rx, ry, rz, ok) (C++ Do_RegularPush)."""
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

    qa = a3 * sx * sy
    qb = a1 * sx + a2 * sy - sz
    qc = a0

    h = _quadratic_root(qa, qb, qc, avg_len)
    if np.isnan(h) or h <= 0.0:
        return 0.0, 0.0, 0.0, False

    # Surface normal via bilinear gradient; reflection: ray - 2*(ray·n)*n
    tempvar1 = a1 + a3 * h * sy
    tempvar2 = a2 + a3 * h * sx
    norm_sq = tempvar1 * tempvar1 + tempvar2 * tempvar2 + 1.0
    tempvar3 = 2.0 * (a0 / h - a3 * h * sx * sy) / norm_sq

    rx = sx + tempvar1 * tempvar3
    ry = sy + tempvar2 * tempvar3
    rz = sz - tempvar3

    return rx, ry, rz, True


# ---------------------------------------------------------------------------
# Batch push-forward
# ---------------------------------------------------------------------------

def ray_trace(push_cloud: np.ndarray, Ref_regular: np.ndarray,
              Regular_side: np.ndarray, P_func) -> np.ndarray:
    """Ray-trace push-forward; returns (K,3) array of (u_south, v_south, density)."""
    k_res = len(Regular_side)
    densities = P_func(push_cloud)

    rows = []
    for ind in range(len(push_cloud)):
        sx, sy, sz = push_cloud[ind]
        rho = float(densities[ind])

        # Match C++ Pushforward_Ref_regular: unsupported source points do not
        # contribute a ray or a projected density row.
        if rho == 0.0:
            continue

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
