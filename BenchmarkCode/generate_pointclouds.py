"""
generate_pointclouds.py
=======================
Regenerates the 3 Monte-Carlo / quasi-Monte-Carlo point-cloud header files
at any desired size.

Method
------
The original code uses **inverse stereographic projection** from a square
[-0.6, 0.6]^2 -> upper hemisphere (source x) / lower hemisphere (target y).

  Given (X, Y) in R^2:
      N2 = X^2 + Y^2
      sphere = ( 2X/(1+N2),  2Y/(1+N2),  ±(1-N2)/(1+N2) )

We sample (X, Y) quasi-randomly inside the disk that maps to the sphere cap
using the **Halton sequence** (base 2 and 3), which is low-discrepancy and
gives much better coverage than pure Monte Carlo.

Output files
------------
  QuasiMonteCarlo/MonteCarlo_Pointcloud_3D_128.h   -> NK   points  (x + y arrays + embedded small arrays)
  SmallGrid/3D_MonteCarlo_Pointcloud_small.h        -> NK_small points
  PushForward/PushForward_Cloud_128.h               -> Push_Cloud_Size points  (same as NK, positive z only)

Usage
-----
  python generate_pointclouds.py [NK] [NK_small]

  python generate_pointclouds.py          # defaults: NK=1600, NK_small=200
  python generate_pointclouds.py 1600 200 # same as defaults (fast, original sizes)
  python generate_pointclouds.py 16488 381 # large (slow, original data sizes)
"""

import sys
import math
import os

# ── Configuration ──────────────────────────────────────────────────────────────

NK        = int(sys.argv[1]) if len(sys.argv) > 1 else 1600
NK_small  = int(sys.argv[2]) if len(sys.argv) > 2 else 200
SQUARE_HALF = 0.6          # domain is [-0.6, 0.6]^2 (matches the C++ code)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

OUT_LARGE  = os.path.join(SCRIPT_DIR, "QuasiMonteCarlo", "MonteCarlo_Pointcloud_3D_128.h")
OUT_SMALL  = os.path.join(SCRIPT_DIR, "SmallGrid",       "3D_MonteCarlo_Pointcloud_small.h")
OUT_PUSH   = os.path.join(SCRIPT_DIR, "PushForward",     "PushForward_Cloud_128.h")

FMT = "{:.21e}"   # 21-digit scientific notation, matching original data

# ── Halton sequence ────────────────────────────────────────────────────────────

def halton(index: int, base: int) -> float:
    """Return the index-th term of the Halton sequence in the given base."""
    result = 0.0
    f = 1.0
    i = index
    while i > 0:
        f /= base
        result += f * (i % base)
        i //= base
    return result


def sphere_point(X: float, Y: float, upper: bool) -> tuple:
    """Inverse stereographic projection of (X,Y) -> unit sphere point."""
    N2 = X * X + Y * Y
    denom = 1.0 + N2
    z_sign = 1.0 if upper else -1.0
    return (
        2.0 * X / denom,
        2.0 * Y / denom,
        z_sign * (1.0 - N2) / denom,
    )


def generate_points(n: int, upper: bool, skip: int = 0) -> list:
    """
    Generate n quasi-random sphere points using the Halton sequence.
    Points are drawn from a square [-SQUARE_HALF, SQUARE_HALF]^2
    (some will map near the equator; all are unit-sphere points).
    """
    points = []
    idx = skip  # starting index in the Halton sequence
    while len(points) < n:
        # Map [0,1)^2 -> [-SQUARE_HALF, SQUARE_HALF]^2
        u = halton(idx, 2)
        v = halton(idx, 3)
        X = (u - 0.5) * 2.0 * SQUARE_HALF
        Y = (v - 0.5) * 2.0 * SQUARE_HALF
        points.append(sphere_point(X, Y, upper))
        idx += 1
    return points


def fmt_row(pt: tuple) -> str:
    """Format a (x,y,z) triple as a C++ initializer row."""
    return ", ".join(FMT.format(v) for v in pt) + ", "


# ── Writers ───────────────────────────────────────────────────────────────────

def write_large(path: str, nk: int, nk_small: int):
    """Write MonteCarlo_Pointcloud_3D_128.h"""
    x_pts = generate_points(nk,       upper=True,  skip=0)
    y_pts = generate_points(nk,       upper=False, skip=0)       # antipodal
    xs_pts = generate_points(nk_small, upper=True,  skip=nk)     # offset to avoid overlap
    ys_pts = generate_points(nk_small, upper=False, skip=nk)

    with open(path, "w", newline="") as f:
        f.write(f"#ifndef MonteCarlo_Pointcloud_3D_{nk}\r\n")
        f.write(f"#define MonteCarlo_Pointcloud_3D_{nk}\r\n")
        f.write("\r\n")
        f.write(f"const int NK={nk};\r\n")
        f.write("const int dim=3;\r\n")
        f.write("\r\n")

        # x array
        f.write("double x[NK][dim]=\r\n{\r\n")
        for pt in x_pts:
            f.write(fmt_row(pt) + "\r\n")
        f.write("};\r\n\r\n\r\n")

        # y array
        f.write("double y[NK][dim]=\r\n{\r\n")
        for pt in y_pts:
            f.write(fmt_row(pt) + "\r\n")
        f.write("};\r\n\r\n\r\n")

        # embedded small arrays (also used via SmallGrid include)
        f.write(f"const int NK_small={nk_small};\r\n")
        f.write("\r\n\r\n")
        f.write("double x_small[NK_small][dim]=\r\n{\r\n")
        for pt in xs_pts:
            f.write(fmt_row(pt) + "\r\n")
        f.write("};\r\n\r\n\r\n")

        f.write("double y_small[NK_small][dim]=\r\n{\r\n")
        for pt in ys_pts:
            f.write(fmt_row(pt) + "\r\n")
        f.write("};\r\n\r\n\r\n")

        f.write("#endif\r\n")

    print(f"  Written {path}  (NK={nk}, NK_small={nk_small})")


def write_small(path: str, nk_small: int):
    """Write 3D_MonteCarlo_Pointcloud_small.h"""
    # Use same skip offset as write_large so points are different from large cloud
    xs_pts = generate_points(nk_small, upper=True,  skip=NK)
    ys_pts = generate_points(nk_small, upper=False, skip=NK)

    with open(path, "w", newline="") as f:
        f.write(f"const int NK_small={nk_small};\r\n")
        f.write("\r\n")
        f.write("//monte-carlo cloud for both x and y. where reflector cost is applied, last values of y should change sign. \r\n")
        f.write("\r\n\r\n")

        f.write("double x_small[NK_small][3]=\r\n{ \r\n")
        for pt in xs_pts:
            f.write(fmt_row(pt) + "\r\n")
        f.write("};\r\n")
        f.write("\r\n\r\n\r\n\r\n\r\n\r\n\r\n")

        f.write("double y_small[NK_small][3]=\r\n{ \r\n")
        for pt in ys_pts:
            f.write(fmt_row(pt) + "\r\n")
        f.write("};\r\n\r\n\r\n")

    print(f"  Written {path}  (NK_small={nk_small})")


def write_push(path: str, nk: int):
    """Write PushForward_Cloud_128.h  (positive-z source hemisphere only)"""
    pts = generate_points(nk, upper=True, skip=0)

    with open(path, "w", newline="") as f:
        f.write("#ifndef PushForward_Cloud_128\r\n")
        f.write("#define PushForward_Cloud_128\r\n")
        f.write("\r\n")
        f.write(f"const int Push_Cloud_Size={nk};\r\n")
        f.write("\r\n")
        f.write("double Push_Cloud[Push_Cloud_Size][3]=\r\n{\r\n")
        for pt in pts:
            f.write(fmt_row(pt) + "\r\n")
        f.write("};\r\n")
        f.write("\r\n\r\n")
        f.write("#endif\r\n")

    print(f"  Written {path}  (Push_Cloud_Size={nk})")


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print(f"Generating point clouds: NK={NK}, NK_small={NK_small}")
    print(f"  (Halton quasi-random sequence, inverse stereographic projection)")
    print()

    write_large(OUT_LARGE, NK, NK_small)
    write_small(OUT_SMALL, NK_small)
    write_push(OUT_PUSH,  NK)

    print()
    print("Done. Rebuild with: make -C BenchmarkCode")
