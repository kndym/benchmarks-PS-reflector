"""
generate_pointclouds.py
=======================
Regenerates the 3 Monte-Carlo / quasi-Monte-Carlo point-cloud header files
at any desired size.

Method
------
Source and target patches are defined by spherical-coordinate ranges
(theta = polar angle from +z, phi = azimuthal angle, both in degrees).
Points are sampled quasi-randomly inside each patch using the **Halton
sequence** (base 2 and 3).  Sampling is uniform in solid angle by drawing
cos(theta) uniformly and phi uniformly, then converting to Cartesian:

  x = sin(theta)*cos(phi),  y = sin(theta)*sin(phi),  z = cos(theta)

Output files
------------
  QuasiMonteCarlo/MonteCarlo_Pointcloud_3D_128.h   -> NK   points  (x + y arrays + embedded small arrays)
  SmallGrid/3D_MonteCarlo_Pointcloud_small.h        -> NK_small points
  PushForward/PushForward_Cloud_128.h               -> Push_Cloud_Size points  (source patch only)

Usage
-----
  python generate_pointclouds.py [NK] [NK_small] [options]

  python generate_pointclouds.py          # defaults: NK=1600, NK_small=200
  python generate_pointclouds.py 1600 200 # same as defaults (fast, original sizes)
  python generate_pointclouds.py 16488 381 # large (slow, original data sizes)

  # Refraction example (both patches upper hemisphere, Figure 4):
  python generate_pointclouds.py 1600 200 \
      --src-theta-min 15 --src-theta-max 60 --src-phi-min 15 --src-phi-max 45 \
      --tgt-theta-min 18 --tgt-theta-max 36 --tgt-phi-min 18 --tgt-phi-max 36
"""

import argparse
import math
import os

# ── Argument parsing ───────────────────────────────────────────────────────────

_p = argparse.ArgumentParser(description='Generate QMC point-cloud header files.')
_p.add_argument('NK',       type=int, nargs='?', default=1600,
                help='Main point-cloud size (default: 1600)')
_p.add_argument('NK_small', type=int, nargs='?', default=200,
                help='Warm-start point-cloud size (default: 200)')
_p.add_argument('--src-theta-min', type=float, default=0.0,   dest='src_theta_min',
                metavar='DEG', help='Source patch min polar angle in degrees (default: 0)')
_p.add_argument('--src-theta-max', type=float, default=60.0,  dest='src_theta_max',
                metavar='DEG', help='Source patch max polar angle in degrees (default: 60)')
_p.add_argument('--src-phi-min',   type=float, default=0.0,   dest='src_phi_min',
                metavar='DEG', help='Source patch min azimuthal angle in degrees (default: 0)')
_p.add_argument('--src-phi-max',   type=float, default=360.0, dest='src_phi_max',
                metavar='DEG', help='Source patch max azimuthal angle in degrees (default: 360)')
_p.add_argument('--tgt-theta-min', type=float, default=120.0, dest='tgt_theta_min',
                metavar='DEG', help='Target patch min polar angle in degrees (default: 120)')
_p.add_argument('--tgt-theta-max', type=float, default=180.0, dest='tgt_theta_max',
                metavar='DEG', help='Target patch max polar angle in degrees (default: 180)')
_p.add_argument('--tgt-phi-min',   type=float, default=0.0,   dest='tgt_phi_min',
                metavar='DEG', help='Target patch min azimuthal angle in degrees (default: 0)')
_p.add_argument('--tgt-phi-max',   type=float, default=360.0, dest='tgt_phi_max',
                metavar='DEG', help='Target patch max azimuthal angle in degrees (default: 360)')
_args = _p.parse_args()

# ── Point-cloud sizes ──────────────────────────────────────────────────────────

NK       = _args.NK
NK_small = _args.NK_small

# ── Patch configuration (spherical coordinates, degrees) ──────────────────────
#
#   theta : polar angle measured from the +z axis  (0 = north pole, 180 = south)
#   phi   : azimuthal angle in the x-y plane       (0 = +x axis, 360 = full circle)
#
SRC_THETA_MIN = _args.src_theta_min
SRC_THETA_MAX = _args.src_theta_max
SRC_PHI_MIN   = _args.src_phi_min
SRC_PHI_MAX   = _args.src_phi_max

TGT_THETA_MIN = _args.tgt_theta_min
TGT_THETA_MAX = _args.tgt_theta_max
TGT_PHI_MIN   = _args.tgt_phi_min
TGT_PHI_MAX   = _args.tgt_phi_max

# ── Output paths ───────────────────────────────────────────────────────────────

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


# ── Sphere sampling ────────────────────────────────────────────────────────────

def sphere_point(theta: float, phi: float) -> tuple:
    """
    Convert spherical coordinates to a unit-sphere Cartesian point.

    theta : polar angle from +z axis  [radians]
    phi   : azimuthal angle           [radians]
    """
    st = math.sin(theta)
    return (st * math.cos(phi), st * math.sin(phi), math.cos(theta))


def generate_points(n: int,
                    theta_min_deg: float, theta_max_deg: float,
                    phi_min_deg: float,   phi_max_deg: float,
                    skip: int = 0) -> list:
    """
    Generate n quasi-random sphere points uniform in solid angle within
    the patch [theta_min, theta_max] x [phi_min, phi_max] (degrees).

    Uniform-in-area sampling:
      - cos(theta) sampled uniformly in [cos(theta_max), cos(theta_min)]
      - phi        sampled uniformly in [phi_min, phi_max]
    Uses the Halton sequence (base-2 for phi, base-3 for cos(theta)).
    """
    theta_min = math.radians(theta_min_deg)
    theta_max = math.radians(theta_max_deg)
    phi_min   = math.radians(phi_min_deg)
    phi_max   = math.radians(phi_max_deg)

    # cos(theta) decreases as theta increases
    cos_lo = math.cos(theta_max)   # smaller value (theta_max > theta_min)
    cos_hi = math.cos(theta_min)   # larger  value

    points = []
    idx = skip
    while len(points) < n:
        u = halton(idx, 2)
        v = halton(idx, 3)
        cos_theta = cos_lo + v * (cos_hi - cos_lo)
        theta = math.acos(cos_theta)
        phi   = phi_min + u * (phi_max - phi_min)
        points.append(sphere_point(theta, phi))
        idx += 1
    return points


def fmt_row(pt: tuple) -> str:
    """Format a (x,y,z) triple as a C++ initializer row."""
    return ", ".join(FMT.format(v) for v in pt) + ", "


# ── Writers ───────────────────────────────────────────────────────────────────

def write_large(path: str, nk: int, nk_small: int):
    """Write MonteCarlo_Pointcloud_3D_128.h"""
    x_pts  = generate_points(nk,       SRC_THETA_MIN, SRC_THETA_MAX, SRC_PHI_MIN, SRC_PHI_MAX, skip=0)
    y_pts  = generate_points(nk,       TGT_THETA_MIN, TGT_THETA_MAX, TGT_PHI_MIN, TGT_PHI_MAX, skip=0)
    xs_pts = generate_points(nk_small, SRC_THETA_MIN, SRC_THETA_MAX, SRC_PHI_MIN, SRC_PHI_MAX, skip=nk)
    ys_pts = generate_points(nk_small, TGT_THETA_MIN, TGT_THETA_MAX, TGT_PHI_MIN, TGT_PHI_MAX, skip=nk)

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
    xs_pts = generate_points(nk_small, SRC_THETA_MIN, SRC_THETA_MAX, SRC_PHI_MIN, SRC_PHI_MAX, skip=NK)
    ys_pts = generate_points(nk_small, TGT_THETA_MIN, TGT_THETA_MAX, TGT_PHI_MIN, TGT_PHI_MAX, skip=NK)

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
    """Write PushForward_Cloud_128.h  (source patch only)"""
    pts = generate_points(nk, SRC_THETA_MIN, SRC_THETA_MAX, SRC_PHI_MIN, SRC_PHI_MAX, skip=0)

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
    print(f"  Source patch: theta=[{SRC_THETA_MIN}, {SRC_THETA_MAX}] deg,  phi=[{SRC_PHI_MIN}, {SRC_PHI_MAX}] deg")
    print(f"  Target patch: theta=[{TGT_THETA_MIN}, {TGT_THETA_MAX}] deg,  phi=[{TGT_PHI_MIN}, {TGT_PHI_MAX}] deg")
    print()

    write_large(OUT_LARGE, NK, NK_small)
    write_small(OUT_SMALL, NK_small)
    write_push(OUT_PUSH,  NK)

    print()
    print("Done. Rebuild with: make -C BenchmarkCode")
