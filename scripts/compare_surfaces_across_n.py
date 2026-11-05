"""Compare entropic refractor surfaces across N on one shared set of test rays.

For each grid cell, the entropic surface R(x) = exp(softmin_eps(g_j - c(x, y_j)))
is evaluated at the same test rays (the N=6400 source points that every cell's
source support covers). For each epsilon, every N is compared with the largest
N solved at that epsilon:
  radius_max / radius_mean   |dlog R - mean dlog R| (relative radius change)
  ray_max / ray_mean         angle between refracted rays, degrees
For scale, each surface is also compared with its neighbour one epsilon step
smaller at the same N.

  python scripts/compare_surfaces_across_n.py RUN_DIR
"""
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from refracter.optics import entropic_surface_normals, refract_using_normals  # noqa: E402

KAPPA = 0.8


def solution(folder):
    path = next((folder / "entropic").rglob("solve_*.npz"))
    with np.load(path) as data:
        return {key: data[key] for key in ("x", "y", "p", "q", "g_raw")}


def evaluate(sol, test, epsilon):
    normals, log_radius = entropic_surface_normals(test, sol["y"], sol["q"], sol["g_raw"], epsilon, KAPPA)
    return log_radius, refract_using_normals(test, normals, KAPPA)


def angle_deg(a, b):
    a = a / np.linalg.norm(a, axis=1, keepdims=True)
    b = b / np.linalg.norm(b, axis=1, keepdims=True)
    return np.degrees(np.arctan2(np.linalg.norm(np.cross(a, b), axis=1), np.sum(a * b, axis=1)))


def difference(one, other):
    d = other[0] - one[0]
    d = np.abs(d - d.mean())
    rays = angle_deg(one[1], other[1])
    return dict(radius_max=float(d.max()), radius_mean=float(d.mean()),
                ray_max=float(rays.max()), ray_mean=float(rays.mean()))


def main(run):
    cells = {}
    for path in (Path(run) / "grid").rglob("cell.json"):
        cell = json.loads(path.read_text())
        cells[cell["nk"], round(1 / (8 * cell["epsilon"]))] = (cell, path.parent)
    big = max(n for n, _ in cells)
    reference_cell = next(folder for (n, _), (_, folder) in cells.items() if n == big)
    ref_sol = solution(reference_cell)
    test = ref_sol["x"][ref_sol["p"] > 0]
    # Keep test rays inside every source support: drop rays outside the convex
    # cone of the coarsest source grid by requiring a close source point there.
    coarse = solution(next(folder for (n, _), (_, folder) in cells.items() if n == min(n for n, _ in cells)))
    spacing = np.median(np.sort(1 - coarse["x"] @ coarse["x"].T, axis=1)[:, 1])
    inside = (1 - test @ coarse["x"].T).min(axis=1) <= 2 * spacing
    test = test[inside]
    print(f"{len(test)} shared test rays")
    surfaces = {}
    for key, (cell, folder) in sorted(cells.items()):
        surfaces[key] = evaluate(solution(folder), test, cell["epsilon"])
    rows = []
    m_eps = sorted({m for _, m in cells})
    for m in m_eps:
        ns = sorted(n for n, mm in cells if mm == m)
        for n in ns:
            row = dict(nk=n, epsilon=f"1/{8 * m}", reference_nk=ns[-1],
                       **{f"vs_bigN_{k}": v for k, v in difference(surfaces[ns[-1], m], surfaces[n, m]).items()})
            smaller = [mm for mm in m_eps if mm > m and (n, mm) in surfaces]
            if smaller:
                row.update({f"vs_next_eps_{k}": v for k, v in
                            difference(surfaces[n, smaller[0]], surfaces[n, m]).items()})
            rows.append(row)
    out = Path(run) / "surface_across_n.json"
    out.write_text(json.dumps(rows, indent=2))
    print(f"{'eps':>8} {'N':>5} {'refN':>5} {'dR max':>9} {'dR mean':>9} {'ray max':>8} {'ray mean':>8} "
          f"| next-eps: {'dR mean':>9} {'ray mean':>8}")
    for r in rows:
        print(f"{r['epsilon']:>8} {r['nk']:5d} {r['reference_nk']:5d} {r['vs_bigN_radius_max']:9.2e} "
              f"{r['vs_bigN_radius_mean']:9.2e} {r['vs_bigN_ray_max']:8.3f} {r['vs_bigN_ray_mean']:8.3f} | "
              f"{r.get('vs_next_eps_radius_mean', np.nan):19.2e} {r.get('vs_next_eps_ray_mean', np.nan):8.3f}")


if __name__ == "__main__":
    main(sys.argv[1])
