"""Measure how much the refractor itself moves when the stopping rule is relaxed.

The refractor is R(x) = exp(softmin_j(g_j - c(x, y_j))) on source rays x, fixed only
up to an overall scale. Metrics against the strict reference cell:
  radius_rel_max   max_x |dlog R - weighted mean dlog R|  (~ max relative radius change)
  radius_rel_mean  source-weighted mean of the same
  normal_deg_max   max angle between unit surface normals (degrees)
  ray_deg_max      max angle between refracted (entropic-map) directions (degrees)
Also reports relative change of raw EOT cost and debiased divergence for reference.

  python scripts/compare_refractor_surfaces.py BENCHMARK_DIR [--threshold 1e-3]
"""
import argparse
import json
from pathlib import Path

import numpy as np

STRICT = (1e-7, 1e-5, 1e-9)


def surface(cell_dir):
    path = next((cell_dir / "entropic").rglob("solve_*.npz"))
    with np.load(path) as data:
        p = data["p"]
        normals = data["surface_normals"]
        return dict(p=p[p > 0] / p.sum(), log_radius=data["log_radius"][p > 0],
                    normals=(normals / np.linalg.norm(normals, axis=1, keepdims=True))[p > 0],
                    rays=data["y_push_3d"])


def angle_deg(a, b):
    return np.degrees(np.arctan2(np.linalg.norm(np.cross(a, b), axis=1), np.sum(a * b, axis=1)))


def compare(reference, other):
    d = other["log_radius"] - reference["log_radius"]
    d = d - np.sum(reference["p"] * d)
    return dict(radius_rel_max=float(np.abs(d).max()), radius_rel_mean=float(np.sum(reference["p"] * np.abs(d))),
                normal_deg_max=float(angle_deg(reference["normals"], other["normals"]).max()),
                ray_deg_max=float(angle_deg(reference["rays"], other["rays"]).max()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--threshold", type=float, default=1e-3)
    args = parser.parse_args()
    cells = {}
    for path in args.root.rglob("cell.json"):
        cell = json.loads(path.read_text())
        setting = (cell["marginal_tolerance"], cell["refractor_tolerance"], cell["pf_tolerance"])
        cells[cell["nk"], round(cell["epsilon"], 12), setting] = (cell, path.parent)
    rows = []
    for (nk, eps, setting), (cell, folder) in sorted(cells.items()):
        ref = cells.get((nk, eps, STRICT))
        if not ref:
            continue
        row = dict(nk=nk, epsilon=eps, marginal=setting[0], change=setting[1], pf=setting[2],
                   minutes=cell["total_seconds"] / 60, speedup=ref[0]["total_seconds"] / cell["total_seconds"])
        try:
            row.update(compare(surface(ref[1]), surface(folder)))
        except StopIteration:
            row["missing_surface"] = True
        for field in ("hard_cost", "entropic_cost", "hard_debiased_divergence", "entropic_debiased_divergence"):
            if field in cell:
                row[field + "_rel"] = abs(cell[field] - ref[0][field]) / abs(ref[0][field])
        row["passes"] = row.get("radius_rel_max", np.inf) <= args.threshold
        rows.append(row)
    (args.root / "surface_comparison.json").write_text(json.dumps(rows, indent=2))
    print(f"{'NK':>5} {'eps':>9} {'marg':>6} {'chg':>6} {'pf':>6} {'min':>7} {'speed':>6} "
          f"{'dR max':>9} {'dR mean':>9} {'nrm deg':>8} {'ray deg':>8} {'hC':>8} {'hD':>8} pass")
    for r in rows:
        print(f"{r['nk']:5d} {r['epsilon']:9.3g} {r['marginal']:6g} {r['change']:6g} {r['pf']:6g} "
              f"{r['minutes']:7.1f} {r['speedup']:6.2f} {r.get('radius_rel_max', np.nan):9.2e} "
              f"{r.get('radius_rel_mean', np.nan):9.2e} {r.get('normal_deg_max', np.nan):8.2e} "
              f"{r.get('ray_deg_max', np.nan):8.2e} {r.get('hard_cost_rel', np.nan):8.1e} "
              f"{r.get('hard_debiased_divergence_rel', np.nan):8.1e} {r['passes']}")


if __name__ == "__main__":
    main()
