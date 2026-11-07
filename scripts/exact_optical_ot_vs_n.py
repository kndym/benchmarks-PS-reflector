"""Exact (unregularized) optical OT cost between N source and N target points.

Source and target are the same Halton patches as the equal-count grid
(kappa = 0.8), cost c(x, y) = -log(1 - kappa x.y). With equal counts and
uniform weights the optimal plan is a permutation, so the exact discrete OT
problem is the linear assignment problem, solved with SciPy's
linear_sum_assignment (no Sinkhorn, no epsilon).

For comparison the script also reads the Sinkhorn dual value <p, f> + <q, g>
of each saved entropic grid solve (the entropic OT cost at that refractor
epsilon), when a run directory is given.

  python scripts/exact_optical_ot_vs_n.py OUT.json [RUN_DIR]
"""
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from refracter.qmc import gen_spherical_patch  # noqa: E402

KAPPA = 0.8
SOURCE = np.deg2rad([15, 60, 15, 45])   # theta_min, theta_max, phi_min, phi_max
TARGET = np.deg2rad([18, 36, 18, 36])
SIZES = [100, 196, 400, 784, 1600, 2500, 3249, 4900, 6400, 8100, 10000]


def exact_cost(n):
    x = gen_spherical_patch(n, *SOURCE, skip=0)
    y = gen_spherical_patch(n, *TARGET, skip=0)
    cost = -np.log1p(-KAPPA * (x @ y.T))
    start = time.time()
    rows, cols = linear_sum_assignment(cost)
    seconds = time.time() - start
    return dict(nk=n, exact_cost=float(cost[rows, cols].mean()), seconds=seconds)


def sinkhorn_costs(run):
    out = []
    for path in (Path(run) / "grid").rglob("cell.json"):
        cell = json.loads(path.read_text())
        npz = next((path.parent / "entropic").rglob("solve_*.npz"))
        with np.load(npz) as d:
            value = float(np.sum(d["p"] * d["f_raw"]) + np.sum(d["q"] * d["g_raw"]))
        out.append(dict(nk=cell["nk"], epsilon=cell["epsilon"], sinkhorn_cost=value))
    return sorted(out, key=lambda r: (r["epsilon"], r["nk"]))


def main(out, run=None):
    out = Path(out)
    result = json.loads(out.read_text()) if out.exists() else {"exact": []}
    done = {r["nk"] for r in result["exact"]}
    for n in SIZES:
        if n in done:
            continue
        row = exact_cost(n)
        print(row, flush=True)
        result["exact"].append(row)
        out.write_text(json.dumps(result, indent=2))
    if run:
        result["sinkhorn"] = sinkhorn_costs(run)
        out.write_text(json.dumps(result, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:])
