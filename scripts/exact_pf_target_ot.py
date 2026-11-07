"""Exact squared-Euclidean OT cost between each grid cell's pushforward and the target.

Pushforward and target both have N points with uniform weights, so the exact
discrete OT plan is a permutation: solved with SciPy's linear_sum_assignment
(no Sinkhorn, no evaluation epsilon). Cost is ||T(x) - y||^2, the same cost
the Sinkhorn evaluation uses.

  python scripts/exact_pf_target_ot.py RUN_DIR [WORKERS]
Writes RUN_DIR/exact_pf_target_ot.json (resumable).
"""
import json
from multiprocessing import Pool
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from refracter.cost import cost_matrix_l2_chunk  # noqa: E402


def solve(job):
    nk, epsilon, mode, npz = job
    with np.load(npz) as d:
        pf, y, p, q = d["y_push_3d"], d["y"], d["p"], d["q"]
    target = y[q > 0]
    if len(pf) != len(target) or np.ptp(p[p > 0]) > 1e-15 or np.ptp(q[q > 0]) > 1e-15:
        raise ValueError(f"not an equal-count uniform problem: {npz}")
    start = time.time()
    cost = cost_matrix_l2_chunk(pf, target)
    rows, cols = linear_sum_assignment(cost)
    return dict(nk=nk, epsilon=epsilon, mode=mode, exact_cost=float(cost[rows, cols].mean()),
                seconds=time.time() - start)


def main(run, workers=3):
    run = Path(run)
    out = run / "exact_pf_target_ot.json"
    rows = json.loads(out.read_text()) if out.exists() else []
    done = {(r["nk"], r["epsilon"], r["mode"]) for r in rows}
    jobs = []
    for path in (run / "grid").rglob("cell.json"):
        cell = json.loads(path.read_text())
        for mode in ("hard", "entropic"):
            if (cell["nk"], cell["epsilon"], mode) not in done:
                npz = next((path.parent / mode).rglob("solve_*.npz"))
                jobs.append((cell["nk"], cell["epsilon"], mode, str(npz)))
    jobs.sort(key=lambda j: -j[0])  # largest first so the pool ends with quick jobs
    print(f"{len(jobs)} solves to run", flush=True)
    with Pool(int(workers)) as pool:
        for row in pool.imap_unordered(solve, jobs):
            rows.append(row)
            out.write_text(json.dumps(rows, indent=2))
            print(row, flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:])
