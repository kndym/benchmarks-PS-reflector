"""Plan the equal-count epsilon x N grid and split it between this PC and Colab.

Axes use m: N = m^2, epsilon = 1/(8m). The joint-movement schedule
epsilon = 1/(8 floor(sqrt N)) is then the grid diagonal, and the fixed
epsilon = 1/320 N-sweep is the m_eps = 40 row. Only cells predicted to take
under --max-minutes are kept. The runtime model was fitted on N=400 cells with
marginal tolerance 1e-3 and potential-change tolerance 1e-3, scaled by (N/400)^2.

  python scripts/plan_equal_count_grid.py RUN_DIR
"""
import argparse
import json
from pathlib import Path

M_N = [10, 14, 20, 28, 40, 57, 80]
M_EPS = [5, 7, 10, 14, 20, 28, 40, 57, 80, 113, 160, 226, 320, 453, 640, 905, 1280]
SETTINGS = dict(marginal_tolerance=1e-3, refractor_tolerance=1e-3, pf_tolerance=1e-9)


def predicted_minutes(nk, epsilon):
    return max(0.3, 1.3 * (1 / (1280 * epsilon)) ** 0.85) * (nk / 400) ** 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--max-minutes", type=float, default=180)
    parser.add_argument("--colab-sessions", nargs="*", default=["refr1", "refr2"])
    parser.add_argument("--colab-share", type=float, default=0.25)
    args = parser.parse_args()
    cells = []
    for mn in M_N:
        for me in M_EPS:
            nk, eps = mn * mn, 1 / (8 * me)
            minutes = predicted_minutes(nk, eps)
            if minutes < args.max_minutes:
                cells.append(dict(nk=nk, epsilon=eps, m_n=mn, m_eps=me, predicted_minutes=minutes, **SETTINGS))
    cells.sort(key=lambda c: -c["predicted_minutes"])
    total = sum(c["predicted_minutes"] for c in cells)
    # Colab gets mid-sized cells (a lost VM loses little); this PC keeps the largest.
    budget = {s: args.colab_share * total / max(1, len(args.colab_sessions)) for s in args.colab_sessions}
    for cell in cells:
        cell["host"] = "local"
        if cell["predicted_minutes"] <= 120:
            session = max(budget, key=budget.get, default=None)
            if session and budget[session] >= cell["predicted_minutes"]:
                cell["host"] = session
                budget[session] -= cell["predicted_minutes"]
    args.run.mkdir(parents=True, exist_ok=True)
    plan = dict(m_n=M_N, m_eps=sorted({c["m_eps"] for c in cells}), settings=SETTINGS,
                max_minutes=args.max_minutes, total_predicted_minutes=total, cells=cells,
                note="N = m_n^2 (source = target), epsilon = 1/(8 m_eps), kappa = 0.8")
    (args.run / "grid_plan.json").write_text(json.dumps(plan, indent=2))
    for host in ["local", *args.colab_sessions]:
        mine = [c for c in cells if c["host"] == host]
        print(f"{host}: {len(mine)} cells, {sum(c['predicted_minutes'] for c in mine):.0f} predicted minutes")


if __name__ == "__main__":
    main()
