"""Plot exact optical OT cost against N, with Sinkhorn costs per epsilon.

  python scripts/plot_exact_optical_ot.py RUN_DIR/exact_optical_ot.json
"""
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

IMAGES = Path(__file__).resolve().parents[1] / "geometric_optics_progress" / "images"


def main(path):
    data = json.loads(Path(path).read_text())
    exact = sorted(data["exact"], key=lambda r: r["nk"])
    n = np.array([r["nk"] for r in exact])
    cost = np.array([r["exact_cost"] for r in exact])
    with plt.rc_context({"font.size": 13, "legend.fontsize": 10}):
        fig, axes = plt.subplots(1, 2, figsize=(14, 6), layout="constrained")
        ax = axes[0]
        sink = data.get("sinkhorn", [])
        eps = sorted({r["epsilon"] for r in sink})
        if eps:
            cmap, norm = plt.get_cmap("viridis"), LogNorm(min(eps), max(eps))
            for e in eps:
                rows = sorted((r for r in sink if r["epsilon"] == e), key=lambda r: r["nk"])
                ax.plot([r["nk"] for r in rows], [r["sinkhorn_cost"] for r in rows], "o-",
                        color=cmap(norm(e)), markersize=3, linewidth=1, alpha=.8)
            bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, pad=0.01)
            bar.set_label(r"Sinkhorn refractor $\varepsilon$")
            ticks = eps[::2]
            bar.set_ticks(ticks, labels=[f"1/{round(1 / t)}" for t in ticks])
            bar.minorticks_off()
        ax.plot(n, cost, "k*-", markersize=9, linewidth=2.2, label="Exact OT (linear assignment)")
        ax.set(xscale="log", xlabel=r"$N$ (source = target)",
               ylabel=r"Optical OT cost, $c=-\log(1-\kappa\,x\cdot y)$",
               title=r"Exact vs Sinkhorn optical OT cost, $\kappa=0.8$")
        ax.legend(loc="best")
        ax.grid(alpha=.25, which="both")

        ax = axes[1]
        ref = cost[-1]
        ax.plot(n[:-1], np.abs(cost[:-1] - ref), "k*-", markersize=9, linewidth=2)
        ax.set(xscale="log", yscale="log", xlabel=r"$N$ (source = target)",
               ylabel=rf"|exact cost(N) $-$ exact cost(N={n[-1]})|",
               title="Exact cost change against the largest N")
        ax.grid(alpha=.25, which="both")
        fig.savefig(IMAGES / "exact_optical_ot_vs_n.png", dpi=200)
    print(np.column_stack([n, cost]))


if __name__ == "__main__":
    main(sys.argv[1])
