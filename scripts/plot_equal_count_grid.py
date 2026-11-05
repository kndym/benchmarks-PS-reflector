"""Plot the equal-count (source N = target N) epsilon x N grid for the write-up.

Grid axes use m: N = m^2 and epsilon = 1/(8m), so the joint-movement schedule
epsilon = 1/(8 floor(sqrt N)) is exactly the grid diagonal and the fixed
epsilon = 1/320 N-sweep is the m_eps = 40 row.

  python scripts/plot_equal_count_grid.py RUN_DIR
"""
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
IMAGES = ROOT / "geometric_optics_progress" / "images"
SERIES = [("hard_cost", "Hard map, raw EOT", "o-", "C0"),
          ("entropic_cost", "Entropic map, raw EOT", "s-", "C1"),
          ("hard_debiased_divergence", "Hard map, debiased", "o--", "C0"),
          ("entropic_debiased_divergence", "Entropic map, debiased", "s--", "C1")]


def load(run):
    cells = {}
    for path in (Path(run) / "grid").rglob("cell.json"):
        cell = json.loads(path.read_text())
        cells[cell["nk"], round(1 / (8 * cell["epsilon"]))] = cell
    return cells


def line_panel(ax, rows, x, title, xlabel):
    for field, label, style, color in SERIES:
        ax.plot(x, [r[field] for r in rows], style, color=color, label=label, markersize=4, linewidth=1.6)
    ax.set(xscale="log", yscale="log", title=title, xlabel=xlabel, ylabel=r"Cost at evaluation $\varepsilon=1/320$")
    ax.grid(alpha=.25, which="both")


def main(run):
    cells = load(run)
    params = json.loads((Path(run) / "grid_plan.json").read_text())
    ms_n, ms_e = params["m_n"], params["m_eps"]
    with plt.rc_context({"font.size": 14, "axes.titlesize": 15, "legend.fontsize": 11}):
        fig, axes = plt.subplots(2, 1, figsize=(7.5, 9.5), layout="constrained")
        joint = [cells[m * m, m] for m in ms_n if (m * m, m) in cells]
        line_panel(axes[0], joint, [c["nk"] for c in joint],
                   f"1. Joint movement, $\\varepsilon=1/(8\\sqrt{{N}})$ ({len(joint)} cases)", r"$N$ (source = target)")
        fixed = sorted((c for (n, me), c in cells.items() if me == 40), key=lambda c: c["nk"])
        line_panel(axes[1], fixed, [c["nk"] for c in fixed],
                   f"2. $N$ movement, fixed $\\varepsilon=1/320$ ({len(fixed)} cases)", r"$N$ (source = target)")
        axes[0].legend(loc="best")
        fig.savefig(IMAGES / "equal_count_joint_and_n.png", dpi=220)
        plt.close(fig)

        # One line per refractor epsilon: error (debiased divergence) against N.
        fields = [("hard_debiased_divergence", "Hard map"), ("entropic_debiased_divergence", "Entropic map")]
        cmap = plt.get_cmap("viridis")
        norm = LogNorm(1 / (8 * max(ms_e)), 1 / (8 * min(ms_e)))
        fig, axes = plt.subplots(1, 2, figsize=(14, 6.5), sharey=True, layout="constrained")
        for ax, (field, title) in zip(axes, fields):
            for me in ms_e:
                row = sorted((c for (n, m), c in cells.items() if m == me), key=lambda c: c["nk"])
                if row:
                    ax.plot([c["nk"] for c in row], [c[field] for c in row], "o-", color=cmap(norm(1 / (8 * me))),
                            markersize=3.5, linewidth=1.4, label=f"1/{8 * me}")
            ax.set(xscale="log", yscale="log", title=f"{title}, debiased divergence",
                   xlabel=r"$N$ (source = target)")
            ax.grid(alpha=.25, which="both")
        axes[0].set_ylabel(r"Error: debiased divergence at evaluation $\varepsilon=1/320$")
        bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes, pad=0.01)
        bar.set_label(r"Refractor $\varepsilon$")
        ticks = [1 / (8 * m) for m in ms_e[::2]]
        bar.set_ticks(ticks, labels=[f"1/{8 * m}" for m in ms_e[::2]])
        bar.minorticks_off()
        fig.savefig(IMAGES / "equal_count_eps_lines.png", dpi=200)
        plt.close(fig)
    print(f"{len(cells)} cells plotted")


if __name__ == "__main__":
    main(sys.argv[1])
