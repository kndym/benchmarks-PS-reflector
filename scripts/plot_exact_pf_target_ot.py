"""Plot exact pushforward-to-target squared-L2 OT cost over the equal-count grid.

  python scripts/plot_exact_pf_target_ot.py RUN_DIR
Writes geometric_optics_progress/images/exact_pf_{joint_and_n,eps_lines}.png.
"""
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

IMAGES = Path(__file__).resolve().parents[1] / "geometric_optics_progress" / "images"
YLABEL = r"Exact OT cost, $\|T(x)-y\|^2$ (pushforward $\to$ target)"
MODES = [("hard", "Hard map", "o-", "C0"), ("entropic", "Entropic map", "s-", "C1")]


def main(run):
    rows = json.loads((Path(run) / "exact_pf_target_ot.json").read_text())
    cost = {(r["nk"], round(1 / (8 * r["epsilon"])), r["mode"]): r["exact_cost"] for r in rows}
    ms_e = sorted({m for _, m, _ in cost})
    with plt.rc_context({"font.size": 14, "axes.titlesize": 15, "legend.fontsize": 11}):
        fig, axes = plt.subplots(2, 1, figsize=(7.5, 9.5), layout="constrained")
        joint = sorted(n for n, m, mode in cost if mode == "hard" and n == m * m)
        fixed = sorted(n for n, m, mode in cost if mode == "hard" and m == 40)
        panels = [(axes[0], [(n, int(round(n ** .5))) for n in joint],
                   rf"1. Joint movement, $\varepsilon=1/(8\sqrt{{N}})$ ({len(joint)} cases)"),
                  (axes[1], [(n, 40) for n in fixed],
                   rf"2. $N$ movement, fixed $\varepsilon=1/320$ ({len(fixed)} cases)")]
        for ax, keys, title in panels:
            for mode, label, style, color in MODES:
                ax.plot([n for n, _ in keys], [cost[n, m, mode] for n, m in keys], style, color=color,
                        label=label, markersize=4, linewidth=1.6)
            ax.set(xscale="log", yscale="log", title=title, xlabel=r"$N$ (source = target)", ylabel=YLABEL)
            ax.grid(alpha=.25, which="both")
        axes[0].legend(loc="best")
        fig.savefig(IMAGES / "exact_pf_joint_and_n.png", dpi=220)
        plt.close(fig)

        cmap = plt.get_cmap("viridis")
        norm = LogNorm(1 / (8 * max(ms_e)), 1 / (8 * min(ms_e)))
        fig, axes = plt.subplots(1, 2, figsize=(14, 6.5), layout="constrained", sharey=True)
        for ax, (mode, title, _, _) in zip(axes, MODES):
            for me in ms_e:
                ns = sorted(n for n, m, md in cost if m == me and md == mode)
                if ns:
                    ax.plot(ns, [cost[n, me, mode] for n in ns], "o-", color=cmap(norm(1 / (8 * me))),
                            markersize=3.5, linewidth=1.4)
            ax.set(xscale="log", yscale="log", title=f"{title}: exact pushforward-target OT",
                   xlabel=r"$N$ (source = target)")
            ax.grid(alpha=.25, which="both")
        axes[0].set_ylabel(YLABEL)
        bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes, pad=0.01)
        bar.set_label(r"Refractor $\varepsilon$")
        bar.set_ticks([1 / (8 * m) for m in ms_e[::2]], labels=[f"1/{8 * m}" for m in ms_e[::2]])
        bar.minorticks_off()
        fig.savefig(IMAGES / "exact_pf_eps_lines.png", dpi=200)
        plt.close(fig)
    print(f"{len(rows)} solves plotted")


if __name__ == "__main__":
    main(sys.argv[1])
