"""Render the fixed October 26 write-up figures with legible stacked panels."""
try:
    from result_paths import artifact_path, artifact_glob, resolve_path, dated_run_dir, ensure_readme, self_cost_dir
except ModuleNotFoundError:
    from scripts.result_paths import artifact_path, artifact_glob, resolve_path, dated_run_dir, ensure_readme, self_cost_dir
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
SNAPSHOTS = ROOT / "results" / "snapshots"
OUTPUT = ROOT / "geometric_optics_progress" / "images"


def read_rows(path):
    path = artifact_path(Path(path).parent, Path(path).name)
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def render(filename, series, ylabel, log_y=False):
    with plt.rc_context({"font.size": 16, "axes.titlesize": 17,
                         "axes.labelsize": 16, "xtick.labelsize": 14,
                         "ytick.labelsize": 14, "legend.fontsize": 13}):
        fig, axes = plt.subplots(3, 1, figsize=(6.4, 8.7), layout="constrained")
        titles = ["1. Joint movement", "2. Fixed refractor epsilon = 1/320",
                  "3. Fixed source count = 1600"]
        for experiment, ax in enumerate(axes, 1):
            field = "refractor_epsilon" if experiment == 3 else "nk"
            for rows, value, label, style in series:
                selected = sorted((r for r in rows if int(r["experiment"]) == experiment),
                                  key=lambda r: float(r[field]))
                ax.plot([float(r[field]) for r in selected],
                        [float(r[value]) for r in selected], style,
                        label=label, markersize=3.5, linewidth=1.8)
            ax.set_xscale("log")
            if log_y:
                ax.set_yscale("log")
            ax.set_title(titles[experiment - 1], pad=8)
            ax.set_xlabel(r"Refractor $\varepsilon$" if experiment == 3
                          else r"Source count $N_K$")
            ax.set_ylabel(ylabel)
            ax.grid(alpha=0.25)
            if len(series) > 1:
                ax.legend(loc="best", framealpha=0.9)
        target = resolve_path(OUTPUT / filename)
        ensure_readme(target.parent)
        fig.savefig(target, dpi=240)
        plt.close(fig)


def main(comparison=None):
    if comparison is None:
        capture = resolve_path(ROOT / "geometric_optics_progress/metadata/results_capture.json")
        comparison = (resolve_path(json.loads(capture.read_text(encoding="utf-8"))["comparison"])
                      if capture.exists() else
                      resolve_path(SNAPSHOTS / "debiased_resources_20261026" / "comparison"))
    comparison = resolve_path(comparison)
    hard = read_rows(comparison / "hard" / "costs.csv")
    entropic = read_rows(comparison / "entropic" / "costs.csv")
    assert len(hard) == 120 and 0 < len(entropic) <= 120
    render("latest_sinkhorn_comparison.png", [
        (hard, "divergence", "Hard map", "o-"),
        (entropic, "divergence", "Entropic map + Snell", "s-")],
        "Sinkhorn divergence", log_y=True)
    mass = read_rows(SNAPSHOTS / "hard_pushforward_max_error_20261026" /
                     "max_mass_error.csv")
    assert len(mass) == 120
    render("latest_hard_mass_error.png", [
        (mass, "max_abs_mass_difference", "Hard map", "o-")],
        r"Mass error $E_\infty$")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path)
    main(parser.parse_args().comparison)
