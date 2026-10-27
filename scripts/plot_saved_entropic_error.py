"""Plot saved cross-EOT objectives without self subtraction or new solves."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/snapshots/entropic_error_complete_20261027"
FIELDS = ["map", "experiment", "nk", "refractor_epsilon", "pf_epsilon",
          "eot", "converged", "iterations"]


def read(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def validate(rows):
    keys = {}
    for mode in ("hard", "entropic"):
        selected = [r for r in rows if r["map"] == mode]
        if len(selected) != 120:
            raise ValueError(f"Expected 120 {mode} rows")
        keys[mode] = {(int(r["experiment"]), int(r["nk"]),
                       float(r["refractor_epsilon"])) for r in selected}
        if len(keys[mode]) != 120:
            raise ValueError("Duplicate experiment cases")
        for exp in (1, 2, 3):
            if sum(int(r["experiment"]) == exp for r in selected) != 40:
                raise ValueError("Expected 40 rows per experiment")
        for r in selected:
            if (r["converged"] != "True" or float(r["pf_epsilon"]) != 1/320
                    or not math.isfinite(float(r["eot"])) or float(r["eot"]) <= 0):
                raise ValueError("Invalid EOT objective or evaluation settings")
    if len(rows) != 240 or keys["hard"] != keys["entropic"]:
        raise ValueError("Unmatched map cases")


def draw(ax, rows, experiment, log_y):
    field = "refractor_epsilon" if experiment == 3 else "nk"
    for mode, label, style in [("hard", "Hard map", "o-"),
                               ("entropic", "Entropic surface + Snell", "s-")]:
        selected = sorted((r for r in rows if r["map"] == mode
                           and int(r["experiment"]) == experiment),
                          key=lambda r: float(r[field]))
        ax.plot([float(r[field]) for r in selected],
                [float(r["eot"]) for r in selected], style,
                label=label, markersize=3.5, linewidth=1.8)
    ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    titles = ["1. Joint movement", "2. Fixed refractor epsilon = 1/320",
              "3. Fixed source count = 1600"]
    ax.set_title(titles[experiment-1], fontsize=16, pad=9)
    ax.set_xlabel(r"Refractor $\varepsilon$" if experiment == 3 else r"Source count $N_K$")
    ax.set_ylabel(r"Cross-EOT$_{1/320}$")
    ax.grid(alpha=.25)
    ax.legend(fontsize=11, loc="best")


def render(rows, output):
    with plt.rc_context({"font.size": 14, "xtick.labelsize": 12,
                         "ytick.labelsize": 12}):
        for log_y, suffix in [(False, ""), (True, "_log")]:
            fig, axes = plt.subplots(3, 1, figsize=(7, 10), layout="constrained")
            for exp, ax in enumerate(axes, 1):
                draw(ax, rows, exp, log_y)
                single, panel = plt.subplots(figsize=(8, 5.5), layout="constrained")
                draw(panel, rows, exp, log_y)
                single.savefig(output/f"experiment_{exp}{suffix}.png", dpi=220)
                plt.close(single)
            fig.savefig(output/f"experiments_1_2_3{suffix}.png", dpi=220)
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--hard", type=Path, help="Saved cross-EOT or debiased-table CSV")
    parser.add_argument("--entropic", type=Path, help="Matching saved CSV")
    args = parser.parse_args()
    if bool(args.hard) != bool(args.entropic):
        parser.error("Provide both source CSVs, or neither to replot the bundled table")
    args.output.mkdir(parents=True, exist_ok=True)
    source_hashes = {}
    if args.hard:
        rows = []
        for mode, path in [("hard", args.hard), ("entropic", args.entropic)]:
            source_hashes[mode] = {"path": str(path.as_posix()),
                                  "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for row in read(path):
                rows.append(dict(map=mode, experiment=row["experiment"], nk=row["nk"],
                    refractor_epsilon=row["refractor_epsilon"], pf_epsilon=row["pf_epsilon"],
                    eot=row["eot"] if "eot" in row else row["cost"],
                    converged=row["converged"],
                    iterations=row["cross_iterations"] if "cross_iterations" in row else row["iterations"]))
    else:
        rows = read(args.output/"costs.csv")
    validate(rows)
    if args.hard:
        with (args.output/"costs.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)
    render(rows, args.output)
    if args.hard:
        report = {"metric": "cross_EOT_KL_relative_to_product", "evaluation_epsilon": 1/320,
                  "ground_cost": "squared_euclidean_3d", "self_subtraction": False,
                  "square_root": False, "target_points": 1600, "cases_per_map": 120,
                  "sources": source_hashes, "all_evaluations_converged": True,
                  "costs_sha256": hashlib.sha256((args.output/"costs.csv").read_bytes()).hexdigest()}
        (args.output/"validation.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(f"Rendered 8 plots from 240 saved EOT objectives in {args.output}")


if __name__ == "__main__":
    main()
