"""Plot saved cross-EOT objectives without self subtraction or new solves."""
try:
    from result_paths import artifact_path, artifact_glob, resolve_path, dated_run_dir, ensure_readme, self_cost_dir
except ModuleNotFoundError:
    from scripts.result_paths import artifact_path, artifact_glob, resolve_path, dated_run_dir, ensure_readme, self_cost_dir
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
DEFAULT_OUTPUT = ROOT / "results/published/regular_entropic_error_oct_27_26"
FIELDS = ["map", "experiment", "nk", "refractor_epsilon", "pf_epsilon",
          "eot", "converged", "iterations"]


def read(path):
    path = artifact_path(Path(path).parent, Path(path).name)
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
                single.savefig(artifact_path(output, f"experiment_{exp}{suffix}.png", for_write=True), dpi=220)
                plt.close(single)
            fig.savefig(artifact_path(output, f"experiments_1_2_3{suffix}.png", for_write=True), dpi=220)
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--hard", type=Path, help="Saved cross-EOT or debiased-table CSV")
    parser.add_argument("--entropic", type=Path, help="Matching saved CSV")
    args = parser.parse_args()
    args.output = resolve_path(args.output)
    if args.hard:
        args.hard = artifact_path(args.hard.parent, args.hard.name)
    if args.entropic:
        args.entropic = artifact_path(args.entropic.parent, args.entropic.name)
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
        rows = read(artifact_path(args.output, "costs.csv"))
    validate(rows)
    if args.hard:
        with artifact_path(args.output, "costs.csv", for_write=True).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)
    render(rows, args.output)
    if args.hard:
        report = {"metric": "cross_EOT_KL_relative_to_product", "evaluation_epsilon": 1/320,
                  "ground_cost": "squared_euclidean_3d", "self_subtraction": False,
                  "square_root": False, "target_points": 1600, "cases_per_map": 120,
                  "sources": source_hashes, "all_evaluations_converged": True,
                  "costs_sha256": hashlib.sha256((artifact_path(args.output, "costs.csv")).read_bytes()).hexdigest()}
        artifact_path(args.output, "validation.json", for_write=True).write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(f"Rendered 8 plots from 240 saved EOT objectives in {args.output}")


if __name__ == "__main__":
    main()
