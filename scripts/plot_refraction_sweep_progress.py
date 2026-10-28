"""Plot saved main-sweep costs without interrupting numerical workers."""
try:
    from result_paths import artifact_path, artifact_glob, resolve_path, dated_run_dir, ensure_readme, self_cost_dir
except ModuleNotFoundError:
    from scripts.result_paths import artifact_path, artifact_glob, resolve_path, dated_run_dir, ensure_readme, self_cost_dir
import argparse
import csv
import json
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_refraction_pushforward_comparison import comparison_experiments, comparison_case_keys


def render(output):
    parameters = json.loads(artifact_path(output, "parameters.json").read_text())
    with artifact_path(output, "costs.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    complete = True
    expected_keys = comparison_case_keys(parameters)
    for experiment in comparison_experiments(parameters):
        values = [row for row in rows if row["experiment"] == str(experiment)]
        cases = {(experiment, int(row['nk']), float(row['refractor_epsilon'])) for row in values}
        experiment_keys = {key for key in expected_keys if key[0] == experiment}
        expected = len(experiment_keys)
        complete &= experiment_keys <= cases
        if not values:
            continue
        field = "refractor_epsilon" if experiment == 3 else "nk"
        values.sort(key=lambda row: float(row[field]))
        title = {1: "Joint movement", 2: f"Fixed epsilon = {parameters['fixed_epsilon']:.6g}",
                 3: f"Fixed NK = {parameters['fixed_nk']}"}[experiment]
        if len(values) < expected:
            title += f" (partial: {len(values)}/{expected})"
        if parameters.get("target_nk") is not None:
            title += f"\nTarget fixed at {parameters['target_nk']} points"
        if (parameters.get('refractor_max_iter', 17), parameters.get('identity_max_iter', 17)) != (17, 17):
            title += f"\nTransport cap: {parameters['refractor_max_iter'] or 'none'}; self cap: {parameters['identity_max_iter'] or 'none'}"
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot([float(row[field]) for row in values], [float(row["cost"]) for row in values], "o-")
        ax.set(title=f"{experiment}. {title}", xlabel=r"Refractor $\epsilon$" if experiment == 3 else r"Source $N_k$",
               ylabel=rf"$EOT_{{{parameters['cost_epsilon']:.6g}}}(PF,T)$ (squared Euclidean, 3D)" if "cost_epsilon" in parameters else "Original solver cost (approx.)", xscale="log")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        temporary = artifact_path(output, f"experiment_{experiment}.progress.png", for_write=True)
        fig.savefig(temporary, dpi=180)
        plt.close(fig)
        temporary.replace(artifact_path(output, f"experiment_{experiment}.png", for_write=True))
    return complete


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--watch", action="store_true")
    args = parser.parse_args()
    deadline = time.monotonic() + 24 * 3600
    while True:
        ready = all(artifact_path(args.output, name).exists() for name in ("parameters.json", "costs.csv"))
        complete = render(args.output) if ready else False
        if complete or not args.watch:
            break
        if time.monotonic() >= deadline:
            raise TimeoutError("Plot monitoring exceeded 24 hours; saved partial plots remain available")
        time.sleep(30)
