"""Plot completed/partial entropic Snell sweeps against matching hard-map data."""
import argparse
import csv
import json
import math
import time
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


from atomic_io import atomic_replace


def read_rows(folder):
    path = folder / "costs.csv"
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def render(output, baseline, *, value_field="cost",
           ylabel="PF-to-target EOT (squared Euclidean, 3D)",
           progress_label="entropic evaluations", yscales=("linear", "log")):
    params = json.loads((output / "parameters.json").read_text(encoding="utf-8"))
    hard_params = json.loads((baseline / "parameters.json").read_text(encoding="utf-8"))
    for key in ("nk", "epsilon", "fixed_nk", "fixed_epsilon", "target_nk", "cost_epsilon",
                "refractor_max_iter", "identity_max_iter", "refractor_tolerance", "pf_tolerance",
                "self_cost_mode"):
        if params[key] != hard_params[key]:
            raise ValueError(f"Comparison parameter mismatch: {key}")
    if params["pushforward"] != "entropic":
        raise ValueError("Expected entropic Snell output")
    hard, soft = read_rows(baseline), read_rows(output)
    def key(r):
        return int(r["experiment"]), int(r["nk"]), float(r["refractor_epsilon"])
    hard_keys = {key(r) for r in hard}
    if any(key(r) not in hard_keys or r["converged"] != "True" or not math.isfinite(float(r[value_field]))
           or float(r["pf_epsilon"]) != params["cost_epsilon"] for r in soft):
        raise ValueError("Unmatched or unconverged entropic evaluation")
    expected = sum(len(params['epsilon'] if e == 3 else params['nk']) for e in (1, 2, 3))
    complete = len({key(r) for r in soft}) == expected
    context = (f"Target fixed at {params['target_nk']} points; "
               f"transport/self caps {params['refractor_max_iter']}/{params['identity_max_iter']}\n"
               f"Evaluation epsilon {params['cost_epsilon']:.6g}; "
               f"{progress_label} {len(soft)}/{expected}" + (" (complete)" if complete else " (partial)"))
    nonpositive = sum(float(r[value_field]) <= 0 for r in hard + soft)

    def save_variants(figure, panels, stem, fontsize):
        for scale in yscales:
            for panel in panels:
                if scale == 'log':
                    panel.set_yscale('log', nonpositive='mask')
                else:
                    panel.set_yscale('linear')
            note = '' if scale == 'linear' else '; logarithmic y axis'
            if scale == 'log' and nonpositive:
                note += f'; {nonpositive} nonpositive values omitted'
            figure.suptitle(context + note, fontsize=fontsize)
            figure.tight_layout(rect=(0, 0, 1, .90))
            suffix = '' if scale == 'linear' else '_log'
            temporary = output / f'{stem}{suffix}.tmp.png'
            figure.savefig(temporary, dpi=180)
            atomic_replace(temporary, output / f'{stem}{suffix}.png')

    combined, axes = plt.subplots(1, 3, figsize=(17, 5.5))
    for e, ax in zip((1, 2, 3), axes):
        single, one = plt.subplots(figsize=(8, 6))
        field = 'refractor_epsilon' if e == 3 else 'nk'
        title = {1:'1. Joint movement', 2:f"2. Fixed refractor epsilon = {params['fixed_epsilon']:.6g}",
                 3:f"3. Fixed source NK = {params['fixed_nk']}"}[e]
        for data, label, style in ((hard, 'Hard c-transform', 'o-'),
                                    (soft, 'Entropic surface + Snell', 's-')):
            values = sorted((r for r in data if int(r['experiment']) == e), key=lambda r:float(r[field]))
            xx, yy = [float(r[field]) for r in values], [float(r[value_field]) for r in values]
            for target in (ax, one):
                target.plot(xx, yy, style, label=label, markersize=3)
        for target in (ax, one):
            target.set(title=title, xlabel=r'Refractor $\epsilon$' if e == 3 else r'Source $N_K$',
                       ylabel=ylabel, xscale='log')
            target.grid(alpha=.3)
            target.legend(fontsize=8)
        save_variants(single, [one], f'comparison_experiment_{e}', 10)
        plt.close(single)
    save_variants(combined, axes, 'comparison_experiments_1_2_3', 12)
    plt.close(combined)
    return complete


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--watch', action='store_true')
    parser.add_argument('--metric', choices=('eot', 'divergence'), default='eot')
    parser.add_argument('--log-only', action='store_true', help='Write only log-y plots alongside an existing linear plot worker')
    args = parser.parse_args()
    deadline = time.monotonic() + 7*24*3600
    while True:
        options = dict(yscales=('log',)) if args.log_only else {}
        if args.metric == 'divergence':
            options.update(value_field='divergence',
                           ylabel=r'$S_{1/320}(PF,T)$ (squared Euclidean, 3D)',
                           progress_label='debiased entropic cases')
        complete = render(args.output, args.baseline, **options)
        if complete or not args.watch:
            break
        if time.monotonic() > deadline:
            raise TimeoutError('Comparison plot watch exceeded seven days')
        time.sleep(30)
