"""Plot max_j |(T_hard # p)_j - q_j| on the full saved target grid."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def target_mass_error(points, weights, target, target_weights):
    points, weights, target, target_weights = map(np.asarray, (points, weights, target, target_weights))
    if len(points) != len(weights) or len(target) != len(target_weights) or not len(target):
        raise ValueError('Point and weight lengths must match')
    for value in (points, weights, target, target_weights):
        if not np.isfinite(value).all():
            raise ValueError('Measures must be finite')
    for value in (weights, target_weights):
        if value.ndim != 1 or np.any(value < 0) or not np.isclose(value.sum(), 1, rtol=0, atol=1e-12):
            raise ValueError('Weights must be nonnegative probability masses')
    lookup = {tuple(point): j for j, point in enumerate(target)}
    if len(lookup) != len(target):
        raise ValueError('Target locations must be distinct')
    try:
        indices = [lookup[tuple(point)] for point in points]
    except KeyError as exc:
        raise ValueError('Hard push-forward points must lie on the saved target grid') from exc
    mass = np.bincount(indices, weights=weights, minlength=len(target))
    difference = mass-target_weights
    worst = int(np.argmax(np.abs(difference)))
    if not np.isclose(mass.sum(), 1, rtol=0, atol=1e-12):
        raise ValueError('Push-forward mass is not conserved')
    return dict(max_abs_mass_difference=float(abs(difference[worst])),
                max_mass_excess=float(max(0, difference.max())),
                max_mass_deficit=float(max(0, -difference.min())),
                worst_target_index=worst, pushforward_mass_at_worst=float(mass[worst]),
                target_mass_at_worst=float(target_weights[worst]),
                empty_target_locations=int((mass == 0).sum()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    params = json.loads((args.source/'parameters.json').read_text())
    if params.get('pushforward', 'hard') != 'hard':
        raise ValueError('Expected hard-transform checkpoints')
    with (args.source/'costs.csv').open(newline='') as stream:
        source_rows = list(csv.DictReader(stream))
    rows, provenance, seen = [], {}, set()
    target_min, target_max = np.inf, -np.inf
    for row in source_rows:
        experiment, nk, eps = int(row['experiment']), int(row['nk']), float(row['refractor_epsilon'])
        if experiment not in (1,2,3):
            continue
        key = experiment, nk, eps
        if key in seen:
            raise ValueError(f'Duplicate case: {key}')
        seen.add(key)
        matches = list(args.source.glob(f'solve_NK{nk}_*eps{eps:.17g}.npz'))
        if len(matches) != 1:
            raise ValueError(f'Missing/ambiguous checkpoint: {key}')
        path = matches[0]
        with np.load(path) as data:
            stats = target_mass_error(data['y_push_3d'], data['p'], data['y'], data['q'])
            if len(data['p']) != nk or len(data['q']) != params['target_nk']:
                raise ValueError('Checkpoint sizes disagree with parameters')
            target_min = min(target_min, data['q'].min())
            target_max = max(target_max, data['q'].max())
        rows.append(dict(experiment=experiment, nk=nk, refractor_epsilon=eps, **stats, checkpoint=path.name))
        provenance[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output/'max_mass_error.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    context = (f'Hard c-transform; target fixed at {params["target_nk"]}; '
               f'transport/self caps {params["refractor_max_iter"]}/{params["identity_max_iter"]}\n'
               'Maximum absolute difference of probability masses on the full target grid; empty locations included')
    combined, axes = plt.subplots(1, 3, figsize=(17, 5.5))

    def save_variants(figure, panels, stem, size):
        for scale in ('linear','log'):
            for panel in panels:
                if scale == 'log':
                    panel.set_yscale('log', nonpositive='mask'); panel.autoscale(axis='y')
            suffix = '' if scale == 'linear' else '_log'
            note = '' if scale == 'linear' else '\nLogarithmic y axis; exact zero errors omitted'
            figure.suptitle(context+note, fontsize=size)
            figure.tight_layout(rect=(0,0,1,.86))
            figure.savefig(args.output/f'{stem}{suffix}.png', dpi=180)

    for experiment, panel in zip((1,2,3), axes):
        single, axis = plt.subplots(figsize=(8,6))
        field = 'refractor_epsilon' if experiment == 3 else 'nk'
        cases = sorted((r for r in rows if r['experiment'] == experiment), key=lambda r:r[field])
        title = {1:'1. Joint movement',2:f'2. Fixed refractor epsilon = {params["fixed_epsilon"]:.6g}',
                 3:f'3. Fixed source NK = {params["fixed_nk"]}'}[experiment]
        for ax in (panel,axis):
            ax.plot([r[field] for r in cases], [r['max_abs_mass_difference'] for r in cases], 'o-', markersize=3)
            ax.set(title=title, xscale='log', xlabel=r'Refractor $\epsilon$' if experiment == 3 else r'Source $N_K$',
                   ylabel=r'$\max_j |\mu_{\mathrm{PF}}(y_j)-q_j|$')
            ax.grid(alpha=.3); ax.set_ylim(bottom=0)
        save_variants(single,[axis],f'experiment_{experiment}',9)
        plt.close(single)
    save_variants(combined,axes,'experiments_1_2_3',11)
    plt.close(combined)
    report = dict(source=args.source.as_posix(), cases=len(rows),
                  counts={str(e):sum(r['experiment'] == e for r in rows) for e in (1,2,3)},
                  metric='max_j abs(sum_{i: T_hard(x_i)=y_j} p_i - q_j)',
                  empty_target_locations_included=True, all_masses_conserved=True,
                  target_mass_range=[float(target_min),float(target_max)],
                  max_abs_error_range=[min(r['max_abs_mass_difference'] for r in rows),max(r['max_abs_mass_difference'] for r in rows)],
                  source_checkpoint_sha256=provenance)
    report['files_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(args.output.iterdir())}
    (args.output/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    (args.output/'README.md').write_text(
        '# Maximum absolute hard push-forward mass error\n\n'
        f'All {len(rows)} saved hard-transform cases; no optical solves were repeated.\n\n'
        '`max_j |mu_PF(y_j) - q_j|`, where `mu_PF(y_j)` sums source masses mapped to target location '
        '`y_j`. Empty target locations have push-forward mass zero and are included. The reference '
        'is the saved target measure `q`, which here has mass **1/1600** per target location; '
        'source masses are **1/NK**. Values are absolute probability-mass differences, not density '
        'errors or a spatial transport distance.\n\n'
        '![Maximum mass error](experiments_1_2_3.png)\n\n'
        '![Log-log maximum mass error](experiments_1_2_3_log.png)\n\n'
        '[Per-case CSV](max_mass_error.csv) | [Validation and checkpoint hashes](validation.json)\n\n'
        'The CSV also records maximum excess/deficit and the saved zero-based target index where '
        'the absolute error is largest (first index if tied).\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if not k.endswith('sha256')},indent=2))


if __name__ == '__main__':
    main()
