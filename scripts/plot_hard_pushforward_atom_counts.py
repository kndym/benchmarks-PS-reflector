"""Count occupied hard-map atoms whose aggregated mass differs from 1/NK."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def count_atoms(points, weights, target):
    points, weights, target = map(np.asarray, (points, weights, target))
    nk = len(weights)
    if len(points) != nk or not np.isfinite(points).all() or not np.isfinite(weights).all():
        raise ValueError('Expected finite pushed points and matching source weights')
    if not np.allclose(weights, 1/nk, rtol=1e-12, atol=0):
        raise ValueError('This diagnostic requires uniform source weights 1/NK')
    lookup = {tuple(point): j for j, point in enumerate(target)}
    if len(lookup) != len(target):
        raise ValueError('Target locations must be distinct')
    try:
        indices = [lookup[tuple(point)] for point in points]
    except KeyError as exc:
        raise ValueError('Hard-map output is not an exact saved target location') from exc
    multiplicities = np.bincount(indices, minlength=len(target))
    masses = np.bincount(indices, weights=weights, minlength=len(target))
    occupied = multiplicities > 0
    unequal = occupied & ~np.isclose(masses, 1/nk, rtol=1e-10, atol=1e-14)
    if not np.array_equal(unequal, multiplicities > 1):
        raise ValueError('Aggregated masses disagree with integer collision counts')
    if not np.isclose(masses.sum(), 1, rtol=0, atol=1e-12):
        raise ValueError('Push-forward mass is not conserved')
    return dict(occupied_atoms=int(occupied.sum()), unequal_weight_atoms=int(unequal.sum()),
                singleton_atoms=int((multiplicities == 1).sum()),
                empty_target_locations=int((multiplicities == 0).sum()),
                source_points_in_merged_atoms=int(multiplicities[unequal].sum()),
                max_multiplicity=int(multiplicities.max()),
                unequal_fraction_of_occupied=float(unequal.sum()/occupied.sum()),
                total_mass=float(masses.sum()))


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
    for row in source_rows:
        experiment, nk, eps = int(row['experiment']), int(row['nk']), float(row['refractor_epsilon'])
        if experiment not in (1, 2, 3):
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
            stats = count_atoms(data['y_push_3d'], data['p'], data['y'])
            if len(data['p']) != nk or len(data['y']) != params['target_nk']:
                raise ValueError('Checkpoint sizes disagree with parameters')
        rows.append(dict(experiment=experiment, nk=nk, refractor_epsilon=eps,
                         reference_atom_mass=1/nk, **stats, checkpoint=path.name))
        provenance[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output/'atom_counts.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    context = (f'Hard c-transform; target fixed at {params["target_nk"]}; '
               f'transport/self caps {params["refractor_max_iter"]}/{params["identity_max_iter"]}\n'
               'Count distinct occupied locations after summing coincident source weights; empty locations excluded')
    combined, axes = plt.subplots(1, 3, figsize=(17, 5.5))
    for experiment, panel in zip((1, 2, 3), axes):
        single, axis = plt.subplots(figsize=(8, 6))
        field = 'refractor_epsilon' if experiment == 3 else 'nk'
        cases = sorted((r for r in rows if r['experiment'] == experiment), key=lambda r:r[field])
        title = {1:'1. Joint movement', 2:f'2. Fixed refractor epsilon = {params["fixed_epsilon"]:.6g}',
                 3:f'3. Fixed source NK = {params["fixed_nk"]}'}[experiment]
        for ax in (panel, axis):
            ax.plot([r[field] for r in cases], [r['unequal_weight_atoms'] for r in cases], 'o-', markersize=3)
            ax.set(title=title, xscale='log', xlabel=r'Refractor $\epsilon$' if experiment == 3 else r'Source $N_K$',
                   ylabel=r'Occupied atoms with mass $\ne 1/N_K$')
            ax.grid(alpha=.3)
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            ax.set_ylim(bottom=0)
        for scale in ('linear', 'log'):
            suffix = '' if scale == 'linear' else '_log'
            note = '' if scale == 'linear' else '\nLogarithmic y axis; zero counts omitted'
            if scale == 'log':
                axis.set_yscale('log', nonpositive='mask'); axis.autoscale(axis='y')
            single.suptitle(context+note, fontsize=9)
            single.tight_layout(rect=(0, 0, 1, .86))
            single.savefig(args.output/f'experiment_{experiment}{suffix}.png', dpi=180)
        plt.close(single)
    for scale in ('linear', 'log'):
        suffix = '' if scale == 'linear' else '_log'
        note = '' if scale == 'linear' else '\nLogarithmic y axis; zero counts omitted'
        if scale == 'log':
            for axis in axes:
                axis.set_yscale('log', nonpositive='mask'); axis.autoscale(axis='y')
        combined.suptitle(context+note, fontsize=11)
        combined.tight_layout(rect=(0, 0, 1, .86))
        combined.savefig(args.output/f'experiments_1_2_3{suffix}.png', dpi=180)
    plt.close(combined)
    summary = dict(source=args.source.as_posix(), cases=len(rows),
                   counts={str(e):sum(r['experiment'] == e for r in rows) for e in (1,2,3)},
                   definition='Distinct occupied atoms with aggregated mass != 1/NK; empty target locations excluded',
                   mass_comparison=dict(rtol=1e-10, atol=1e-14),
                   all_source_weights_uniform=True, all_masses_conserved=True,
                   unequal_atom_range=[min(r['unequal_weight_atoms'] for r in rows),max(r['unequal_weight_atoms'] for r in rows)],
                   source_checkpoint_sha256=provenance)
    summary['files_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(args.output.iterdir())}
    (args.output/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    (args.output/'README.md').write_text(
        '# Hard push-forward atom weights\n\n'
        f'{len(rows)} saved cases from the hard-transform 2000/2000-cap sweep; no optical solves were repeated.\n\n'
        'At each distinct occupied push-forward location, sum the weights of all source points mapped there. '
        'Count the resulting atoms whose mass differs from `1/NK` (source NK). For these uniform-source runs, '
        'this is exactly the number of locations receiving at least two source points. '
        'Empty target locations carry zero mass and are excluded from this count.\n\n'
        '![Atom counts](experiments_1_2_3.png)\n\n![Log-log atom counts](experiments_1_2_3_log.png)\n\n'
        '[Per-case counts](atom_counts.csv) include occupied/singleton counts, empty target locations, '
        'source points in merged atoms, maximum multiplicity, and unequal fraction of occupied support. '
        '[Validation](validation.json) records uniform-source and mass-conservation checks and checkpoint hashes.\n',
        encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if not k.endswith('sha256')},indent=2))


if __name__ == '__main__':
    main()
