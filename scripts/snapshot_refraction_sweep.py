"""Export only CSV-completed cases from a running fixed-target sweep.

Example: python scripts/snapshot_refraction_sweep.py results/live results/snapshots/name
The worker is never interrupted. No process logs, PIDs, or incomplete cases are copied.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np


def snapshot(source, destination, baseline=None):
    if destination.exists():
        raise FileExistsError(f"Snapshot already exists: {destination}")
    parameters = json.loads((source / 'parameters.json').read_text(encoding='utf-8-sig'))
    with (source / 'costs.csv').open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        fields, rows = reader.fieldnames, list(reader)
    if not rows:
        raise ValueError('No completed evaluations to snapshot')
    selected = {'parameters.json'}
    keys, residuals = set(), {}
    allowed = set()
    for experiment in parameters['experiments']:
        if experiment == 1:
            pairs = [(n, 1/(8*int(np.floor(np.sqrt(n))))) for n in parameters['nk']]
        elif experiment == 2:
            pairs = [(n, parameters['fixed_epsilon']) for n in parameters['nk']]
        elif experiment == 3:
            pairs = [(parameters['fixed_nk'], e) for e in parameters['epsilon']]
        else:
            raise ValueError('Snapshot exporter supports experiments 1–3 only')
        allowed.update((experiment, n, e, parameters['cost_epsilon']) for n, e in pairs)
    for row in rows:
        nk, eps, eta = int(row['nk']), float(row['refractor_epsilon']), float(row['pf_epsilon'])
        key = (int(row['experiment']), nk, eps, eta)
        if key not in allowed or key in keys or row['converged'] != 'True' or not np.isfinite(float(row['cost'])):
            raise ValueError(f'Duplicate, nonfinite or unconverged evaluation: {key}')
        keys.add(key)
        matches = list(source.glob(f'solve_NK{nk}_*eps{eps:.17g}.npz'))
        if len(matches) != 1:
            raise ValueError(f'Ambiguous/missing solve checkpoint: {key}')
        solve = matches[0]
        suffix = solve.stem.removeprefix('solve_')
        pf = source / f'pf_{suffix}_eval{eta:.17g}.json'
        residual = source / f'residual_{suffix}.json'
        cost, info = json.loads(pf.read_text(encoding='utf-8'))
        if not info['converged'] or not np.isclose(cost, float(row['cost']), atol=1e-15, rtol=0):
            raise ValueError(f'PF checkpoint does not match CSV: {key}')
        diagnostic = json.loads(residual.read_text(encoding='utf-8'))
        if not all(np.isfinite(diagnostic[k]) for k in ('source_marginal_l1', 'target_marginal_l1')):
            raise ValueError('Nonfinite marginal diagnostic')
        residuals[solve.name] = diagnostic
        with np.load(solve) as data:
            if not all(np.isfinite(data[k]).all() for k in data.files):
                raise ValueError(f'Nonfinite checkpoint arrays: {solve}')
            if len(data['y']) != len(data['q']) or len(data['y']) != parameters['target_nk']:
                raise ValueError('Target length mismatch')
            if len(data['p']) != nk or len(data['y_push_3d']) != nk:
                raise ValueError('Source length mismatch')
            if parameters.get('pushforward') == 'entropic':
                if not np.allclose(np.sum(data['x']*data['surface_normals'], axis=1), 1, atol=1e-10, rtol=0):
                    raise ValueError('Invalid entropic normals')
                if not np.allclose(np.linalg.norm(data['y_push_3d'], axis=1), 1, atol=1e-10, rtol=0):
                    raise ValueError('Nonunit transmitted rays')
        selected.update((solve.name, pf.name, residual.name))
    expected = len(allowed)
    destination.mkdir(parents=True)
    for name in selected:
        shutil.copy2(source / name, destination / name)
    with (destination / 'costs.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    if baseline:
        from plot_refraction_pushforward_comparison import render
        render(destination, baseline)
    report = dict(source=str(source.as_posix()), captured_at=datetime.now(timezone.utc).isoformat(),
                  status='complete' if len(rows) == expected else 'partial', rows=len(rows),
                  expected_rows=expected, unique_solves=len(residuals),
                  counts={str(e):sum(int(r['experiment']) == e for r in rows) for e in (1, 2, 3)},
                  all_pf_converged=True,
                  marginal_ranges={k:[min(r[k] for r in residuals.values()), max(r[k] for r in residuals.values())]
                                   for k in ('source_marginal_l1', 'target_marginal_l1')})
    report['files_sha256'] = {p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in sorted(destination.iterdir()) if p.is_file()}
    (destination / 'snapshot.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    (destination / 'README.md').write_text(
        f"# Entropic-surface sweep snapshot\n\n**{report['status'].upper()}: {len(rows)}/{expected} evaluations.**\n\n"
        "This immutable snapshot contains only CSV-completed cases and their validated solve, "
        "PF, and residual checkpoints. The live worker continues in the source directory named "
        "in `snapshot.json`. Comparison images were rendered from this snapshot's CSV, not copied "
        "from the changing live output. Target fixed at 1600 points; transport/self caps 2000/2000; "
        "PF-to-target EOT uses squared Euclidean 3D cost and epsilon 1/320.\n\n"
        "The hard comparison is `../../refraction_sweeps_extended_high_iter/`. "
        "The snapshot records hashes and validation counts. Incomplete experiments are not final results.\n",
        encoding='utf-8')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--baseline', type=Path)
    args = parser.parse_args()
    result = snapshot(args.source, args.destination, args.baseline)
    print(f"Saved {result['status']} snapshot: {result['rows']}/{result['expected_rows']} rows")
