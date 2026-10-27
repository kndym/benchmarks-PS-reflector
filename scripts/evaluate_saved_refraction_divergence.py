"""Debiased hard-vs-entropic comparison from saved maps; never runs optics.

Reuses each saved cross EOT and caches only missing self EOT terms. The target
self term is shared. With --watch, incorporates completed source CSV rows as
the independent original worker advances. Run this worker at low priority.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from refracter.sinkhorn import entropic_ot_cost, sinkhorn_divergence_from_terms
from plot_refraction_pushforward_comparison import read_rows, render
from atomic_io import atomic_replace
from sinkhorn_resources import enable_threaded_reductions


def atomic_json(path, value):
    temporary = path.with_suffix('.tmp.json')
    temporary.write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8')
    atomic_replace(temporary, path)


def row_key(row):
    return int(row['experiment']), int(row['nk']), float(row['refractor_epsilon'])


def measure_hash(points, weights):
    weights = weights / weights.sum()
    h = hashlib.sha256()
    for value in (points, weights):
        value = np.ascontiguousarray(value, dtype='<f8')
        h.update(str(value.shape).encode())
        h.update(value.tobytes())
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hard', type=Path, required=True)
    parser.add_argument('--entropic', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--watch', action='store_true')
    parser.add_argument('--chunk-size', type=int, default=256)
    parser.add_argument('--cache-entries', type=int, default=2_000_000)
    parser.add_argument('--workers', type=int, default=1, help='Threads for independent SciPy reductions')
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('Workers must be positive')
    enable_threaded_reductions(args.workers)
    params = json.loads((args.entropic/'parameters.json').read_text(encoding='utf-8'))
    hard_params = json.loads((args.hard/'parameters.json').read_text(encoding='utf-8'))
    for name in ('nk', 'epsilon', 'fixed_nk', 'fixed_epsilon', 'target_nk',
                 'cost_epsilon', 'refractor_max_iter', 'identity_max_iter',
                 'refractor_tolerance', 'pf_tolerance', 'self_cost_mode'):
        if params[name] != hard_params[name]:
            raise ValueError(f'Unmatched hard/entropic parameters: {name}')
    eta, tolerance, max_iter = params['cost_epsilon'], params['pf_tolerance'], params['pf_max_iter']
    if hard_params['pf_max_iter'] != max_iter or params['pushforward'] != 'entropic':
        raise ValueError('Unmatched evaluation settings or map type')
    settings = dict(epsilon=eta, tolerance=tolerance, max_iter=max_iter,
                    ground_cost='squared_euclidean_3d', objective='KL_relative_to_product')
    args.output.mkdir(parents=True, exist_ok=True)
    cache = args.output/'self_cache'; cache.mkdir(exist_ok=True)
    settings_path = args.output/'evaluation_settings.json'
    if settings_path.exists() and json.loads(settings_path.read_text()) != settings:
        raise ValueError('Existing self cache has different numerical settings')
    atomic_json(settings_path, settings)
    atomic_json(args.output/'resource_settings.json', dict(workers=args.workers,
        chunk_size=args.chunk_size,cache_entries=args.cache_entries,
        note='Execution settings only; objective, tolerance and cache numerical settings are unchanged'))
    tables, completed = {}, {}
    fields = ['experiment', 'nk', 'refractor_epsilon', 'pf_epsilon', 'eot',
              'source_self', 'target_self', 'divergence', 'raw_divergence',
              'converged', 'cross_iterations', 'source_self_iterations', 'target_self_iterations']
    for mode in ('hard', 'entropic'):
        folder = args.output/mode; folder.mkdir(exist_ok=True)
        metadata = dict(params, pushforward=mode,
                        cost_metric='debiased_sinkhorn_divergence_squared_euclidean_3d',
                        experiments=[1, 2, 3])
        atomic_json(folder/'parameters.json', metadata)
        tables[mode] = read_rows(folder)
        completed[mode] = {row_key(r) for r in tables[mode]}

    def save_tables():
        for mode, rows in tables.items():
            temporary = args.output/mode/'costs.tmp.csv'
            with temporary.open('w', newline='', encoding='utf-8') as stream:
                writer = csv.DictWriter(stream, fieldnames=fields)
                writer.writeheader(); writer.writerows(rows)
            atomic_replace(temporary, args.output/mode/'costs.csv')

    def plot():
        render(args.output/'entropic', args.output/'hard', value_field='divergence',
               ylabel=rf'$S_{{{eta:.6g}}}(PF,T)$ (squared Euclidean, 3D)',
               progress_label='debiased entropic cases')

    def self_term(points, weights):
        identity = measure_hash(points, weights)
        path = cache/f'{identity}.json'
        if path.exists():
            payload = json.loads(path.read_text())
            if payload['measure_sha256'] != identity or payload['settings'] != settings:
                raise ValueError('Self checkpoint metadata mismatch')
            value, info = payload['value'], payload['info']
        else:
            print(f'Self EOT: {len(points)} points, {identity[:12]}', flush=True)
            value, info = entropic_ot_cost(
                points, points, weights, weights, eta,
                chunk_size=args.chunk_size, max_iter=max_iter, tolerance=tolerance,
                return_info=True, max_cache_entries=args.cache_entries)
            if info['converged'] and np.isfinite(value):
                atomic_json(path, dict(measure_sha256=identity, settings=settings, value=value, info=info))
        if not info['converged'] or not np.isfinite(value):
            raise RuntimeError(f'Self EOT failed for {identity}: {info}')
        return value, info

    def evaluate(mode, source_row):
        folder = args.hard if mode == 'hard' else args.entropic
        experiment, nk, epsilon = row_key(source_row)
        if source_row['converged'] != 'True' or float(source_row['pf_epsilon']) != eta:
            raise ValueError('Invalid saved cross evaluation')
        matches = list(folder.glob(f'solve_NK{nk}_*eps{epsilon:.17g}.npz'))
        if len(matches) != 1:
            raise ValueError(f'Missing or ambiguous saved map: {mode}, NK={nk}, epsilon={epsilon}')
        solve = matches[0]
        suffix = solve.stem.removeprefix('solve_')
        cross_path = folder/f'pf_{suffix}_eval{eta:.17g}.json'
        cross = json.loads(cross_path.read_text())
        if not cross[1]['converged'] or not np.isclose(cross[0], float(source_row['cost']), atol=1e-15, rtol=0):
            raise ValueError('Saved cross checkpoint does not match source CSV')
        with np.load(solve) as data:
            points, weights = data['y_push_3d'], data['p']
            points_weights = weights[weights > 0]
            y, q = data['y'], data['q']
            if len(points) != len(points_weights) or len(y) != params['target_nk']:
                raise ValueError('Saved point/weight lengths do not match')
            source_self = self_term(points, points_weights)
            target_self = self_term(y[q > 0], q[q > 0])
        value, info = sinkhorn_divergence_from_terms(cross, source_self, target_self)
        if not info['all_converged'] or info['raw_divergence'] < -1e-7:
            raise RuntimeError(f'Invalid debiased evaluation: {mode} {nk} {info}')
        result = dict(experiment=experiment, nk=nk, refractor_epsilon=epsilon, pf_epsilon=eta,
                      eot=cross[0], source_self=source_self[0], target_self=target_self[0],
                      divergence=value, raw_divergence=info['raw_divergence'], converged=True,
                      **{k:info[k] for k in ('cross_iterations', 'source_self_iterations', 'target_self_iterations')})
        atomic_json(args.output/mode/f'divergence_NK{nk}_eps{epsilon:.17g}.json',
                    dict(result=result, diagnostics=info, cross_checkpoint=str(cross_path),
                         cross_checkpoint_sha256=hashlib.sha256(cross_path.read_bytes()).hexdigest(),
                         saved_map=str(solve), saved_map_sha256=hashlib.sha256(solve.read_bytes()).hexdigest()))
        tables[mode].append(result); completed[mode].add(row_key(result))
        print(f'{mode} experiment={experiment}, NK={nk}, S={value:.10g}', flush=True)

    save_tables(); plot()
    deadline = time.monotonic()+7*24*3600
    expected = sum(len(params['epsilon'] if e == 3 else params['nk']) for e in (1,2,3))
    while True:
        hard_rows = [r for r in read_rows(args.hard) if int(r['experiment']) in (1, 2, 3)]
        soft_rows = [r for r in read_rows(args.entropic) if int(r['experiment']) in (1, 2, 3)]
        hard_by_key = {row_key(r):r for r in hard_rows}
        # Compute paired points first, then remaining available hard-only cases.
        for soft in soft_rows:
            key = row_key(soft)
            changed = False
            for mode, row in (('hard', hard_by_key[key]), ('entropic', soft)):
                if key not in completed[mode]:
                    evaluate(mode, row)
                    changed = True
            if changed:
                save_tables(); plot()
                atomic_json(args.output/'progress.json', dict(
                    hard_completed=len(completed['hard']), entropic_completed=len(completed['entropic']),
                    entropic_available=len(soft_rows), expected_per_map=expected))
        for row in hard_rows:
            if row_key(row) not in completed['hard']:
                evaluate('hard', row); save_tables(); plot()
        if len(completed['entropic']) == expected and len(completed['hard']) == expected:
            print('Completed debiased comparison', flush=True); return
        if not args.watch:
            return
        if time.monotonic() > deadline:
            raise TimeoutError('Debiased evaluation watch exceeded seven days; saved results preserved')
        time.sleep(30)


if __name__ == '__main__':
    main()
