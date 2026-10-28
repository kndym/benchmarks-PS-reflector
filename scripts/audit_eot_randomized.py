"""Reproduce the 72 small cached/chunked numerical EOT audit evaluations.

Run ``python -B scripts/audit_eot_randomized.py --output /absolute/audit.json``.
Nonconvergence at the iteration cap is recorded, not counted as a success.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from refracter.sinkhorn import entropic_ot_cost


def run_audit():
    seed = 628
    rng = np.random.default_rng(seed)
    records = []
    converged = 0
    for trial in range(36):
        n, m = rng.integers(1, 7, size=2)
        x = rng.normal(size=(n, 2)) * .2
        y = rng.normal(size=(m, 2)) * .2
        p = np.exp(rng.uniform(-15, 0, size=n))
        q = np.exp(rng.uniform(-15, 0, size=m))
        p, q = p/p.sum(), q/q.sum()
        epsilon = (.0001, .01, .2)[trial % 3]
        results = []
        costs = np.sum((x[:, None, :] - y[None, :, :])**2, axis=2)
        for cache in (0, 1000):
            value, info = entropic_ot_cost(x, y, p, q, epsilon, max_iter=500,
                chunk_size=2, max_cache_entries=cache, return_info=True,
                return_potentials=True, tolerance=1e-9, marginal_tolerance=1e-7)
            f, g = info.pop('potentials')
            plan = np.exp(np.log(p[:, None]) + np.log(q[None, :])
                          + (f[:, None] + g[None, :] - costs)/epsilon)
            source_error = float(np.abs(plan.sum(axis=1) - p).sum())
            target_error = float(np.abs(plan.sum(axis=0) - q).sum())
            if not np.isfinite(value):
                raise AssertionError(f'Nonfinite objective: trial={trial}, cache={cache}')
            if info['converged']:
                if max(source_error, target_error) > 1e-7:
                    raise AssertionError(f'False convergence: trial={trial}, cache={cache}')
                converged += 1
            results.append(dict(cache_entries=cache, value=value, diagnostics=info,
                independently_measured_source_l1=source_error,
                independently_measured_target_l1=target_error))
        if all(result['diagnostics']['converged'] for result in results):
            if abs(results[0]['value'] - results[1]['value']) >= 1e-10:
                raise AssertionError(f'Cached/chunked disagreement: trial={trial}')
        records.append(dict(trial=trial, epsilon=epsilon, x=x.tolist(), y=y.tolist(),
                            p=p.tolist(), q=q.tolist(), results=results))
    return dict(schema_version=1, seed=seed, numpy_version=np.__version__,
        configuration=dict(trials=36, evaluations_per_trial=2, point_dimension=2,
            point_count_min=1, point_count_max=6, coordinates='Normal(0, 0.2^2)',
            weights='exp(Uniform(-15, 0)), then normalized separately',
            epsilon_cycle=[.0001, .01, .2], cache_entries=[0, 1000], chunk_size=2,
            max_iter=500, potential_tolerance=1e-9, marginal_tolerance=1e-7,
            objective_agreement_tolerance=1e-10, ground_cost='squared_euclidean'),
        summary=dict(evaluations=72, converged_and_feasible=converged,
            reported_nonconverged=72-converged, false_successes=0,
            cached_chunked_agreement='passed'),
        nonconverged_cases=[dict(trial=case['trial'], epsilon=case['epsilon'],
            cache_entries=result['cache_entries']) for case in records
            for result in case['results'] if not result['diagnostics']['converged']],
        cases=records)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = run_audit()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({key: value for key, value in report.items() if key != 'cases'}, indent=2))
