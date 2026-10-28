"""Validate and snapshot saved debiased evaluations without recomputing OT."""
try:
    from result_paths import artifact_path, artifact_glob, resolve_path, dated_run_dir, ensure_readme, self_cost_dir
except ModuleNotFoundError:
    from scripts.result_paths import artifact_path, artifact_glob, resolve_path, dated_run_dir, ensure_readme, self_cost_dir
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
from plot_refraction_pushforward_comparison import read_rows, render, comparison_case_keys


def snapshot(source, destination):
    source, destination = resolve_path(source), resolve_path(destination)
    if destination.exists():
        raise FileExistsError(destination)
    tables = {mode: read_rows(source / mode) for mode in ('hard', 'entropic')}
    selected = {artifact_path(source, 'evaluation_settings.json')}
    for mode, rows in tables.items():
        folder = source / mode
        selected.add(artifact_path(folder, 'parameters.json'))
        keys = set()
        for row in rows:
            nk, eps = int(row['nk']), float(row['refractor_epsilon'])
            key = int(row['experiment']), nk, eps
            if key in keys or row['converged'] != 'True':
                raise ValueError(f'Duplicate or unconverged row: {mode} {key}')
            keys.add(key)
            path = artifact_path(folder, f'divergence_NK{nk}_eps{eps:.17g}.json')
            record = json.loads(path.read_text())
            diagnostics = record['diagnostics']
            if not diagnostics['all_converged']:
                raise ValueError(f'Unconverged terms: {path}')
            for field in ('eot', 'source_self', 'target_self', 'divergence', 'raw_divergence'):
                value = float(row[field])
                if not np.isfinite(value) or value != record['result'][field]:
                    raise ValueError(f'Invalid/mismatched value: {path} {field}')
            raw = float(row['eot']) - .5*(float(row['source_self'])+float(row['target_self']))
            if abs(raw-float(row['raw_divergence'])) > 1e-15 or raw < -1e-7:
                raise ValueError(f'Invalid divergence formula: {path}')
            if float(row['divergence']) != max(0., raw):
                raise ValueError(f'Invalid plotted divergence: {path}')
            for name, hash_key in (('saved_map', 'saved_map_sha256'),
                                   ('cross_checkpoint', 'cross_checkpoint_sha256')):
                if hashlib.sha256(resolve_path(record[name]).read_bytes()).hexdigest() != record[hash_key]:
                    raise ValueError(f'Changed source checkpoint: {path}')
            selected.add(path)
    destination.mkdir(parents=True)
    for path in sorted(selected):
        target = artifact_path(destination / path.relative_to(source).parts[0], path.name, for_write=True) if path.relative_to(source).parts[0] in ('hard', 'entropic') else artifact_path(destination, path.name, for_write=True)
        shutil.copy2(path, target)
    for mode, rows in tables.items():
        with artifact_path(destination/mode, 'costs.csv', for_write=True).open('w', newline='', encoding='utf-8') as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
    render(destination / 'entropic', destination / 'hard', value_field='divergence',
           ylabel=r'$S_{1/320}(PF,T)$ (squared Euclidean, 3D)',
           progress_label='debiased entropic cases')
    params = json.loads(artifact_path(source/'entropic', 'parameters.json').read_text())
    expected = len(comparison_case_keys(params))
    report = dict(captured_at=datetime.now(timezone.utc).isoformat(),
                  expected_per_map=expected, all_terms_converged=True,
                  formula='EOT(PF,T) - 0.5 * (EOT(PF,PF) + EOT(T,T))',
                  counts={mode: {str(e): sum(int(r['experiment']) == e for r in rows)
                                 for e in (1, 2, 3)} for mode, rows in tables.items()},
                  completed={mode: len(rows) for mode, rows in tables.items()})
    hard = {(r['experiment'], r['nk'], r['refractor_epsilon']):r for r in tables['hard']}
    pairs = [(hard[(r['experiment'], r['nk'], r['refractor_epsilon'])], r)
             for r in tables['entropic']]
    report['entropic_lower_divergence_pairs'] = sum(float(s['divergence']) < float(h['divergence'])
                                                      for h, s in pairs)
    report['paired_cases'] = len(pairs)
    report['hash_base'] = '..'
    report['files_sha256'] = {str(p.relative_to(destination).as_posix()):
                              hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in sorted(destination.rglob('*')) if p.is_file() and p.name != 'README.md'}
    artifact_path(destination, 'snapshot.json', for_write=True).write_text(json.dumps(report, indent=2)+'\n')
    (destination/'README.md').write_text(
        '# Debiased hard versus entropic comparison\n\n'
        f"Hard: **{len(tables['hard'])}/{expected}**. Entropic: **{len(tables['entropic'])}/{expected}**. " +
        ('Complete snapshot; both maps cover every experiment grid.\n\n'
         if all(len(rows) == expected for rows in tables.values()) else
         'Partial snapshot; missing entropic cases are still being computed.\n\n') +
        'All three EOT terms converged for every published row. Evaluation uses squared '
        'Euclidean distance in 3D, epsilon 1/320, and the KL-relative-to-product convention. '
        'The graph shows Sinkhorn divergence itself, without a square root. Optical potential '
        'stopping and marginal feasibility are separate; see the optical snapshot residuals.\n\n'
        '![Hard versus entropic debiased comparison](entropic/comparison_experiments_1_2_3.png)\n\n'
        '![Log-log debiased comparison](entropic/comparison_experiments_1_2_3_log.png)\n\n'
        '[Hard data](hard/costs.csv) | [Entropic data](entropic/costs.csv) | '
        '[Validation and counts](metadata/snapshot.json)\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k != 'files_sha256'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    snapshot(args.source, args.destination)
