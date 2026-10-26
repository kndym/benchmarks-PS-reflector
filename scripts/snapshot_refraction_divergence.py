"""Validate and snapshot saved debiased evaluations without recomputing OT."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
from plot_refraction_pushforward_comparison import read_rows, render


def snapshot(source, destination):
    if destination.exists():
        raise FileExistsError(destination)
    tables = {mode: read_rows(source / mode) for mode in ('hard', 'entropic')}
    selected = {source / 'evaluation_settings.json'}
    for mode, rows in tables.items():
        folder = source / mode
        selected.add(folder / 'parameters.json')
        keys = set()
        for row in rows:
            nk, eps = int(row['nk']), float(row['refractor_epsilon'])
            key = int(row['experiment']), nk, eps
            if key in keys or row['converged'] != 'True':
                raise ValueError(f'Duplicate or unconverged row: {mode} {key}')
            keys.add(key)
            path = folder / f'divergence_NK{nk}_eps{eps:.17g}.json'
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
                if hashlib.sha256(Path(record[name]).read_bytes()).hexdigest() != record[hash_key]:
                    raise ValueError(f'Changed source checkpoint: {path}')
            selected.add(path)
    destination.mkdir(parents=True)
    for path in sorted(selected):
        target = destination / path.relative_to(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    for mode, rows in tables.items():
        with (destination/mode/'costs.csv').open('w', newline='', encoding='utf-8') as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
    render(destination / 'entropic', destination / 'hard', value_field='divergence',
           ylabel=r'$S_{1/320}(PF,T)$ (squared Euclidean, 3D)',
           progress_label='debiased entropic cases')
    params = json.loads((source/'entropic/parameters.json').read_text())
    expected = 2*len(params['nk'])+len(params['epsilon'])
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
    report['files_sha256'] = {str(p.relative_to(destination).as_posix()):
                              hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in sorted(destination.rglob('*')) if p.is_file()}
    (destination/'snapshot.json').write_text(json.dumps(report, indent=2)+'\n')
    (destination/'README.md').write_text(
        '# Debiased hard versus entropic comparison\n\n'
        f"Hard: **{len(tables['hard'])}/{expected}**. Entropic: **{len(tables['entropic'])}/{expected}**. "
        'Partial snapshot; missing entropic cases are still being computed.\n\n'
        'All three EOT terms converged for every published row. Evaluation uses squared '
        'Euclidean distance in 3D, epsilon 1/320, and the KL-relative-to-product convention. '
        'The graph shows Sinkhorn divergence itself, without a square root. Optical potential '
        'stopping and marginal feasibility are separate; see the optical snapshot residuals.\n\n'
        '![Hard versus entropic debiased comparison](entropic/comparison_experiments_1_2_3.png)\n\n'
        '![Log-log debiased comparison](entropic/comparison_experiments_1_2_3_log.png)\n\n'
        '[Hard data](hard/costs.csv) | [Entropic data](entropic/costs.csv) | '
        '[Validation and counts](snapshot.json)\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k != 'files_sha256'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    snapshot(args.source, args.destination)
