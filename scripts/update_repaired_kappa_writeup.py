"""Capture completed repaired-kappa cases and refresh the two write-up figures."""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from result_paths import artifact_path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / 'geometric_optics_progress'
RUN = ROOT / 'results/runs/repaired_refractor_kappa_oct_28_26'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def key(row):
    return int(row['experiment']), int(row['nk']), float(row['refractor_epsilon'])


def render(name, series, ylabel, log_y=False):
    with plt.rc_context({'font.size': 16, 'axes.titlesize': 17,
                         'axes.labelsize': 16, 'xtick.labelsize': 14,
                         'ytick.labelsize': 14, 'legend.fontsize': 13}):
        fig, axes = plt.subplots(3, 1, figsize=(6.4, 8.7), layout='constrained')
        for experiment, ax in enumerate(axes, 1):
            field = 'refractor_epsilon' if experiment == 3 else 'nk'
            counts = []
            for rows, value, label, style in series:
                selected = sorted((r for r in rows if key(r)[0] == experiment),
                                  key=lambda r: float(r[field]))
                counts.append(len(selected))
                if selected:
                    ax.plot([float(r[field]) for r in selected],
                            [float(r[value]) for r in selected], style,
                            label=label, markersize=3.5, linewidth=1.8)
            ax.set_xscale('log')
            if log_y:
                ax.set_yscale('log')
            titles = ['1. Joint movement', '2. Fixed refractor epsilon = 1/320',
                      '3. Fixed source count = 1600']
            ax.set_title(titles[experiment-1] + f' ({min(counts)}/40)', pad=8)
            ax.set_xlabel(r'Refractor $\varepsilon$' if experiment == 3 else r'Source count $N_K$')
            ax.set_ylabel(ylabel)
            ax.grid(alpha=.25)
            if not any(counts):
                ax.text(.5, .5, 'Pending: no completed cases', ha='center', va='center',
                        transform=ax.transAxes, fontsize=13,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=.9))
            elif len(series) > 1:
                ax.legend(loc='best')
        fig.savefig(DOC / 'images' / name, dpi=240)
        plt.close(fig)


def main():
    merged = {m: {} for m in ('hard', 'entropic')}
    provenance = []
    reference = None
    fields = ['kappa', 'source_theta', 'source_phi', 'target_theta', 'target_phi',
              'fixed_nk', 'target_nk', 'fixed_epsilon', 'cost_epsilon',
              'solver_version', 'self_cost_mode', 'optical_marginal_tolerance']
    for group in ('kappa_only', 'parallel_epsilon/kappa_only'):
        for mode in merged:
            folder = RUN / group / mode
            params = json.loads(artifact_path(folder, 'parameters.json').read_text(encoding='utf-8'))
            settings = {k: params[k] for k in fields}
            if reference is None:
                reference = settings
            assert settings == reference, (group, mode, 'incompatible physics or tolerances')
            table = artifact_path(RUN / group / 'debiased' / mode, 'costs.csv')
            rows = list(csv.DictReader(table.open(encoding='utf-8', newline='')))
            for row in rows:
                case = key(row)
                record_path = artifact_path(table.parent.parent,
                    f'divergence_NK{case[1]}_eps{case[2]:.17g}.json')
                record = json.loads(record_path.read_text(encoding='utf-8'))
                assert row['converged'] == 'True' and record['diagnostics']['all_converged']
                for field in ('eot', 'source_self', 'target_self', 'divergence', 'raw_divergence'):
                    assert np.isfinite(float(row[field])) and float(row[field]) == record['result'][field]
                raw = float(row['eot']) - .5*(float(row['source_self'])+float(row['target_self']))
                assert abs(raw-float(row['raw_divergence'])) < 1e-15
                assert float(row['divergence']) == max(0, raw) and raw >= -1e-7
                for name, hashfield in [('saved_map', 'saved_map_sha256'),
                                         ('cross_checkpoint', 'cross_checkpoint_sha256')]:
                    assert digest(Path(record[name])) == record[hashfield]
                if case in merged[mode]:
                    old = merged[mode][case][0]
                    assert all(np.isclose(float(old[f]), float(row[f]), rtol=1e-7, atol=1e-12)
                               for f in ('eot', 'source_self', 'target_self', 'divergence'))
                    continue
                merged[mode][case] = (row, record)
                provenance.append(dict(group=group, mode=mode, case=case,
                                       record=str(record_path.relative_to(ROOT)),
                                       record_sha256=digest(record_path),
                                       saved_map_sha256=record['saved_map_sha256']))
    common = set(merged['hard']) & set(merged['entropic'])
    tables = {m: [merged[m][k][0] for k in sorted(common)] for m in merged}
    mass = []
    for case in sorted(common):
        row, record = merged['hard'][case]
        with np.load(record['saved_map']) as data:
            y, rays, p, q = (data[n] for n in ('y', 'y_push_3d', 'p', 'q'))
            assert len(p) == case[1] and len(q) == 1600
            assert np.isclose(p.sum(), 1) and np.isclose(q.sum(), 1)
            lookup = {tuple(v): i for i, v in enumerate(y)}
            assert len(lookup) == len(y)
            assigned = np.array([lookup[tuple(v)] for v in rays])
            pushed = np.bincount(assigned, weights=p, minlength=len(y))
            mass.append(dict(experiment=case[0], nk=case[1], refractor_epsilon=case[2],
                             max_abs_mass_difference=float(np.max(np.abs(pushed-q)))))
    output = DOC / 'metadata/repaired_kappa_capture'
    output.mkdir(exist_ok=True)
    for name, rows in [*tables.items(), ('hard_mass_error', mass)]:
        with (output/f'{name}.csv').open('w', encoding='utf-8', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)
    counts = [sum(k[0] == e for k in common) for e in (1, 2, 3)]
    now = datetime.now(timezone.utc)
    report = dict(captured_at=now.isoformat(), settings=reference, counts=counts,
                  paired_completed=len(common), expected_per_map=120, provenance=provenance)
    report['data_sha256'] = {p.name: digest(p) for p in output.glob('*.csv')}
    (output/'capture.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    render('repaired_kappa_sinkhorn_comparison.png', [
        (tables['hard'], 'divergence', 'Hard map', 'o-'),
        (tables['entropic'], 'divergence', 'Entropic map + Snell', 's-')],
        'Sinkhorn divergence', log_y=True)
    render('repaired_kappa_hard_mass_error.png', [
        (mass, 'max_abs_mass_difference', 'Hard map', 'o-')], r'Mass error $E_\infty$')
    tex = r'''\clearpage
\subsection{Repaired kappa run: current progress}
\label{sec:repaired-kappa-progress}
The following figures use the repaired solver with $\kappa=0.8$ and the
same source and target geometry as the original $\kappa=0.6$ sweep above.
They are a separate, partial run. Captured CAPTUREDATE (UTC), there are
TOTAL paired completed cases per map out of 120: COUNTS in experiments
1, 2, and 3, respectively (40 planned per experiment).
The independent epsilon worker contributes completed experiment-3 cases;
matching case keys are counted once. Pending cases are omitted.
The optical marginal tolerance is $10^{-7}$. Each plotted divergence has
three converged independent EOT terms, evaluated at $\varepsilon=1/320$.
The target contains 1600 atoms. These partial curves do not establish a
convergence rate.

\begin{figure}[p]
 \centering
 \includegraphics[width=\textwidth]{repaired_kappa_sinkhorn_comparison.png}
 \caption{\large Repaired $\kappa=0.8$ run: debiased Sinkhorn divergence,
 with logarithmic axes and the same three-panel layout as the original run.
 Panel titles show paired completed cases out of 40. Experiment 2 is pending.}
 \label{fig:repaired-kappa-divergence}
\end{figure}

\begin{figure}[p]
 \centering
 \includegraphics[width=\textwidth]{repaired_kappa_hard_mass_error.png}
 \caption{\large Repaired $\kappa=0.8$ run: maximum absolute hard-map
 target-mass error $E_\infty=\max_j|(T_\#p)_j-q_j|$, computed directly
 from the saved hard-ray assignments and source weights. All target atoms,
 including empty locations, are included. No original-run mass errors are reused.}
 \label{fig:repaired-kappa-mass}
\end{figure}
\clearpage
'''
    tex = tex.replace('CAPTUREDATE', now.strftime('%B %d, %Y at %H:%M')).replace('TOTAL', str(len(common)))
    tex = tex.replace('COUNTS', ', '.join(map(str, counts)))
    if counts[1]:
        tex = tex.replace('Experiment 2 is pending.', 'Unfinished cases are omitted.')
    (DOC/'repaired_kappa_results.tex').write_text(tex, encoding='utf-8')
    included = DOC/'new_developments.tex'
    text = included.read_text(encoding='utf-8')
    if '\\input{repaired_kappa_results}' not in text:
        included.write_text(text+'\n\\input{repaired_kappa_results}\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k != 'provenance'}, indent=2))


if __name__ == '__main__':
    main()
