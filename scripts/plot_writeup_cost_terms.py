"""Plot EOT terms and certified hard-map OT bounds from saved write-up captures."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from threadpoolctl import threadpool_limits

from result_paths import artifact_path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT/'geometric_optics_progress'
OUTPUT = DOC/'metadata/cost_terms'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hard_ot_bounds(points, pushed_mass, target_mass, costs=None):
    """Keep common mass on the diagonal, then couple excess to deficits.

    Bounds are on unregularized squared-distance OT, not cross EOT. The
    product residual coupling is a feasible witness, not an optimal solve.
    """
    y = np.asarray(points, dtype=float)
    m, q = np.asarray(pushed_mass, dtype=float), np.asarray(target_mass, dtype=float)
    if (y.ndim != 2 or m.shape != (len(y),) or q.shape != m.shape
            or not np.isfinite(y).all() or not np.isfinite(m).all()
            or not np.isfinite(q).all() or (m < 0).any() or (q < 0).any()
            or not np.isclose(m.sum(), 1, rtol=0, atol=1e-12)
            or not np.isclose(q.sum(), 1, rtol=0, atol=1e-12)):
        raise ValueError('Expected finite points and nonnegative probability masses')
    if costs is None:
        costs = np.sum((y[:, None, :]-y[None, :, :])**2, axis=2)
    costs = np.asarray(costs)
    if costs.shape != (len(y), len(y)) or not np.isfinite(costs).all() or (costs < 0).any():
        raise ValueError('Invalid cost matrix')
    if np.max(np.abs(np.diag(costs))) > 1e-15:
        raise ValueError('Diagonal costs must vanish for retained mass')
    stay = np.minimum(m, q)
    excess, deficit = m-stay, q-stay
    moved = float(excess.sum())
    if moved == 0:
        return dict(excess_mass=0., retained_mass=float(stay.sum()), worst_all_bound=0.,
                    worst_deficit_bound=0., feasible_plan_cost=0., witness_marginal_linf=float(np.max(abs(m-q))))
    if not np.isclose(moved, deficit.sum(), rtol=0, atol=2e-12) or not (deficit > 0).any():
        raise ValueError('Unbalanced residual mass')
    # pi_ij = min(m_i,q_i) delta_ij + excess_i deficit_j / moved.
    row_mass = stay + excess*(deficit.sum()/moved)
    col_mass = stay + deficit*(excess.sum()/moved)
    residual = float(max(np.max(abs(row_mass-m)), np.max(abs(col_mass-q))))
    witness = float(excess @ costs @ deficit / moved)
    restricted = float(excess @ costs[:, deficit > 0].max(axis=1))
    worst = float(excess @ costs.max(axis=1))
    if residual > 2e-12 or not (-1e-15 <= witness <= restricted+1e-14 <= worst+2e-14):
        raise ValueError('Feasible coupling / bound verification failed')
    return dict(excess_mass=moved, retained_mass=float(stay.sum()), worst_all_bound=worst,
                worst_deficit_bound=restricted, feasible_plan_cost=witness,
                witness_marginal_linf=residual)


def render(filename, rows, lines, ylabel):
    with plt.rc_context({'font.size': 14, 'axes.titlesize': 16,
                         'axes.labelsize': 15, 'legend.fontsize': 10.5}):
        fig, axes = plt.subplots(3, 1, figsize=(7.3, 9.4), layout='constrained')
        for experiment, ax in enumerate(axes, 1):
            field = 'refractor_epsilon' if experiment == 3 else 'nk'
            selected = sorted((r for r in rows if int(r['experiment']) == experiment),
                              key=lambda r: float(r[field]))
            titles = ['1. Joint movement', '2. Fixed refractor epsilon = 1/320',
                      '3. Fixed source count = 1600']
            ax.set_title(titles[experiment-1]+f' ({len(selected)}/40)')
            for value, label, style in lines:
                if selected:
                    ax.plot([float(r[field]) for r in selected],
                            [float(r[value]) for r in selected], style,
                            label=label, markersize=3, linewidth=1.7)
            ax.set_xscale('log')
            ax.set_xlabel(r'Refractor $\varepsilon$' if experiment == 3 else r'Source count $N_K$')
            ax.set_ylabel(ylabel)
            ax.grid(alpha=.25)
            if selected:
                ax.legend(loc='best', framealpha=.9, ncol=2 if len(lines) > 3 else 1)
            else:
                ax.set_xticks([]); ax.set_yticks([])
                ax.text(.5, .5, 'Pending: no completed cases', transform=ax.transAxes,
                        ha='center', va='center', fontsize=13)
        fig.savefig(DOC/'images'/filename, dpi=220)
        plt.close(fig)


def read_csv(path):
    return list(csv.DictReader(path.open(encoding='utf-8', newline='')))


def write_csv(name, rows):
    with (OUTPUT/name).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    OUTPUT.mkdir(exist_ok=True)
    original = json.loads((DOC/'metadata/results_capture.json').read_text(encoding='utf-8'))
    repaired = json.loads((DOC/'metadata/repaired_kappa_capture/capture.json').read_text(encoding='utf-8'))
    repaired_records = {(int(p['case'][0]), int(p['case'][1]), float(p['case'][2])): ROOT/p['record']
                        for p in repaired['provenance'] if p['mode'] == 'hard'}
    manifest = dict(cost='squared Euclidean distance in 3D', original_capture=original['captured_at'],
                    repaired_capture=repaired['captured_at'], sources=[], maps=[], figures=[])
    matrix_cache = {}
    for run in ('original', 'repaired_kappa'):
        base = ROOT/original['comparison'] if run == 'original' else DOC/'metadata/repaired_kappa_capture'
        for mode in ('hard', 'entropic'):
            path = artifact_path(base/mode, 'costs.csv') if run == 'original' else base/f'{mode}.csv'
            expected = original['csv_sha256'][mode] if run == 'original' else repaired['data_sha256'][path.name]
            assert sha(path) == expected
            rows = read_csv(path)
            manifest['sources'].append(dict(path=path.relative_to(ROOT).as_posix(), sha256=sha(path)))
            for row in rows:
                a, b, cross = (float(row[k]) for k in ('source_self', 'target_self', 'eot'))
                assert row['converged'] == 'True' and np.isfinite([a,b,cross]).all()
                assert abs(cross-.5*(a+b)-float(row['raw_divergence'])) < 1e-15
                row['self_sum'], row['half_self_sum'] = a+b, .5*(a+b)
            filename = f'{run}_{mode}_cost_terms.png'
            render(filename, rows, [('eot', 'Raw cross EOT', 'o-'),
                ('source_self', 'Push-forward self EOT', 's-'),
                ('target_self', 'Target self EOT', '^-'),
                ('self_sum', 'Sum of self EOTs', 'd-'),
                ('half_self_sum', 'Half-sum (debiasing term)', '--')], 'Regularized cost')
            write_csv(f'{run}_{mode}_cost_terms.csv', rows)
            manifest['figures'].append(filename)
            if mode != 'hard':
                continue
            bounds = []
            for row in rows:
                nk, eps = int(row['nk']), float(row['refractor_epsilon'])
                record_path = (artifact_path(base/mode, f'divergence_NK{nk}_eps{eps:.17g}.json')
                               if run == 'original' else repaired_records[(int(row['experiment']),nk,eps)])
                record = json.loads(record_path.read_text(encoding='utf-8'))
                path = ROOT/record['saved_map']
                assert sha(path) == record['saved_map_sha256']
                with np.load(path, allow_pickle=False) as data:
                    y, rays, p, q = (data[k] for k in ('y','y_push_3d','p','q'))
                    lookup = {tuple(v): i for i,v in enumerate(y)}
                    assert len(lookup) == len(y) and len(rays) == len(p) == nk
                    assigned = [lookup[tuple(v)] for v in rays]
                    mass = np.bincount(assigned, weights=p, minlength=len(y))
                    identity = hashlib.sha256(y.tobytes()).hexdigest()
                    if identity not in matrix_cache:
                        matrix_cache[identity] = np.sum((y[:,None,:]-y[None,:,:])**2, axis=2)
                    result = hard_ot_bounds(y,mass,q,matrix_cache[identity])
                bounds.append(dict(experiment=row['experiment'], nk=nk, refractor_epsilon=eps,
                                   raw_eot=float(row['eot']), **result))
                manifest['maps'].append(dict(run=run, case=[int(row['experiment']),nk,eps],
                    record=record_path.relative_to(ROOT).as_posix(), record_sha256=sha(record_path),
                    map=path.relative_to(ROOT).as_posix(), map_sha256=record['saved_map_sha256']))
            write_csv(f'{run}_hard_ot_bounds.csv', bounds)
            filename = f'{run}_hard_ot_bounds.png'
            render(filename, bounds, [('worst_all_bound', 'Worst cost over all targets', 'o-'),
                ('worst_deficit_bound', 'Worst cost over deficit targets', 's-'),
                ('feasible_plan_cost', 'Explicit feasible-plan cost', '^-')], 'Upper bound on OT cost')
            manifest['figures'].append(filename)
    tex = r'''\clearpage
\subsection{EOT components and hard-map OT upper bounds}
\label{sec:cost-terms-and-bounds}
For each run, map type, and parameter experiment, we plot the raw cross-EOT,
$A=\mathrm{EOT}_{1/320}(T_\#p,T_\#p)$, $B=\mathrm{EOT}_{1/320}(q,q)$,
their full sum $A+B$, and the debiasing term $(A+B)/2$.
Thus the divergence plotted earlier is the raw cross-EOT minus the half-sum,
not minus the full sum. ``Push-forward self'' refers to the mapped measure,
not the optical source-self correction. These EOT terms use squared
Euclidean distance and KL relative to the product measure. The component
plots use a linear cost axis and the same completed-case captures as above.

For a hard map, all pushed points are target atoms $y_j$. Aggregate their
source weights into $m_j=(T_\#p)(\{y_j\})$. Put
\[
 s_j=\min(m_j,q_j),\qquad e_j=(m_j-q_j)_+,\qquad
 d_j=(q_j-m_j)_+,\qquad M=\sum_j e_j=\sum_j d_j.
\]
Retain $s_j$ at $y_j$ at zero cost. For $M>0$, the residual coupling
\[
 \pi_{ij}=s_i\delta_{ij}+\frac{e_i d_j}{M}
\]
has row masses $m_i$ and column masses $q_j$. With
$c_{ij}=\|y_i-y_j\|^2$, this gives
\[
 \mathrm{OT}(T_\#p,q)
 \le C_\pi=\frac1M\sum_{i,j}e_i d_jc_{ij}
 \le U_D=\sum_i e_i\max_{j:d_j>0}c_{ij}
 \le U=\sum_i e_i\max_jc_{ij}.
\]
The last expression is the requested per-pushed-point worst-cost bound.
The deficit-only maximum is a tighter valid variant. If $M=0$, all three
upper bounds are zero. We verify the witness marginals and bound ordering
numerically for every plotted case. These bounds apply to unregularized
OT, not to regularized EOT; EOT need not lie below these bounds.
Keeping common mass on target atoms at zero cost is specific here to the
hard maps; we do not apply this construction to the smooth Snell maps.
'''
    for run, kappa in [('original','0.6'),('repaired_kappa','0.8')]:
        for mode in ('hard','entropic'):
            tex += '\n'+rf'''\begin{{figure}}[p]
 \centering
 \includegraphics[width=\textwidth]{{{run}_{mode}_cost_terms.png}}
 \caption{{\large $\kappa={kappa}$, {mode} map: raw EOT and self-EOT
 components for all three experiments. The full self sum and its half-sum
 are shown separately. Panel titles give completed cases out of 40.}}
 \label{{fig:{run}-{mode}-cost-terms}}
\end{{figure}}
'''
        tex += '\n'+rf'''\begin{{figure}}[p]
 \centering
 \includegraphics[width=\textwidth]{{{run}_hard_ot_bounds.png}}
 \caption{{\large $\kappa={kappa}$ hard map: worst-cost OT upper bounds
 after retaining common mass, together with the cost of an explicit
 feasible residual coupling. All three curves bound unregularized OT.
 Pending cases are omitted.}}
 \label{{fig:{run}-hard-ot-bounds}}
\end{{figure}}
'''
    (DOC/'cost_terms_results.tex').write_text(tex+'\n\\clearpage\n',encoding='utf-8')
    include = DOC/'new_developments.tex'
    text = include.read_text(encoding='utf-8')
    if '\\input{cost_terms_results}' not in text:
        include.write_text(text+'\n\\input{cost_terms_results}\n', encoding='utf-8')
    manifest['outputs_sha256'] = {p.name:sha(p) for p in OUTPUT.glob('*.csv')}
    manifest['figure_sha256'] = {name:sha(DOC/'images'/name) for name in manifest['figures']}
    manifest['max_witness_marginal_linf'] = max(float(r['witness_marginal_linf'])
        for p in OUTPUT.glob('*_hard_ot_bounds.csv') for r in read_csv(p))
    (OUTPUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(figures=manifest['figures'],
                         max_witness_marginal_linf=manifest['max_witness_marginal_linf'])))


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
