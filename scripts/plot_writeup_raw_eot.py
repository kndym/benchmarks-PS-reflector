"""Add raw cross-EOT figures for both existing write-up captures."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from plot_writeup_results import render as render_original
from update_repaired_kappa_writeup import render as render_repaired
from result_paths import artifact_path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / 'geometric_optics_progress'


def main():
    original = json.loads((DOC/'metadata/results_capture.json').read_text(encoding='utf-8'))
    repaired = json.loads((DOC/'metadata/repaired_kappa_capture/capture.json').read_text(encoding='utf-8'))
    sources = {}
    for name, renderer, base in [
        ('original', render_original, ROOT/original['comparison']),
        ('repaired_kappa', render_repaired, DOC/'metadata/repaired_kappa_capture')]:
        tables = {}
        for mode in ('hard', 'entropic'):
            path = (artifact_path(base/mode, 'costs.csv') if name == 'original'
                    else base/f'{mode}.csv')
            sha = hashlib.sha256(path.read_bytes()).hexdigest()
            expected = (original['csv_sha256'][mode] if name == 'original'
                        else repaired['data_sha256'][f'{mode}.csv'])
            assert sha == expected, f'Capture hash mismatch: {path}'
            rows = list(csv.DictReader(path.open(encoding='utf-8', newline='')))
            assert all(r['converged'] == 'True' for r in rows)
            for r in rows:
                eot = float(r['eot'])
                assert np.isfinite(eot) and eot > 0
                raw = eot - .5*(float(r['source_self'])+float(r['target_self']))
                assert abs(raw-float(r['raw_divergence'])) < 1e-15
            tables[mode] = rows
            sources[f'{name}/{mode}'] = dict(path=path.relative_to(ROOT).as_posix(), sha256=sha, cases=len(rows))
        renderer(f'{name}_raw_eot_comparison.png', [
            (tables['hard'], 'eot', 'Hard map', 'o-'),
            (tables['entropic'], 'eot', 'Entropic map + Snell', 's-')],
            'Raw EOT cost', log_y=True)
    text = r'''\clearpage
\subsection{Comparison without debiasing: raw EOT}
\label{sec:raw-eot-comparison}
For both runs we also plot the uncorrected cross-EOT objective
$\mathrm{EOT}_{1/320}(T_\#p,q)$, using squared Euclidean distance in 3D
and the KL-relative-to-product convention. No source-self or target-self
term is subtracted. Thus these are raw regularized costs; their ordering
can differ from the debiased Sinkhorn divergences. The hard and entropic
curves use exactly the corresponding captured cases in the figures above.
Hard-map mass error is a separate diagnostic and has no debiasing step.

\begin{figure}[p]
 \centering
 \includegraphics[width=\textwidth]{original_raw_eot_comparison.png}
 \caption{\large Original $\kappa=0.6$ run without debiasing:
 raw cross-EOT at evaluation $\varepsilon=1/320$ for the hard and entropic
 maps, with logarithmic axes. All 120 cases per map are included, using
 the same completed October 26 capture as the original divergence figure.}
 \label{fig:original-raw-eot}
\end{figure}

\begin{figure}[p]
 \centering
 \includegraphics[width=\textwidth]{repaired_kappa_raw_eot_comparison.png}
 \caption{\large Repaired $\kappa=0.8$ run without debiasing:
 raw cross-EOT at evaluation $\varepsilon=1/320$ for the hard and entropic
 maps, with logarithmic axes. This uses the same TOTAL paired cases as
 the repaired divergence figure: COUNTS in experiments 1, 2, and 3.
 Pending cases are omitted.}
 \label{fig:repaired-kappa-raw-eot}
\end{figure}
\clearpage
'''
    text = text.replace('TOTAL', str(repaired['paired_completed'])).replace(
        'COUNTS', ', '.join(map(str, repaired['counts'])))
    (DOC/'raw_eot_results.tex').write_text(text, encoding='utf-8')
    include = DOC/'new_developments.tex'
    text = include.read_text(encoding='utf-8')
    if '\\input{raw_eot_results}' not in text:
        include.write_text(text+'\n\\input{raw_eot_results}\n', encoding='utf-8')
    (DOC/'metadata/raw_eot_capture.json').write_text(json.dumps(dict(
        metric='Raw cross EOT; no self-term subtraction', evaluation_epsilon=1/320,
        sources=sources), indent=2)+'\n', encoding='utf-8')
    print(json.dumps(sources, indent=2))


if __name__ == '__main__':
    main()
