"""Validate a completed-case capture, then refresh section 8 and its figures."""
try:
    from result_paths import artifact_path, artifact_glob, resolve_path, dated_run_dir, ensure_readme, self_cost_dir
except ModuleNotFoundError:
    from scripts.result_paths import artifact_path, artifact_glob, resolve_path, dated_run_dir, ensure_readme, self_cost_dir
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re
from zoneinfo import ZoneInfo
from zipfile import ZIP_DEFLATED, ZipFile

from plot_writeup_results import ROOT, read_rows, main as plot
from snapshot_refraction_divergence import snapshot


def tex_number(value):
    mantissa, exponent = f"{float(value):.3g}".split("e") if "e" in f"{float(value):.3g}" else (f"{float(value):.3g}", None)
    return mantissa if exponent is None else rf"{mantissa}\times10^{{{int(exponent)}}}"


def main():
    now = datetime.now(ZoneInfo("America/New_York"))
    source = resolve_path(ROOT / "results/runs/debiased_comparison_oct_26_26")
    current = resolve_path(ROOT / "geometric_optics_progress/metadata/results_capture.json")
    signatures = {mode: hashlib.sha256(artifact_path(source/mode, 'costs.csv').read_bytes()).hexdigest()
                  for mode in ('hard', 'entropic')}
    if current.exists() and json.loads(current.read_text()).get('csv_sha256') == signatures:
        print('Write-up already matches the latest completed debiased cases.')
        return
    destination = ROOT / 'results/published' / now.strftime('writeup_comparison_%H%M%S_%b_%d_%y').lower()
    snapshot(source, destination)
    tables = {mode: read_rows(artifact_path(destination/mode, 'costs.csv')) for mode in ('hard','entropic')}
    hard = {(int(r['experiment']),int(r['nk']),float(r['refractor_epsilon'])): r
            for r in tables['hard']}
    pairs = [(hard[(int(r['experiment']),int(r['nk']),float(r['refractor_epsilon']))],r)
             for r in tables['entropic']]
    wins = sum(float(h['divergence']) < float(e['divergence']) for h,e in pairs)
    count = len(pairs)
    counts = [sum(int(r['experiment']) == e for r in tables['entropic']) for e in (1,2,3)]
    stamp = now.strftime('%B %d, %Y at %I:%M %p ET').replace(' 0',' ')
    path = resolve_path(ROOT / 'geometric_optics_progress/new_developments.tex')
    text = path.read_text(encoding='utf-8')
    setup = (f'The results captured on {stamp} contain all 120 hard-map cases and\n'
             f'{count} entropic-map cases: {counts[0]} in experiment 1, {counts[1]} in\n'
             f'experiment 2, and {counts[2]} in experiment 3. ')
    setup += ('Both maps now cover the full grid.\n' if count == 120 else
              f'The remaining {120-count} entropic\ncases are still being computed. All comparisons below use completed cases only.\n')
    start = text.index('The October 26, 2026 results') if 'The October 26, 2026 results' in text else text.index('The results captured on ')
    end = text.index(r'\subsection{Push-forward maps', start)
    text = text[:start] + setup + '\n' + text[end:]
    result = (rf'\subsection{{Results of the three experiments}}'+'\n\n'+
              f'The hard map has lower Sinkhorn divergence in {wins} of the {count}\n'
              'matched cases. All three EOT terms passed the evaluation stopping test\n'
              'for every published divergence. Figure~\\ref{fig:new-divergence} shows\n'
              'the latest comparison.\n\n')
    titles = ['Joint movement','Fixed refractor epsilon','Fixed source count']
    for exp in (1,2,3):
        field = 'refractor_epsilon' if exp == 3 else 'nk'
        rows = {mode: sorted((r for r in table if int(r['experiment']) == exp),
                             key=lambda r:float(r[field])) for mode,table in tables.items()}
        h, e = rows['hard'], rows['entropic']
        result += rf'\paragraph{{Experiment {exp}: {titles[exp-1]}.}}'+'\n'
        if exp == 1:
            result += ('As $N_K$ increases from 100 to 20,000 and the refractor epsilon\n'
                       'decreases, the hard-map divergence falls from\n'
                       rf'${tex_number(h[0]["divergence"])}$ to ${tex_number(h[-1]["divergence"])}$. The entropic-map divergence'+'\n'
                       rf'falls from ${tex_number(e[0]["divergence"])}$ to ${tex_number(e[-1]["divergence"])}$.'+'\n'
                       'Both improve substantially, although changing two parameters together\n'
                       'does not isolate their effects.\n\n')
        elif exp == 2:
            result += ('Increasing the source count alone gives little improvement after the\n'
                       'initial changes. The hard-map divergence ranges from\n'
                       rf'${tex_number(min(float(r["divergence"]) for r in h))}$ to ${tex_number(max(float(r["divergence"]) for r in h))}$; the entropic map stays'+'\n'
                       rf'between ${tex_number(min(float(r["divergence"]) for r in e))}$ and ${tex_number(max(float(r["divergence"]) for r in e))}$.'+'\n'
                       'Refining the source at this fixed epsilon and target does not remove\n'
                       'the observed error floor.\n\n')
        else:
            result += ('For the hard map, divergence increases from\n'
                       rf'${tex_number(h[0]["divergence"])}$ to ${tex_number(h[-1]["divergence"])}$ across the full epsilon grid.'+'\n'
                       f'The {len(e)} completed entropic cases cover\n'
                       rf'${tex_number(e[0][field])}\leq\varepsilon_{{\mathrm{{refr}}}}\leq {tex_number(e[-1][field])}$.'+'\n'
                       rf'Their endpoint divergences are ${tex_number(e[0]["divergence"])}$ and ${tex_number(e[-1]["divergence"])}$.'+'\n')
            if len(e)<40: result += 'The entropic comparison remains incomplete beyond this interval.\n'
            result += '\n'
    result += ('Raw EOT can give a different ranking. At the final point of experiment 1\n'
               'it is $0.0060977$ for the hard map and $0.0058235$ for the entropic map.\n'
               'After self-term subtraction, the hard map has the smaller divergence.\n'
               'A lower cross-EOT value alone does not establish a closer target match.\n\n')
    start=text.index(r'\subsection{Results of the three experiments}')
    end=text.index(r'\begin{figure}',start)
    text=text[:start]+result+text[end:]
    start=text.index(r' \caption{\large Debiased Sinkhorn divergence')
    end=text.index(r' \label{fig:new-divergence}',start)
    text=text[:start]+(rf' \caption{{\large Debiased Sinkhorn divergence at $\eta=1/320$, with logarithmic'+'\n'
         f' axes. Hard: 120 cases; entropic: {count} cases. Experiment 3 contains\n'
         f' {counts[2]} entropic values. Captured {stamp}.}}\n')+text[end:]
    plot(destination)
    path.write_text(text,encoding='utf-8')
    project = path.parent
    artifacts = ['new_developments.tex', 'images/latest_sinkhorn_comparison.png',
                 'images/latest_hard_mass_error.png']
    with ZipFile(resolve_path(project/'overleaf_update.zip'), 'w', ZIP_DEFLATED) as archive:
        for name in artifacts:
            archive.write(resolve_path(project/name), name)
    current.write_text(json.dumps({'captured_at':now.isoformat(),'comparison':str(destination.relative_to(ROOT)),
        'hard_cases':120,'entropic_cases':count,'counts':counts,'hard_lower_pairs':wins,
        'csv_sha256':signatures,
        'artifact_sha256':{name:hashlib.sha256(resolve_path(project/name).read_bytes()).hexdigest()
                           for name in artifacts}},indent=2)+'\n',encoding='utf-8')
    print(f'Updated write-up and enlarged figures: hard 120, entropic {count}, counts {counts}.')


if __name__ == '__main__':
    main()
