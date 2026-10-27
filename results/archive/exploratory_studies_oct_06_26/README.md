# Archived laptop exploratory refraction research

Completed exploratory output from `.research/refraction/`, originally collected on 6 October 2026 using active pipeline commit `6905ba2`. Packaged on 22 October 2026 at Kevin's request. **No numerical experiments were started or restarted, no numerical code was changed, and no solver defaults were changed.** The original local research folder remains Git-ignored; this is its results-only archive.

Start with the [GitHub-readable illustrated overview](report/README.md) or the
[full HTML report](report/interactive/index.html). The HTML contains six scientific figures
in PNG and SVG, interpretations, exact scope, metric definitions, failed-case
qualifications and links to the archived CSVs. Its archived copy has two
packaging-only adjustments: its original-folder description and the link to the
omitted local report generator. All other copied payloads preserve their original
bytes. `manifest.json` records original locations and working-copy SHA-256 hashes
(Git may normalize text line endings). The desktop added the Markdown overview
and `desktop_import_validation.json`; the original payload hashes remain unchanged.

## Completed studies

| Folder | Authoritative table | Scope and qualifications |
| --- | --- | --- |
| `density/` | `density_costs.csv` | 48 source problems: uniform source, four target densities, NK=100/225/400, four source epsilons; 192 PF evaluations, all passing existing PF tests. All 48 inherited source marginal errors exceed 1e-6. `refinement.csv` contains original/refined comparisons for 12 controls reaching source residual below 1e-8. `cold_probe.csv` intentionally retains failed cold-start diagnostics; these are not solved transport. |
| `convergence/` | `diagnostics.csv`, also JSON | 24 uniform/tilted cases, NK=100/225/400, four epsilons from 0.000125 to 0.0125. Six remain above marginal L1=1e-6 after 2,000 additional iterations. These are research continuation controls, not changed production defaults. |
| `pf/` | `data.csv` | 72 evaluated rows; 64 pass cross and both self stopping tests, eight excluded from divergence rankings. Original uniform kappa=0.6; NK=100/225/400, three source epsilons, four PF epsilons, native and fixed independent 225-point target grids. |
| `report/` | `summary.json`, `density_winners.csv` | Consolidated interpretation, six PNG/SVG figure pairs, HTML report. |

The individual `report.md` / `REPORT.md` files retain the detailed study notes. Their references to Python scripts describe the original local investigation; scripts are intentionally not packaged here. `pf/eval_*.json` preserves raw evaluator diagnostics, including failures, cached aliases and probes not promoted into the completed CSV (including partial NK=900 probes). `pf/solve_*.npz` preserves available source hard-map checkpoints. **Do not count the JSON filenames or NPZ files as additional completed sweep rows.** Configured lower-kappa variants were not completed. The authoritative completed PF sweep remains exactly the 72 rows in `pf/data.csv`.

## Metrics must remain distinct

- Density `source_objective`: reconstructed raw cross dual bound, not the inherited identity-corrected source solver `total_cost`, and not a certified primal cost at infeasible marginals.
- PF `ot` / density `pf_cost`: separately evaluated KL-regularized PF-to-target OT with squared Euclidean ground cost on 3D spherical points; includes the entropic contribution.
- PF `divergence`: PF OT minus half the two separately evaluated self costs. Only rows with all three passing evaluations are used in divergence rankings. This differs from the original refractor solver's inherited identity correction.
- A `converged` flag here means the evaluator passed its recorded potential-change stopping test. It does not certify the proposed full-plan marginal L1 tolerance. Numerical finiteness does not establish source feasibility.
- Native target sampling varies with NK. Only the PF study's `common` rows hold an independent empirical 225-point target fixed (Halton skip=5000). These exploratory results are distinct from the newer desktop fixed-target sweeps and entropic-normal/Snell comparison.

The main conclusions are conditional: NK rankings depend on density and both epsilons; source refinement can reverse a ranking; raw entropic OT and debiased divergence can rank the same maps differently. Small grids and inherited approximate source solutions do not establish universal asymptotic behavior.

## Archive inventory

205 copied research output files, plus this README and `manifest.json`: **207 archive files**. Ten additional existing surface PNGs are committed in `results/examples/uniform_sensitivity_oct_02_26/`; each corresponds by filename to an already-tracked NPZ. Their paths, paired NPZs and hashes are listed in `manifest.json`. They are visualizations of existing results, not new numerical experiments.

Excluded: process logs/PIDs, helper and numerical Python scripts, compiler files, extracted paper text, smoke tests and duplicate preview folders. The original ignored folder retains those local provenance/development materials. The HTML report's local links were checked after packaging.

<!-- directory-index:start -->

## Directory guide

[Parent directory](../README.md)

| Folder | Contents |
| --- | --- |
| [convergence/](convergence/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |
| [density/](density/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |
| [metadata/](metadata/README.md) | Run settings, validation reports, and provenance records. Hash manifests specify their relative hash base where needed. |
| [pf/](pf/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |
| [report/](report/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |

<!-- directory-index:end -->
