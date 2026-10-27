# Results catalog

Start with **published** for completed plots, **runs** for reusable numerical
work, or **examples** for notebook NPZ bundles. Every collection has a README;
tables, figures, settings, and checkpoints live in separate subfolders.

## Latest completed comparisons

The extended fixed-target comparison completed October 26, 2026 at about
11 p.m. ET. It contains **120 hard-map cases and 120 smooth-map cases**, with
40 cases in each of three experiments. Target size is fixed at 1600 and
evaluation epsilon at `1/320`; the ground cost is squared Euclidean distance
between three-dimensional unit directions.

| Collection | What to inspect | Settings and validation |
| --- | --- | --- |
| [Regular entropic error](published/regular_entropic_error_oct_27_26/) | [Large combined plot](published/regular_entropic_error_oct_27_26/plots/experiments_1_2_3.png), individual plots, and [240 saved objectives](published/regular_entropic_error_oct_27_26/data/costs.csv). Cross-EOT, without debiasing or a square root. | [Metric and source hashes](published/regular_entropic_error_oct_27_26/metadata/validation.json) |
| [Debiased comparison](published/debiased_comparison_complete_oct_26_26/) | [Log-log comparison](published/debiased_comparison_complete_oct_26_26/entropic/plots/comparison_experiments_1_2_3_log.png); hard and smooth tables are in their own `data/` folders. | [Complete capture validation](published/debiased_comparison_complete_oct_26_26/metadata/snapshot.json) |
| [Complete smooth-refraction capture](published/smooth_refraction_complete_oct_26_26/) | Optical solution NPZs, cross-EOT checkpoints, plots, and marginal diagnostics. | [Parameters](published/smooth_refraction_complete_oct_26_26/metadata/parameters.json) and [120-case validation](published/smooth_refraction_complete_oct_26_26/metadata/snapshot.json) |
| [Hard-map mass error](published/hard_map_mass_error_oct_26_26/) | [Per-case maximum mass errors](published/hard_map_mass_error_oct_26_26/data/max_mass_error.csv), including empty target locations. | [Validation](published/hard_map_mass_error_oct_26_26/metadata/validation.json) |
| [Hard-map occupancy](published/hard_map_occupancy_oct_26_26/) | [Empty, singleton, and merged-location counts](published/hard_map_occupancy_oct_26_26/data/atom_counts.csv). | [Validation](published/hard_map_occupancy_oct_26_26/metadata/validation.json) |

The hard map has lower **debiased divergence** in all 120 paired cases. All
three EOT terms passed their recorded stopping tests. Optical marginal
feasibility is a separate diagnostic; these results do not assert zero optical
marginal residuals.

## Choose a category

| Folder | Contents |
| --- | --- |
| [published/](published/) | Validated completed captures and plot galleries listed above. |
| [runs/](runs/) | Canonical sweep outputs, reusable solve checkpoints, and the debiased evaluator's tables/caches. |
| [examples/](examples/) | Uniform refraction, sixteen density pairs, random surface examples, and sensitivity controls. |
| [reference_runs/](reference_runs/) | Dated default output from standalone C++/Python reflector benchmarks. |
| [Notebook exports](../notebooks/README.md) | New saves create `results/notebook_outputs/` with dated plot and interactive-HTML folders. |
| [archive/](archive/) | Historical varying-target studies, incomplete progress captures, desktop shards, smoke checks, and earlier reference arrays. |
| [local_runtime/](local_runtime/) | Machine-local worker supervision state and logs, excluded from GitHub. |

## Canonical working runs

| Run | Numerical scope | Parameters |
| --- | --- | --- |
| [Standard fixed-target](runs/standard_fixed_target_oct_22_26/) | Smaller hard-map sweep; 114 evaluations across experiments 1–4. | [Settings](runs/standard_fixed_target_oct_22_26/metadata/parameters.json) |
| [Hard map, 17 updates](runs/hard_map_17_updates_oct_08_26/) | Extended 120-case baseline with inherited final-loop limits. | [Settings](runs/hard_map_17_updates_oct_08_26/metadata/parameters.json) |
| [Hard map, 2000 updates](runs/hard_map_2000_updates_oct_08_26/) | The matching extended hard-map run with higher transport/self caps. | [Settings](runs/hard_map_2000_updates_oct_08_26/metadata/parameters.json) |
| [Smooth refraction, 2000 updates](runs/smooth_refraction_2000_updates_oct_26_26/) | Softmin normals and vector Snell refraction; full matching 120-case grid. | [Settings](runs/smooth_refraction_2000_updates_oct_26_26/metadata/parameters.json) |
| [Debiased evaluator](runs/debiased_comparison_oct_26_26/) | Independent hard/smooth tables and content-hashed self-EOT cache. | [Evaluation settings](runs/debiased_comparison_oct_26_26/metadata/evaluation_settings.json) |

## Keep the metrics separate

- **Cross-EOT** is the regularized transport objective from pushed points to the
  target, including its entropy penalty. It can remain positive for equal measures.
- **Sinkhorn divergence** subtracts half of each measure's self-EOT objective.
  The completed comparison plots divergence itself, without a square root.
- **Target-mass error and occupancy** compare masses/counts at target atoms;
  they do not measure spatial transport distance.
- **Optical solver cost and marginal residuals** diagnose the source solve.
  They are not interchangeable with the independent push-forward evaluations.

Older studies may use varying target sizes, different costs, evaluation epsilons,
or internal solver metrics. Their archived READMEs retain those distinctions.

## Names, locations, and regeneration

Run/update folders end with a readable date such as `oct_27_26`. Recorded
computation/capture dates are retained. If no original run date was recorded,
the suffix identifies its saved/publication or organization date, as explained
in that folder's README; it does not imply a rerun.

Use [scripts/result_paths.py](../scripts/result_paths.py) for readers and writers.
The [path catalog](path_catalog.json) maps old locations to their canonical
destinations and records migration provenance. See the
[layout guide](../docs/repository_layout.md) for typed run folders and the rule
that a directory with subdirectories keeps at most five regular files beside them.

From the repository root:

```bash
# Replot the published regular-EOT gallery from its saved CSV
python scripts/plot_saved_entropic_error.py

# Capture new completed debiased rows and refresh the write-up and large figures
python scripts/update_writeup_results.py
```

The active [LaTeX project](../geometric_optics_progress/) records its selected
capture in [results_capture.json](../geometric_optics_progress/metadata/results_capture.json).
The full previous catalog, including run history and old reproduction commands,
is retained in [catalog history](archive/catalog_history_oct_27_26/source/results_catalog_before_cleanup.md).

<!-- directory-index:start -->

## Directory guide

[Parent directory](../README.md)

| Folder | Contents |
| --- | --- |
| [archive/](archive/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |
| [examples/](examples/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |
| [local_runtime/](local_runtime/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |
| [published/](published/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |
| [reference_runs/](reference_runs/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |
| [runs/](runs/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |

<!-- directory-index:end -->
