# Superseded extended varying-target attempt

This run stopped with 39 experiment-1 PF-to-target EOT evaluations before the
target-size correction. The target point count varies with source NK.
It is preserved for provenance, not resumed or relabeled as a fixed-target run.

The completed replacement is `../../refraction_sweeps_extended/`, with target
fixed at 1600 points and 120 evaluations. The 2000-cap comparison is
`../../refraction_sweeps_extended_high_iter/`.

<!-- directory-index:start -->

## Directory guide

[Parent directory](../README.md)

| Folder | Contents |
| --- | --- |
| [checkpoints/](checkpoints/README.md) | Saved numerical state, separated by solve type. These files allow evaluation or resumption without recomputing completed optical solves. |
| [data/](data/README.md) | Numerical arrays and tabular results for the parent run. NPZ/NPY files retain their original array names; CSV files retain their numerical columns. |
| [diagnostics/](diagnostics/README.md) | Numerical diagnostics, separated from plotted transport costs and optical solution arrays. |
| [logs/](logs/README.md) | Local process logs and runtime identifiers. These are operational records and are excluded from GitHub. |
| [metadata/](metadata/README.md) | Run settings, validation reports, and provenance records. Hash manifests specify their relative hash base where needed. |
| [plots/](plots/README.md) | Static figures for the parent run or archived export. See the parent README and parameter metadata for the metric, units, and experiment settings. |
| [previous_refractor_metric/](previous_refractor_metric/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |

<!-- directory-index:end -->
