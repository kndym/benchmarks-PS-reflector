# Standard fixed-target sweep

Completed corrected small-grid sweep: experiments 1/2/3 contain 5/5/4
PF-to-target EOT evaluations, and experiment 4 contains 100 PF-epsilon
evaluations. Target fixed at 1600 points, kappa=0.6, hard c-transform map,
and inherited transport/self final-loop caps 17/17.

For experiments 1–3, evaluation epsilon is fixed at 1/320. Experiment 4 varies
it independently as recorded in `costs.csv`. Every row is finite and passes
the existing PF convergence test. Source feasibility is separately described
by `residual_*.json`; no optical marginal-convergence guarantee follows from
the PF convergence flag.

This is not the older laptop pilot with the same historical folder name.
That distinct varying-target/internal-cost dataset is preserved at
`../archive/laptop_pilot/`.

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

<!-- directory-index:end -->
