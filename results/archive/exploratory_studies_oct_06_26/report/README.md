# Exploratory refraction research: illustrated overview

This is a GitHub-readable entry point to the archived laptop research.
For the full narrative, open [`index.html`](interactive/index.html) locally after cloning
the repository; its local figure and table links are self-contained.
The [archive README](../README.md) defines the scopes, metrics and limitations.

These are historical small-grid hard-map investigations, not the newer
fixed-1600-target or entropic-surface/Snell production sweeps. A passed
potential-change stopping test does not certify optical marginal feasibility.
Raw PF EOT and debiased divergence are different quantities.

| Study | Saved scope | Detailed notes and data |
| --- | --- | --- |
| Density controls | 48 optical source cases and 192 PF evaluations; all PF tests passed, but inherited optical marginal errors exceeded 1e-6 | [Notes](../density/metadata/report.md), [costs](../density/data/density_costs.csv) |
| Source refinement | 12 paired controls; maximum refined residual about 9.82e-9 | [Comparisons](../density/data/refinement.csv) |
| Continuation probes | 24 cases; six still exceeded marginal L1 1e-6 after 2000 additional iterations | [Notes](../convergence/metadata/REPORT.md), [diagnostics](../convergence/data/diagnostics.csv) |
| PF bias/reference grid | 72 evaluated rows; 64 passed all three cross/self tests; eight remain excluded from divergence rankings | [Notes](../pf/metadata/REPORT.md), [data](../pf/data/data.csv) |

Native target grids vary with NK; only the PF study's `common` rows use an
independent fixed target with 225 points. Rankings and comparisons are
conditional on these settings, not universal convergence claims. Failed
cold-start diagnostics and partial probes remain archived as diagnostics.

## Density costs and rankings

![Density PF cost comparisons](plots/density_costs.png)

![Density-dependent NK rankings](plots/density_rankings.png)

## Source convergence and refinement

![Source convergence diagnostics](plots/convergence.png)

![Effect of source refinement](plots/refinement.png)

## PF bias and evaluation status

![Raw entropic OT versus debiased divergence](plots/entropy_bias.png)

![PF evaluation status, including failed cases](plots/pf_status.png)

Each figure also has a same-named SVG. Numerical summary counts are retained
in [`summary.json`](metadata/summary.json). File hashes and original locations are in
[`../manifest.json`](../metadata/manifest.json); desktop checks are in
[`../desktop_import_validation.json`](../metadata/desktop_import_validation.json).

<!-- directory-index:start -->

## Directory guide

[Parent directory](../README.md)

| Folder | Contents |
| --- | --- |
| [data/](data/README.md) | Numerical arrays and tabular results for the parent run. NPZ/NPY files retain their original array names; CSV files retain their numerical columns. |
| [interactive/](interactive/README.md) | Interactive HTML visualizations. Open the HTML file in a browser; local assets are linked relative to this directory. |
| [metadata/](metadata/README.md) | Run settings, validation reports, and provenance records. Hash manifests specify their relative hash base where needed. |
| [plots/](plots/README.md) | Static figures for the parent run or archived export. See the parent README and parameter metadata for the metric, units, and experiment settings. |

<!-- directory-index:end -->
