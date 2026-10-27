> Consolidated from laptop commit `9d4cd21`; missing desktop raw checkpoints and laptop NK10000 PF checkpoints are included here. See `desktop_consolidation.json`. The original laptop assembly inputs and manifest below are provenance, not independent results.

# Historical laptop expanded sweeps (completed 6 October 2026)

These are the completed original-pipeline sweeps from source commit `6905ba2`, preserved without rerunning or correcting the numerical experiments. **They are distinct from the newer fixed-independent-target experiments on the desktop branch.**

## Metrics and scope

`costs.csv` is the canonical merged table: **460 rows**.

| Experiment | Rows | Swept variables | Meaning of `cost` |
| --- | ---: | --- | --- |
| 1 | 20 | NK and coupled source epsilon | Inherited identity-corrected source solver `total_cost`; NOT PF-to-target EOT |
| 2 | 20 | NK, source epsilon fixed at 0.003125 | Same internal source metric; NOT PF-to-target EOT |
| 3 | 20 | Source epsilon, NK fixed at 1600 | Same internal source metric; NOT PF-to-target EOT |
| 4 | 400 | 20 NK x 20 PF epsilon; coupled source epsilon | Separately solved KL-regularized PF-to-target EOT, squared Euclidean ground cost on 3D spherical points |

All experiments use the original uniform refraction patches, kappa=0.6, normalized weights and Halton sampling. **The target grid varies with NK**: it has NK points from the original target sampler. It is not the later fixed-target benchmark. Experiment 4 varies source epsilon only through `1 / (8 * floor(sqrt(NK)))`; it is not an independent 20 x 20 x 20 scan.

NK: `[100,127,162,207,264,336,428,546,695,886,1129,1438,1833,2336,2976,3793,4833,6158,7848,10000]`.

Independent source-epsilon grid (experiment 3) and PF-epsilon grid (experiment 4): 20 geometrically spaced values from 0.000125 to 0.0125 inclusive. Exact floating-point values and rendering arguments are in `parameters.json`; exact per-row source/PF epsilons are in `costs.csv`.

## Solver settings and qualifications

The source solver retains the inherited warm-start/continuation and screening behavior: cap_iter=16 (17 final loop iterations under the inherited `i > cap_iter` condition), cap_thr=1e-5, and original identity corrections. No stopping defaults were corrected for these runs. Source residual diagnostics are in `source_residuals.csv`: none of the 60 source cases meets maximum marginal L1 <= 1e-6; the largest is approximately 0.890247. These are approximate source solutions and hard maps, not certified converged source OT costs.

PF evaluation uses the active `entropic_ot_cost`, potential-change tolerance 1e-9, max_iter=10000, chunk_size=512, descending-epsilon warm starts within live workers and bounded caching. All 400 PF rows passed that existing potential-change stopping test; maximum recorded iterations=1121. This is not certification under the proposed marginal L1 test. Final rendering imported completed rows and ran no additional solves.

## Shards and provenance

- Laptop root checkpoints: 333 PF rows, including the first 13 PF epsilons for NK=7848.
- Separate laptop NK=10000 shard: 20 PF rows, folded into this directory during desktop consolidation (no duplicate shard table).
- Desktop imports: `pc_pf_costs.csv` contains 40 PF rows (NK=3793 and 6158); `pc_pf_7848_costs.csv` contains the remaining 7 low-epsilon PF rows for NK=7848.
- Source experiments 1/2: 40 laptop cases with saved source NPZs. Experiment 3: 20 desktop rows in `pc_e3_costs.csv`, with residuals in `pc_e3_residuals.json`; the original laptop bundle lacked raw experiment-3 NPZs; desktop consolidation supplied them here.
- Desktop used NumPy 2.2.1 / SciPy 1.16 / Python 3.13. Sampler reconstruction can differ by about one ulp across environments; imported costs are preserved as recorded, not recomputed.

`combined_costs.csv` and `source_costs.csv` are assembly intermediates retained for provenance. The final `costs.csv` contains all 460 rows and supersedes intermediate tables for analysis. Imported rows do not all have local PF JSON checkpoints: completeness must be checked from the merged table.

`validation.json` records 9 distinct NK rankings and 14 pairwise reversals across PF epsilon. `ordering.json` gives full rankings. The figures for experiments 1–3 explicitly label the internal source metric; do not relabel them as PF EOT.

`laptop_bundle_manifest.json` records SHA-256 hashes of the laptop working-copy bytes for root-level artifacts in all three historical laptop sweep folders. Git may normalize text line endings on checkout; binary hashes remain unchanged. Duplicate preview subfolders are omitted from the bundle. Original run logs are retained; parameters.json describes the final assembly, not every worker command.

## Other laptop-only research

The separate `.research/refraction/` folder remains Git-ignored at the user's request. It contains completed density (192 PF evaluations), convergence (24 probes), dense refinement (12 controls), and PF bias/reference-grid (72 evaluations, 64 passing existing tests) studies; source CSVs, scripts and the illustrated `report/index.html` with six PNG/SVG figures remain local. These are additional relevant result bundles, **not included in this tracked bundle**. Their scope and failures are documented in their local reports. No numerical work was restarted during packaging.

An audit against desktop commit `26b5b49` found no distinct numerical data in `results/examples/uniform_sensitivity_oct_02_26/` or `results/examples/refractor_surfaces_sep_24_26/`; apparent text differences were line endings only, and NPZ files matched Git blob hashes. Laptop-only rendered artifacts are the ten `results/examples/uniform_sensitivity_oct_02_26/*_surface.png` files and `results/examples/refractor_surfaces_sep_24_26/main (2).pdf`. The latter also has local compiler `.aux`/`.log` files. These rendered/temporary files were reported but not added to this sweep-results bundle, as requested by desktop coordination.

<!-- directory-index:start -->

## Directory guide

[Parent directory](../README.md)

| Folder | Contents |
| --- | --- |
| [checkpoints/](checkpoints/README.md) | Saved numerical state, separated by solve type. These files allow evaluation or resumption without recomputing completed optical solves. |
| [data/](data/README.md) | Numerical arrays and tabular results for the parent run. NPZ/NPY files retain their original array names; CSV files retain their numerical columns. |
| [diagnostics/](diagnostics/README.md) | Numerical diagnostics, separated from plotted transport costs and optical solution arrays. |
| [metadata/](metadata/README.md) | Run settings, validation reports, and provenance records. Hash manifests specify their relative hash base where needed. |
| [plots/](plots/README.md) | Static figures for the parent run or archived export. See the parent README and parameter metadata for the metric, units, and experiment settings. |

<!-- directory-index:end -->
