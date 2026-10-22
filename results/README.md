# Results catalog

Start here to choose a dataset. **Internal optical solver cost, PF-to-target
EOT, projected-plane Wasserstein distance, and marginal residuals are different
quantities. Do not concatenate these families into a single unlabeled table.**

## Current fixed-target experiments

All three main sweeps use kappa=0.6, Halton source points, a target fixed at
1600 points, and PF-to-target EOT with squared Euclidean 3D cost and evaluation
epsilon 1/320. Exact grids and tolerances are in each `parameters.json`.

| Dataset | Status | Maps / final-loop caps | Main outputs |
| --- | --- | --- | --- |
| [Standard sweep](refraction_sweeps/) | Complete: 114 evaluations across experiments 1–4 | Hard c-transform; 17/17 | [costs](refraction_sweeps/costs.csv), [main plots](refraction_sweeps/experiments_1_2_3.png), [PF-epsilon sweep](refraction_sweeps/experiment_4_pushforward.png) |
| [Extended baseline](refraction_sweeps_extended/) | Complete: 120 evaluations, 40 per experiment 1/2/3 | Hard c-transform; 17/17 | [costs](refraction_sweeps_extended/costs.csv), [plots](refraction_sweeps_extended/experiments_1_2_3.png), [validation](refraction_sweeps_extended/validation_summary.json) |
| [Extended higher iteration](refraction_sweeps_extended_high_iter/) | Complete: 120 evaluations, 40 per experiment 1/2/3 | Hard c-transform; 2000/2000 | [costs](refraction_sweeps_extended_high_iter/costs.csv), [plots](refraction_sweeps_extended_high_iter/experiments_1_2_3.png), [validation](refraction_sweeps_extended_high_iter/validation_summary.json) |
| [Entropic-surface snapshot](snapshots/entropic_high_iter_20261022/) | Partial; exact capture count in snapshot manifest | Softmin normals + Snell refraction; 2000/2000 | [comparison](snapshots/entropic_high_iter_20261022/comparison_experiments_1_2_3.png), [costs](snapshots/entropic_high_iter_20261022/costs.csv), [manifest and hashes](snapshots/entropic_high_iter_20261022/snapshot.json) |

Extended grids use 40 source counts from 100 to 20,000 and 40 refractor epsilons
from 0.0000625 to 0.025. Experiment 1 couples source count and refractor epsilon;
experiment 2 fixes refractor epsilon at 1/320; experiment 3 fixes source count at
1600. The separate evaluation epsilon stays fixed. The baseline parameter file
predates explicit cap fields and inherits 17 final updates per loop.

PF evaluation convergence is checked independently from the optical transport's
potential changes and marginal feasibility. All complete fixed-target cost rows
pass the recorded PF stopping test. The 2000-cap hard run meets potential-change
tolerance in all 120 transport and 240 self loops without hitting a cap; its
worst source/target marginal L1 residuals are approximately 0.01472/0.01472.

The live entropic worker remains in the Git-ignored local directory
`refraction_sweeps_extended_entropic_high_iter/`. Published snapshots contain
only completed CSV cases, their matching finite checkpoints, and comparison
plots rendered from that exact captured CSV. This avoids committing a changing
mix of plot and data versions. No active worker directory was moved or restarted.

## Historical sweeps: different targets or metrics

| Archive | Status | Interpretation |
| --- | --- | --- |
| [Laptop expanded master](archive/laptop_expanded/) | Complete: 460 rows | Experiments 1–3: 60 **internal solver cost** rows. Experiment 4: 400 PF EOT rows with **NK-dependent target sizes**, not target fixed at 1600. |
| [Laptop pilot](archive/laptop_pilot/) | Complete: 114 rows | Older internal-cost experiments 1–3 and varying-target PF experiment 4. Distinct from the corrected fixed-target standard sweep. |
| [Superseded extended attempt](archive/extended_varying_target_20261006/) | Incomplete: 39 experiment-1 rows | PF EOT with target varying alongside source; superseded before the corrected fixed-target sweep. Never use as a completed 120-case run. |

The laptop expanded master incorporates the completed desktop experiment-3 and
PF shards. Missing raw desktop checkpoints and the completed laptop NK=10000 PF
checkpoint shard have been folded into the master. We retain one canonical
`costs.csv` rather than publishing duplicate shard result directories.
`desktop_consolidation.json` records the import commit and consolidation choices.
The laptop's original per-shard CSV assembly inputs and manifests remain as
provenance; they are not additional independent experiments.

## Other existing benchmark bundles

The [laptop exploratory research archive](archive/laptop_research/) is now
included, with a [GitHub-readable illustrated overview](archive/laptop_research/report/README.md)
and the original [HTML report](archive/laptop_research/report/index.html).
It contains density controls (192 PF evaluations), convergence probes (24),
refinement controls (12), and PF bias/reference-grid evaluations (72, including
eight that did not pass the existing tests). Failed evaluations and incomplete
diagnostic probes are preserved and identified; they are not additional solved
sweep cases. These exploratory studies differ from the fixed-target production
sweeps and the currently running entropic-normal comparison.

The remaining laptop transfer is complete: all 205 payload hashes, 34 local
HTML links and the ten surface image mappings were checked. See the
[desktop import validation](archive/laptop_research/desktop_import_validation.json)
and the [original provenance manifest](archive/laptop_research/manifest.json).

- [Uniform sensitivity](uniform_sensitivity/): ten kappa/source-angle controls;
  [numeric summary and surface gallery](uniform_sensitivity/README.md). The current summary
  uses debiased Sinkhorn divergence with squared Euclidean 3D cost and epsilon
  0.005, its square root, and a separate target-atom mass L2 metric. These differ
  from the non-debiased sweep EOT; the square root is not exact unregularized W2.
- [Representative refractor examples](sinkhorn_refracter_examples/): saved
  surfaces, metadata, figures and original LaTeX source.
- `results_refraction_*_NK1600.npz`: the source/target density-pair examples
  generated by `scripts/helper_all_pairs.py`, already tracked at this level.
- `results_NK*.npz`, `output_cpp_*`, `output_py_*`, `output_fast`,
  `output_python` and `output_refraction`: earlier C++/Python reference outputs.
- [Figures](../figures/) and [portable LaTeX project](../geometric_optics_progress/):
  visual outputs and the current research write-up, with supporting images.

## Reproduction and maintenance

The numerical pipeline remains shared: `refracter/sinkhorn.py` performs the
optical and PF solves, `scripts/generate_results.py` builds the existing hard-map
artifacts, and `refracter/optics.py` supplies the entropic-surface normal and
vector Snell functions. `scripts/run_refraction_sweeps.py` selects the map;
all experiments use its one PF evaluation/cache path. Comparison plotting
does not implement or rerun a solver.

Run commands are in the [repository README](../README.md). To publish a later
entropic snapshot without interrupting its worker:

```bash
python -B scripts/snapshot_refraction_sweep.py results/refraction_sweeps_extended_entropic_high_iter results/snapshots/NEW_CAPTURE_NAME --baseline results/refraction_sweeps_extended_high_iter
```

Use a new destination name; existing snapshots are not overwritten. Disposable
smoke tests, process logs/PIDs, compiler files, and duplicate ZIP previews are
kept out of Git. Superseded results are preserved in archives rather than deleted.
