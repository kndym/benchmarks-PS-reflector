# Results catalog

Start here to choose a dataset. **Internal optical solver cost, PF-to-target
EOT, projected-plane Wasserstein distance, and marginal residuals are different
quantities. Do not concatenate these families into a single unlabeled table.**

## Current fixed-target experiments

**Complete as of October 26, about 23:00 ET:** hard and entropic maps each have
120/120 cross-EOT and 120/120 debiased evaluations (40 per experiment).
All three EOT terms passed the recorded stopping test in every debiased row.
The hard map has lower divergence in all 120 pairs. The
[final optical capture](../../../published/smooth_refraction_complete_oct_26_26) validates
finite arrays, unit transmitted rays, normal identities, and matching checkpoints.
The final debiased capture and artifact hashes are indexed in
`../geometric_optics_progress/metadata/results_capture.json`.

**Overnight run resumed October 26 at 22:18 ET.** The user's overnight request
superseded the earlier cross-EOT pause. Sweep and plot pause markers were
archived with a `pause.resumed_20261026_2220.json` suffix, and their existing
scheduled tasks restarted. Both computation workers and the plot worker then
completed without a restart. The write-up's
current validated capture is indexed in
`../geometric_optics_progress/metadata/results_capture.json`; run
`python scripts/update_writeup_results.py` to refresh its text and large figures.
The overnight heartbeat is no longer needed after final validation and the
local/Overleaf update. The 89-case resource snapshot below is historical.

**[Comparison and resource benchmark before overnight resumption](../../progress_captures/debiased_resources_oct_26_26):
89/120 entropic cases, all 89 available maps processed.** Experiments 1 and 2
are complete at 40 cases each; experiment 3 has 9 completed cases. Cross-EOT
generation was paused at this capture, so the remaining 31 maps were unavailable. The
debiased task used four threads, normal CPU priority, and a
32 MB cost cache. The prior pause was lifted for the overnight run above.

**[Earlier comparison images: October 26, 18:39 ET](../../progress_captures/comparison_images_1839_oct_26_26).**
Entropic cross-EOT: 79/120 cases (40 experiment 1, 39 experiment 2).
Debiased entropic: 78/120 (40 experiment 1, 38 experiment 2). Hard: 120/120.
Experiment 3 entropic cases remain pending. The gallery includes both metrics,
linear-y and log-log variants, individual experiments, and captured CSVs.
The scheduler recovered memory-allocation interruptions in the largest cases;
both compute workers were running again at this check, with one restart each.

### Background workers independent of Codex

Windows Task Scheduler runs `Reflector-EntropicSweep`,
`Reflector-DebiasedEvaluation`, and `Reflector-ComparisonPlots` under the current
Windows user. Their supervisors are children of the Windows scheduling service,
so closing Codex or its SSH/CLI connection does not terminate them. Workers retry
after 30 seconds if they exit before all expected cases are complete; Task
Scheduler retries a failed supervisor after one minute. Tasks have no execution
time limit and resume at the user's next Windows login after a reboot or logout.
These are interactive-user tasks: they require the Windows user to be logged in.
The computer must remain powered on and awake for computation to continue.

Cross-EOT was paused at the user's request earlier on October 26, then resumed
under the overnight request above. During the pause, the sweep and its
comparison plot task were stopped, and persistent `sweep.pause.json` and
`plots.pause.json` markers in the local supervision folder prevented relaunch
after login or an automatic restart. Debiased evaluation remains active on
already-saved cross results. Removing a pause marker is an explicit resume
operation; starting or reinstalling tasks alone does not resume a paused role.

The debiased supervisor now assigns four threads to independent SciPy
logsumexp outputs, normal CPU priority, and up to four million cached cost
entries (32 MB for a float64 cost matrix). The measured 1600-point self solve
took 13.52 seconds with one thread and 8.19 seconds with four (1.65x faster),
with exactly matching objective, iteration count, convergence, and diagnostics.
That is a sample benchmark, not a guaranteed speedup for every case. Numerical
settings and completed checkpoints are preserved. The thread pool does not
parallelize separate large self-cost solves, limiting peak memory use.

Local live state and append-only logs are in
`results/local_runtime/worker_supervision/` (Git-ignored). Read `sweep.json`,
`debiased.json`, and `plots.json` there for current worker PIDs, start times,
restart counts, and completion state. Reinstall/start the tasks with
`powershell -File scripts/install_refraction_worker_tasks.ps1`; inspect them with
`Get-ScheduledTask -TaskName 'Reflector-*'`. Tasks stop doing work once their
expected completed-case sets are present. Parameter JSON publication is atomic
so the evaluator cannot observe a partially rewritten file during resumption.

[Maximum absolute hard push-forward mass error](../../../published/hard_map_mass_error_oct_26_26)
plots `max_j |mu_PF(y_j) - q_j|` for all 120 cases, including empty target
locations. The target reference is the saved mass `q_j = 1/1600`; source
weights are `1/NK`. Linear-y and log-log plots and per-case CSVs are included.

[Hard push-forward atom counts](../../../published/hard_map_occupancy_oct_26_26)
show how many target locations have aggregated mass different from
`1/NK`, across all 120 hard-map cases. Linear-y and log-log plots are included.
Coincident source weights are summed first; empty target locations have zero
mass and are included in the total. Plots also show the empty and merged-weight
components separately.

### October 26 progress and numerical fixes

[Debiased comparison and plots](../../progress_captures/debiased_progress_1559_oct_26_26):
**hard 120/120, entropic 75/120**, with all three EOT terms converged for each
published row. The metric is `EOT(PF,T) - (EOT(PF,PF) + EOT(T,T))/2` at epsilon
1/320, using squared Euclidean 3D cost. Hard interpolation has lower divergence
in all 75 currently paired cases. This is a partial comparison: experiment 1
has 40 pairs, experiment 2 has 35, and experiment 3 has no entropic cases yet.
At the October 26, 15:59 ET check, optical and debiased entropic evaluations
had both reached 75/120. The workers were no longer running, with no failure
reported in either error log; the cause of termination is unknown. The sweep,
debiased evaluator, and linear/log plot worker were restarted from saved
checkpoints. The next optical case is experiment 2, NK=11615.

Earlier progress check:
At the October 26, 13:48 ET progress check, the live entropic optical sweep
had 72/120 completed cases (40 in experiment 1, 32 in experiment 2). The optical
worker was solving experiment 2 at NK=7727; debiasing was evaluating the saved
NK=6746 push-forward. Both workers were active with no new errors in their logs.

[Entropic optical snapshot](../../progress_captures/smooth_refraction_partial_oct_26_26) contains
63 validated optical/PF checkpoints and the separate cross-EOT plots.
The sweep had stopped on a Windows CSV sharing violation. CSV and plot
publication now retry temporary sharing failures. Small-grid kernel smoothing
and scaling updates now use logsumexp to avoid underflow and reciprocal
overflow. Tests check agreement with the original updates at moderate
regularization and finite outputs at inverse regularization 16000.
Existing validated checkpoints are retained; the restarted optical worker
uses log-domain warm starts. These changes preserve the experiment
grids, objectives, caps and stopping tolerances. Potential-change stopping
still does not guarantee marginal feasibility; the published residuals remain
part of the interpretation of these results.

The canonical local checkout is `C:\Users\Kevin\Github\benchmarks-PS-reflector`.
The former `benchmarks-PS-reflector-v2` path is a Windows junction to the same
directory. Both Git histories were combined, including the latest main-branch
figure uploads. All 3353 copied files were verified byte-for-byte before the
workers restarted from saved checkpoints in the canonical folder. A recovery
copy is stored outside the GitHub directory at
`C:\Users\Kevin\.codex\backups\reflector-consolidation-20261026`.

### Cross-EOT datasets

All three main sweeps use kappa=0.6, Halton source points, a target fixed at
1600 points, and PF-to-target EOT with squared Euclidean 3D cost and evaluation
epsilon 1/320. Exact grids and tolerances are in each `parameters.json`.

| Dataset | Status | Maps / final-loop caps | Main outputs |
| --- | --- | --- | --- |
| [Standard sweep](../../../runs/standard_fixed_target_oct_22_26) | Complete: 114 evaluations across experiments 1–4 | Hard c-transform; 17/17 | [costs](../../../runs/standard_fixed_target_oct_22_26/data/costs.csv), [main plots](../../../runs/standard_fixed_target_oct_22_26/plots/experiments_1_2_3.png), [PF-epsilon sweep](../../../runs/standard_fixed_target_oct_22_26/plots/experiment_4_pushforward.png) |
| [Extended baseline](../../../runs/hard_map_17_updates_oct_08_26) | Complete: 120 evaluations, 40 per experiment 1/2/3 | Hard c-transform; 17/17 | [costs](../../../runs/hard_map_17_updates_oct_08_26/data/costs.csv), [plots](../../../runs/hard_map_17_updates_oct_08_26/plots/experiments_1_2_3.png), [validation](../../../runs/hard_map_17_updates_oct_08_26/metadata/validation_summary.json) |
| [Extended higher iteration](../../../runs/hard_map_2000_updates_oct_08_26) | Complete: 120 evaluations, 40 per experiment 1/2/3 | Hard c-transform; 2000/2000 | [costs](../../../runs/hard_map_2000_updates_oct_08_26/data/costs.csv), [plots](../../../runs/hard_map_2000_updates_oct_08_26/plots/experiments_1_2_3.png), [validation](../../../runs/hard_map_2000_updates_oct_08_26/metadata/validation_summary.json) |
| [Entropic-surface snapshot](../../progress_captures/smooth_refraction_partial_oct_22_26) | Partial; exact capture count in snapshot manifest | Softmin normals + Snell refraction; 2000/2000 | [comparison](../../progress_captures/smooth_refraction_partial_oct_22_26/plots/comparison_experiments_1_2_3.png), [costs](../../progress_captures/smooth_refraction_partial_oct_22_26/data/costs.csv), [manifest and hashes](../../progress_captures/smooth_refraction_partial_oct_22_26/metadata/snapshot.json) |

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

The completed entropic run remains in the Git-ignored local directory
`refraction_sweeps_extended_entropic_high_iter/`. Published snapshots contain
only completed CSV cases, their matching finite checkpoints, and comparison
plots rendered from that exact captured CSV. This avoids committing a changing
mix of plot and data versions. Workers resumed from checkpoints during the
October 26 checkout consolidation described above.

## Historical sweeps: different targets or metrics

| Archive | Status | Interpretation |
| --- | --- | --- |
| [Laptop expanded master](../../historical_runs/expanded_varying_target_oct_06_26) | Complete: 460 rows | Experiments 1–3: 60 **internal solver cost** rows. Experiment 4: 400 PF EOT rows with **NK-dependent target sizes**, not target fixed at 1600. |
| [Laptop pilot](../../historical_runs/pilot_varying_target_oct_22_26) | Complete: 114 rows | Older internal-cost experiments 1–3 and varying-target PF experiment 4. Distinct from the corrected fixed-target standard sweep. |
| [Superseded extended attempt](../../superseded_runs/varying_target_attempt_oct_06_26) | Incomplete: 39 experiment-1 rows | PF EOT with target varying alongside source; superseded before the corrected fixed-target sweep. Never use as a completed 120-case run. |

The laptop expanded master incorporates the completed desktop experiment-3 and
PF shards. Missing raw desktop checkpoints and the completed laptop NK=10000 PF
checkpoint shard have been folded into the master. We retain one canonical
`costs.csv` rather than publishing duplicate shard result directories.
`desktop_consolidation.json` records the import commit and consolidation choices.
The laptop's original per-shard CSV assembly inputs and manifests remain as
provenance; they are not additional independent experiments.

## Other existing benchmark bundles

The [laptop exploratory research archive](../../exploratory_studies_oct_06_26) is now
included, with a [GitHub-readable illustrated overview](../../exploratory_studies_oct_06_26/report/README.md)
and the original [HTML report](../../exploratory_studies_oct_06_26/report/interactive/index.html).
It contains density controls (192 PF evaluations), convergence probes (24),
refinement controls (12), and PF bias/reference-grid evaluations (72, including
eight that did not pass the existing tests). Failed evaluations and incomplete
diagnostic probes are preserved and identified; they are not additional solved
sweep cases. These exploratory studies differ from the fixed-target production
sweeps and the currently running entropic-normal comparison.

The remaining laptop transfer is complete: all 205 payload hashes, 34 local
HTML links and the ten surface image mappings were checked. See the
[desktop import validation](../../exploratory_studies_oct_06_26/metadata/desktop_import_validation.json)
and the [original provenance manifest](../../exploratory_studies_oct_06_26/metadata/manifest.json).

- [Uniform sensitivity](../../../examples/uniform_sensitivity_oct_02_26): ten kappa/source-angle controls;
  [numeric summary and surface gallery](../../../examples/uniform_sensitivity_oct_02_26/README.md). The current summary
  uses debiased Sinkhorn divergence with squared Euclidean 3D cost and epsilon
  0.005, its square root, and a separate target-atom mass L2 metric. These differ
  from the non-debiased sweep EOT; the square root is not exact unregularized W2.
- [Representative refractor examples](../../../examples/refractor_surfaces_sep_24_26): saved
  surfaces, metadata, figures and original LaTeX source.
- `results_refraction_*_NK1600.npz`: the source/target density-pair examples
  generated by `scripts/helper_all_pairs.py`, already tracked at this level.
- `results_NK*.npz`, `output_cpp_*`, `output_py_*`, `output_fast`,
  `output_python` and `output_refraction`: earlier C++/Python reference outputs.
- [Figures](../../../../archive/legacy_figures_oct_27_26) and [portable LaTeX project](../../../../geometric_optics_progress):
  visual outputs and the current research write-up, with supporting images.

## Reproduction and maintenance

The numerical pipeline remains shared: `refracter/sinkhorn.py` performs the
optical and PF solves, `scripts/generate_results.py` builds the existing hard-map
artifacts, and `refracter/optics.py` supplies the entropic-surface normal and
vector Snell functions. `scripts/run_refraction_sweeps.py` selects the map;
all experiments use its one PF evaluation/cache path. Comparison plotting
does not implement or rerun a solver.

Run commands are in the [repository README](../../../../README.md). To publish a later
entropic snapshot without interrupting its worker:

```bash
python -B scripts/snapshot_refraction_sweep.py results/runs/smooth_refraction_2000_updates_oct_26_26 results/published/NEW_CAPTURE_NAME --baseline results/runs/hard_map_2000_updates_oct_08_26
```

Use a new destination name; existing snapshots are not overwritten. Disposable
smoke tests, process logs/PIDs, compiler files, and duplicate ZIP previews are
kept out of Git. Superseded results are preserved in archives rather than deleted.

<!-- directory-index:start -->

## Directory guide

[Parent directory](../../../../README.md)

| Folder | Contents |
| --- | --- |
| [archive/](../../README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |
| [examples/](../../../examples/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |
| [local_runtime/](../../../local_runtime/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |
| [published/](../../../published/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |
| [reference_runs/](../../../reference_runs/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |
| [runs/](../../../runs/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |

<!-- directory-index:end -->
