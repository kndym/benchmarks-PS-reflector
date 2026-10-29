# Repaired kappa=0.8 progress

[October 29 write-up figures: 60 paired completed cases](../../../geometric_optics_progress/README.md) | [Overleaf](https://www.overleaf.com/project/69695773fcd59c00e6008659). The repaired-run figures are appended on pages 32–33.

**Bounded cost caching is implemented:** [256 MiB budget, benchmarks and rollout](../../runs/repaired_refractor_kappa_oct_28_26/#bounded-cost-cache). Existing cases finish before cached workers take over.

**October 29 execution upgrade:** [Additional parallel epsilon results](../parallel_kappa_epsilon_oct_29_26_0657/). The large-N workers switch to four threads after their current saved case; optical solutions are now shared. [Settings and benchmarks](../../runs/repaired_refractor_kappa_oct_28_26/#execution-acceleration).

**Latest capture: [October 29, 06:32 AM ET](../repaired_kappa_progress_oct_29_26_0632/).**

Captured October 28, 2026 at 08:59 AM Eastern Daylight Time. Original source and target runs are paused. Only the repaired kappa run is active.

Hard optical: **8/120**. Entropic optical: **8/120**. Debiased hard/entropic: **7/7**. Missing panels and points are pending work.

The optical cutoff is removed. Optical transport, source self solve, target self solve, and pushforward evaluation must satisfy marginal feasibility. Old invalid checkpoints are not reused. This is an early partial snapshot; it does not yet reproduce the full repaired joint-movement curve at large N.

[Audit and rerun settings](../../runs/repaired_refractor_kappa_oct_28_26/) | [Cutoff diagnosis](../../runs/kappa_cutoff_diagnostic_oct_28_26/)

## Cross-EOT

![Repaired cross-EOT](entropic/plots/comparison_experiments_1_2_3_log.png)

## Debiased divergence

![Repaired divergence](debiased/entropic/plots/comparison_experiments_1_2_3_log.png)

[All linear/log cross-EOT plots](entropic/plots/) | [All linear/log divergence plots](debiased/entropic/plots/)
