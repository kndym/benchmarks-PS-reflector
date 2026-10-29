# Repaired kappa-only rerun and numerical audit — October 28, 2026

[Published progress plots](../../published/repaired_kappa_progress_oct_28_26/)

All three original isolated sweeps are paused and their scheduled tasks disabled. Only kappa=0.8 is being rerun, in `kappa_only/`, using solver version `untruncated_marginal_v2`. The original data remains separate.

The rerun repeats the original 40-size and 40-epsilon grids for experiments 1–3, with 1,600 target points, hard and entropic maps, cross-EOT and debiased divergence diagnostics. It has 120 cases per map. Optical transport and both optical self solves require marginal L1 error <=1e-7. Iteration caps are 10,000 instead of 2,000; reaching the cap without feasibility raises an error. Pushforward EOT also requires marginal feasibility, independently of its potential-change threshold. Old solver checkpoints are rejected by provenance guards and versioned names.

## Verification

The complete suite passed **79 tests**. The randomized audit ran **72 evaluations**: 68 feasible convergences, four explicit nonconvergences, and **zero false successes**. Cached and chunked objectives agree within 1e-10. The randomized report records seed 628, all point clouds, weights, epsilons, and diagnostics.

The repaired end-to-end smoke run completed all three experiments, both maps, debiased evaluation, and plots. Production's first N=100 case has source/target optical L1 errors 5.8e-15 and 4.6e-8; source and target self-solve errors are below 8.4e-8. Those checks used to be omitted from acceptance.

| Confirmed numerical issue | Repair |
| --- | --- |
| Optical cutoff can freeze entire rows and report zero change | Removed cutoff from main and self updates; supported nonfinite reductions raise |
| Tiny epsilon can make an EOT potential update look converged while marginals are wrong | Require source and target marginal feasibility |
| Huge finite weights overflow sums; large potential gauges lose precision | Stable normalization and gauge centering |
| Dot-product squared distances lose nearby-ray distances; tiny-kappa costs round to zero | Coordinate differences and log1p |
| Extreme normal scales overflow/underflow; near-0/1 kappa amplifies cancellation | Scale-safe normalization and stable normal/tangent formulas |
| Resume/comparison can mix physics or solver settings | Provenance checks, separate outputs, and versioned checkpoints |
| Partial experiment selections can wait forever | Completion uses the selected experiment keys |
| Polar endpoint can be excluded by undefined azimuth | Pole-aware patch membership |
| Failed supervisor can leave sibling workers running | Failure cleanup terminates owned siblings, preserving checkpoints |

## Remaining limits

Small-epsilon problems can be ill-conditioned. The audit includes cases that properly report nonconvergence instead of accepting inaccurate output. This repair does not promise convergence for arbitrary parameters under a finite iteration limit.

Finite cost is also weaker than physical same-side Snell reachability: the latter needs x dot y >= kappa. An incident north-pointing ray with a target 30 degrees away and kappa=.9 illustrates this domain restriction. The current kappa=.8 source and target patches satisfy it; the general domain policy is unchanged.

[Audit summary](metadata/audit_report.json) | [Reproducible randomized cases](metadata/randomized_eot_audit.json) | [Randomized audit script](../../../scripts/audit_eot_randomized.py) | [Cutoff diagnosis and large-case controls](../kappa_cutoff_diagnostic_oct_28_26/)

<a id="execution-acceleration"></a>
## Execution acceleration - October 29

Two additional optical processes compute the original epsilon sweep (experiment 3) in `parallel_epsilon/kappa_only/`, with two reduction threads each and 1,024-row chunks. The epsilon values, physics, iteration caps and tolerances are unchanged; larger epsilon values execute first. Source and target parameter variations remain paused.

The large-N workers are adopted by `accelerate_kappa_workers.py`. Each saves its current case before restarting with four reduction threads and identical numerical settings. Handoff state and targets are recorded in `kappa_only/metadata/worker_status.json` and `kappa_only/logs/acceleration.log`.

Hard and entropic consumers reuse a shared optical solution instead of repeating the same optical transport and self solves. The key hashes point clouds, weights, kappa, epsilon, method version, iteration caps, chunk size and tolerances. Only complete finite solutions are atomically published; consumers derive their own maps. Failed solves are never cached.

Validated epsilon-case optical/map checkpoints and cross-EOT evaluations are prefetched into the main output without modifying its live CSV. The main loop reuses that work when it reaches those cases. The parallel sweep has its own CSV and plots until then; count the union of case keys instead of adding overlapping rows.

A 4,096-point identity update measured 1.557 s with one thread, 1.071 s with two, and 0.921 s with four (about 1.69x faster), with zero measured output difference. These are update timings rather than whole-sweep ETA claims.

[Execution settings](metadata/acceleration_plan.json) | [Thread benchmark](metadata/parallel_benchmark.json) | [Benchmark script](../../../scripts/benchmark_refraction_resources.py)

The operational supervisor and resource benchmark use psutil, installed on this run machine. The numerical solver, shared-solution cache and cost-block cache do not require it.

<a id="bounded-cost-cache"></a>
## Bounded cost cache - October 29

The optical cost cache is implemented and configured at **256 MiB of retained blocks per active optical solve**. Both worker groups switch each worker after its current case is saved, preserving in-flight work. Status includes `cost_cache_mib_per_optical`; producing workers emit `Cost cache NK=...` and save per-case `cost_cache_NK*.json` statistics when a fresh optical solve finishes.

The cache owns private read-only point-cloud snapshots. Keys distinguish kappa, source-self, target-self and cross blocks, and row ranges. Cached costs are read-only. Fixed-prefix retention prevents sequential-scan cache thrashing. Budgets allocate 80% to source-self costs and 10% each to target-self and cross costs; these bound retained blocks, not every solver working array. Allocations are demand-driven and released after each solve.

Warm-start self-smoothing now computes the same full row reductions in tiles. At N=20,000, one 256-row source-self cost block needs about 41 MB instead of the full 3.2 GB source-self cost matrix. The cross warm-start matrix and usual reduction buffers still exist. Point counts, equations, iteration schedules, caps and tolerances are unchanged.

| Four-thread source-self update | Cache off | 256 MiB cache | Speedup |
| --- | ---: | ---: | ---: |
| 4,096 points, warm median | 1.217 s | 0.562 s | 2.17x |
| 8,192 points, warm median | 4.738 s | 3.903 s | 1.21x |

The larger matrix only partly fits the cache. Measured output differences were zero and retained blocks stayed within budget. These are timings under production load, not whole-sweep speedup claims. The full suite passed **97 tests**, including dense-versus-tiled warm-start equivalence, cached/uncached optical solves and CLI outputs, memory bounds, input immutability, kappa changes and cleanup after exceptions.

[Cache rollout settings](metadata/cost_cache_plan.json) | [Raw cache benchmark](metadata/cost_cache_benchmark.json) | [Cache benchmark script](../../../scripts/benchmark_cost_cache.py)
