# Repaired kappa-only rerun and numerical audit â€” October 28, 2026

[Published progress plots](../../published/repaired_kappa_progress_oct_28_26/)

All three original isolated sweeps are paused and their scheduled tasks disabled. Only kappa=0.8 is being rerun, in `kappa_only/`, using solver version `untruncated_marginal_v2`. The original data remains separate.

The rerun repeats the original 40-size and 40-epsilon grids for experiments 1â€“3, with 1,600 target points, hard and entropic maps, cross-EOT and debiased divergence diagnostics. It has 120 cases per map. Optical transport and both optical self solves require marginal L1 error <=1e-7. Iteration caps are 10,000 instead of 2,000; reaching the cap without feasibility raises an error. Pushforward EOT also requires marginal feasibility, independently of its potential-change threshold. Old solver checkpoints are rejected by provenance guards and versioned names.

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

## Execution acceleration — October 29

Two additional optical processes now compute the original epsilon sweep (experiment 3) in `parallel_epsilon/kappa_only/`, with two reduction threads each and 1,024-row chunks. The epsilon values, physics, iteration caps and tolerances are unchanged; larger epsilon values execute first. Source and target parameter variations remain paused.

The original large-N workers are adopted by `accelerate_kappa_workers.py`. Each finishes its current N=10,140, epsilon=.00125 case before restarting with four reduction threads and identical numerical settings. This preserves the in-flight case. Handoff state and targets are recorded in `kappa_only/metadata/worker_status.json` and `kappa_only/logs/acceleration.log`.

Hard and entropic consumers now reuse a shared optical solution instead of repeating the same optical transport and self solves. The cache key hashes point clouds, weights, kappa, epsilon, method version, iteration caps, chunk size and tolerances. Only complete finite solutions are atomically published; consumers derive their own maps. Failed solves are never cached.

Validated epsilon-case optical/map checkpoints and cross-EOT evaluations are prefetched into the main output without modifying its live CSV. When the main loop reaches those cases it reuses the completed work. The parallel sweep has its own CSV and plots until then; count the union of case keys, rather than adding rows that already overlap.

An actual 4,096-point identity update measured 1.557 s with one thread, 1.071 s with two, and 0.921 s with four (about 1.69x faster). Output differences were zero in this benchmark. Pure reduction benchmarks were about 2.4x faster at four threads; that is not an end-to-end speedup claim. The acceleration code and numerical regression suite pass 87 tests.

Full cost-matrix caching is not enabled: one 20,000-by-20,000 float64 matrix needs 3.2 GB, comparable to currently free RAM. A bounded cache or SSD-backed cache would need its own measurements before enabling it.

[Execution settings](metadata/acceleration_plan.json) | [Full-update benchmark](metadata/parallel_benchmark.json) | [Benchmark script](../../../scripts/benchmark_refraction_resources.py)

The execution supervisor and resource benchmark use `psutil`, already installed on this run machine. The numerical solver and shared-solution cache do not require it.
