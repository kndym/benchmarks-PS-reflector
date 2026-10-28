# Three isolated refractor problems — October 27, 2026

Three new problems change one physical parameter relative to the previous uniform-patch sweep. Angles below are degrees; theta is azimuth and phi is polar.

| Problem | Source theta | Source phi | Target theta | Target phi | Kappa |
| --- | --- | --- | --- | --- | --- |
| source_only | 15–60 | **15–55** | 18–36 | 18–36 | 0.6 |
| target_only | 15–60 | 15–45 | **18–45** | 18–36 | 0.6 |
| kappa_only | 15–60 | 15–45 | 18–36 | 18–36 | **0.8** |

Each problem repeats experiments 1–3 from the prior extended sweep: 40 source sizes from 100 to 20,000 with coupled epsilon, those same sizes with epsilon 1/320, and 40 epsilons from 0.0000625 to 0.025 at 1,600 source points. The target stays at 1,600 points. There are 120 cases per map, 240 per problem, and 720 across all three problems.

Hard c-transform and entropic-surface/Snell maps run independently. Optical transport and identity iteration caps are 2,000 with tolerance 1e-5. Pushforward EOT and debiased Sinkhorn divergence use squared Euclidean 3D cost, epsilon 1/320, tolerance 1e-9, and a 2,000 iteration cap. Sampling and uniform patch densities follow the previous sweep. Previous images are not reused because these isolated parameter changes do not match the old random cases.

Diagnostics include source/target optical marginal L1 residuals, final iteration counts and potential changes, entropic normal and transmitted-direction checks, cross-EOT convergence, both self-EOT terms, raw and reported divergence, and checkpoint hashes. Cross-EOT and divergence comparison plots refresh in linear and logarithmic scales. These metrics match the previous comparison; the divergence is not the older planar/angle diagnostic from the random surface examples.

Each problem contains `hard/`, `entropic/`, and `debiased/` outputs, `logs/`, and `metadata/worker_status.json` for live worker status. Numerical settings are recorded in each mode's `metadata/parameters.json`. One reduction thread per worker, BLAS threads limited to one, and chunk size 256 limit concurrent memory use; numerical tolerances and grids match the prior sweep.

Start or resume a problem with `python -B scripts/run_three_refractor_problems.py --case source_only` (or `target_only`, `kappa_only`). Checkpoints are reused. Windows tasks named `Refractor-Isolated-source_only`, `Refractor-Isolated-target_only`, and `Refractor-Isolated-kappa_only` run independently of the chat. Earlier paused sweep tasks are left paused.

On October 28, the target and kappa diagnostic workers were recovered after an asynchronous pairing bug: an entropic CSV row could arrive before its hard counterpart. The evaluator now waits for the matching row, and comparison plotting accepts independent progress on each curve. All six original optical workers continued without interruption. The three Windows tasks now run with `--diagnostics-only`; status records retain the existing optical worker IDs in `optical_workers`. Diagnostic logs retain historical tracebacks from before this repair.

<!-- directory-index:start -->

## Directory guide

[Parent directory](../README.md)

| Folder | Contents |
| --- | --- |
| [kappa_only/](kappa_only/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |
| [source_only/](source_only/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |
| [target_only/](target_only/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |

<!-- directory-index:end -->
