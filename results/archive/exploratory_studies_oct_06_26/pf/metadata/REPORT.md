# Isolated pushforward evaluation research

Completed 72 evaluated rows across 1 problems. 64 cross evaluations converged; 64 had converged cross and both self evaluations. Rankings and plots exclude failures.

All refractor solutions and evaluation objectives import the existing active pipeline. No tracked code was edited. Fixed source geometry and target geometry follow the original benchmark. Completed control study: original kappa=0.6 uniform. Lower-kappa and skew-density variants were configured but not run to preserve the time budget; other agents cover problem variants.

Completed grid: NK=100,225,400 (900 configured but omitted from completed rows); refractor epsilon=0.0125,0.00125,0.000125; PF epsilon=0.0125,0.003125,0.00078125,0.000125. Native target changes with NK. Common target is an independent fixed 225-point Halton grid (skip=5000), held fixed for every NK. This controls changing target discretization but remains a finite empirical reference.

PF evaluator uses 1200 iterations and potential-change tolerance 1e-8 for the bounded exploratory runs; some cached evaluations have stricter 1e-9 /2500 settings. Every cached result retains its exact final potential change and iteration count. These are exploratory evaluator settings, not changes to production defaults. Numerical finiteness is not proof of convergence. The original refractor solver retains its inherited iteration cap and threshold, so its costs and hard maps are approximations.

## Findings

Of 22 converged NK comparisons, 6 have nonincreasing raw EOT and 11 have nonincreasing divergence. Raw EOT and divergence disagree on complete NK ranking in 12 comparisons.

The reported raw EOT includes the KL regularization contribution. For fixed measures it is nondecreasing in epsilon, because its objective is transport + epsilon times nonnegative KL. It is not a pure matching error, and even identical measures generally have positive raw EOT. Debiasing subtracts half the two self costs. These self costs depend on support and weights; therefore raw EOT rankings can change even when matching quality improves.

Converged original example, refractor epsilon=0.00125 and PF epsilon=0.0125: NK100→225 worsens native raw EOT (0.01200766→0.01209485), while divergence improves (1.1276e-5→1.0737e-5). With the common reference, raw EOT also worsens (0.01201524→0.01209576), while divergence improves (1.0710e-5→9.901e-6). Changing native targets is therefore not the only cause.

Strict low-epsilon timing probe: original NK100, refractor epsilon=0.0125, PF epsilon=0.000125 native cross converges in 1751 iterations at tolerance1e-9, 2.18 seconds excluding startup, with potential change9.984e-10. The pushed measure self solve at the same epsilon does not converge in2500 iterations (change1.43e-8). Lower refractor epsilon=0.000125 native cross still has change2.215e-7 after1200 iterations. Low epsilon is finite here but not uniformly converged under a fixed small iteration budget.

The transport-only component of the entropic coupling was not separately measured: the active evaluator returns the dual objective and diagnostics, without coupling potentials. Raw EOT and debiased divergence provide the requested bias comparison without duplicating the evaluator.

## Interpretation and scope

There is no universal monotonic NK improvement demonstrated. Hard argmin pushforwards discretize mass and can change abruptly with regularization and target spacing; inherited solver caps add approximation error. Hold PF epsilon, density and geometry fixed when comparing NK. A common high-resolution independent target plus divergence is a more defensible quality comparison than comparing raw EOT across different epsilons. Small research grids are suggestive and do not establish universal laws.

Feydy et al. (2019) analyze debiased Sinkhorn divergences and their geometric properties: [primary paper](https://proceedings.mlr.press/v89/feydy19a.html). Their theory concerns properly computed divergences; it does not imply monotone finite-NK behavior of this hard refractor map.

Artifacts: research.py reproduces the bounded scan; data.csv has all rows including failures; eval_*.json has individual diagnostics; solve_*.npz has hard PF supports; comparison.png compares native/common EOT/divergence; comparison_summary.json and rankings.json record rankings.