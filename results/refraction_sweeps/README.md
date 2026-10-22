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
