# Historical pilot sweeps

Distinct earlier pilot, preserved without numerical reruns. `costs.csv` contains 114 rows: experiments 1/2/3 have 5/5/4 rows and experiment 4 has 100 rows. Exact grids are in `parameters.json`.

Experiments 1–3 report the inherited internal source solver `total_cost`, **not pushforward-to-target EOT**. Experiment 4 separately evaluates PF-to-target KL-regularized OT using squared Euclidean 3D cost and an NK-dependent target grid. It covers the coupled source-epsilon schedule plus four fixed-source-epsilon schedules, each with five NK and four PF-epsilon values.

Original pipeline: uniform refraction patches, kappa=0.6, Halton samples, inherited source stopping cap_iter=16 (17 final iterations), threshold=1e-5. PF evaluator uses its then-default max_iter=2000 and potential-change tolerance=1e-9. See each checkpoint's convergence status; the table and finite costs do not certify source marginal convergence. This pilot predates expanded-run warm-start/caching/threading options.

These are historical varying-target results, not the newer desktop fixed-target EOT results. For the larger 460-row scan, provenance and limitations, see `../laptop_expanded/README.md`.
