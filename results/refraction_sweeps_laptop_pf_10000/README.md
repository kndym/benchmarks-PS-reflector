# Laptop NK=10000 PF shard

Completed 20/20 PF evaluations. All rows are already included in `../refraction_sweeps_expanded/costs.csv`; **do not concatenate them a second time**. This directory preserves the raw checkpoint diagnostics and execution log for provenance.

NK=10000, source epsilon=0.00125; 20 geometrically spaced PF epsilons from 0.0125 down to 0.000125. The source checkpoint was copied from the original laptop refractor solve, not recomputed. Uses the original uniform patches, kappa=0.6, NK-dependent target, squared Euclidean 3D PF cost, max_iter=10000, potential-change tolerance=1e-9, chunk_size=512, cache budget=100000000 entries, four independent SciPy reduction threads, descending-epsilon warm starts. Exact arguments are in `parameters.json`.

All 20 evaluations passed the existing PF stopping test. The saved source map inherits the original capped source approximation. This is not an independent fixed-target or refined-source experiment.
