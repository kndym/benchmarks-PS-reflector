# Entropic-surface sweep snapshot

**PARTIAL: 63/120 evaluations.**

This immutable snapshot contains only CSV-completed cases and their validated solve, PF, and residual checkpoints. The live worker continues in the source directory named in `snapshot.json`. Comparison images were rendered from this snapshot's CSV, not copied from the changing live output. Target fixed at 1600 points; transport/self caps 2000/2000; PF-to-target EOT uses squared Euclidean 3D cost and epsilon 1/320.

The hard comparison is `../../refraction_sweeps_extended_high_iter/`. The snapshot records hashes and validation counts. Incomplete experiments are not final results.
