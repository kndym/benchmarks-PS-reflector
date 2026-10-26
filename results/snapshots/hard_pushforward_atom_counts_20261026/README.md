# Hard push-forward atom weights

120 saved cases from the hard-transform 2000/2000-cap sweep; no optical solves were repeated.

At each distinct occupied push-forward location, sum the weights of all source points mapped there. Count the resulting atoms whose mass differs from `1/NK` (source NK). For these uniform-source runs, this is exactly the number of locations receiving at least two source points. Empty target locations carry zero mass and are excluded from this count.

![Atom counts](experiments_1_2_3.png)

![Log-log atom counts](experiments_1_2_3_log.png)

[Per-case counts](atom_counts.csv) include occupied/singleton counts, empty target locations, source points in merged atoms, maximum multiplicity, and unequal fraction of occupied support. [Validation](validation.json) records uniform-source and mass-conservation checks and checkpoint hashes.
