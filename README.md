# Point Source Far-Field Refractor — Entropic Optimal Transport

Solve the **far-field refractor problem** via Sinkhorn-based entropic optimal transport.
Given source and target intensity distributions on the unit sphere, the solver finds the refractor surface that redirects light from the source to match the target.

**Cost function:** `c(x, y) = -log(1 - κ·(x·y))`, κ = 0.6 (refractor), κ = 1.0 (reflector).

---

## Quickstart

```bash
pip install numpy scipy matplotlib
```

**Generate refractor results (uniform source → uniform target):**

```bash
python scripts/generate_results.py
# Saves: results/results_refraction_NK1600.npz
```

The primary exploration notebook loads this NPZ directly:

```text
notebooks/Benchmark_refraction.ipynb
```

**Generate all 16 source×target density-pair results:**

```bash
python scripts/helper_all_pairs.py
# Densities: uniform, gaussian, donut, cross (4×4 = 16 combos)
# Saves: results/results_refraction_{src}_{tgt}_NK1600.npz
```

**Load and explore results in Python:**

```python
import numpy as np

r = np.load('results/results_refraction_NK1600.npz')

x, y   = r['x'], r['y']          # source / target point clouds (NK, 3)
f, g   = r['f'], r['g']          # Kantorovich potentials
R, Ref = r['R'], r['Ref']        # refractor scale and surface points
Y_Pushed_projected = r['Y_Pushed_projected']  # push-forward (u, v, density)
```

**Primary notebook:** `notebooks/Benchmark_refraction.ipynb`

---

## Repository Structure

```
benchmarks-PS-reflector/
├── refracter/                    # Python library (core solver)
│   ├── sinkhorn.py               # Sinkhorn divergence algorithm
│   ├── cost.py                   # Cost function c(x,y) = -log(1 - κ·(x·y))
│   ├── distributions.py          # Density functions and stereographic projections
│   ├── build.py                  # Surface construction and c-transforms
│   ├── pushforward.py            # Ray tracing and push-forward
│   └── qmc.py                    # Halton sampling and QMC point cloud loading
│
├── scripts/
│   ├── generate_results.py       # Refractor benchmark, NK=1600 (default)
│   ├── helper_all_pairs.py       # Optional 16 density-pair NPZ generation
│   └── archived/                 # Reflector and notebook-specific legacy utilities
│
├── notebooks/
│   ├── Benchmark_refraction.ipynb     # Primary refractor notebook
│   ├── Benchmark_Python.ipynb         # Reflector exploration
│   ├── Benchmark_fast.ipynb           # Fast reflector benchmarks
│   └── archive/                       # Legacy notebooks
│
├── results/                      # Computed output bundles (.npz)
├── figures/                      # Generated plots
├── BenchmarkCode/                # C++ reference implementation
├── benchmark_reflector.py        # Full reflector pipeline (CLI)
└── Dockerfile                    # Docker environment (includes Intel oneAPI)
```

---

## Refractor Pipeline

The primary pipeline is `refracter/` → `scripts/generate_results.py` →
`notebooks/Benchmark_refraction.ipynb`. Run the generator first; it writes the
common notebook-compatible result schema to `results/results_refraction_NK1600.npz`.
The all-pairs helper reuses the cold-start solver from `refracter/` and the
shared sampling, artifact construction, and NPZ writer from `generate_results.py`.
The notebook measures the push-forward against the target with the entropic OT
objective from equation 2.1, using the refraction cost and ε = 1/k_final.

Run the four uniform refractor parameter experiments:

```bash
python -B scripts/run_refraction_sweeps.py
```

This reuses the original multi-scale solver and hard pushforward builder.
Experiments 1–3 plot pushforward-to-target EOT with squared Euclidean cost,
holding evaluation epsilon at 1/320 (override with `--cost-epsilon`). They vary
the refractor parameters: joint ε = 1/(8 floor(sqrt(NK))), fixed ε with varying
NK, and fixed NK with varying ε. The refractor's approximate `total_cost` is
retained only in solve checkpoints and is not the plotted EOT objective.
Experiment 4 evaluates equation 2.1 with the notebook's current default,
squared Euclidean cost on 3D pushed points, independently varying the PF
evaluation ε. Each panel fixes the refractor ε (plus a joint-schedule panel);
each line fixes the PF ε. This is entropic OT, with no self-cost subtraction,
and is not the paper's projected-plane ray-tracing Wasserstein error.
Defaults: NK = 100, 225, 400, 900, 1600; fixed NK = 1600; fixed ε = 1/320;
refractor and PF ε grids = 1/80, 1/160, 1/320, 1/640.
Plots, `costs.csv`, `parameters.json`, and `ordering.json` are written to
`results/refraction_sweeps/`. Use `--help` to change these grids.
The original refractor iteration caps are retained; PF cost evaluations must
pass the existing convergence check. The smaller density display grid used by
the sweep runner does not affect the solver or pushed measure.

Target point count defaults to 1600. Source grids (`--nk`), refractor grids
(`--epsilon`), evaluation epsilon (`--cost-epsilon`), transport and self-solve
limits (`--refractor-max-iter`, `--identity-max-iter`), stopping tolerances
(`--refractor-tolerance`, `--pf-tolerance`), and `--chunk-size` can all be
set from the command line. Iteration limits default to 17 final updates;
0 removes a final-loop limit. Small-grid and continuation schedules are unchanged.
Residual JSON files record final iteration counts and potential changes as
well as marginal errors. A potential-change stopping tolerance alone does
not guarantee small marginal error. Changed refractor limits/tolerances use
distinct checkpoint names to prevent reusing capped solves.

For the comparable extended sweep with much higher limits:

```bash
python -B scripts/run_refraction_sweeps.py --extended --experiments 1 2 3 --target-nk 1600 --refractor-max-iter 2000 --identity-max-iter 2000 --workers 4 --output results/refraction_sweeps_extended_high_iter
```

To rerun the original non-expanded sweeps with only the source NK varying and
the target held at 1600 points:

```bash
python -B scripts/run_refraction_sweeps.py --target-nk 1600 --workers 4 --output results/refraction_sweeps
```

Fixed-target checkpoint names include the target count, so old varying-target
checkpoints are not reused. Rectangular transports use each marginal's own
support for warm-start smoothing and self-cost corrections; the inherited
equal-size solver behavior and stopping rules are preserved.

For denser, wider main experiments 1–3 (40 values per axis), run:

```bash
python -B scripts/run_refraction_sweeps.py --extended --experiments 1 2 3 --workers 4 --output results/refraction_sweeps_extended
```

NK spans 100–20,000; epsilon spans 0.0000625–0.025. Experiment 2 keeps
epsilon = 1/320, and experiment 3 keeps NK = 1600. The original stopping
rules are retained. Inspect the saved marginal residuals before interpreting
the pushed measures as coming from converged refractor solves. The CSV is saved after every
case, and completed solve checkpoints can be reused when resuming the command.
This is a main experiment output directory, separate from autonomous research.

Completed fixed-target extended results are saved for both final-loop caps:

| Transport / self cap | Results | Evaluations | Source marginal L1 range | Target marginal L1 range |
| --- | --- | --- | --- | --- |
| 17 / 17 | [Baseline plots and data](results/refraction_sweeps_extended/) | 120, all PF-converged | 2.51e-7–0.17678 | 0.000174–0.54792 |
| 2000 / 2000 | [Higher-iteration plots and data](results/refraction_sweeps_extended_high_iter/) | 120, all PF-converged | 5.24e-9–0.01472 | 0.000174–0.01472 |

Each directory includes the three individual plots, a combined plot, `costs.csv`,
parameters, solve/PF checkpoints, residual diagnostics, and a validation summary.
Both runs use the same 40-value grids, target of 1600 points, and squared-Euclidean
3D PF-to-target EOT evaluation at epsilon 1/320. The baseline predates explicit
cap fields in its parameters file and uses the inherited 17-update limits.
In the 2000-cap run, all 120 transport and both sets of 120 self solves met the
1e-5 potential-change tolerance; none hit the cap. The longest final loops used
659 transport, 551 source-self, and 15 target-self updates. PF convergence and
potential-change stopping remain distinct from marginal feasibility.

The expanded bounded experiment uses 20 logarithmic NK values from 100 to
10,000 and 20 epsilon values from 0.000125 to 0.0125 (100x ranges):

```bash
python -B scripts/run_refraction_sweeps.py --expanded --joint-only --pf-warm-start --pf-max-iter 10000 --workers 4 --output results/refraction_sweeps_expanded
```

Experiment 4 then has 400 evaluations on the original coupled refractor
schedule, with 20 PF-epsilon lines. PF warm starts reuse the preceding larger
epsilon's potentials; the objective and stopping tolerance are unchanged.
Optional reduction workers split independent row/column outputs and evaluate
each using the same SciPy `logsumexp`, rather than changing the solver.
Raw refractor marginal residuals are saved separately as `residual_*.json`;
they diagnose the inherited iteration caps without changing them.
Autonomous side investigations are isolated under the locally Git-excluded
`.research/refraction/` folder.

The solver runs on **spherical patches** (upper hemisphere) with κ = 0.6:

1. **Generate QMC cloud** — Halton quasi-Monte Carlo points on the source and target patches
2. **Evaluate densities** — P(x) on source patch, Q(y) on target patch
3. **Sinkhorn solve** — Iterative log-domain solver; produces Kantorovich potentials f, g
4. **Build refractor surface** — R = exp(f), Ref = 2·x·R
5. **C-transforms** — Verify optimality: gc ≈ f, fc ≈ g
6. **Push-forward** — Argmin OT map x_i → y_{j*(i)}; project via north-pole stereographic

Default patch geometry (matching C++ reference):
- Source: θ ∈ [π/12, π/3], φ ∈ [π/12, π/4]
- Target: θ ∈ [π/10, π/5], φ ∈ [π/10, π/5]

---

## Density Shapes

The `helper_all_pairs.py` script runs all 16 combinations of:

| Name | Description |
|------|-------------|
| `uniform` | Flat weight across the entire patch |
| `gaussian` | Isotropic Gaussian centred on the patch |
| `donut` | Soft annulus peaked at 0.5·d_max from centre |
| `cross` | Four Gaussians at N/S/E/W of patch centre |

---

## Using the Library

```python
from refracter.cost import set_kappa
set_kappa(0.6)   # must be called before other imports that cache cost computations

from refracter.sinkhorn import sinkhorn_step, sinkhorn_identity_f_step, sinkhorn_identity_g_step
from refracter.build import c_transform_gc, c_transform_fc
from refracter.distributions import stereo_north, make_patch_gaussian
```

**κ = 0.6** (refractor, default), **κ = 1.0** (reflector).

---

## Reflector Benchmark (C++ replication)

For the reflector case (κ = 1.0, SquareToCircle / SquareToTwoGaussSide):

```bash
python benchmark_reflector.py --benchmark SquareToCircle --output_dir my_output
```

Compares against the C++ reference in `BenchmarkCode/`. Requires Intel MKL for the C++ code.
