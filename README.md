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
Experiments 1–3 plot its approximate `total_cost`: joint
ε = 1/(8 floor(sqrt(NK))), fixed ε with varying NK, and fixed NK with varying ε.
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
