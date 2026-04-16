<img src="images/romsoclogo-logo.png" alt="ROMSOC logo" width="150"/>

# Point Source Far-Field Refractor — Entropic Optimal Transport

Solve the **far-field refractor problem** via Sinkhorn divergence (entropic optimal transport).  
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

**Generate all 16 source×target density-pair results:**

```bash
python scripts/generate_results_all_pairs.py
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
│   └── qmc.py                    # Quasi-Monte Carlo point cloud loading
│
├── scripts/
│   ├── generate_results.py       # Refractor benchmark, NK=1600 (default)
│   ├── generate_results_all_pairs.py  # All 16 density-pair combinations
│   ├── refractor_sinkhorn.py     # Standalone refractor Sinkhorn (self-contained)
│   ├── visualize.py              # Plot utilities
│   │
│   ├── generate_results_reflector.py  # Reflector benchmark (κ=1.0, SquareToCircle)
│   ├── compare_reflector.py      # Compare C++ vs Python reflector outputs
│   ├── compare_nk_reflector.py   # Per-NK resolution comparison
│   ├── run_compare_reflector.py  # Runs Python reflector for C++ comparison
│   └── run_fast_reflector.py     # Fast reflector runner (NK=381)
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

The solver runs on **spherical patches** (upper hemisphere) with κ = 0.6:

1. **Generate QMC cloud** — Halton quasi-Monte Carlo points on the source and target patches
2. **Evaluate densities** — P(x) on source patch, Q(y) on target patch
3. **Sinkhorn divergence** — Iterative log-domain solver; produces Kantorovich potentials f, g
4. **Build refractor surface** — R = exp(f), Ref = 2·x·R
5. **C-transforms** — Verify optimality: gc ≈ f, fc ≈ g
6. **Push-forward** — Argmin OT map x_i → y_{j*(i)}; project via north-pole stereographic

Default patch geometry (matching C++ reference):
- Source: θ ∈ [π/12, π/3], φ ∈ [π/12, π/4]
- Target: θ ∈ [π/3, 5π/12], φ ∈ [π/10, π/5]

---

## Density Shapes

The `generate_results_all_pairs.py` script runs all 16 combinations of:

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
