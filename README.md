<img src="images/romsoclogo-logo.png" alt="ROMSOC logo"  width="150"/>

# Benchmarks for Point Source to Far Field Reflector Problem
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5171811.svg)](https://doi.org/10.5281/zenodo.5171811) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ROMSOC/benchmarks-ps-reflector/HEAD?labpath=notebooks/Benchmark.ipynb)

This repository contains benchmarks for the **Point Source to Far Field reflector problem**, solved via **Entropic Optimal Transport** (Sinkhorn divergence). It includes both a C++ reference implementation and a full Python rewrite.

**Reference:** Benamou, Ijzerman, Rukhaia — *An Entropic Optimal Transport Numerical Approach to the Reflector Problem* (2020)

---

## Repository Structure

```
benchmarks-ps-reflector/
├── reflector/              # Python module (core library)
│   ├── sinkhorn.py         # Sinkhorn divergence algorithm
│   ├── build.py            # Reflector surface construction and c-transforms
│   ├── distributions.py    # Benchmark definitions (SquareToCircle, etc.)
│   ├── pushforward.py      # Ray tracing and push-forward computation
│   ├── qmc.py              # Quasi-Monte Carlo point cloud loading
│   └── cost.py             # Cost matrix computation
├── BenchmarkCode/          # C++ reference implementation
│   ├── main.cpp            # Full Sinkhorn pipeline
│   ├── main_fast.cpp       # Optimised variant
│   ├── Benchmarks/         # Test case headers
│   ├── PushForward/        # Reflection computation and point clouds
│   ├── QuasiMonteCarlo/    # QMC grid containers
│   └── SmallGrid/          # Coarse grid for warm-start initialisation
├── notebooks/              # Jupyter notebooks
│   ├── Benchmark.ipynb     # Main walkthrough (mirrors C++ pipeline)
│   ├── Benchmark_Python.ipynb  # Python-focused exploration
│   ├── Benchmark_fast.ipynb    # Fast implementation benchmarks
│   └── Benchmark_old.ipynb     # Legacy notebook
├── scripts/                # Utility scripts
│   ├── compare.py          # Compare C++ vs Python results
│   ├── compare_nk.py       # Per-resolution (NK) comparison
│   ├── run_compare.py      # Driver for comparisons
│   ├── run_fast.py         # Fast implementation runner
│   ├── generate_results.py # Batch result generation
│   └── visualize.py        # Visualisation utilities
├── figures/                # Generated plots and visualisations
├── results/                # Computed output data (NPZ, TXT)
├── images/                 # Repository images
├── benchmark.py            # Main Python entry point
├── Dockerfile              # Docker environment
└── paper.txt               # Reference paper notes
```

---

## Python Implementation

### Quickstart

```python
from benchmark import run_benchmark

results = run_benchmark(benchmark='SquareToCircle', output_dir='my_output')
```

Or from the command line:

```bash
python benchmark.py --benchmark SquareToCircle --output_dir my_output
```

Available benchmarks: `SquareToCircle`, `SquareToTwoGaussSide`

### Pipeline

The Python rewrite (`benchmark.py` + `reflector/` module) replicates the full C++ pipeline:

1. Load QMC point clouds
2. Evaluate source/target density functions P and Q
3. Run the Sinkhorn divergence algorithm (with warm-start from coarse grid)
4. Build the reflector surface (R, Ref)
5. Compute c-transforms (gc, fc)
6. Build a regular grid and compute the reflector on it
7. Ray-trace the push-forward
8. Save all outputs (compatible with C++ output format)

### Python Dependencies

```bash
pip install numpy scipy
```

---

## C++ Implementation

The C++ code requires **Intel's MKL** (Math Kernel Library), available via [Intel oneAPI](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html).

Edit the `Makefile` in `BenchmarkCode/` to adjust MKL paths for your installation. The C++ code is Linux-native (uses system file/folder creation commands).

---

## Running via Jupyter Notebook

### Binder (cloud, no setup required)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ROMSOC/benchmarks-ps-reflector/HEAD?labpath=notebooks/Benchmark.ipynb)

> Note: mybinder cloud computations are limited to 2 GB of RAM.

### Docker (local, includes Intel oneAPI)

Build the container:
```bash
docker build -t benchmarks-ps-reflector .
```

Run with Jupyter Lab:
```bash
docker run -it --rm -p 8888:8888 benchmarks-ps-reflector \
    jupyter-lab --ip=0.0.0.0 --port=8888 --allow-root
```

---

## Comparing C++ and Python Results

After generating outputs from both implementations, use the comparison scripts:

```bash
# Compare results for a specific NK resolution
python scripts/compare_nk.py

# Full comparison
python scripts/compare.py
```

---

## Disclaimer

In downloading this SOFTWARE you are deemed to have read and agreed to the following terms:
This SOFTWARE has been designed with an exclusive focus on civil applications. It is not to be used for any illegal, deceptive, misleading or unethical purpose or in any military applications. This includes ANY APPLICATION WHERE THE USE OF THE SOFTWARE MAY RESULT IN DEATH, PERSONAL INJURY OR SEVERE PHYSICAL OR ENVIRONMENTAL DAMAGE. Any redistribution of the software must retain this disclaimer. BY INSTALLING, COPYING, OR OTHERWISE USING THE SOFTWARE, YOU AGREE TO THE TERMS ABOVE. IF YOU DO NOT AGREE TO THESE TERMS, DO NOT INSTALL OR USE THE SOFTWARE.

---

## Acknowledgments

<img src="images/EU_Flag.png" alt="EU Flag"  width="150" height="100" />

The ROMSOC project has received funding from the European Union's Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie Grant Agreement No. 765374.
This repository reflects the views of the author(s) and does not necessarily reflect the views or policy of the European Commission. The REA cannot be held responsible for any use that may be made of the information this repository contains.
