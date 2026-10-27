# Notebooks

Start with [Benchmark_refraction.ipynb](Benchmark_refraction.ipynb) for the active refractor workflow.

| Notebook | Purpose |
| --- | --- |
| [Benchmark_refraction.ipynb](Benchmark_refraction.ipynb) | Load saved refractor NPZ data, inspect surfaces, and evaluate push-forward behavior. |
| [Benchmark_Python.ipynb](Benchmark_Python.ipynb) | Run and inspect the Python reflector reference pipeline. |
| [Benchmark_fast.ipynb](Benchmark_fast.ipynb) | Visualize a saved fast reflector NPZ archive. |
| [Benchmark.ipynb](Benchmark.ipynb) | Build and inspect the original C++/MKL benchmark. |

Each notebook locates the repository from its current working directory. Static plots and interactive HTML are saved under `results/notebook_outputs/<notebook_name>_<mon_dd_yy>/plots/` and `interactive/`, with README indexes refreshed on every save. They are not written beside these notebooks.

The refraction notebook reads the [uniform refraction example](../results/examples/uniform_refraction_oct_27_26/README.md); the fast reflector notebook reads the [archived reflector NPZ bundle](../results/archive/reference_outputs/reflector_npz_bundles_oct_27_26/README.md). Both loaders use the central path resolver. The C++ notebook selects the newest matching dated run under `results/reference_runs/`, then reads its `data/` leaf.

[Archive](archive/) contains historical notebooks and the loose figure exports collected during the October 27, 2026 cleanup. Archived artifacts are reference material, not fresh numerical runs.

<!-- directory-index:start -->

## Directory guide

[Parent directory](../README.md)

| Folder | Contents |
| --- | --- |
| [archive/](archive/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |

<!-- directory-index:end -->
