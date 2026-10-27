# Point-source far-field refractor

Design a refractor that redirects a source light distribution toward a prescribed
target on the unit sphere, using entropic optimal transport and Sinkhorn solves.
The optical cost is `c(x,y) = -log(1 - kappa * (x dot y))`, with `kappa=0.6`
for the refractor and `kappa=1` for the reflector reference problem.

## Start here

| What you need | Open |
| --- | --- |
| Regular entropic error, without debiasing | [Complete cross-EOT plots and numerics](results/published/regular_entropic_error_oct_27_26/) |
| Debiased hard-versus-smooth comparison | [Complete Sinkhorn-divergence capture](results/published/debiased_comparison_complete_oct_26_26/) |
| Complete optical run and validation | [Smooth-refraction capture](results/published/smooth_refraction_complete_oct_26_26/) |
| All datasets and metric definitions | [Results catalog](results/README.md) |
| Explore a saved refractor | [Primary notebook](notebooks/Benchmark_refraction.ipynb) |
| Current research write-up | [LaTeX project](geometric_optics_progress/) |
| Folder naming and output conventions | [Repository layout guide](docs/repository_layout.md) |

The completed extended comparison has **120 cases per map**, 40 per experiment.
The hard map has lower debiased Sinkhorn divergence in all 120 paired cases.
Regular cross-EOT is a different metric; the galleries label it separately.

## Quickstart

Run commands from the repository root:

```bash
pip install numpy scipy matplotlib
python scripts/generate_results.py
```

The generator saves the uniform refractor NPZ under
`results/examples/uniform_refraction_<mon_dd_yy>/data/`. Open
[Benchmark_refraction.ipynb](notebooks/Benchmark_refraction.ipynb) to load a saved
bundle and inspect its surface and push-forward. Notebook visualization also
uses Jupyter and Plotly. Existing examples are indexed in
[results/examples](results/examples/).

```bash
# All sixteen uniform/Gaussian/donut/cross source-target combinations
python scripts/helper_all_pairs.py

# Regenerate the published regular-EOT plots from their bundled CSV
python scripts/plot_saved_entropic_error.py
```

New density-pair bundles go to
`results/examples/density_pairs_<mon_dd_yy>/data/`. Plot generation above reads
saved numerics; it does not repeat the optical solves.

## Repository map

```text
refracter/                 Shared numerical and optical library
scripts/                   Runners, plotting, validation, and maintenance
notebooks/                 Notebook source; generated exports live in results
BenchmarkCode/             C++ reference source and headers
geometric_optics_progress/ Active LaTeX project and its required images
results/
  published/               Validated captures and readable plot galleries
  runs/                    Canonical runs and reusable checkpoints
  examples/                NPZ examples and sensitivity controls
  reference_runs/          New standalone C++/Python reflector output
  notebook_outputs/        Dated notebook exports; created on save
  local_runtime/           Ignored worker state, logs, and process IDs
  archive/                 Historical, superseded, and partial numerics
docs/                      Layout guide, research notes, references, and credits
packaging/                 Docker/build environment
archive/                   Inactive figures, branding, binaries, and exports
tests/                     Existing numerical and workflow checks
```

Folders explain their contents in READMEs. Numerical runs separate `data/`,
`plots/`, `metadata/`, `checkpoints/`, and `diagnostics/`; files are not written
beside notebook source or scattered through the repository root. See the
[layout guide](docs/repository_layout.md) for the date suffix and five-file
parent-directory conventions.

## Core pipeline

`refracter/` supplies the shared solver, sampling, geometry, and push-forward
operations. [generate_results.py](scripts/generate_results.py) constructs the
notebook-compatible NPZ schema. [run_refraction_sweeps.py](scripts/run_refraction_sweeps.py)
varies source discretization and refractor regularization; independent transport
evaluations measure the pushed distribution. Use each runner's `--help` for its
grids, caps, and tolerances.

[scripts/result_paths.py](scripts/result_paths.py) resolves migrated paths and
routes new artifacts into their typed subfolders. The
[migration catalog](results/path_catalog.json) records old and new locations.

The standalone reflector entry point is
[benchmark_reflector.py](benchmark_reflector.py). Its default output is a dated
run under `results/reference_runs/`; an explicit `--output_dir` is honored.
Container instructions are in [packaging](packaging/).

Detailed historical commands, API examples, and numerical background are retained
in the [pre-organization guide](docs/research_notes/pre_organization_repository_guide_oct_27_26.md).
