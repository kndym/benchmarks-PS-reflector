# Repository layout

The repository separates code, current numerical work, published captures,
examples, and inactive exports. Read a folder's README before using its files.

| Location | Purpose |
| --- | --- |
| `refracter/` | Shared numerical/optical implementation; unchanged by the organization pass. |
| `tests/` | Existing checks; unchanged by the organization pass. |
| `scripts/` | Runners, plotting, validation, migration, and maintenance tools. |
| `BenchmarkCode/` | C++ reference source and headers. Generated run output belongs in results. |
| `notebooks/` | Notebook entry points, with generated exports separated from notebook source. |
| `geometric_optics_progress/` | Active LaTeX source and the images needed to compile it. |
| `results/runs/` | Canonical working numerical runs and checkpoints. |
| `results/published/` | Validated complete captures and readable plot galleries. |
| `results/examples/` | Curated NPZ examples, density pairs, and sensitivity controls. |
| `results/reference_runs/` | New standalone C++/Python reflector benchmark output. |
| `results/local_runtime/` | Ignored local supervision state and process logs. |
| `results/archive/` | Historical studies, partial captures, shards, and superseded numerical output. |
| `archive/` | Inactive repository-level images, exports, scratch files, and old binaries. |
| `docs/` | References, research notes, credits, and this guide. |
| `packaging/` | Container/build environment files. |

## Descriptive names and dates

Use lowercase descriptive names, followed by a date suffix such as
`regular_entropic_error_oct_27_26`. The label should identify the map, metric,
or experiment without opening the folder. A capture time can precede the date
when several updates occur on one day, as in
`writeup_comparison_223420_oct_26_26`.

Recorded computation or capture dates are retained. Some older bundles have
only a saved/publication date. Undated files organized during this pass end in
`oct_27_26`, and their README identifies that as the organization date. A new
folder name does not imply that the numerical computation was repeated.

## Inside a numerical run

```text
descriptive_run_oct_27_26/
  README.md
  data/                    # CSV summaries, NPZ/NPY bundles, numerical TXT
  plots/                   # PNG/PDF/SVG figures
  metadata/                # Parameters, validation, provenance, hash manifests
  checkpoints/
    solutions/             # Optical solve NPZs
    pushforward/           # Cross-EOT evaluations
    self_costs/            # Reusable self-EOT cache, when applicable
  diagnostics/
    marginals/             # Optical marginal residuals
    divergence/            # Debiased evaluations and their three EOT terms
  logs/                    # Local logs/errors/PIDs, excluded from GitHub
  interactive/             # HTML visualizations, when present
```

Only the subfolders a dataset needs are created. Hard/entropic comparison
bundles may have separate `hard/` and `entropic/` branches with this same layout.
An optical solution, a regularized cross-EOT objective, and a debiased
divergence remain distinct files and quantities.

## Keep parent folders small

Every folder should explain its contents in a README. A directory containing
subdirectories should have no more than five regular files beside them.
Numerical tables, arrays, figures, metadata, and logs belong in their typed
subfolders rather than accumulating beside category or run directories.
Script source directories are the exception: many source files are acceptable.
Leaf directories may contain the naturally related files of a run.

Do not archive assets that are still required to build an active artifact.
For example, `geometric_optics_progress/images/` remains active even though
the old root `images/` collection moved to the archive.

## Read and write through shared paths

[`scripts/result_paths.py`](../scripts/result_paths.py) centralizes locations:

- `resolve_path(path)` translates known old directory/file names through
  [`results/path_catalog.json`](../results/path_catalog.json).
- `artifact_path(bundle, name)` finds a file in its typed folder, with legacy
  flat-bundle read fallback. Set `for_write=True` to create its canonical
  destination and missing explanatory READMEs.
- `artifact_glob(bundle, pattern)` searches typed checkpoint folders as well as
  older flat fixtures.
- `self_cost_dir(bundle)` locates the dedicated self-EOT cache.
- `dated_run_dir(category, label)` creates a descriptive dated destination.

New writers should use these helpers rather than concatenating a run root
with `solve_*.npz`, `costs.csv`, or figure filenames. For example:

```python
from scripts.result_paths import artifact_path, dated_run_dir

run = dated_run_dir("runs", "smooth_refraction")
table = artifact_path(run, "costs.csv", for_write=True)
figure = artifact_path(run, "comparison.png", for_write=True)
```

The path catalog records the migration and original hashes. Numerical archives
and images were moved byte-for-byte; embedded paths and manifest entries were
updated to their new locations. Preserve original provenance and use each
manifest's documented `hash_base` when verifying its hashes.

When adding results, update the dataset README and the relevant category's
navigation. Put partial captures in the archive, reserve published galleries
for clearly described validated captures, and keep process state out of both.

On Windows, enable long checkpoint paths before checkout with
`git config --global core.longpaths true`. This checkout already enables it
locally. Artifact folders use `.gitattributes` to preserve exact bytes across
platforms, so recorded SHA256 values remain valid after cloning.
