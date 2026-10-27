# Scripts

Python and PowerShell entry points for generation, plotting, validation, and maintenance. Shared output routing is implemented in result_paths.py.

| Task | Entry point |
| --- | --- |
| Uniform refraction example | `generate_results.py` |
| Sixteen density pairs | `helper_all_pairs.py` |
| Parameter sweeps | `run_refraction_sweeps.py` |
| Debiased evaluation of saved maps | `evaluate_saved_refraction_divergence.py` |
| Regular entropic-error plots | `plot_saved_entropic_error.py` |
| Update the paper and figures | `update_writeup_results.py` |
| Optical/debiased capture validation | `snapshot_refraction_sweep.py`, `snapshot_refraction_divergence.py` |
| Refresh folder READMEs and check links/layout | `check_repository_layout.py --refresh-readmes` |
| Verify relocated artifact hashes | `repair_result_manifests.py` |
| Regenerate C++ sampling headers | `generate_cpp_pointclouds.py` |

Run these from the repository root with `python scripts/NAME.py`. Writers use
dated run directories and typed subfolders through `result_paths.py`. Existing
legacy path arguments resolve through `results/path_catalog.json`; new output
folders receive READMEs and update their parent directory guides automatically.
Use a manifest repair with `--apply` only after an intentional, documented
relocation or documentation change; ordinary verification never edits arrays.

<!-- directory-index:start -->

## Directory guide

[Parent directory](../README.md)

| Folder | Contents |
| --- | --- |
| [archived/](archived/README.md) | Python and PowerShell entry points for generation, plotting, validation, and maintenance. Shared output routing is implemented in result_paths.py. |

<!-- directory-index:end -->
