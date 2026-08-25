# Generate the standard uniform-to-uniform refraction result.
# Usage: python scripts/generate_results.py [NK]

from __future__ import annotations

import os
import sys


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from refracter.distributions import P_refraction_patch, Q_refraction_patch
from refracter.refraction_pipeline import (
    run_refraction_case,
    save_refraction_result,
    save_surface_figure,
)


def main(argv=None) -> None:
    nk = int(argv[0]) if argv else int(sys.argv[1]) if len(sys.argv) > 1 else 1600
    print(f"=== Refraction benchmark (kappa=0.6, NK={nk}) ===")

    result = run_refraction_case(
        P_refraction_patch,
        Q_refraction_patch,
        nk=nk,
        use_identity=False,
        verbose=True,
    )

    results_dir = os.path.join(REPO_ROOT, "results")
    output_path = os.path.join(results_dir, f"results_refraction_NK{nk}.npz")
    save_refraction_result(result, output_path)
    print(f"Saved {output_path}")

    figure_path = os.path.join(REPO_ROOT, "figures", f"fig_refractor_3d_NK{nk}.png")
    save_surface_figure(result, figure_path)
    print(f"Saved {figure_path}")


if __name__ == "__main__":
    main()
