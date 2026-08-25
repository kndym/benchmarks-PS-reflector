# Generate results for every source/target pair in DENSITY_NAMES.
# Usage: python scripts/generate_results_all_pairs.py [NK]

from __future__ import annotations

import os
import sys


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from refracter.refraction_pipeline import (
    DENSITY_NAMES,
    SRC_PHI,
    SRC_THETA,
    TGT_PHI,
    TGT_THETA,
    make_patch_density,
    run_refraction_case,
    save_refraction_result,
)


def main(argv=None) -> None:
    nk = int(argv[0]) if argv else int(sys.argv[1]) if len(sys.argv) > 1 else 1600
    density_names = DENSITY_NAMES
    results_dir = os.path.join(REPO_ROOT, "results")
    total = len(density_names) ** 2

    print(f"=== Refraction density sweep (NK={nk}, {total} pairs) ===")
    print(f"  Source theta/phi: {SRC_THETA} / {SRC_PHI}")
    print(f"  Target theta/phi: {TGT_THETA} / {TGT_PHI}")

    for source_index, source_name in enumerate(density_names):
        for target_index, target_name in enumerate(density_names):
            pair_number = source_index * len(density_names) + target_index + 1
            print(f"\n[{pair_number:2d}/{total}] {source_name} -> {target_name}")
            result = run_refraction_case(
                make_patch_density(source_name, source=True),
                make_patch_density(target_name, source=False),
                nk=nk,
                use_identity=True,
                verbose=True,
            )
            output_path = os.path.join(
                results_dir,
                f"results_refraction_{source_name}_{target_name}_NK{nk}.npz",
            )
            save_refraction_result(
                result,
                output_path,
                source_density=source_name,
                target_density=target_name,
            )
            print(f"  Saved {output_path}")

    print(f"\nDone. {total} result files written to {results_dir}/")


if __name__ == "__main__":
    main()
