"""Run one equal-count (source N = target N) grid cell: hard + entropic maps and debiased terms.

Each cell gets its own folder so cells can run on any machine (local or Colab)
and be merged afterwards. Writes cell.json with costs, timings and settings.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def cell_name(nk, epsilon, marginal_tolerance, refractor_tolerance):
    return f"NK{nk}_eps{epsilon:.6g}_m{marginal_tolerance:g}_t{refractor_tolerance:g}"


def read_cost(folder):
    path = next(folder.rglob("costs.csv"))
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return rows[-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nk", type=int, required=True)
    parser.add_argument("--epsilon", type=float, required=True)
    parser.add_argument("--kappa", type=float, default=0.8)
    parser.add_argument("--marginal-tolerance", type=float, default=1e-7)
    parser.add_argument("--refractor-tolerance", type=float, default=1e-5)
    parser.add_argument("--pf-tolerance", type=float, default=1e-9)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--chunk-size", type=int, default=512)
    parser.add_argument("--cost-cache-mib", type=int, default=0)
    parser.add_argument("--no-divergence", action="store_true")
    args = parser.parse_args()
    folder = args.output_root.resolve() / cell_name(args.nk, args.epsilon, args.marginal_tolerance,
                                                    args.refractor_tolerance)
    done = folder / "cell.json"
    if done.exists():
        print(f"Already complete: {done}")
        return
    folder.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    common = ["--extended", "--experiments", "3", "--nk", str(args.nk), "--fixed-nk", str(args.nk),
              "--target-nk", str(args.nk), "--epsilon", repr(args.epsilon), "--kappa", str(args.kappa),
              "--refractor-max-iter", "10000", "--identity-max-iter", "10000",
              "--identity-refinement-budget", "190000", "--transport-refinement-budget", "190000",
              "--optical-marginal-tolerance", repr(args.marginal_tolerance),
              "--refractor-tolerance", repr(args.refractor_tolerance),
              "--pf-tolerance", repr(args.pf_tolerance),
              "--workers", str(args.workers), "--chunk-size", str(args.chunk_size),
              "--cost-cache-mib", str(args.cost_cache_mib),
              "--shared-optical-cache", str(folder / "shared_optical")]
    timings = {}
    for mode in ("hard", "entropic"):
        start = time.time()
        with (folder / f"{mode}.log").open("w", encoding="utf-8") as log:
            subprocess.run([sys.executable, "-u", "-B", "scripts/run_refraction_sweeps.py", *common,
                            "--pushforward", mode, "--output", str(folder / mode)],
                           cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        timings[mode + "_seconds"] = time.time() - start
    result = dict(nk=args.nk, target_nk=args.nk, epsilon=args.epsilon, kappa=args.kappa,
                  marginal_tolerance=args.marginal_tolerance,
                  refractor_tolerance=args.refractor_tolerance, pf_tolerance=args.pf_tolerance,
                  hard_cost=float(read_cost(folder / "hard")["cost"]),
                  entropic_cost=float(read_cost(folder / "entropic")["cost"]))
    residual = json.loads(next((folder / "hard").rglob("residual_*.json")).read_text())
    result.update({k: residual[k] for k in ("final_iterations", "source_identity_iterations",
                                            "target_identity_iterations")})
    if not args.no_divergence:
        start = time.time()
        with (folder / "debiased.log").open("w", encoding="utf-8") as log:
            subprocess.run([sys.executable, "-u", "-B", "scripts/evaluate_saved_refraction_divergence.py",
                            "--hard", str(folder / "hard"), "--entropic", str(folder / "entropic"),
                            "--output", str(folder / "debiased"), "--workers", str(args.workers)],
                           cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        timings["debiased_seconds"] = time.time() - start
        for mode in ("hard", "entropic"):
            row = read_cost(folder / "debiased" / mode)
            for name, value in row.items():
                try:
                    result[f"{mode}_debiased_{name}"] = float(value)
                except (TypeError, ValueError):
                    pass
    result.update(timings)
    result["total_seconds"] = sum(timings.values())
    done.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
