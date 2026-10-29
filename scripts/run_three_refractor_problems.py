"""Run three isolated geometry/kappa variations with the previous sweep diagnostics."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from atomic_io import atomic_write_json
from result_paths import artifact_path, ensure_readme

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/runs/isolated_refractor_problems_oct_27_26"
CASES = {
    "source_only": ["--source-phi", "15", "55"],
    "target_only": ["--target-theta", "18", "45"],
    "kappa_only": ["--kappa", "0.8"],
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CASES, required=True)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument('--output-root', type=Path, default=OUTPUT)
    parser.add_argument('--experiments', nargs='+', type=int, choices=[1,2,3], default=[1,2,3])
    parser.add_argument('--nk', nargs='+', type=int)
    parser.add_argument('--epsilon',nargs='+',type=float)
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--chunk-size',type=int,default=256)
    parser.add_argument('--shared-optical-cache',type=Path)
    parser.add_argument('--refractor-max-iter', type=int, default=10000)
    parser.add_argument('--identity-max-iter', type=int, default=10000)
    parser.add_argument("--diagnostics-only", action="store_true",
                        help="Recover diagnostics while existing optical workers continue")
    args = parser.parse_args()
    folder = (ROOT / "results/archive/smoke_checks/repaired_problems_oct_28_26" if args.smoke else args.output_root.resolve()) / args.case
    ensure_readme(folder)
    logs = folder / "logs"
    ensure_readme(logs)
    common = ["--extended", "--experiments", *map(str,args.experiments), "--target-nk", "1600",
              "--refractor-max-iter", str(args.refractor_max_iter), "--identity-max-iter", str(args.identity_max_iter),
              "--workers", str(args.workers), "--chunk-size", str(args.chunk_size), *CASES[args.case]]
    if args.nk:
        common += ['--nk', *map(str,args.nk)]
    if args.epsilon:
        common += ['--epsilon',*map(str,args.epsilon)]
    if args.shared_optical_cache:
        common += ['--shared-optical-cache',str(args.shared_optical_cache.resolve())]
    if args.smoke:
        common += ["--nk", "16", "--fixed-nk", "16", "--target-nk", "16", "--epsilon", "0.01"]
    children, streams = {}, []
    state_path = artifact_path(folder, "worker_status.json", for_write=True)
    optical_workers = {}
    if args.diagnostics_only and state_path.exists():
        previous = json.loads(state_path.read_text())
        optical_workers = previous.get("optical_workers") or {
            mode: previous["children"][mode] for mode in ("hard", "entropic")}

    def status(state):
        atomic_write_json(state_path, dict(state=state, supervisor_pid=os.getpid(),
            updated_at=datetime.now(timezone.utc).isoformat(), case=args.case,
            optical_workers=optical_workers,
            children={role: dict(pid=child.pid, exit_code=child.poll()) for role, child in children.items()}))

    def start(role, command):
        stream = (logs / (role + ".log")).open("a", encoding="utf-8")
        streams.append(stream)
        env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
        children[role] = subprocess.Popen([sys.executable, "-u", "-B", *command], cwd=ROOT,
            stdin=subprocess.DEVNULL, stdout=stream, stderr=subprocess.STDOUT, env=env,
            creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        status("running")

    def fail(message):
        # Stop this supervisor's own children, preserving their checkpoints.
        for child in children.values():
            if child.poll() is None:
                child.terminate()
        for child in children.values():
            try:
                child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                child.kill()
                child.wait()
        status('failed')
        for stream in streams:
            stream.close()
        raise RuntimeError(message)

    for mode in (() if args.diagnostics_only else ("hard", "entropic")):
        start(mode, ["scripts/run_refraction_sweeps.py", *common, "--pushforward", mode,
                     "--output", str(folder / mode)])
    while not all(artifact_path(folder / mode, "parameters.json").exists() for mode in ("hard", "entropic")):
        if any(child.poll() is not None and child.returncode != 0 for child in children.values()):
            fail("Sweep failed before writing parameters; inspect logs")
        time.sleep(1)
    start("debiased", ["scripts/evaluate_saved_refraction_divergence.py", "--hard", str(folder / "hard"),
        "--entropic", str(folder / "entropic"), "--output", str(folder / "debiased"),
        "--workers", "1", "--cache-entries", "4000000", "--watch"])
    start("plots", ["scripts/plot_refraction_pushforward_comparison.py", str(folder / "entropic"),
                    "--baseline", str(folder / "hard"), "--watch"])
    # The evaluator refreshes comparison plots after every completed pair.
    while any(child.poll() is None for child in children.values()):
        if any(child.poll() is not None and child.returncode != 0 for child in children.values()):
            fail("Worker failed; checkpoints preserved; inspect logs and resume this case")
        status("running")
        time.sleep(2 if args.smoke else 30)
    if any(child.returncode != 0 for child in children.values()):
        fail('Worker failed; inspect logs before resuming')
    status('complete')
    for stream in streams:
        stream.close()


if __name__ == "__main__":
    main()
