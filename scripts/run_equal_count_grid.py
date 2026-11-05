"""Run this PC's share of the equal-count grid, a few cells at a time.

A new cell starts only when fewer than --slots solver processes are running on
this machine (counting any other refraction jobs), so it shares the CPU politely.

  python scripts/run_equal_count_grid.py RUN_DIR [--slots 5] [--extra EXTRA_JOBS.json]
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import threading
import time

import psutil

ROOT = Path(__file__).resolve().parents[1]


def busy():
    count = 0
    for proc in psutil.process_iter(["cmdline"]):
        line = " ".join(proc.info["cmdline"] or [])
        if "run_refraction_sweeps.py" in line or "evaluate_saved_refraction_divergence.py" in line:
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--slots", type=int, default=5)
    parser.add_argument("--extra", type=Path, help="Extra jobs (with output_root) to run first")
    args = parser.parse_args()
    run = args.run.resolve()
    plan = json.loads((run / "grid_plan.json").read_text())
    jobs = json.loads(args.extra.read_text()) if args.extra else []
    jobs += [dict(c, output_root=str(run / "grid" / "local")) for c in plan["cells"] if c["host"] == "local"]
    log = run / "local_runner.jsonl"
    threads = []

    def run_job(job):
        start = time.time()
        code = subprocess.run([sys.executable, "-B", "scripts/run_equal_count_cell.py", "--nk", str(job["nk"]),
                               "--epsilon", repr(job["epsilon"]), "--output-root", job["output_root"],
                               "--marginal-tolerance", repr(job["marginal_tolerance"]),
                               "--refractor-tolerance", repr(job["refractor_tolerance"]),
                               "--pf-tolerance", repr(job["pf_tolerance"])],
                              cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT).returncode
        with log.open("a") as stream:
            stream.write(json.dumps(dict(nk=job["nk"], epsilon=job["epsilon"], exit_code=code,
                                         minutes=(time.time() - start) / 60)) + "\n")

    for job in jobs:
        while busy() >= args.slots:
            time.sleep(20)
        thread = threading.Thread(target=run_job, args=(job,))
        thread.start()
        threads.append(thread)
        time.sleep(15)  # let the new cell's solver process appear before counting again
    for thread in threads:
        thread.join()
    with log.open("a") as stream:
        stream.write('"done"\n')


if __name__ == "__main__":
    main()
