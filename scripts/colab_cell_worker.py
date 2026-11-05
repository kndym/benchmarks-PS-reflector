"""Run a list of equal-count cells sequentially on a remote (Colab) VM.

Usage: python scripts/colab_cell_worker.py JOBS.json OUTPUT_ROOT
JOBS.json holds [{"nk":..,"epsilon":..,"marginal_tolerance":..,"refractor_tolerance":..,
"pf_tolerance":..}, ...]. Finished cells are skipped, so the worker can be restarted.
"""
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    jobs = json.loads(Path(sys.argv[1]).read_text())
    output = Path(sys.argv[2])
    output.mkdir(parents=True, exist_ok=True)
    status = output / "worker_status.jsonl"
    for job in jobs:
        root = output / job.get("subdir", "")
        command = [sys.executable, "-B", "scripts/run_equal_count_cell.py", "--nk", str(job["nk"]),
                   "--epsilon", repr(job["epsilon"]), "--output-root", str(root),
                   "--marginal-tolerance", repr(job.get("marginal_tolerance", 1e-7)),
                   "--refractor-tolerance", repr(job.get("refractor_tolerance", 1e-5)),
                   "--pf-tolerance", repr(job.get("pf_tolerance", 1e-9))]
        code = subprocess.run(command, cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT).returncode
        with status.open("a") as stream:
            stream.write(json.dumps(dict(job, exit_code=code)) + "\n")
    with status.open("a") as stream:
        stream.write('"done"\n')


if __name__ == "__main__":
    main()
