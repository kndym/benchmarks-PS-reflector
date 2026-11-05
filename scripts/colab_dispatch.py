"""Send equal-count cells to Colab CPU sessions and pull finished results back.

  python scripts/colab_dispatch.py launch SESSION JOBS.json REMOTE_NAME
  python scripts/colab_dispatch.py pull SESSION REMOTE_NAME LOCAL_DIR
  python scripts/colab_dispatch.py status SESSION REMOTE_NAME

Requires the official Colab CLI (`pip install google-colab-cli`) and a signed-in
account. The repository is cloned on the VM; the new cell scripts are copied
over in case they are not yet pushed.
"""
import base64
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ["scripts/run_equal_count_cell.py", "scripts/colab_cell_worker.py"]
BRANCH = "codex/expanded-refraction-sweeps"


def remote(session, code, timeout=600):
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    result = subprocess.run(["colab", "exec", "-s", session, "--timeout", str(timeout)], input=code,
                            text=True, capture_output=True, env=env, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(result.stdout[-2000:] + result.stderr[-2000:])
    return result.stdout


def launch(session, jobs_path, name):
    jobs = json.loads(Path(jobs_path).read_text())
    lines = ["import subprocess, base64, os",
             "subprocess.run('cd /content && (test -d repo || git clone -q --depth 1 -b %s "
             "https://github.com/kndym/benchmarks-PS-reflector repo)', shell=True, check=True)" % BRANCH]
    for path in SCRIPTS:
        data = base64.b64encode((ROOT / path).read_bytes()).decode()
        lines.append(f"open('/content/repo/{path}','wb').write(base64.b64decode('{data}'))")
    lines.append(f"open('/content/{name}.jobs.json','w').write({json.dumps(json.dumps(jobs))})")
    lines.append("env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')")
    lines.append(f"subprocess.Popen(['python','-B','scripts/colab_cell_worker.py','/content/{name}.jobs.json',"
                 f"'/content/results/{name}'], cwd='/content/repo', env=env, start_new_session=True, "
                 f"stdout=open('/content/{name}.log','a'), stderr=subprocess.STDOUT)")
    lines.append(f"print('launched {len(jobs)} jobs')")
    print(remote(session, "\n".join(lines)))


def status(session, name):
    code = (f"import subprocess;print(subprocess.run('cat /content/results/{name}/worker_status.jsonl 2>/dev/null;"
            f" tail -3 /content/{name}.log; pgrep -fc run_refraction_sweeps', shell=True, capture_output=True,"
            f" text=True).stdout)")
    print(remote(session, code))


def pull(session, name, local_dir):
    archive = f"/content/{name}_cells.tar.gz"
    # Summaries, CSVs, residuals, settings, logs, and the entropic map checkpoint (refractor
    # surface: log radius, normals, refracted rays). Shared optical potentials are left behind.
    code = (f"import subprocess;subprocess.run(\"cd /content/results && find {name} -type f \\\\( -name '*.json' "
            f"-o -name '*.csv' -o -name '*.log' -o -name '*.jsonl' -o \\\\( -path '*/entropic/*' "
            f"-name 'solve_*.npz' \\\\) \\\\) -not -path '*shared_optical*' "
            f"| tar czf {archive} -T -\", shell=True, check=True);print('ok')")
    remote(session, code)
    local_dir = Path(local_dir)
    local_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp) / "cells.tar.gz"
        env = dict(os.environ, PYTHONIOENCODING="utf-8")
        subprocess.run(["colab", "download", "-s", session, archive, str(target)], check=True, env=env,
                       capture_output=True)
        with tarfile.open(target) as tar:
            tar.extractall(local_dir, filter="data")
    print(f"pulled into {local_dir / name}")


if __name__ == "__main__":
    command, *rest = sys.argv[1:]
    dict(launch=launch, pull=pull, status=status)[command](*rest)
