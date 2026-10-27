"""Checkpoint-resuming workers launched by Windows Task Scheduler, outside Codex."""
try:
    from result_paths import artifact_path, artifact_glob, resolve_path, dated_run_dir, ensure_readme, self_cost_dir
except ModuleNotFoundError:
    from scripts.result_paths import artifact_path, artifact_glob, resolve_path, dated_run_dir, ensure_readme, self_cost_dir
import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from atomic_io import atomic_write_json

ROOT = Path(__file__).resolve().parents[1]
HARD = 'results/runs/hard_map_2000_updates_oct_08_26'
SOFT = 'results/runs/smooth_refraction_2000_updates_oct_26_26'
DEBIASED = 'results/runs/debiased_comparison_oct_26_26'
COMMANDS = {
    'sweep': ['scripts/run_refraction_sweeps.py', '--extended', '--experiments', '1','2','3',
              '--target-nk','1600','--pushforward','entropic','--refractor-max-iter','2000',
              '--identity-max-iter','2000','--workers','4','--output',SOFT],
    'debiased': ['scripts/evaluate_saved_refraction_divergence.py','--hard',HARD,
                 '--entropic',SOFT,'--output',DEBIASED,'--workers','4',
                 '--cache-entries','4000000','--watch'],
    'plots': ['scripts/plot_refraction_pushforward_comparison.py',SOFT,'--baseline',HARD,'--watch'],
}


def case_keys(folder):
    try:
        with artifact_path(ROOT/folder, 'costs.csv').open(newline='',encoding='utf-8') as stream:
            return {(int(r['experiment']),int(r['nk']),float(r['refractor_epsilon']))
                    for r in csv.DictReader(stream) if r['converged'] == 'True'}
    except (OSError,ValueError,KeyError):
        return set()


def expected_keys():
    params=json.loads(artifact_path(ROOT/SOFT, 'parameters.json').read_text(encoding='utf-8'))
    import math
    return ({(1,n,1/(8*math.floor(math.sqrt(n)))) for n in params['nk']}
            | {(2,n,params['fixed_epsilon']) for n in params['nk']}
            | {(3,params['fixed_nk'],e) for e in params['epsilon']})


def complete(role):
    try:
        expected=expected_keys()
    except (OSError,ValueError,KeyError):
        return False
    folders = [DEBIASED+'/hard',DEBIASED+'/entropic'] if role == 'debiased' else [SOFT]
    return all(expected <= case_keys(folder) for folder in folders)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--role',choices=COMMANDS,required=True)
    args=parser.parse_args()
    os.chdir(ROOT)
    state_dir=resolve_path(ROOT/'results/refraction_worker_supervision')
    state_dir.mkdir(exist_ok=True)
    legacy_fixture = not state_dir.resolve().is_relative_to(Path(__file__).resolve().parents[1])
    status=(state_dir/f'{args.role}.json' if legacy_fixture else artifact_path(state_dir, f'{args.role}.json', for_write=True))
    python=Path(sys.executable).with_name('python.exe')
    command=[str(python),'-u','-B',*COMMANDS[args.role]]
    restarts=0
    while not complete(args.role):
        pause_file=(state_dir/f'{args.role}.pause.json' if legacy_fixture else artifact_path(state_dir, f'{args.role}.pause.json'))
        if pause_file.exists():
            atomic_write_json(status,dict(state='paused',role=args.role,supervisor_pid=os.getpid(),
                worker_pid=None,restarts=restarts,pause_file=str(pause_file),
                updated_at=datetime.now(timezone.utc).isoformat()))
            time.sleep(30)
            continue
        with artifact_path(state_dir, f'{args.role}.worker.log', for_write=True).open('a',encoding='utf-8') as log:
            child=subprocess.Popen(command,cwd=ROOT,stdin=subprocess.DEVNULL,stdout=log,
                                   stderr=subprocess.STDOUT,
                                   creationflags=(subprocess.CREATE_NO_WINDOW |
                                       (subprocess.NORMAL_PRIORITY_CLASS if args.role == 'debiased' else 0))
                                       if os.name == 'nt' else 0)
            atomic_write_json(status,dict(state='running',role=args.role,supervisor_pid=os.getpid(),
                worker_pid=child.pid,started_at=datetime.now(timezone.utc).isoformat(),
                restarts=restarts,command=command))
            code=child.wait()
        if complete(args.role):
            break
        restarts+=1
        atomic_write_json(status,dict(state='retrying',role=args.role,supervisor_pid=os.getpid(),
            last_exit_code=code,restarts=restarts,retry_seconds=30,
            updated_at=datetime.now(timezone.utc).isoformat()))
        time.sleep(30)
    atomic_write_json(status,dict(state='complete',role=args.role,supervisor_pid=os.getpid(),
        restarts=restarts,completed_at=datetime.now(timezone.utc).isoformat()))


if __name__ == '__main__':
    main()
