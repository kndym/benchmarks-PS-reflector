"""Adopt the live kappa workers and upgrade each at a completed-case boundary."""
import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import numpy as np
from atomic_io import atomic_write_json, atomic_replace
from result_paths import artifact_path, artifact_glob
from plot_refraction_pushforward_comparison import COMPARISON_SETTINGS, comparison_case_keys, read_rows

ROOT=Path(__file__).resolve().parents[1]


def current_case(log):
    if not Path(log).exists(): return None
    matches=re.findall(r'Refractor NK=(\d+), epsilon=([0-9.eE+-]+)',Path(log).read_text(errors='replace'))
    return (int(matches[-1][0]),float(matches[-1][1])) if matches else None


def case_saved(folder, case):
    return case is not None and any(int(r['nk'])==case[0]
        and math.isclose(float(r['refractor_epsilon']),case[1],rel_tol=1e-7,abs_tol=1e-15)
        and r['converged']=='True' for r in read_rows(folder))


def worker_command(original, workers, shared_cache, cost_cache_mib=None):
    command=list(original)
    command[command.index('--workers')+1]=str(workers)
    if '--shared-optical-cache' not in command:
        command += ['--shared-optical-cache',str(shared_cache)]
    if cost_cache_mib is not None:
        if '--cost-cache-mib' in command:
            command[command.index('--cost-cache-mib')+1]=str(cost_cache_mib)
        else:
            command += ['--cost-cache-mib',str(cost_cache_mib)]
    return command


def matching_supervisor(command,folder):
    if any(Path(part).name=='run_three_refractor_problems.py' for part in command):
        return ('--output-root' in command and '--case' in command
            and Path(command[command.index('--output-root')+1]).resolve()
                /command[command.index('--case')+1]==folder)
    if any(Path(part).name=='accelerate_kappa_workers.py' for part in command):
        return '--case-root' in command and Path(command[command.index('--case-root')+1]).resolve()==folder
    return False


def mirror_completed(source, destination):
    """Prefetch immutable solved cases; never modify the live destination CSV."""
    total=0
    for mode in ('hard','entropic'):
        src,dst=source/mode,destination/mode
        source_parameters=artifact_path(src,'parameters.json')
        if not source_parameters.exists():
            continue
        sp=json.loads(source_parameters.read_text());dp=json.loads(artifact_path(dst,'parameters.json').read_text())
        for key in COMPARISON_SETTINGS:
            equal = set(sp[key])==set(dp[key]) if key in ('epsilon','nk') else sp.get(key)==dp.get(key)
            if not equal:
                raise ValueError(f'Prefetch numerical parameter mismatch: {key}')
        # Once the destination itself starts experiment 3, leave its writes alone.
        if any(int(row['experiment'])==3 for row in read_rows(dst)):
            continue
        active=current_case(destination/'logs'/f'{mode}.log')
        for row in read_rows(src):
            nk,eps=int(row['nk']),float(row['refractor_epsilon'])
            if active and nk==active[0] and math.isclose(eps,active[1],rel_tol=1e-7):
                continue
            solves=artifact_glob(src,f'solve_NK{nk}_*eps{eps:.17g}.npz')
            if len(solves)!=1:
                raise ValueError('Missing or ambiguous completed prefetch solve')
            solve=solves[0];suffix=solve.stem.removeprefix('solve_')
            paths=[solve,artifact_path(src,f'residual_{suffix}.json'),
                   artifact_path(src,f'pf_{suffix}_eval{float(row["pf_epsilon"]):.17g}.json')]
            if artifact_path(dst,solve.name).exists():
                continue
            diagnostic=json.loads(paths[1].read_text())
            errors=[diagnostic[key] for key in ('source_marginal_l1','target_marginal_l1',
                'source_identity_marginal_l1','target_identity_marginal_l1')]
            cross,info=json.loads(paths[2].read_text())
            if (not np.isfinite(errors).all() or max(errors)>dp['optical_marginal_tolerance']
                    or row['converged']!='True' or not info['converged']
                    or not np.isclose(cross,float(row['cost']),atol=1e-15,rtol=0)):
                raise ValueError('Invalid completed prefetch evaluation')
            # Publish the optical NPZ last, so a resumed reader sees a complete trio.
            for path in paths[1:]+paths[:1]:
                target=artifact_path(dst,path.name,for_write=True)
                if not target.exists():
                    temporary=target.with_name(target.name+'.prefetch.tmp')
                    shutil.copyfile(path,temporary);atomic_replace(temporary,target)
            total+=1
    return total


def main():
    import psutil  # Operational Windows supervisor dependency, not needed by the solver.
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-root',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--shared-cache',type=Path,required=True)
    parser.add_argument('--prefetch-case-root',type=Path)
    parser.add_argument('--cost-cache-mib',type=int,default=0)
    args=parser.parse_args()
    if args.workers<1: parser.error('workers must be positive')
    if args.cost_cache_mib<0: parser.error('cost-cache-mib must be nonnegative')
    folder=args.case_root.resolve();state_path=artifact_path(folder,'worker_status.json')
    logs=folder/'logs'
    log=(logs/'acceleration.log').open('a',encoding='utf-8',buffering=1)
    sys.stdout=sys.stderr=log
    old=json.loads(state_path.read_text())
    children={role:psutil.Process(data['pid']) for role,data in old['children'].items()}
    commands={}
    targets={}
    for mode in ('hard','entropic'):
        command=children[mode].cmdline()
        if Path(command[command.index('--output')+1]).resolve()!=folder/mode:
            raise RuntimeError('Refusing to adopt unrelated optical worker')
        commands[mode]=command;targets[mode]=current_case(logs/f'{mode}.log')
    supervisor=psutil.Process(old['supervisor_pid'])
    if not matching_supervisor(supervisor.cmdline(),folder):
        raise RuntimeError('Refusing to replace an unrelated supervisor')
    supervisor.terminate();supervisor.wait(timeout=10)
    upgraded=set();streams=[];prefetched=0
    def alive(proc):
        return proc.poll() is None if isinstance(proc,subprocess.Popen) else proc.is_running()
    def status(state):
        atomic_write_json(state_path,dict(state=state,supervisor_pid=os.getpid(),
            updated_at=datetime.now(timezone.utc).isoformat(),case='kappa_only',
            handoff_targets=targets,upgraded_workers=sorted(upgraded),
            workers_per_upgraded_optical=args.workers,prefetched_cases=prefetched,
            cost_cache_mib_per_optical=args.cost_cache_mib,
            children={role:dict(pid=proc.pid,exit_code=None if alive(proc) else 0)
                      for role,proc in children.items()}))
    print('Adopted live workers; checkpoint-boundary targets',targets,flush=True)
    while True:
        for mode in ('hard','entropic'):
            if mode in upgraded: continue
            if case_saved(folder/mode,targets[mode]):
                proc=children[mode]
                if alive(proc): proc.terminate();proc.wait(timeout=10)
                stream=(logs/f'{mode}.log').open('a',encoding='utf-8');streams.append(stream)
                command=worker_command(commands[mode],args.workers,args.shared_cache.resolve(),args.cost_cache_mib)
                children[mode]=subprocess.Popen(command,cwd=ROOT,stdin=subprocess.DEVNULL,
                    stdout=stream,stderr=subprocess.STDOUT,
                    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'),
                    creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
                upgraded.add(mode)
                print('Upgraded',mode,'to',args.workers,'threads; saved target',targets[mode],flush=True)
        if args.prefetch_case_root and args.prefetch_case_root.exists():
            prefetched+=mirror_completed(args.prefetch_case_root.resolve(),folder)
        status('running')
        if not all(alive(proc) for proc in children.values()):
            params=json.loads(artifact_path(folder/'entropic','parameters.json').read_text())
            expected=comparison_case_keys(params)
            def complete(path):
                return expected <= {(int(r['experiment']),int(r['nk']),float(r['refractor_epsilon']))
                                    for r in read_rows(path) if r['converged']=='True'}
            if all(complete(folder/path) for path in ('hard','entropic','debiased/hard','debiased/entropic')):
                status('complete');return
            if any(not alive(children[mode]) for mode in ('hard','entropic','debiased')):
                status('failed');raise RuntimeError('Worker stopped before completion; inspect logs')
        time.sleep(2)


if __name__=='__main__': main()
