"""Benchmark execution-only reduction parallelism against identical outputs."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.special import logsumexp
import psutil
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import refracter.sinkhorn as sinkhorn
from refracter.cost import set_kappa
from refracter.qmc import gen_spherical_patch


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    records=[]
    original=sinkhorn.logsumexp
    x=gen_spherical_patch(4096,np.pi/12,np.pi/3,np.pi/12,np.pi/4)
    logp=np.full(len(x),-np.log(len(x)));f=np.zeros(len(x));set_kappa(.8)
    reference=None
    for workers in (1,2,4):
        with ThreadPoolExecutor(workers) as pool:
            def reduction(a,axis=None,**kwargs):
                if workers==1 or a.ndim!=2 or axis not in (0,1) or a.size<1_000_000 or kwargs:
                    return original(a,axis=axis,**kwargs)
                pieces=np.array_split(a,min(workers,a.shape[1-axis]),axis=1-axis)
                return np.concatenate(list(pool.map(lambda p:original(p,axis=axis),pieces)))
            sinkhorn.logsumexp=reduction
            samples=[]
            for repeat in range(3):
                start=time.perf_counter()
                result,change=sinkhorn.sinkhorn_identity_f_step(x,x,logp,f,600,256)
                samples.append(time.perf_counter()-start)
            if reference is None: reference=result
            error=float(np.max(np.abs(result-reference)))
            assert error<=1e-14
            record=dict(workers=workers,seconds=float(np.median(samples)),
                        max_output_error=error,operation='one_4096_point_source_identity_update')
            records.append(record);print(record,flush=True)
    sinkhorn.logsumexp=original
    report=dict(created_at=datetime.now(timezone.utc).isoformat(),
        physical_cores=psutil.cpu_count(logical=False),logical_processors=psutil.cpu_count(),
        free_memory_bytes=psutil.virtual_memory().available,
        note='Measured with production work active; same inputs, algorithm and reductions.',results=records)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__': main()
