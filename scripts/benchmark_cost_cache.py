"""Measure bounded-cost-cache benefit without changing the update equations."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import sys
import time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from refracter.cost import set_kappa
from refracter.cost_cache import bounded_cost_cache
from refracter.qmc import gen_spherical_patch
from refracter.sinkhorn import sinkhorn_identity_f_step
from sinkhorn_resources import enable_threaded_reductions


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();enable_threaded_reductions(4);set_kappa(.8)
    records=[]
    for n in (4096,8192):
        x=gen_spherical_patch(n,np.pi/12,np.pi/3,np.pi/12,np.pi/4)
        y=gen_spherical_patch(1600,np.pi/10,np.pi/5,np.pi/10,np.pi/5)
        lp=np.full(n,-np.log(n));f=np.zeros(n);reference=None
        for mib in (0,256):
            with bounded_cost_cache(x,y,mib*1024**2,256) as cache:
                samples=[]
                for repeat in range(4):
                    start=time.perf_counter()
                    result,_=sinkhorn_identity_f_step(cache.x,cache.x,lp,f,600,256)
                    samples.append(time.perf_counter()-start)
                if reference is None: reference=result
                error=float(np.max(np.abs(result-reference)))
                assert error<=1e-14
                record=dict(n=n,cache_mib=mib,cold_seconds=samples[0],
                    warm_median_seconds=float(np.median(samples[1:])),max_output_error=error,
                    cache=cache.stats())
                records.append(record);print(json.dumps(record),flush=True)
    report=dict(created_at=datetime.now(timezone.utc).isoformat(),workers=4,
        note='Same four-thread optical update; only retained cost blocks differ. Production work active.',results=records)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__': main()
