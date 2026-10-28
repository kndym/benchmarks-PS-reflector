"""Compare saved kappa=0.8 maps to an untruncated optical OT solve."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.special import logsumexp
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from refracter.cost import set_kappa, cost_matrix_chunk, cost_matrix_l2_chunk
from refracter.sinkhorn import entropic_ot_cost, _coalesce_weighted_points
from refracter.optics import entropic_surface_normals, refract_using_normals
from result_paths import artifact_path
from atomic_io import atomic_write_json

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'results/runs/isolated_refractor_problems_oct_27_26/kappa_only/entropic'
OUTPUT = ROOT/'results/runs/kappa_cutoff_diagnostic_oct_28_26'


def kernel_control(costs,p,q,g,epsilon):
    """Independent untruncated scaling, for kernels safely above underflow."""
    offset = float(costs.min())
    exponents = -(costs-offset)/epsilon
    if exponents.min() < -650:
        raise ValueError('Use log-domain control for this kernel range')
    kernel = np.exp(exponents)
    b = np.exp((g-g.max())/epsilon)
    for iteration in range(1,10001):
        a = 1/(kernel@(q*b))
        b = 1/(kernel.T@(p*a))
        scale = b.max()
        b /= scale
        a *= scale
        if iteration % 10 == 0:
            source_l1 = float(np.abs(p*a*(kernel@(q*b))-p).sum())
            if source_l1 < 1e-8:
                return offset+epsilon*np.log(a),epsilon*np.log(b),dict(
                    iterations=iteration,converged=True,source_l1=source_l1,
                    kernel_exponent_min=float(exponents.min()),method='untruncated_kernel_scaling')
    raise RuntimeError('Kernel control did not converge')


def kernel_eot(x,y,p,q):
    x,p = _coalesce_weighted_points(x,p)
    y,q = _coalesce_weighted_points(y,q)
    p,q = p/p.sum(),q/q.sum()
    costs = cost_matrix_l2_chunk(x,y)
    f,g,info = kernel_control(costs,p,q,np.zeros(len(y)),1/320)
    return float(p@f+q@g),info


def compare_divergence(nk):
    with np.load(next(SOURCE.rglob(f'solve_NK{nk}_*.npz'))) as data:
        x,y,p,q = [data[k] for k in ('x','y','p','q')]
        maps = dict(original_entropic=data['y_push_3d'],original_hard=data['hard_y_push_3d'])
    with np.load(artifact_path(OUTPUT,f'corrected_NK{nk}.npz')) as data:
        maps['corrected_entropic'] = data['y_push_3d']
        set_kappa(.8)
        indices = np.argmin(cost_matrix_chunk(x,y)-data['g_raw'][None,:],axis=1)
        maps['corrected_hard'] = y[indices]
    target_self,_ = kernel_eot(y,y,q,q)
    results = {}
    for label,points in maps.items():
        cross,cross_info = kernel_eot(points,y,p,q)
        self_cost,self_info = kernel_eot(points,points,p,p)
        results[label] = dict(cross_eot=cross,source_self=self_cost,target_self=target_self,
                             divergence=cross-.5*(self_cost+target_self),
                             cross_info=cross_info,self_info=self_info)
        print('DIVERGENCE',nk,label,results[label],flush=True)
        atomic_write_json(artifact_path(OUTPUT,f'divergence_control_NK{nk}.json',for_write=True),results)


def diagnose(nk,kernel=False):
    start = time.monotonic()
    with np.load(next(SOURCE.rglob(f'solve_NK{nk}_*.npz'))) as data:
        x,y,p,q,f,g = [data[k] for k in ('x','y','p','q','f_raw','g_raw')]
        old_push = data['y_push_3d'].copy()
    epsilon = 1/(8*int(np.sqrt(nk)))
    set_kappa(.8)
    costs = cost_matrix_chunk(x,y)

    def marginals(ff,gg):
        plan = np.exp(np.log(p[:,None])+np.log(q[None,:])+(ff[:,None]+gg[None,:]-costs)/epsilon)
        rows, cols = plan.sum(axis=1), plan.sum(axis=0)
        return dict(mass=float(rows.sum()), source_l1=float(abs(rows-p).sum()),
                    target_l1=float(abs(cols-q).sum()))

    original = marginals(f,g)
    print('ORIGINAL',nk,original,flush=True)
    if kernel:
        ff,gg,info = kernel_control(costs,p,q,g,epsilon)
    else:
        _,info = entropic_ot_cost(x,y,p,q,epsilon,cost_fn=cost_matrix_chunk,
            initial_potentials=(f,g),max_iter=10000,tolerance=1e-10,
            max_cache_entries=30_000_000,return_info=True,return_potentials=True)
        ff,gg = info.pop('potentials')
    corrected = marginals(ff,gg)
    print('CORRECTED',nk,corrected,info,flush=True)
    normals,_ = entropic_surface_normals(x,y,q,gg,epsilon,.8,256)
    new_push = refract_using_normals(x,normals,.8)
    np.savez(artifact_path(OUTPUT,f'corrected_NK{nk}.npz',for_write=True),
             x=x,y=y,p=p,q=q,f_raw=ff,g_raw=gg,y_push_3d=new_push)
    metrics = {}
    for label,points in (('original',old_push),('untruncated',new_push)):
        if kernel:
            value,diagnostics = kernel_eot(points,y,p,q)
        else:
            value,diagnostics = entropic_ot_cost(points,y,p,q,1/320,return_info=True,
                 max_iter=2000,tolerance=1e-9,max_cache_entries=30_000_000)
        metrics[label] = dict(cross_eot=value,diagnostics=diagnostics)
        print('METRIC',nk,label,metrics[label],flush=True)
    report = dict(nk=nk,epsilon=epsilon,kappa=.8,original=original,corrected=corrected,
                  solver=info,metrics=metrics,elapsed_seconds=time.monotonic()-start)
    atomic_write_json(artifact_path(OUTPUT,f'cutoff_NK{nk}.json',for_write=True),report)
    print('SAVED',report,flush=True)


if __name__ == '__main__':
    kernel = '--kernel' in sys.argv
    for nk in map(int,[arg for arg in sys.argv[1:] if not arg.startswith('--')] or [2275,6746]):
        if '--divergence' in sys.argv:
            compare_divergence(nk)
        else:
            diagnose(nk,kernel=kernel)
