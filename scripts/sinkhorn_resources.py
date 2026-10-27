"""Parallelize independent SciPy reductions while preserving each reduction."""
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import refracter.sinkhorn as sinkhorn_module


def enable_threaded_reductions(workers, min_entries=1_000_000):
    if workers < 1:
        raise ValueError('Reduction workers must be positive')
    if workers == 1:
        return None
    pool=ThreadPoolExecutor(max_workers=workers)
    original=sinkhorn_module.logsumexp

    def threaded(a, axis=None, **kwargs):
        if a.ndim != 2 or axis not in (0,1) or a.size < min_entries or kwargs:
            return original(a,axis=axis,**kwargs)
        pieces=np.array_split(a,min(workers,a.shape[1-axis]),axis=1-axis)
        return np.concatenate(list(pool.map(lambda piece:original(piece,axis=axis),pieces)))

    sinkhorn_module.logsumexp=threaded
    return pool
