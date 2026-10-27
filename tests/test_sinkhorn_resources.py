import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from sinkhorn_resources import enable_threaded_reductions
import refracter.sinkhorn as sinkhorn


class ParallelReductionTests(unittest.TestCase):
    def test_parallel_solver_preserves_objective_and_convergence(self):
        rng=np.random.default_rng(8)
        x,y=rng.normal(size=(31,3)),rng.normal(size=(43,3))
        p,q=np.full(31,1/31),np.full(43,1/43)
        options=dict(epsilon=.5,chunk_size=11,max_iter=1000,tolerance=1e-9,
                     return_info=True,max_cache_entries=0)
        reference=sinkhorn.entropic_ot_cost(x,y,p,q,**options)
        original=sinkhorn.logsumexp
        pool=enable_threaded_reductions(4,min_entries=0)
        try:
            actual=sinkhorn.entropic_ot_cost(x,y,p,q,**options)
            self.assertEqual(actual,reference)
            self.assertTrue(actual[1]['converged'])
        finally:
            pool.shutdown()
            sinkhorn.logsumexp=original
