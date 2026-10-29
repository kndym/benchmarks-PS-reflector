import numpy as np
import pytest
from scipy.special import logsumexp
from unittest.mock import patch
import refracter.sinkhorn as sinkhorn
from refracter.cost import cost_matrix_chunk,set_kappa,get_kappa
from refracter.cost_cache import bounded_cost_cache,ACTIVE_COST_CACHE
from refracter.qmc import gen_spherical_patch


@pytest.fixture(autouse=True)
def restore_kappa():
    old=get_kappa();set_kappa(.8)
    yield
    set_kappa(old)


def clouds(n=16,m=8):
    return gen_spherical_patch(n,.2,.8,.2,.5),gen_spherical_patch(m,.3,.6,.25,.45)


def test_partial_cache_keeps_prefix_and_never_exceeds_budget():
    x,y=clouds()
    reference=cost_matrix_chunk(x[:4],x)
    with bounded_cost_cache(x,y,1024,4) as cache:
        first=cost_matrix_chunk(cache.x[:4],cache.x)
        for start in range(4,16,4): cost_matrix_chunk(cache.x[start:start+4],cache.x)
        assert cost_matrix_chunk(cache.x[:4],cache.x) is first
        np.testing.assert_array_equal(first,reference)
        assert cache.stats()['peak_retained_bytes']<=1024
        assert cache.stats()['admission_rejections']>0
        assert not first.flags.writeable
        with pytest.raises(ValueError): first[0,0]=0
    assert not cache.blocks and sum(cache.used.values())==0


def test_snapshots_prevent_stale_costs_when_original_inputs_change():
    x,y=clouds()
    with bounded_cost_cache(x,y,1024**2,4) as cache:
        original=cost_matrix_chunk(cache.x[:4],cache.x).copy()
        x[:]=np.roll(x,1,axis=0)
        np.testing.assert_array_equal(cost_matrix_chunk(cache.x[:4],cache.x),original)
        np.testing.assert_array_equal(cost_matrix_chunk(x[:4],x),cost_matrix_chunk(x[:4].copy(),x.copy()))
        assert not cache.x.flags.writeable


def test_kappa_is_part_of_cache_identity():
    x,y=clouds()
    with bounded_cost_cache(x,y,1024**2,4) as cache:
        first=cost_matrix_chunk(cache.x[:4],cache.x)
        set_kappa(.6);other=cost_matrix_chunk(cache.x[:4],cache.x)
        assert not np.array_equal(first,other)
        set_kappa(.8)
        assert cost_matrix_chunk(cache.x[:4],cache.x) is first


def test_context_restored_on_exception_and_nested_scopes():
    x,y=clouds()
    with bounded_cost_cache(x,y,1024**2,4) as outer:
        first=cost_matrix_chunk(outer.x[:4],outer.x)
        with pytest.raises(RuntimeError):
            with bounded_cost_cache(x,y,1024**2,4) as inner:
                cost_matrix_chunk(inner.x[:4],inner.x)
                raise RuntimeError('failure')
        assert ACTIVE_COST_CACHE.get() is outer
        assert not inner.blocks
        assert cost_matrix_chunk(outer.x[:4],outer.x) is first
    assert ACTIVE_COST_CACHE.get() is None


def test_tiled_warmstart_matches_dense_reference_without_self_matrix():
    x,y=clouds(37,19);p=np.arange(1,38,dtype=float);q=np.arange(19,0,-1,dtype=float)
    p/=p.sum();q/=q.sum();k_reg=16
    c=cost_matrix_chunk(x,y)
    lp=logsumexp(-k_reg*cost_matrix_chunk(x,x)+np.log(p)[None,:],axis=1);lp-=logsumexp(lp)
    lq=logsumexp(-k_reg*cost_matrix_chunk(y,y)+np.log(q)[None,:],axis=1);lq-=logsumexp(lq)
    f,g=np.zeros(len(x)),np.zeros(len(y))
    for k in range(1,k_reg,int(np.floor(k_reg**(1/3)))):
        matrix=-k*(c-f[:,None]-g[None,:])
        lf=-logsumexp(matrix+lq[None,:],axis=1)
        lg=-logsumexp(matrix+(lf+lp)[:,None],axis=0)
        f+=lf/k;g+=lg/k
    calls=[]
    def cost(a,b):
        calls.append((len(a),len(b)))
        return cost_matrix_chunk(a,b)
    with patch.object(sinkhorn,'cost_matrix_chunk',side_effect=cost):
        actual=sinkhorn.run_small_sinkhorn(x,y,p,q,k_reg,separate_marginals=True,chunk_size=7)
    np.testing.assert_allclose(actual[0],f,atol=1e-14,rtol=1e-14)
    np.testing.assert_allclose(actual[1],g,atol=1e-14,rtol=1e-14)
    assert (37,37) not in calls and (19,19) not in calls


def test_cached_full_solve_preserves_potentials_and_feasibility():
    x,y=clouds(19,13);p=np.ones(19)/19;q=np.ones(13)/13
    kwargs=dict(epsilon=.0015,separate_marginals=True,verbose=False,chunk_size=7,
        max_iter=10000,identity_max_iter=10000,marginal_tolerance=1e-7)
    reference=sinkhorn.run_sinkhorn_divergence(x,y,p,q,**kwargs)
    with bounded_cost_cache(x,y,64*1024,7) as cache:
        actual=sinkhorn.run_sinkhorn_divergence(cache.x,cache.y,p,q,**kwargs)
        assert cache.hits>0 and cache.peak_bytes<=64*1024
    for key in reference:
        np.testing.assert_allclose(actual[key],reference[key],atol=1e-13,rtol=1e-13)


def test_invalid_budget_is_rejected():
    x,y=clouds()
    with pytest.raises(ValueError):
        with bounded_cost_cache(x,y,-1): pass


def test_sweep_cli_enables_cache_without_changing_saved_map(tmp_path):
    import json
    from test_repaired_run_supervisor import call_sweep
    for name,budget in (('reference',0),('cached',8)):
        result=call_sweep('--kappa','.8','--cost-cache-mib',budget,'--output',tmp_path/name)
        assert result.returncode==0,result.stderr
    with np.load(next((tmp_path/'reference').rglob('solve_*.npz'))) as ref, \
         np.load(next((tmp_path/'cached').rglob('solve_*.npz'))) as actual:
        np.testing.assert_array_equal(ref['y_push_3d'],actual['y_push_3d'])
        np.testing.assert_allclose(ref['total_cost'],actual['total_cost'],atol=1e-13,rtol=1e-13)
    report=json.loads(next((tmp_path/'cached').rglob('cost_cache_NK*.json')).read_text())
    assert report['hits']>0 and report['peak_retained_bytes']<=8*1024**2
