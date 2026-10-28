"""Regression and small-parameter checks for untruncated optical transport."""
import numpy as np
import pytest
from refracter.cost import set_kappa, get_kappa, cost_matrix_chunk
from refracter.qmc import gen_spherical_patch
from refracter.sinkhorn import (sinkhorn_step, sinkhorn_identity_f_step,
    sinkhorn_identity_g_step, run_sinkhorn_divergence, optical_marginal_errors,
    entropic_ot_cost)


@pytest.fixture(autouse=True)
def restore_kappa():
    old = get_kappa()
    yield
    set_kappa(old)


def test_dead_cutoff_rows_cannot_report_zero_change():
    set_kappa(.8)
    x = np.array([[0.,0.,1.],[.1,0.,np.sqrt(.99)]])
    logw = np.log([.5,.5])
    f,g,change = sinkhorn_step(x,x,logw,logw,np.full(2,-100.),np.zeros(2),100.)
    assert change > 1
    assert np.isfinite(f).all() and np.isfinite(g).all()
    for step in (sinkhorn_identity_f_step,sinkhorn_identity_g_step):
        potential,change = step(x,x,logw,np.full(2,-100.),100.)
        assert change > 1
        assert np.isfinite(potential).all()


@pytest.mark.parametrize('kappa,epsilon',[(.2,.01),(.6,.003125),(.8,.0015),(.95,.005)])
def test_full_optical_matches_untruncated_reference(kappa,epsilon):
    set_kappa(kappa)
    x = gen_spherical_patch(19,.2,.8,.2,.5)
    y = gen_spherical_patch(13,.3,.6,.25,.45)
    p,q = np.ones(19)/19,np.ones(13)/13
    result = run_sinkhorn_divergence(x,y,p,q,epsilon=epsilon,separate_marginals=True,
        max_iter=10000,identity_max_iter=10000,tolerance=1e-5,
        marginal_tolerance=1e-7,verbose=False,chunk_size=7)
    assert max(optical_marginal_errors(x,y,p,q,result['f_raw'],result['g_raw'],epsilon,7)) <= 1e-7
    _,info = entropic_ot_cost(x,y,p,q,epsilon,cost_fn=cost_matrix_chunk,
        max_iter=10000,tolerance=1e-10,return_info=True,return_potentials=True)
    assert info['converged']
    ref_g = info['potentials'][1]
    difference = result['g_raw']-ref_g
    assert np.max(np.abs(difference-difference.mean())) < 2e-6


def test_capped_infeasible_optics_raise_instead_of_writing_results():
    set_kappa(.8)
    x=gen_spherical_patch(15,.2,.8,.2,.5)
    y=gen_spherical_patch(11,.3,.6,.25,.45)
    with pytest.raises(RuntimeError,match='marginal feasibility'):
        run_sinkhorn_divergence(x,y,np.ones(15),np.ones(11),epsilon=.002,
            separate_marginals=True,max_iter=1,identity_max_iter=1,
            marginal_tolerance=1e-10,verbose=False)
