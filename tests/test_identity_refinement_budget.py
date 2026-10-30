import numpy as np
import pytest

from refracter.cost import get_kappa, set_kappa
from refracter.sinkhorn import run_sinkhorn_divergence


def test_refinement_reaches_same_strict_solution_as_larger_cap():
    old = get_kappa()
    set_kappa(.8)
    try:
        a = np.linspace(.15, .8, 12)
        x = np.array([np.sin(a), np.zeros(12), np.cos(a)]).T
        a = np.linspace(.2, .4, 9)
        y = np.array([np.sin(a)*.7, np.sin(a)*np.sqrt(.51), np.cos(a)]).T
        args = dict(epsilon=.0003949243034490774, separate_marginals=True,
                    max_iter=20000, marginal_tolerance=1e-7, verbose=False)
        weights = (np.ones(12)/12, np.ones(9)/9)
        with pytest.raises(RuntimeError, match='Source self solve failed marginal feasibility'):
            run_sinkhorn_divergence(x, y, *weights, identity_max_iter=1, **args)
        extended = run_sinkhorn_divergence(x, y, *weights, identity_max_iter=1,
                                          identity_refinement_budget=19999, **args)
        reference = run_sinkhorn_divergence(x, y, *weights, identity_max_iter=20000, **args)
        assert extended['source_identity_iterations'] > 1
        for field in extended:
            np.testing.assert_array_equal(extended[field], reference[field])
        assert extended['source_identity_marginal_l1'] <= 1e-7
        assert extended['target_identity_marginal_l1'] <= 1e-7
        with pytest.raises(RuntimeError, match='Source self solve failed marginal feasibility'):
            run_sinkhorn_divergence(x, y, *weights, identity_max_iter=1,
                                   identity_refinement_budget=1, **args)
    finally:
        set_kappa(old)


@pytest.mark.parametrize('budget', [-1, .5, None])
def test_invalid_refinement_budget(budget):
    points = np.array([[0., 0., 1.]])
    with pytest.raises(ValueError, match='identity_refinement_budget'):
        run_sinkhorn_divergence(points, points, [1.], [1.],
                               identity_refinement_budget=budget, verbose=False)


def test_transport_refinement_preserves_strict_solution_and_rejects_short_budget():
    old = get_kappa()
    set_kappa(.8)
    try:
        a = np.linspace(.15, .8, 12)
        x = np.array([np.sin(a), np.zeros(12), np.cos(a)]).T
        a = np.linspace(.2, .4, 9)
        y = np.array([np.sin(a)*.7, np.sin(a)*np.sqrt(.51), np.cos(a)]).T
        args = dict(epsilon=.0003949243034490774, separate_marginals=True,
                    identity_max_iter=20000, marginal_tolerance=1e-7, verbose=False)
        weights = (np.ones(12)/12, np.ones(9)/9)
        with pytest.raises(RuntimeError, match='Optical transport failed marginal feasibility'):
            run_sinkhorn_divergence(x, y, *weights, max_iter=1,
                                   transport_refinement_budget=1, **args)
        extended = run_sinkhorn_divergence(x, y, *weights, max_iter=1,
                                          transport_refinement_budget=19999, **args)
        reference = run_sinkhorn_divergence(x, y, *weights, max_iter=20000, **args)
        assert extended['final_iterations'] > 1
        for field in extended:
            np.testing.assert_array_equal(extended[field], reference[field])
    finally:
        set_kappa(old)


@pytest.mark.parametrize('budget', [-1, .5, None])
def test_invalid_transport_refinement_budget(budget):
    points = np.array([[0., 0., 1.]])
    with pytest.raises(ValueError, match='transport_refinement_budget'):
        run_sinkhorn_divergence(points, points, [1.], [1.],
                               transport_refinement_budget=budget, verbose=False)
