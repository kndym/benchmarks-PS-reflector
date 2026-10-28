"""Small independent audits of EOT feasibility and parameter edge cases."""
import numpy as np
import pytest

from refracter.sinkhorn import entropic_ot_cost, entropic_sinkhorn_divergence


def plan_marginals(x, y, p, q, epsilon, potentials):
    f, g = potentials
    costs = np.sum((x[:, None, :] - y[None, :, :]) ** 2, axis=2)
    plan = np.exp(np.log(p[:, None]) + np.log(q[None, :])
                  + (f[:, None] + g[None, :] - costs) / epsilon)
    return plan.sum(axis=1), plan.sum(axis=0)


@pytest.mark.parametrize('cache_entries', [0, 1000])
def test_tiny_epsilon_cannot_claim_convergence_with_wrong_marginals(cache_entries):
    x = np.array([[0.], [1.]])
    p, q = np.array([.9, .1]), np.array([.1, .9])
    _, info = entropic_ot_cost(x, x, p, q, 1e-10, max_iter=5,
                              max_cache_entries=cache_entries,
                              return_info=True, return_potentials=True)
    source, target = plan_marginals(x, x, p, q, 1e-10, info['potentials'])
    assert abs(source.sum() - 1) < 1e-12
    assert abs(source - p).sum() > 1
    assert not info['converged']
    assert info['source_marginal_l1'] == pytest.approx(abs(source - p).sum())
    assert info['target_marginal_l1'] == pytest.approx(abs(target - q).sum(), abs=1e-14)


@pytest.mark.parametrize('epsilon', [.005, .1, 1.])
def test_rectangular_nonuniform_problem_matches_marginals_and_chunking(epsilon):
    x = np.array([[0.], [.3]])
    y = np.array([[.1], [.5], [.8]])
    p, q = np.array([.3, .7]), np.array([.2, .4, .4])
    values = []
    for cache_entries in (0, 1000):
        value, info = entropic_ot_cost(x, y, p, q, epsilon, chunk_size=1,
            max_cache_entries=cache_entries, return_info=True, return_potentials=True)
        assert info['converged']
        source, target = plan_marginals(x, y, p, q, epsilon, info['potentials'])
        np.testing.assert_allclose(source, p, atol=1e-7, rtol=0)
        np.testing.assert_allclose(target, q, atol=1e-7, rtol=0)
        values.append(value)
    assert values[0] == pytest.approx(values[1], abs=1e-12)


@pytest.mark.parametrize('scale', [1e-320, 1., 1e308])
def test_common_weight_scale_does_not_change_objective(scale):
    x, y = np.array([[0.], [.3]]), np.array([[.1], [.5]])
    expected = entropic_ot_cost(x, y, [1., 1.], [1., 1.], .1)
    value, info = entropic_ot_cost(x, y, [scale, scale], [scale, scale], .1,
                                  return_info=True)
    assert info['converged']
    assert value == pytest.approx(expected, abs=1e-12)


@pytest.mark.parametrize('gauge', [1e6, 1e12, 1e100])
def test_large_warm_start_gauge_is_removed_before_iteration(gauge):
    x, y = np.array([[0.], [.3]]), np.array([[.1], [.5], [.8]])
    p, q = np.array([.3, .7]), np.array([.2, .4, .4])
    expected = entropic_ot_cost(x, y, p, q, .1)
    value, info = entropic_ot_cost(x, y, p, q, .1, return_info=True,
        initial_potentials=(np.full(2, gauge), np.full(3, -gauge)))
    assert info['converged']
    assert value == pytest.approx(expected, abs=1e-12)
    assert max(info['source_marginal_l1'], info['target_marginal_l1']) < 1e-7


def test_duplicate_atoms_and_zero_weights_preserve_discrete_measure():
    x = np.array([[.3], [0.], [.3], [10.]])
    y = np.array([[.2], [.6], [.2], [-10.]])
    value = entropic_ot_cost(x, y, [.2, .3, .5, 0.], [.1, .6, .3, 0.], .1)
    expected = entropic_ot_cost([[0.], [.3]], [[.2], [.6]], [.3, .7], [.4, .6], .1)
    assert value == pytest.approx(expected, abs=1e-12)
    divergence, info = entropic_sinkhorn_divergence(x, x, [.2, .3, .5, 0.],
                                                 [.2, .3, .5, 0.], .1)
    assert info['all_converged']
    assert divergence == pytest.approx(0., abs=1e-12)


@pytest.mark.parametrize('tolerance', [np.nan, np.inf, 0., -1.])
def test_invalid_stopping_tolerance_is_rejected(tolerance):
    with pytest.raises(ValueError):
        entropic_ot_cost([[0.]], [[1.]], [1.], [1.], .1, tolerance=tolerance)


def test_inverse_regularization_must_be_representable():
    with pytest.raises(ValueError, match='epsilon'):
        entropic_ot_cost([[0.]], [[1.]], [1.], [1.], 1e-320)
