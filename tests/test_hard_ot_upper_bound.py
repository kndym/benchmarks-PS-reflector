import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import linprog

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from plot_writeup_cost_terms import hard_ot_bounds


def test_nonuniform_hard_mass_bounds_enclose_exact_ot():
    y = np.array([[0.,0.,0.], [1.,0.,0.], [3.,0.,0.], [4.,0.,0.]])
    m, q = np.array([.4,.1,.1,.4]), np.array([.2,.4,.3,.1])
    c = np.sum((y[:,None,:]-y[None,:,:])**2,axis=2)
    result = hard_ot_bounds(y,m,q)
    # Independent linear programming optimum, not the residual witness formula.
    constraints = np.vstack([np.kron(np.eye(4),np.ones((1,4))),
                             np.kron(np.ones((1,4)),np.eye(4))])
    optimum = linprog(c.ravel(), A_eq=constraints, b_eq=np.r_[m,q],
                      bounds=(0,None), method='highs')
    assert optimum.success
    assert optimum.fun <= result['feasible_plan_cost']+1e-14
    assert result['feasible_plan_cost'] <= result['worst_deficit_bound']+1e-14
    assert result['worst_deficit_bound'] < result['worst_all_bound']
    assert result['witness_marginal_linf'] < 1e-14
    assert result['excess_mass'] == pytest.approx(.5)
    # Common mass must be kept, not charged the worst cost.
    assert result['retained_mass'] == pytest.approx(.5)


def test_exact_target_mass_has_zero_bound_despite_nonzero_distances():
    y = np.array([[0.,0.,0.], [2.,0.,0.]])
    result = hard_ot_bounds(y,[.25,.75],[.25,.75])
    assert result['worst_all_bound'] == result['feasible_plan_cost'] == 0


def test_all_mass_at_one_atom_and_empty_targets():
    y = np.array([[0.,0.,0.], [1.,0.,0.], [2.,0.,0.]])
    result = hard_ot_bounds(y,[1.,0.,0.],[1/3]*3)
    assert result['feasible_plan_cost'] == pytest.approx(5/3)
    assert result['worst_all_bound'] == pytest.approx(8/3)


@pytest.mark.parametrize('m,q', [([.2,.2],[.5,.5]), ([-.1,1.1],[.5,.5]),
                               ([np.nan,0.],[.5,.5])])
def test_invalid_mass_is_rejected(m,q):
    with pytest.raises(ValueError,match='probability masses'):
        hard_ot_bounds(np.array([[0.,0.,0.],[1.,0.,0.]]),m,q)
