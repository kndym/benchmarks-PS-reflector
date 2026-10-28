"""Small parameter/scale probes independent of production sweep outputs."""
import numpy as np
import pytest

from refracter.cost import cost_matrix_l2_chunk, cost_matrix_refraction_chunk
from refracter.optics import entropic_surface_normals, refract_using_normals


@pytest.mark.parametrize('offset', [0., 1e9])
def test_l2_preserves_distance_with_common_offset(offset):
    x = np.array([[offset, offset, offset]])
    y = np.array([[offset+1, offset+2, offset+3]])
    np.testing.assert_array_equal(cost_matrix_l2_chunk(x,y), [[14.]])


def test_l2_resolves_nearby_unit_directions():
    x = np.array([[1., 1e-10, 0.]])
    y = np.array([[1., 0., 0.]])
    np.testing.assert_allclose(cost_matrix_l2_chunk(x,y), [[1e-20]], rtol=1e-15, atol=0)


@pytest.mark.parametrize('kappa', [1e-18, 1e-12, .8, 1-1e-12])
def test_refraction_cost_agrees_with_log1p_at_parameter_edges(kappa):
    x = np.array([[1.,0.,0.]])
    expected = -np.log1p(-kappa)
    np.testing.assert_allclose(cost_matrix_refraction_chunk(x,x,kappa), [[expected]], rtol=1e-15, atol=0)


@pytest.mark.parametrize('scale', [1e-300, 1., 1e300])
def test_snell_is_invariant_to_extreme_finite_normal_scale(scale):
    x = np.array([[0.,0.,1.], [.1,0.,np.sqrt(.99)]])
    normal = np.array([[0.,0.,1.], [0.,0.,1.]])
    expected = refract_using_normals(x,normal,.8)
    with np.errstate(over='raise', invalid='raise', divide='raise'):
        actual = refract_using_normals(x,normal*scale,.8)
    np.testing.assert_allclose(actual,expected,atol=1e-14)


def test_normal_weights_are_invariant_to_large_common_scale():
    x=np.array([[0.,0.,1.]])
    y=np.array([[.1,0.,np.sqrt(.99)], [0.,.2,np.sqrt(.96)]])
    expected=entropic_surface_normals(x,y,[1.,1.],[0.,0.],.01,.8)
    with np.errstate(over='raise',invalid='raise',divide='raise'):
        actual=entropic_surface_normals(x,y,[1e308,1e308],[0.,0.],.01,.8)
    for result,reference in zip(actual,expected):
        np.testing.assert_allclose(result,reference,atol=1e-14)


def test_normals_avoid_large_gauge_cancellation():
    x=np.array([[0.,0.,1.]])
    y=np.array([[.1,0.,np.sqrt(.99)], [0.,.2,np.sqrt(.96)]])
    g=1e10-np.log(1-.8*y[:,2])
    expected,_=entropic_surface_normals(x,y,[.4,.6],g-g.max(),.0001,.8)
    actual,_=entropic_surface_normals(x,y,[.4,.6],g,.0001,.8)
    np.testing.assert_allclose(actual,expected,rtol=1e-12,atol=1e-12)


@pytest.mark.parametrize('kappa', [1e-12, .01, .8, 1-1e-12, 1-1e-15])
def test_coincident_rays_remain_normal_incidence_across_kappa(kappa):
    x=np.array([[.3,.4,np.sqrt(.75)]])
    normal,_=entropic_surface_normals(x,x,[1.],[0.],.001,kappa)
    np.testing.assert_allclose(np.einsum('ij,ij->i',x,normal),1,atol=1e-14)
    np.testing.assert_allclose(refract_using_normals(x,normal,kappa),x,atol=1e-14)


def test_supported_tiny_weight_can_be_compensated_by_potential():
    x=np.array([[0.,0.,1.]])
    y=np.array([[.1,0.,np.sqrt(.99)],[0.,.2,np.sqrt(.96)]])
    q=np.array([1e-300,1e300])
    cost=-np.log1p(-.8*y[:,2])
    g=cost-.01*np.log(q)
    actual,_=entropic_surface_normals(x,y,q,g,.01,.8)
    expected,_=entropic_surface_normals(x,y,[1.,1.],cost,.01,.8)
    np.testing.assert_allclose(actual,expected,atol=1e-12)


def test_small_but_physical_total_internal_reflection_is_rejected():
    tangent=.6*(1+1e-10)
    incident=np.array([[tangent,0.,np.sqrt(1-tangent*tangent)]])
    with pytest.raises(ValueError,match='Total internal reflection'):
        refract_using_normals(incident,[[0.,0.,1.]],.6)
