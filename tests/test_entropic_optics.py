"""Run with python -B -m unittest discover -s tests -p test_entropic_optics.py."""
import unittest
import numpy as np
from refracter.optics import entropic_surface_normals, refract_using_normals


class OpticsTests(unittest.TestCase):
    def setUp(self):
        self.x = np.array([[0., 0., 1.], [.1, 0., np.sqrt(.99)]])
        self.y = np.array([[.2, 0., np.sqrt(.96)], [0., .3, np.sqrt(.91)]])
        self.q = np.array([.3, .7])
        self.g = np.array([.1, -.1])

    def test_single_support_refracts_exactly_to_target(self):
        n, _ = entropic_surface_normals(self.x, self.y[:1], [1], [0], .01, .6)
        t = refract_using_normals(self.x, n, .6)
        np.testing.assert_allclose(t, np.tile(self.y[0], (2, 1)), atol=1e-13)

    def test_normal_matches_implicit_surface_gradient(self):
        n, log_r = entropic_surface_normals(self.x, self.y, self.q, self.g, .1, .6)
        for i in range(len(self.x)):
            p = np.exp(log_r[i])*self.x[i]
            def surface(v):
                d = np.linalg.norm(v) - .6*(self.y @ v)
                return np.sum(self.q*np.exp(self.g/.1)*d**10)-1
            h = 1e-6
            grad = np.array([(surface(p+h*e)-surface(p-h*e))/(2*h) for e in np.eye(3)])
            np.testing.assert_allclose(grad/np.linalg.norm(grad), n[i]/np.linalg.norm(n[i]), atol=1e-8)

    def test_snell_unit_length_and_normal_sign(self):
        n, _ = entropic_surface_normals(self.x, self.y, self.q, self.g, .01, .6)
        np.testing.assert_allclose(np.sum(self.x*n, axis=1), 1, atol=1e-13)
        t = refract_using_normals(self.x, n, .6)
        np.testing.assert_allclose(np.linalg.norm(t, axis=1), 1, atol=1e-13)
        np.testing.assert_allclose(t, refract_using_normals(self.x, -3*n, .6), atol=1e-13)
        unit_n = n / np.linalg.norm(n, axis=1)[:, None]
        tangent = lambda v: v - np.sum(v*unit_n, axis=1)[:, None]*unit_n
        np.testing.assert_allclose(tangent(self.x), .6*tangent(t), atol=1e-13)

    def test_stability_gauge_zero_weights_and_chunks(self):
        n, r = entropic_surface_normals(self.x, self.y, self.q, self.g, .0000625, .6, 1)
        nn, rr = entropic_surface_normals(self.x, self.y, self.q, self.g+100, .0000625, .6, 512)
        np.testing.assert_allclose(n, nn, atol=1e-11)
        np.testing.assert_allclose(r-100, rr, atol=1e-11)
        n, _ = entropic_surface_normals(self.x, self.y, [1, 0], [0, np.nan], .01, .6)
        self.assertTrue(np.isfinite(n).all())

    def test_total_internal_reflection_raises(self):
        with self.assertRaisesRegex(ValueError, 'Total internal reflection'):
            refract_using_normals([[1, 0, 0]], [[0, 0, 1]], .6)


if __name__ == '__main__':
    unittest.main()
