import unittest
import numpy as np
from refracter.cost import cost_matrix_chunk, get_kappa, set_kappa
from refracter.sinkhorn import run_small_sinkhorn


class SmallSinkhornStabilityTests(unittest.TestCase):
    def setUp(self):
        self.kappa = get_kappa()
        set_kappa(.6)
        self.x = np.array([[0., 0., 1.], [.1, 0., np.sqrt(.99)]])
        self.y = np.array([[0., .1, np.sqrt(.99)], [.2, 0., np.sqrt(.96)]])

    def tearDown(self):
        set_kappa(self.kappa)

    def test_matches_kernel_updates_at_moderate_regularization(self):
        p, q, k_reg = np.array([.3, .7]), np.array([.6, .4]), 16
        c = cost_matrix_chunk(self.x, self.y)
        kernel = np.exp(-k_reg*c)
        pp, qq = kernel @ p, kernel @ q
        pp, qq = pp/pp.sum(), qq/qq.sum()
        f, g = np.zeros(2), np.zeros(2)
        for k in range(1, k_reg, int(np.floor(k_reg**(1/3)))):
            matrix = np.exp(-k*(c-f[:, None]-g[None, :]))
            ff = 1/(matrix @ qq)
            gg = 1/(matrix.T @ (ff*pp))
            f += np.log(ff)/k
            g += np.log(gg)/k
        actual = run_small_sinkhorn(self.x, self.y, p, q, k_reg)
        np.testing.assert_allclose(actual, (f, g), atol=2e-14, rtol=2e-14)

    def test_underflowing_smoothing_and_zero_weights_stay_finite(self):
        # At k=16000 all entries of the old same-support kernel are zero.
        for separate in (False, True):
            with np.errstate(over='raise', divide='raise', invalid='raise'):
                f, g = run_small_sinkhorn(self.x, self.y, [1., 0.], [.4, .6],
                                         16000, separate_marginals=separate)
            self.assertTrue(np.isfinite(f).all() and np.isfinite(g).all())


if __name__ == '__main__':
    unittest.main()
