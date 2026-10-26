import unittest
import numpy as np
from refracter.sinkhorn import (
    entropic_ot_cost, entropic_sinkhorn_divergence, sinkhorn_divergence_from_terms,
)


class CachedDivergenceTests(unittest.TestCase):
    def test_cached_terms_match_direct_evaluation(self):
        x = np.array([[0., 0., 1.], [.1, 0., np.sqrt(.99)]])
        y = np.array([[0., .1, np.sqrt(.99)], [.2, 0., np.sqrt(.96)]])
        p = q = np.array([.5, .5])
        terms = [entropic_ot_cost(a, b, p, q, .01, return_info=True)
                 for a, b in [(x, y), (x, x), (y, y)]]
        cached, cached_info = sinkhorn_divergence_from_terms(*terms)
        direct, direct_info = entropic_sinkhorn_divergence(x, y, p, q, .01)
        self.assertEqual(cached, direct)
        self.assertEqual(cached_info, direct_info)
        self.assertTrue(cached_info['all_converged'])

    def test_failed_term_is_not_marked_converged(self):
        ok = (1., dict(iterations=3, converged=True))
        failed = (1., dict(iterations=2000, converged=False))
        _, info = sinkhorn_divergence_from_terms(ok, failed, ok)
        self.assertFalse(info['all_converged'])

    def test_nonfinite_cached_term_rejected(self):
        ok = (1., dict(iterations=3, converged=True))
        with self.assertRaises(ValueError):
            sinkhorn_divergence_from_terms((np.nan, ok[1]), ok, ok)


if __name__ == '__main__':
    unittest.main()
