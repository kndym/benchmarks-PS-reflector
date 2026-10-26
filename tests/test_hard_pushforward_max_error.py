import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from plot_hard_pushforward_max_error import target_mass_error


class MaxMassErrorTests(unittest.TestCase):
    def test_empty_target_can_have_largest_error(self):
        target=np.eye(3)
        result=target_mass_error(target[:2],[.5,.5],target,[.1,.1,.8])
        self.assertAlmostEqual(result['max_abs_mass_difference'],.8)
        self.assertEqual(result['worst_target_index'],2)
        self.assertEqual(result['pushforward_mass_at_worst'],0)

    def test_merged_sources_are_summed(self):
        target=np.eye(3)
        result=target_mass_error(target[[0,0,1]],[1/3]*3,target,[1/3]*3)
        self.assertAlmostEqual(result['max_abs_mass_difference'],1/3)
        self.assertAlmostEqual(result['max_mass_excess'],1/3)
        self.assertAlmostEqual(result['max_mass_deficit'],1/3)

    def test_identical_measures_have_zero_error(self):
        result=target_mass_error(np.eye(3),[.2,.3,.5],np.eye(3),[.2,.3,.5])
        self.assertEqual(result['max_abs_mass_difference'],0)

    def test_rejects_point_outside_target_support(self):
        with self.assertRaises(ValueError):
            target_mass_error([[1,1,1]],[1],np.eye(3),[1/3]*3)
