import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from plot_hard_pushforward_atom_counts import count_atoms


class AtomCountTests(unittest.TestCase):
    def test_merging_and_empty_locations(self):
        target = np.eye(3)
        result = count_atoms(target[[0, 0, 1]], np.full(3, 1/3), target)
        self.assertEqual(result['unequal_weight_atoms'], 1)
        self.assertEqual(result['unequal_target_locations'], 2)
        self.assertEqual(result['unequal_target_locations'],
                         result['unequal_weight_atoms'] + result['empty_target_locations'])
        self.assertEqual(result['target_locations'] - result['singleton_atoms'], 2)
        self.assertEqual(result['occupied_atoms'], 2)
        self.assertEqual(result['singleton_atoms'], 1)
        self.assertEqual(result['empty_target_locations'], 1)
        self.assertEqual(result['source_points_in_merged_atoms'], 2)

    def test_one_to_one_has_no_unequal_atoms(self):
        target = np.eye(3)
        result = count_atoms(target, np.full(3, 1/3), target)
        self.assertEqual(result['unequal_weight_atoms'], 0)
        self.assertEqual(result['unequal_target_locations'], 0)

    def test_empty_locations_count_even_without_merging(self):
        result = count_atoms(np.eye(3)[:2], np.full(2, .5), np.eye(3))
        self.assertEqual(result['unequal_weight_atoms'], 0)
        self.assertEqual(result['unequal_target_locations'], 1)

    def test_rejects_nonuniform_source(self):
        with self.assertRaises(ValueError):
            count_atoms(np.eye(3), np.array([.2,.3,.5]), np.eye(3))
