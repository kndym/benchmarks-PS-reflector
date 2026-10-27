import sys
from pathlib import Path
import tempfile
import unittest
import json
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import supervise_refraction_runs as supervisor


class SupervisorCompletionTests(unittest.TestCase):
    def test_persistent_pause_never_launches_worker(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            state=root/'results/refraction_worker_supervision'
            state.mkdir(parents=True)
            (state/'sweep.pause.json').write_text('{}')
            with patch.object(supervisor,'ROOT',root), patch.object(supervisor,'complete',return_value=False), \
                 patch.object(supervisor.os,'chdir'), \
                 patch.object(supervisor.sys,'argv',['supervisor','--role','sweep']), \
                 patch.object(supervisor.subprocess,'Popen') as launch, \
                 patch.object(supervisor.time,'sleep',side_effect=RuntimeError('end test')):
                with self.assertRaisesRegex(RuntimeError,'end test'):
                    supervisor.main()
                launch.assert_not_called()
                self.assertEqual(json.loads((state/'sweep.json').read_text())['state'],'paused')

    def test_requires_both_debiased_tables(self):
        expected={(1,100,.01)}
        with patch.object(supervisor,'expected_keys',return_value=expected):
            with patch.object(supervisor,'case_keys',side_effect=[expected,set()]):
                self.assertFalse(supervisor.complete('debiased'))
            with patch.object(supervisor,'case_keys',side_effect=[expected,expected]):
                self.assertTrue(supervisor.complete('debiased'))

    def test_partial_or_malformed_csv_is_not_complete(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            (root/'sample').mkdir()
            (root/'sample/costs.csv').write_text('experiment,nk,refractor_epsilon,converged\n1,')
            with patch.object(supervisor,'ROOT',root):
                self.assertEqual(supervisor.case_keys('sample'),set())
