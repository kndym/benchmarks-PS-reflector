"""Live workers may publish an entropic case before its hard counterpart."""
import json
import sys
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import evaluate_saved_refraction_divergence as evaluator
import plot_refraction_pushforward_comparison as plots


def parameters():
    return dict(nk=[100], epsilon=[0.01], fixed_nk=100, fixed_epsilon=0.01,
                target_nk=100, cost_epsilon=0.003125, refractor_max_iter=2000,
                identity_max_iter=2000, refractor_tolerance=1e-5, pf_tolerance=1e-9,
                self_cost_mode='separate_marginals_v1', pf_max_iter=2000,
                pushforward='entropic')


def test_plot_accepts_entropic_ahead_and_waits_for_both(tmp_path):
    params = parameters()
    for mode in ('hard', 'entropic'):
        path = tmp_path / mode
        path.mkdir()
        (path / 'parameters.json').write_text(json.dumps(dict(params, pushforward=mode)))
    row = dict(experiment='1', nk='100', refractor_epsilon='0.0125', converged='True',
               pf_epsilon='0.003125', cost='0.001')
    with patch.object(plots, 'read_rows', side_effect=[[], [row]]):
        assert plots.render(tmp_path/'entropic', tmp_path/'hard') is False


def test_evaluator_leaves_unpaired_entropic_case_pending(tmp_path):
    params = parameters()
    for mode in ('hard', 'entropic'):
        path = tmp_path / mode
        path.mkdir()
        (path / 'parameters.json').write_text(json.dumps(dict(params, pushforward=mode)))
    row = dict(experiment='1', nk='100', refractor_epsilon='0.0125')
    with patch.object(sys, 'argv', ['evaluate', '--hard', str(tmp_path/'hard'),
            '--entropic', str(tmp_path/'entropic'), '--output', str(tmp_path/'debiased')]), \
         patch.object(evaluator, 'render'), \
         patch.object(evaluator, 'read_rows', side_effect=[[], [], [], [row]]):
        evaluator.main()
