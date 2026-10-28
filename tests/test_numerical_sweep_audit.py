"""Numerical provenance and completion checks for partial comparison sweeps."""
import json
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import plot_refraction_pushforward_comparison as plots
import evaluate_saved_refraction_divergence as evaluator
import plot_refraction_sweep_progress as progress


def parameters(**changes):
    return dict(dict(nk=[100], epsilon=[.01], fixed_nk=100, fixed_epsilon=.01,
        target_nk=100, cost_epsilon=.003125, refractor_max_iter=2000,
        identity_max_iter=2000, refractor_tolerance=1e-5, pf_tolerance=1e-9,
        self_cost_mode='separate_marginals_v1', pf_max_iter=2000,
        kappa=.8, source_theta=[15,60], source_phi=[15,45], target_theta=[18,36],
        target_phi=[18,36], experiments=[1], solver_version='untruncated_marginal_v2',
        optical_marginal_tolerance=1e-6), **changes)


def input_folders(tmp_path, params):
    for mode in ('hard','entropic'):
        folder = tmp_path/mode
        folder.mkdir()
        (folder/'parameters.json').write_text(json.dumps(dict(params,pushforward=mode)))


@pytest.mark.parametrize('change', [dict(kappa=.6), dict(source_phi=[15,55]),
    dict(target_theta=[18,45]), dict(solver_version='legacy'),
    dict(optical_marginal_tolerance=1e-4), dict(experiments=[1,2])])
def test_plot_rejects_different_physics_or_solver_policy(tmp_path, change):
    params = parameters()
    input_folders(tmp_path, params)
    (tmp_path/'hard/parameters.json').write_text(json.dumps(dict(params,pushforward='hard',**change)))
    with pytest.raises(ValueError, match='Comparison parameter mismatch'):
        plots.render(tmp_path/'entropic',tmp_path/'hard',yscales=())


def test_joint_only_comparison_completes_with_its_one_expected_pair(tmp_path):
    input_folders(tmp_path, parameters())
    row = dict(experiment='1',nk='100',refractor_epsilon='.0125',
               pf_epsilon='.003125',cost='.001',converged='True')
    with patch.object(plots,'read_rows',side_effect=[[row],[row]]):
        assert plots.render(tmp_path/'entropic',tmp_path/'hard',yscales=()) is True


def test_evaluator_rejects_stale_output_before_overwriting_its_provenance(tmp_path):
    params = parameters()
    input_folders(tmp_path, params)
    old = parameters(kappa=.6)
    previous = tmp_path/'debiased/hard'
    previous.mkdir(parents=True)
    original = json.dumps(dict(old,pushforward='hard'))
    (previous/'parameters.json').write_text(original)
    with patch.object(sys,'argv',['evaluate','--hard',str(tmp_path/'hard'),
            '--entropic',str(tmp_path/'entropic'),'--output',str(tmp_path/'debiased')]):
        with pytest.raises(ValueError,match='kappa'):
            evaluator.main()
    assert (previous/'parameters.json').read_text() == original


def test_legacy_settings_can_match_when_both_lack_new_solver_fields():
    params = parameters()
    for field in ('solver_version','optical_marginal_tolerance'):
        params.pop(field)
    plots.validate_comparison_settings(params,dict(params))


def test_individual_progress_finishes_selected_joint_experiment(tmp_path):
    (tmp_path/'parameters.json').write_text(json.dumps(parameters()))
    (tmp_path/'costs.csv').write_text(
        'experiment,nk,refractor_epsilon,cost\n1,100,.0125,.001\n')
    with patch.object(progress.plt.Figure,'savefig', side_effect=lambda p,**kw: p.write_bytes(b'plot')):
        assert progress.render(tmp_path) is True
