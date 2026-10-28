import json
import os
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import run_three_refractor_problems as runner


def test_failed_worker_stops_siblings_and_preserves_failed_status(tmp_path):
    class Child:
        def __init__(self,pid,code):
            self.pid,self.returncode,self.terminated=pid,code,False
        def poll(self): return self.returncode
        def terminate(self): self.terminated=True;self.returncode=-1
        def wait(self,timeout=None): return self.returncode
    failed,running=Child(1,1),Child(2,None)
    with patch.object(runner,'ROOT',tmp_path), \
         patch.object(sys,'argv',['run','--case','kappa_only','--smoke']), \
         patch.object(runner.subprocess,'Popen',side_effect=[failed,running]):
        with pytest.raises(RuntimeError,match='before writing parameters'):
            runner.main()
    assert running.terminated
    statuses=list(tmp_path.rglob('worker_status.json'))
    assert len(statuses)==1
    assert json.loads(statuses[0].read_text())['state']=='failed'


def call_sweep(*args):
    root=Path(__file__).resolve().parents[1]
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    return subprocess.run([sys.executable,'-B',str(root/'scripts/run_refraction_sweeps.py'),
        '--experiments','1','--nk','8','--target-nk','8','--refractor-max-iter','10000',
        '--identity-max-iter','10000',*map(str,args)],cwd=root,env=env,capture_output=True,text=True)


def test_polar_endpoint_keeps_uniform_source_weight(tmp_path):
    result=call_sweep('--kappa','.8','--source-phi','0','30','--output',tmp_path/'pole')
    assert result.returncode==0,result.stderr
    checkpoint=next((tmp_path/'pole').rglob('solve_*.npz'))
    with np.load(checkpoint) as data:
        assert data['p'][0]>0
        np.testing.assert_allclose(data['p'],np.ones(8)/8)


def test_sweep_refuses_legacy_checkpoint_provenance(tmp_path):
    output=tmp_path/'old'
    output.mkdir()
    (output/'parameters.json').write_text(json.dumps({'kappa':.6}))
    result=call_sweep('--output',output)
    assert result.returncode==2
    assert 'solver_version' in result.stderr
    assert not list(output.rglob('solve_*.npz'))
