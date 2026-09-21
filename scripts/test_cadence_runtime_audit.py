import json
import os
import textwrap
import numpy as np
import pytest
from patch_cadence_runtime_audit import BLOCK, patched, ANCHOR


@pytest.mark.parametrize('cadence,required,fails', [('long','1',True),('short','1',False),('long','0',False)])
def test_runtime_selected_cadence(tmp_path,monkeypatch,cadence,required,fails):
    monkeypatch.setenv('ALDERAAN_REQUIRE_SC',required)
    env=dict(np=np,os=os,quarters=[0],all_dtype=[cadence],all_time=[np.arange(4)],
        all_mask=[np.array([[1,1,0,0],[0,0,1,1]])],texp=[.001],oversample=[1],
        RESULTS_DIR=str(tmp_path),TARGET='K00001',RUN_ID='test')
    if fails:
        with pytest.raises(RuntimeError,match='requires selected'):
            exec(textwrap.dedent(BLOCK),env)
    else:
        exec(textwrap.dedent(BLOCK),env)
    data=json.loads((tmp_path/'K00001_sampler_cadence.json').read_text())
    assert data['quarters'][0]['per_planet_mask_points']==[2,2]


def test_unknown_source_rejected():
    with pytest.raises(ValueError,match='not uniquely'):
        patched('print(1)\n')
