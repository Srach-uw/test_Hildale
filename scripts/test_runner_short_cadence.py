"""Exercise each runner's cadence argument block without invoking ALDERAAN."""
from pathlib import Path
import shutil
import subprocess
import pytest

ROOT = Path(__file__).resolve().parents[1]
BASH = shutil.which('bash') or 'C:/Program Files/Git/bin/bash.exe'


@pytest.mark.parametrize('relative', ['cloud/run_one_target.sh',
    'cloud/ld_validation/run_one_target.sh','cloud/recovery_142/run_one_target.sh'])
@pytest.mark.parametrize('mode,expected', [('both','--use_sc\nTrue\n'),('long','')])
def test_actual_runner_cadence_block(relative,mode,expected):
    if not Path(BASH).exists():
        pytest.skip('Bash unavailable')
    path=ROOT/relative
    source=path.read_text()
    start=source.index('SC_ARGS=()')
    end=source.index('python bin/detrend_and_estimate_ttvs.py',start)
    command=source[end:].splitlines()[0]
    assert '"${SC_ARGS[@]}"' in command
    script='CADENCE_MODE='+mode+'\n'+source[start:end]+'\nfor arg in "${SC_ARGS[@]}"; do printf "%s\\n" "$arg"; done\n'
    result=subprocess.run([BASH,'--noprofile','--norc','-c',script],capture_output=True,text=True,check=True)
    assert result.stdout==expected
    subprocess.run([BASH,'-n',str(path)],check=True,capture_output=True)
