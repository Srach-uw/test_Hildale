"""Check the executed patch, including its idempotence and failure behavior."""

import ast
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/alderaan_combined_confirmation_20260806/provenance/patched_source/alderaan/dynesty_helpers.py'
PATCHES = ['cloud/patch_alderaan_paper_priors.py',
           'cloud/ld_validation/patch_alderaan_paper_priors.py',
           'cloud/recovery_142/patch_alderaan_paper_priors.py']


@pytest.mark.parametrize('script', PATCHES)
def test_patch_replaces_function_not_substring_and_is_idempotent(tmp_path, script):
    folder = tmp_path / 'alderaan'
    folder.mkdir()
    target = folder / 'dynesty_helpers.py'
    original = SOURCE.read_text().replace('0.0, 1.0)', '0.0, 0.1)')
    target.write_text(original)
    command = [sys.executable, str(ROOT/script), str(tmp_path)]
    subprocess.run(command, check=True, capture_output=True)
    patched = target.read_text()
    calls = [node for node in ast.walk(ast.parse(patched)) if isinstance(node, ast.Call)
             and isinstance(node.func, ast.Name) and len(node.args) == 3]
    radius = [node for node in calls if ast.unparse(node.args[1]) == '1e-05'
              and ast.unparse(node.args[2]) == '0.99']
    assert len(radius) == 1
    assert radius[0].func.id == 'uniform_ppf'
    assert 'loguniform_ppf(u_[4 + npl * 5], scit, 3 * durations[npl])' in patched
    assert 'norm_ppf(u_[0 + npl * 5], 0.0, 1.0)' in patched
    assert 'norm_ppf(u_[1 + npl * 5], 0.0, 1.0)' in patched
    subprocess.run(command, check=True, capture_output=True)
    assert target.read_text() == patched


@pytest.mark.parametrize('script', PATCHES)
def test_patch_refuses_unknown_bounds_without_writing(tmp_path, script):
    folder = tmp_path / 'alderaan'
    folder.mkdir()
    target = folder / 'dynesty_helpers.py'
    original = SOURCE.read_text().replace('1e-5, 0.99', '1e-5, 0.98')
    target.write_text(original)
    result = subprocess.run([sys.executable, str(ROOT/script), str(tmp_path)], capture_output=True)
    assert result.returncode != 0
    assert target.read_text() == original
