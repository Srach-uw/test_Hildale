from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np


MODULE = Path(__file__).parents[2] / "scripts" / "summarize_direct_chord_diagnostic.py"
SPEC = importlib.util.spec_from_file_location("summarize_direct", MODULE)
module = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(module)


def test_nested_weights_normalize() -> None:
    weights = module.normalized_nested_weights(np.log([0.2, 0.8]), 0.0)
    assert np.isclose(weights.sum(), 1.0)
    assert np.allclose(weights, [0.2, 0.8])


def test_summary_marks_low_ess_unsuitable(tmp_path: Path) -> None:
    path = tmp_path / "samples.npz"
    samples = np.zeros((3, 8))
    samples[:, 4] = [0.1, 0.2, 0.3]
    np.savez(path, samples=samples, logwt=np.log([0.999, 0.0005, 0.0005]), logz=np.array([0.0]), logl=np.array([-3.0, -2.0, -1.0]))
    result = module.summarize(path)
    assert result["suitable_for_route_comparison"] is False
