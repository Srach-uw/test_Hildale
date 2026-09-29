from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import prepare_diagnostic_injection_background as background  # noqa: E402


def write_input(path: Path, *, include_errors: bool = True) -> Path:
    payload = {
        "time": np.arange(12, dtype=np.float64),
        "flux": np.array([100.0, 101.0, 99.0, 100.0, 200.0, 201.0, 199.0, 200.0, 100.0, 101.0, 99.0, 100.0]),
        "mask_index": np.array([0, 1, 2, 3, 8, 9, 10, 11], dtype=np.int64),
        "source_file_index": np.array([0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0], dtype=np.int32),
        "known_transit_mask": np.array([False, False, False, False, True, True, True, True, False, False, False, False]),
    }
    if include_errors:
        payload["flux_err"] = np.full(12, 0.5, dtype=np.float64)
    np.savez(path, **payload)
    return path


def test_observed_pdcsap_preserves_flux_and_contract(tmp_path):
    source = write_input(tmp_path / "input.npz")
    manifest = background.run_background(source, tmp_path / "observed", mode="observed_pdcsap", seed=8)
    with np.load(source, allow_pickle=False) as expected, np.load(manifest.parent / background.OUTPUT_NAME, allow_pickle=False) as result:
        assert np.array_equal(result["flux"], expected["flux"])
        assert np.array_equal(result["time"], expected["time"])
        assert np.array_equal(result["mask_index"], expected["mask_index"])
    assert json.loads(manifest.read_text(encoding="utf-8"))["background"]["random_noise_added"] is False


@pytest.mark.parametrize("mode", ["quoted_error_gaussian", "variance_matched_iid"])
def test_random_background_is_deterministic_and_preserves_point_arrays(tmp_path, mode):
    source = write_input(tmp_path / "input.npz")
    first = background.run_background(source, tmp_path / "first", mode=mode, seed=23)
    second = background.run_background(source, tmp_path / "second", mode=mode, seed=23)
    with np.load(first.parent / background.OUTPUT_NAME, allow_pickle=False) as one, np.load(second.parent / background.OUTPUT_NAME, allow_pickle=False) as two, np.load(source, allow_pickle=False) as original:
        assert np.array_equal(one["flux"], two["flux"])
        assert not np.array_equal(one["flux"], original["flux"])
        for name in ("time", "mask_index", "source_file_index", "known_transit_mask", "flux_err"):
            assert np.array_equal(one[name], original[name])
    assert len(json.loads(first.read_text(encoding="utf-8"))["background"]["groups"]) == 2


def test_quoted_error_requires_positive_error_array(tmp_path):
    source = write_input(tmp_path / "input.npz", include_errors=False)
    with pytest.raises(background.BackgroundInputError, match="flux_err"):
        background.run_background(source, tmp_path / "output", mode="quoted_error_gaussian", seed=7)


def test_variance_matched_rejects_zero_scatter(tmp_path):
    source = tmp_path / "flat.npz"
    np.savez(source, time=np.arange(5, dtype=np.float64), flux=np.ones(5), mask_index=np.arange(5, dtype=np.int64))
    with pytest.raises(background.BackgroundInputError, match="robust scale"):
        background.run_background(source, tmp_path / "output", mode="variance_matched_iid", seed=7)


def test_background_refuses_to_replace_existing_provenance(tmp_path):
    source = write_input(tmp_path / "input.npz")
    output = tmp_path / "output"
    background.run_background(source, output, mode="observed_pdcsap", seed=1)
    with pytest.raises(background.BackgroundInputError, match="output_dir"):
        background.run_background(source, output, mode="observed_pdcsap", seed=1)
