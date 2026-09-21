from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from compare_alderaan_repeat_posteriors import (  # noqa: E402
    ComparisonFailure,
    normalize_log_weights,
    read_fit,
    validate_pair,
    weighted_quantiles,
)


def write_fit(
    path: Path,
    *,
    target: str = "K00001",
    npl: int = 1,
    log_weights: np.ndarray | None = None,
    include_weights: bool = True,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    log_weights = np.zeros(4) if log_weights is None else np.asarray(log_weights, dtype=float)
    count = len(log_weights)
    columns = []
    for index in range(npl):
        for name, offset in (("ROR", 0.01), ("IMPACT", 0.1), ("DUR14", 0.2)):
            columns.append(
                fits.Column(
                    name=f"{name}_{index}",
                    format="D",
                    array=np.arange(count, dtype=float) + offset + index,
                )
            )
    if include_weights:
        columns.append(fits.Column(name="LN_WT", format="D", array=log_weights))
    primary = fits.PrimaryHDU()
    primary.header["TARGET"] = target
    primary.header["NPL"] = npl
    fits.HDUList([primary, fits.BinTableHDU.from_columns(columns, name="SAMPLES")]).writeto(path)
    return path


def test_normalized_weights_and_weighted_quantiles_preserve_pairing():
    weights = normalize_log_weights(np.log([1.0, 3.0]))
    assert weights.sum() == pytest.approx(1.0)
    assert weights.tolist() == pytest.approx([0.25, 0.75])
    assert 1.0 / np.sum(weights**2) == pytest.approx(1.6)
    assert weighted_quantiles(np.array([0.0, 10.0]), weights).tolist() == pytest.approx(
        [0.0, 10.0 / 3.0, 7.866666666666667]
    )


def test_read_fit_rejects_missing_ln_wt(tmp_path):
    path = write_fit(tmp_path / "K00001-results.fits", include_weights=False)
    with pytest.raises(ComparisonFailure, match="missing LN_WT") as error:
        read_fit(path)
    assert error.value.code == "missing_ln_wt"


@pytest.mark.parametrize("invalid", [np.nan, np.inf, -np.inf])
def test_read_fit_rejects_invalid_ln_wt(tmp_path, invalid):
    path = write_fit(
        tmp_path / "K00001-results.fits", log_weights=np.array([0.0, invalid, -1.0])
    )
    with pytest.raises(ComparisonFailure, match="non-finite") as error:
        read_fit(path)
    assert error.value.code == "invalid_ln_wt"


def test_validate_pair_rejects_mismatched_npl(tmp_path):
    original = read_fit(write_fit(tmp_path / "original.fits", npl=1))
    recovery = read_fit(write_fit(tmp_path / "recovery.fits", npl=2))
    with pytest.raises(ComparisonFailure, match="original NPL=1; recovery NPL=2") as error:
        validate_pair("K00001", original, recovery)
    assert error.value.code == "npl_mismatch"


def test_validate_pair_rejects_target_mismatch(tmp_path):
    original = read_fit(write_fit(tmp_path / "original.fits", target="K00001"))
    recovery = read_fit(write_fit(tmp_path / "recovery.fits", target="K99999"))
    with pytest.raises(ComparisonFailure, match="recovery TARGET=K99999") as error:
        validate_pair("K00001", original, recovery)
    assert error.value.code == "target_mismatch"
