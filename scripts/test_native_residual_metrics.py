import json

import numpy as np
import pytest

from native_residual_metrics import residual_metrics, write_residual_report


def test_largest_one_percent_squared_share_and_robust_scatter():
    metrics = residual_metrics([0, 1, 2, 100], [0, 1, 2, 3], 1.0)

    assert metrics["largest_one_percent_count"] == 1
    assert metrics["largest_one_percent_squared_residual_fraction"] == pytest.approx(10000 / 10005)
    assert metrics["centered_robust_mad"] == pytest.approx(1.4826)


def test_gap_exclusion_and_global_report_do_not_cross_quarters(tmp_path):
    metrics = residual_metrics([1, 2, 3, 4], [0, 1, 10, 11], 1.0)
    assert metrics["lag1_pair_count"] == 2
    assert metrics["lag1_pearson_correlation"] == pytest.approx(1.0)

    _, global_path, report = write_residual_report(
        [(1, [0, 1], [1, 2]), (2, [100, 101], [3, 4])], tmp_path,
        exposure=1.0, sample_index=7, stored_lnlike=-4.0,
    )
    assert report["global_metrics"]["lag1_pair_count"] == 2
    assert json.loads(global_path.read_text())["scope"].endswith("not posterior predictive.")


def test_duplicate_times_within_a_quarter_are_rejected():
    with pytest.raises(ValueError, match="must not contain duplicates"):
        residual_metrics([1, 2], [0, 0], 1.0)


@pytest.mark.parametrize(
    ("residuals", "times", "status"),
    [([1], [0], "undefined_inadequate_pairs"), ([2, 2, 2], [0, 1, 2], "undefined_zero_variance")],
)
def test_short_and_constant_series_have_explicit_undefined_lag1(residuals, times, status):
    assert residual_metrics(residuals, times, 1.0)["lag1_correlation_status"] == status


@pytest.mark.parametrize(
    ("residuals", "times", "exposure", "message"),
    [([1, np.nan], [0, 1], 1.0, "residuals must be finite"), ([1], [0, 1], 1.0, "same length"), ([1], [0], 0.0, "positive and finite")],
)
def test_invalid_inputs_are_rejected(residuals, times, exposure, message):
    with pytest.raises(ValueError, match=message):
        residual_metrics(residuals, times, exposure)
