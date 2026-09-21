import numpy as np
import pytest
from prepare_matched_lc_processing import match_unique_times


def test_exact_alignment_preserves_correspondence():
    raw = np.array([3., 1., 2.])
    processed = np.array([2., 3., 4.])
    common, ri, pi = match_unique_times(raw, processed)
    np.testing.assert_array_equal(common, [2., 3.])
    np.testing.assert_array_equal(raw[ri], processed[pi])


@pytest.mark.parametrize('bad', [[1., 1.], [1., np.nan]])
def test_ambiguous_times_fail(bad):
    with pytest.raises(ValueError):
        match_unique_times(np.array(bad), np.array([1., 2.]))


def test_nearby_is_not_exact():
    with pytest.raises(ValueError):
        match_unique_times(np.array([1.]), np.array([1. + 1e-10]))
