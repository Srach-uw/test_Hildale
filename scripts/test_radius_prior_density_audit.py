import numpy as np
import pytest

from radius_prior_density_audit import prior_weights


def test_joint_radius_ratio_uses_all_paired_planets():
    actual = prior_weights(np.array([.2, .8]), np.array([[.1, .2], [.3, .4]]))
    expected = np.array([.2*.1*.2, .8*.3*.4])
    np.testing.assert_allclose(actual, expected/expected.sum())


def test_loguniform_to_uniform_recovers_uniform_mean():
    r = np.exp(np.linspace(np.log(1e-5), np.log(.99), 100001))
    weights = np.ones(len(r))
    weights[[0, -1]] = .5
    weights /= weights.sum()
    corrected = prior_weights(weights, r[:, None].clip(1e-5, .99))
    assert abs(corrected @ r - (.99+1e-5)/2) < 1e-7


def test_radius_support_is_checked():
    with pytest.raises(ValueError, match='support'):
        prior_weights(np.array([1.]), np.array([[0.]]))
