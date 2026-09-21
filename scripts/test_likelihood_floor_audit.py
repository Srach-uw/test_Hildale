import numpy as np
import pytest
from audit_likelihood_floor import stable_terms
from hierarchical_rayleigh import rayleigh_grid


@pytest.mark.parametrize('mode', ['none', 'manuscript_reciprocal', 'legacy_forward_norm'])
def test_agrees_when_representable(mode):
    grid = np.linspace(0, .95, 300)
    matrix = np.random.default_rng(7).random((3, 300))
    sigmas = np.array([.1, .3])
    rays, norm = rayleigh_grid(grid, sigmas, mode != 'none', mode)
    direct = matrix @ rays
    if mode == 'legacy_forward_norm':
        direct /= norm
    np.testing.assert_allclose(stable_terms(matrix, grid, sigmas, mode), np.log(direct), atol=1e-12)


def test_extreme_point_mass_ratio():
    grid = np.array([.01, .2, .8, .95])
    matrix = np.array([[0., 1., 0., 0.], [0., 0., 1., 0.]])
    terms = stable_terms(matrix, grid, np.array([.01]), 'none')[:, 0]
    expected = np.log(.8/.2)-(.8**2-.2**2)/(2*.01**2)
    assert terms[1] < -3000
    assert terms[1]-terms[0] == pytest.approx(expected)


def test_zero_support_rejected():
    with pytest.raises(ValueError, match='support'):
        stable_terms(np.zeros((1, 3)), np.array([0., .4, .8]), np.array([.1]), 'none')
