import numpy as np
import pytest

from audit_population_seed_leverage import contrast
from replace_likelihood_convergence import replace_row


def test_contrast_ignores_per_planet_normalization():
    grid = np.linspace(0, .95, 300)
    mass = np.ones((2, len(grid)))
    expected = contrast(mass, grid, .06, .11, 'manuscript_reciprocal')
    np.testing.assert_allclose(expected, contrast(mass*np.array([[.001], [10.]]),
                                                grid, .06, .11, 'manuscript_reciprocal'), atol=1e-12)


def test_replacement_preserves_other_rows_and_source():
    planets = np.array(['a', 'b'])
    matrix = np.array([[1., 2.], [3., 4.]])
    result = replace_row(planets, matrix, 'b', np.array([5., 6.]))
    np.testing.assert_array_equal(result, [[1., 2.], [5., 6.]])
    np.testing.assert_array_equal(matrix, [[1., 2.], [3., 4.]])
    with pytest.raises(ValueError):
        replace_row(planets, matrix, 'c', np.array([5., 6.]))


def test_contrast_remains_finite_for_extreme_tail():
    grid = np.linspace(0, .95, 300)
    mass = np.zeros((1, len(grid)))
    mass[0, -1] = 1.
    value = contrast(mass, grid, .0175, .2, 'manuscript_reciprocal')[0]
    assert np.isfinite(value) and value > 1000
