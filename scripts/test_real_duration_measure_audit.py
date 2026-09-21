import numpy as np
import pandas as pd
import pytest

from duration_measure_control import duration_from_density
from extract_eccentricity_posteriors_direct import direct_importance_posterior
from real_duration_measure_audit import compare


def test_likelihood_only_replay_matches_existing_extractor():
    t = float(duration_from_density(1., 0., 0.))
    samples = pd.DataFrame({'LN_WT': [0.], 'DUR14_0': [t],
                            'ROR_0': [.05], 'IMPACT_0': [.5]})
    replay = compare(samples, np.array([12.]), 0, 1., .1, .1, 3, 30000)
    existing = direct_importance_posterior(
        np.array([.05]), np.array([.5]), np.array([t]), np.array([1.]),
        np.array([12.]), 1., .1, .1, n_proposals=30000, e_max=.95,
        e_grid=np.linspace(0, .95, 300),
        omega_grid=np.linspace(-np.pi/2, 3*np.pi/2, 240, endpoint=False),
        density_error_mode='symmetric-average', seed=3,
    )
    np.testing.assert_allclose([replay[0][k] for k in ['e16', 'e50', 'e84']],
                               existing['e_quantiles'], rtol=1e-10)
    assert replay[0]['valid_fraction'] == 1.
    assert replay[1]['ess'] > 0
    assert replay[2]['mode'] == 'uniform_focal_radius_only'
    np.testing.assert_allclose([replay[2][k] for k in ['e16', 'e50', 'e84']],
                               existing['e_quantiles'], rtol=1e-10)


def test_population_masses_preserve_normalization_and_prior_ratio():
    t = float(duration_from_density(1., 0., 0.))
    samples = pd.DataFrame({'LN_WT': [0.], 'DUR14_0': [t],
                            'ROR_0': [.05], 'IMPACT_0': [.5], 'ROR_1': [.1]})
    grid = np.linspace(0, .95, 300)
    result = compare(samples, np.array([12.]), 0, 1., .1, .1, 3, 10000, grid)
    assert result[-1]['mode'] == 'uniform_all_radii'
    for row in result:
        assert np.isclose(row['_masses']['none'].sum(), 1.)
        assert all(np.all(np.isfinite(v) & (v >= 0)) for v in row['_masses'].values())
    for selection in result[0]['_masses']:
        np.testing.assert_allclose(result[0]['_masses'][selection],
                                   result[-1]['_masses'][selection], atol=1e-12)


def test_sobol_preserves_uniform_proposal_and_requires_balanced_count():
    t = float(duration_from_density(1., 0., 0.))
    samples = pd.DataFrame({'LN_WT': [0.], 'DUR14_0': [t],
                            'ROR_0': [.05], 'IMPACT_0': [.5]})
    args = (samples, np.array([12.]), 0, 1., 1e100, 1e100, 3)
    result = compare(*args, 4096, proposal_mode='sobol')
    assert abs(result[0]['e50']-.475) < .95/4096
    assert compare(*args, 4096, proposal_mode='sobol') == result
    with pytest.raises(ValueError, match='power of two'):
        compare(*args, 5000, proposal_mode='sobol')


def test_recorded_rms_and_support_replay_matches_extractor():
    t = float(duration_from_density(1., 0., 0.))
    samples = pd.DataFrame({'LN_WT': [0.], 'DUR14_0': [t],
                            'ROR_0': [.05], 'IMPACT_0': [.5]})
    result = compare(samples, np.array([12.]), 0, 1., .2, .1, 7, 10000,
                     e_max=.6, density_error_mode='symmetric-rms')
    existing = direct_importance_posterior(
        np.array([.05]), np.array([.5]), np.array([t]), np.array([1.]),
        np.array([12.]), 1., .2, .1, n_proposals=10000, e_max=.6,
        e_grid=np.linspace(0, .95, 300),
        omega_grid=np.linspace(-np.pi/2, 3*np.pi/2, 240, endpoint=False),
        density_error_mode='symmetric-rms', seed=7)
    np.testing.assert_allclose([result[0][k] for k in ['e16', 'e50', 'e84']],
                               existing['e_quantiles'], rtol=1e-10)
