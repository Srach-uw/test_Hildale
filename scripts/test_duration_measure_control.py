import numpy as np

from duration_measure_control import density_and_derivative, duration_from_density


def test_density_derivative_matches_finite_difference():
    for e in [0.0, 0.3, 0.8]:
        for w in [0.0, np.pi/2, 3*np.pi/2]:
            t = duration_from_density(1.0, e, w)
            rho, derivative = density_and_derivative(t, e, w)
            eps = t*1e-5
            plus, _ = density_and_derivative(t+eps, e, w)
            minus, _ = density_and_derivative(t-eps, e, w)
            np.testing.assert_allclose(rho, 1.0, rtol=1e-12)
            np.testing.assert_allclose(derivative, -(plus-minus)/(2*eps), rtol=1e-8)
