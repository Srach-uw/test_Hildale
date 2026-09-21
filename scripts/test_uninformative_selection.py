import numpy as np
from audit_uninformative_selection import null_contrasts


def test_null_measurement_preserves_forward_population_prior():
    curves=null_contrasts()
    assert np.ptp(curves['legacy_forward_norm']) < 1e-12
    assert np.ptp(curves['none']) < 1e-12


def test_reciprocal_null_measurement_has_cold_preference():
    curve=null_contrasts()['manuscript_reciprocal']
    assert curve[-1]-curve[0] < -.01
