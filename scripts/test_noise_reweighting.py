import numpy as np
from scipy.special import softmax
from audit_noise_reweighting import reweight


def test_unit_scale_preserves_weights():
    np.testing.assert_allclose(reweight([-2,-3],[-8,-5],[-1,-2],1),softmax([-2,-3]))


def test_limb_penalty_not_tempered():
    # Equal photometric likelihoods must preserve original posterior odds.
    np.testing.assert_allclose(reweight([-2,-3],[-11,-12],[-1,-2],2),softmax([-2,-3]))


def test_analytic_likelihood_ratio():
    w = reweight(np.log([.4,.6]),[-10,-6],[0,0],2)
    np.testing.assert_allclose(w[0]/w[1],(.4/.6)*np.exp(3))
