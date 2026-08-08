import pandas as pd

from sagear_full_postfit_qc_sensitivity import build_qc_ledger


def test_build_qc_ledger_applies_documented_thresholds():
    summary = pd.DataFrame(
        {
            "kepoi_name": ["K1.01", "K2.01"],
            "rho_err_hi_solar": [0.1, 0.2],
            "rho_err_lo_solar": [0.1, 0.1],
        }
    )
    postfit = pd.DataFrame(
        {
            "kepoi_name": ["K1.01", "K2.01"],
            "gilbert_qc_available": [True, True],
            "nested_grazing_fraction": [0.01, 0.06],
            "planet_radius_fractional_uncertainty_approx": [0.1, 0.1],
        }
    )
    viability = pd.DataFrame(
        {
            "kepoi_name": ["K1.01", "K2.01"],
            "intended_guard_raises": [False, False],
            "status": ["ok", "ok"],
        }
    )
    ledger = build_qc_ledger(summary, postfit, viability).set_index("kepoi_name")
    assert not bool(ledger.loc["K1.01", "full_postfit_exclude"])
    assert bool(ledger.loc["K2.01", "exclude_grazing"])
    assert bool(ledger.loc["K2.01", "exclude_density_asymmetry"])
