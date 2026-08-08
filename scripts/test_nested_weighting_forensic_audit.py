import pandas as pd

from nested_weighting_forensic_audit import comparison_rows


def test_default_root_is_repository_root():
    from nested_weighting_forensic_audit import ROOT

    assert (ROOT / "scripts" / "nested_weighting_forensic_audit.py").is_file()


def test_comparison_rows_uses_half_gaussian_sigma():
    rayleigh = pd.DataFrame(
        {
            "population": ["thin_multis"],
            "n": [10],
            "expected_e": [0.03],
            "expected_e_lo": [0.02],
            "expected_e_hi": [0.04],
        }
    )
    shapes = pd.DataFrame(
        {
            "population": ["thin_multis"],
            "n": [10],
            "model": ["half_gaussian"],
            "mean_e": [0.02],
            "sigma": [0.033],
        }
    )
    rows = comparison_rows("test", rayleigh, shapes)
    assert rows[1]["local_value"] == 0.033
