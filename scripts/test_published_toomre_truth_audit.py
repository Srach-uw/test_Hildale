import numpy as np
import pandas as pd

from published_toomre_truth_audit import summarize_coordinates


def test_astropy_coordinate_summary_is_tighter_than_direct_proxy():
    published = pd.DataFrame(
        {
            "vphi_used_kms": [-220.0, -210.0],
            "vr_used_kms": [3.0, -4.0],
            "vz_used_kms": [4.0, 3.0],
            "V_phi": [230.0, 195.0],
            "V_perp": [20.0, 30.0],
            "V_phi_astropy": [-220.0, -210.0],
            "V_perp_astropy": [5.0, 5.0],
            "V_phi_geom": [-219.8, -210.2],
            "V_perp_geom": [5.1, 4.9],
        }
    )
    s = summarize_coordinates(published)
    direct = s[(s.source == "direct_angus") & (s.axis == "vphi")].iloc[0]
    astropy = s[(s.source == "old_astropy") & (s.axis == "vphi")].iloc[0]
    assert direct.median_abs_delta_kms > astropy.median_abs_delta_kms
    assert astropy.median_abs_delta_kms == 0


def test_published_perpendicular_speed_uses_vr_and_vz():
    d = pd.DataFrame(
        {
            "vphi_used_kms": [-220.0],
            "vr_used_kms": [3.0],
            "vz_used_kms": [4.0],
            "V_phi": [220.0],
            "V_perp": [99.0],
            "V_phi_astropy": [-220.0],
            "V_perp_astropy": [5.0],
            "V_phi_geom": [-220.0],
            "V_perp_geom": [5.0],
        }
    )
    s = summarize_coordinates(d)
    row = s[(s.source == "old_astropy") & (s.axis == "vperp")].iloc[0]
    assert np.isclose(row.median_abs_delta_kms, 0.0)
