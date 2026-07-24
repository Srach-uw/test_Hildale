import pandas as pd

from build_missing_recovery_feasibility import build


def test_recovery_classes_preserve_system_and_planet_totals():
    targets = pd.DataFrame(
        {
            "koi_target": ["K00001", "K00002", "K00003", "K00004"],
            "disk": ["thick", "thick", "thin", "thin"],
            "system": ["multi", "single", "multi", "single"],
            "missing_planets": [2, 1, 3, 1],
            "max_snr": [30.0, 10.0, 25.0, 5.0],
            "prior_launch_status": [
                "not_in_previous_592_target_launch",
                "not_in_previous_592_target_launch",
                "previously_launched_no_usable_result",
                "previously_launched_no_usable_result",
            ],
        }
    )
    detail, summary = build(targets)
    assert len(detail) == 4
    assert summary.systems.sum() == 4
    assert summary.planets.sum() == 7
    assert set(detail.recovery_class) == {
        "A_never_attempted_snr_ge_20",
        "B_never_attempted_snr_lt_20",
        "C_prior_failed_snr_ge_20",
        "D_prior_failed_snr_lt_20",
    }
