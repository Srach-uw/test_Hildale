import pandas as pd
import pytest

from visual_qc_contract_audit import (
    rank_unpublished_rejection_candidates,
    validate_contract,
)


def make_manifest() -> pd.DataFrame:
    rows = []
    counts = {
        ("thin", "single"): 1121,
        ("thick", "single"): 275,
        ("thin", "multi"): 883,
        ("thick", "multi"): 212,
    }
    row = 0
    for (disk, system), count in counts.items():
        for index in range(count):
            row += 1
            rows.append(
                {
                    "kepid": row,
                    "kepoi_name": f"K{row:05d}.01",
                    "koi_target": f"K{row:05d}",
                    "disk": disk,
                    "system": system,
                    "posterior_status": (
                        "missing_after_uniform_assembly"
                        if system == "multi" and index < 30
                        else "posterior_available"
                    ),
                    "qc_primary_exclude": False,
                    "koi_impact": 0.5,
                    "koi_model_snr": 20.0,
                    "koi_num_transits": 10,
                }
            )
    frame = pd.DataFrame(rows)
    # Reuse host identifiers so the synthetic inventory has 1,888 hosts and
    # 2,491 planets, matching the published host/planet relationship.
    frame["kepid"] = (frame.index % 1888) + 1
    return frame


def test_contract_reconstructs_published_final_rejections():
    audit = validate_contract(make_manifest())
    removals = {
        (row.disk, row.system): row.implied_visual_or_convergence_rejections
        for row in audit.itertuples()
    }
    assert removals == {
        ("thin", "single"): 0,
        ("thick", "single"): 0,
        ("thin", "multi"): 21,
        ("thick", "multi"): 5,
    }


def test_candidate_ranking_is_sensitivity_only_and_count_constrained():
    frame = make_manifest()
    audit = validate_contract(frame)
    candidates = rank_unpublished_rejection_candidates(frame, audit)
    selected = candidates.loc[candidates["count_constrained_candidate_set"]]
    assert len(selected) == 26
    assert selected.groupby("disk").size().to_dict() == {"thick": 5, "thin": 21}
    assert set(selected["identity_status"]) == {
        "sensitivity_candidate_not_author_confirmed"
    }
    assert selected["missing_posterior"].all()


def test_contract_rejects_wrong_host_or_planet_inventory():
    frame = make_manifest().iloc[:-1].copy()
    with pytest.raises(ValueError):
        validate_contract(frame)
