import pandas as pd

from common import load_config, read_sagear2026_kinematic_hosts
from published_sagear_audit import EXPECTED, relabel_planets, validate_published_hosts


def test_published_table_count_and_threshold_contract() -> None:
    hosts = read_sagear2026_kinematic_hosts(load_config())
    assert validate_published_hosts(hosts) == EXPECTED
    assert hosts["kepid"].is_unique
    assert hosts["p_thick_published"].between(0, 1).all()


def test_published_disk_relabel_preserves_pre_cut_multiplicity(tmp_path) -> None:
    hosts = pd.DataFrame(
        {
            "kepid": [100, 200],
            "disk_published": ["thin", "thick"],
            "p_thick_published": [0.1, 0.9],
        }
    )
    sample = pd.DataFrame(
        {
            "kepid": [100, 200],
            "kepoi_name": ["K00100.01", "K00200.01"],
            "disk": ["thick", "thin"],
            "system": ["multi", "single"],
        }
    )
    sample_path = tmp_path / "sample.csv"
    sample.to_csv(sample_path, index=False)

    relabeled, counts = relabel_planets(hosts, sample_path)

    by_kic = relabeled.set_index("kepid")
    assert by_kic.loc[100, "system_published_overlap"] == "multi"
    assert by_kic.loc[100, "system_overlap_recount"] == "single"
    assert bool(by_kic.loc[100, "multiplicity_recount_disagrees"])
    assert counts["published_relabel_thin_multi_planets"] == 1
    assert counts["published_relabel_multiplicity_recount_disagreements"] == 1
