from common import load_config, read_sagear2026_kinematic_hosts
from published_sagear_audit import EXPECTED, validate_published_hosts


def test_published_table_count_and_threshold_contract() -> None:
    hosts = read_sagear2026_kinematic_hosts(load_config())
    assert validate_published_hosts(hosts) == EXPECTED
    assert hosts["kepid"].is_unique
    assert hosts["p_thick_published"].between(0, 1).all()
