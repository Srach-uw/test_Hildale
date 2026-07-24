from pathlib import Path

from article_data_provenance_audit import EXPECTED_SHA256, audit


def test_journal_host_table_is_canonical() -> None:
    root = Path(__file__).resolve().parents[1]
    canonical = root / "reference" / "data" / "sagear2026_table1_kinematic_hosts_mrt.txt"
    result = audit(canonical, canonical)
    assert result["official_sha256"] == EXPECTED_SHA256
    assert result["byte_identical"] is True
