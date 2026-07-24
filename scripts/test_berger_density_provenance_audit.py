from pathlib import Path

import pytest

from berger_density_provenance_audit import parse_b18, parse_b20


ROOT = Path(__file__).resolve().parents[1]


def test_berger2020_table_has_published_density_fields():
    path = ROOT / "data" / "berger2020_table2.dat.gz"
    if not path.exists():
        pytest.skip("optional Berger 2020 source table is not distributed with the repository")
    table = parse_b20(path)
    assert len(table) == 186301
    assert {"rho_log", "rho_log_upper", "rho_log_lower"}.issubset(table.columns)


def test_berger2018_table_has_radius_but_no_density():
    path = ROOT / "data" / "berger2018_table1_min.tsv"
    if not path.exists():
        pytest.skip("optional Berger 2018 source table is not distributed with the repository")
    table = parse_b18(path)
    assert len(table) > 100000
    assert {"radius_b18", "evol_b18", "bin_b18"}.issubset(table.columns)
    assert "rho_log" not in table.columns


def test_official_berger2018_fixed_width_table_has_radius_but_no_density():
    path = ROOT / "data" / "berger2018_table1.dat.gz"
    if not path.exists():
        pytest.skip("optional Berger 2018 source table is not distributed with the repository")
    table = parse_b18(path)
    assert len(table) == 177911
    assert {"radius_b18", "evol_b18", "bin_b18"}.issubset(table.columns)
    assert "rho_log" not in table.columns
