from source_contract_matrix import build_matrix


def test_matrix_contains_the_high_risk_unresolved_contracts():
    matrix = build_matrix().set_index("contract")
    assert matrix.loc["posterior export weighting", "status"] == "critical unresolved fork"
    assert matrix.loc["posterior coverage", "risk"] == "critical"
    assert matrix.loc["stellar density prior", "status"] == "unresolved source ambiguity"


def test_matrix_separates_published_targets_from_local_coverage():
    matrix = build_matrix()
    assert "2465 planets / 1888 stars" in matrix.loc[matrix["contract"] == "published planet total", "sagear_statement"].iloc[0]
    assert "547 systems" in matrix.loc[matrix["contract"] == "posterior coverage", "local_implementation"].iloc[0]
