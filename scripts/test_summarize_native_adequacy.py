import csv

import pytest

from summarize_native_adequacy import summarize


FIELDS = ["component", "points", "low_lnlike", "high_lnlike", "high_minus_low"]


def _write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def test_derives_photometric_chi2_and_excludes_ld(tmp_path):
    source = tmp_path / "quarter_contrast.csv"
    output = tmp_path / "adequacy.csv"
    _write(
        source,
        [
            {"component": "ld_penalty", "points": 2, "low_lnlike": -9, "high_lnlike": -1},
            {"component": "3", "points": 4, "low_lnlike": -10, "high_lnlike": -8},
            {"component": "9", "points": 2, "low_lnlike": -8, "high_lnlike": -9},
        ],
    )

    rows = summarize(source, output, convention_attested=True)

    assert [row["quarter"] for row in rows] == [3, 9]
    assert rows[0]["low_chi2_per_point"] == pytest.approx(5.0)
    assert rows[0]["high_minus_low_lnlike"] == pytest.approx(2.0)
    assert rows[1]["low_worst_rank"] == 1
    assert rows[1]["high_worst_rank"] == 1
    assert len(list(csv.DictReader(output.open(encoding="utf-8")))) == 2


def test_duplicate_quarter_is_rejected(tmp_path):
    source = tmp_path / "quarter_contrast.csv"
    row = {"component": "3", "points": 4, "low_lnlike": -10, "high_lnlike": -8}
    _write(source, [row, row])

    with pytest.raises(ValueError, match="duplicate photometric quarter"):
        summarize(source, tmp_path / "out.csv", convention_attested=True)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [("points", 0, "points must be positive"), ("low_lnlike", "nan", "likelihood must be finite")],
)
def test_invalid_photometric_rows_are_rejected(tmp_path, field, value, message):
    source = tmp_path / "quarter_contrast.csv"
    row = {"component": "3", "points": 4, "low_lnlike": -10, "high_lnlike": -8}
    row[field] = value
    _write(source, [row])

    with pytest.raises(ValueError, match=message):
        summarize(source, tmp_path / "out.csv", convention_attested=True)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("component", 18, "quarter must be in 0..17"),
        ("low_lnlike", 1, "positive term"),
    ],
)
def test_non_photometric_terms_are_rejected(tmp_path, field, value, message):
    source = tmp_path / "quarter_contrast.csv"
    row = {"component": "3", "points": 4, "low_lnlike": -10, "high_lnlike": -8}
    row[field] = value
    _write(source, [row])

    with pytest.raises(ValueError, match=message):
        summarize(source, tmp_path / "out.csv", convention_attested=True)


def test_requires_gate_or_explicit_convention_attestation(tmp_path):
    source = tmp_path / "quarter_contrast.csv"
    _write(source, [{"component": "3", "points": 4, "low_lnlike": -10, "high_lnlike": -8}])

    with pytest.raises(ValueError, match="requires replay-gated provenance"):
        summarize(source, tmp_path / "out.csv")


def test_rejects_failed_replay_gate(tmp_path):
    source = tmp_path / "quarter_contrast.csv"
    provenance = tmp_path / "provenance.json"
    _write(source, [{"component": "3", "points": 4, "low_lnlike": -10, "high_lnlike": -8}])
    provenance.write_text('{"replay_gate_passed": false}', encoding="utf-8")

    with pytest.raises(ValueError, match="replay_gate_passed"):
        summarize(source, tmp_path / "out.csv", provenance_json=provenance)
