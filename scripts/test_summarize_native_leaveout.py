import csv
import json

import pytest

from summarize_native_leaveout import summarize, validate


QUARTERS = [1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 13, 14, 16, 17]
EXCLUSIONS = {
    "all": "",
    "without_q3": "3",
    "without_q6": "6",
    "without_q10": "10",
    "without_season0": "4,8,12,16",
    "without_season1": "1,5,9,13,17",
    "without_season2": "2,6,10,14",
    "without_season3": "3",
}
FIELDS = [
    "case",
    "excluded_quarters",
    "impact",
    "points",
    "lnlike",
    "success",
    "converged_starts",
    "parameter",
]


def _valid_rows():
    rows = []
    for index, (case, exclusions) in enumerate(EXCLUSIONS.items()):
        points = 4000 - index
        for impact in ("0.0", "0.96"):
            high = impact == "0.96"
            rows.append(
                {
                    "case": case,
                    "excluded_quarters": exclusions,
                    "impact": impact,
                    "points": points,
                    "lnlike": -100.0 + high,
                    "success": "True",
                    "converged_starts": 3,
                    "parameter": 2.0 + high,
                }
            )
    return rows


def _write_inputs(path, rows=None, replay_gate_passed=True):
    path.mkdir()
    (path / "provenance.json").write_text(
        json.dumps({"replay_gate_passed": replay_gate_passed, "quarters": QUARTERS}),
        encoding="utf-8",
    )
    with (path / "leave_out.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(_valid_rows() if rows is None else rows)


def test_writes_scoped_report_and_high_minus_low_differences(tmp_path):
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    _write_inputs(input_dir)

    report_path, differences_path = summarize(input_dir, output_dir)

    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["status"] == "passed"
    assert report["validation"]["case_count"] == 8
    assert "not posterior or evidence" in report["interpretation_limit"]
    with differences_path.open(newline="", encoding="utf-8") as handle:
        differences = list(csv.DictReader(handle))
    assert len(differences) == 8
    assert float(differences[0]["lnlike_high_minus_low"]) == pytest.approx(1.0)
    assert float(differences[0]["parameter_high_minus_low"]) == pytest.approx(1.0)


def test_missing_impact_pair_is_rejected(tmp_path):
    input_dir = tmp_path / "input"
    rows = _valid_rows()
    rows.pop(1)
    _write_inputs(input_dir, rows=rows)

    with pytest.raises(ValueError, match="missing impact pair"):
        validate(input_dir)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [("converged_starts", 2, "converged_starts"), ("success", "False", "success")],
)
def test_failed_starts_are_rejected(tmp_path, field, value, message):
    input_dir = tmp_path / "input"
    rows = _valid_rows()
    rows[0][field] = value
    _write_inputs(input_dir, rows=rows)

    with pytest.raises(ValueError, match=message):
        validate(input_dir)


def test_bad_replay_gate_is_rejected(tmp_path):
    input_dir = tmp_path / "input"
    _write_inputs(input_dir, replay_gate_passed=False)

    with pytest.raises(ValueError, match="replay_gate_passed"):
        validate(input_dir)


def test_duplicate_case_impact_is_rejected(tmp_path):
    input_dir = tmp_path / "input"
    rows = _valid_rows()
    rows.append(rows[0].copy())
    _write_inputs(input_dir, rows=rows)

    with pytest.raises(ValueError, match="duplicate row"):
        validate(input_dir)
