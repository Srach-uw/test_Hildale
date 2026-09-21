"""Validate and summarize the K02712 native leave-out conditional refits."""

import argparse
import csv
import hashlib
import json
import math
from decimal import Decimal, InvalidOperation
from pathlib import Path


EXPECTED_CASES = (
    "all",
    "without_q3",
    "without_q6",
    "without_q10",
    "without_season0",
    "without_season1",
    "without_season2",
    "without_season3",
)
EXPECTED_IMPACTS = (Decimal("0"), Decimal("0.96"))
NON_DIFFERENCE_COLUMNS = {
    "case",
    "excluded_quarters",
    "impact",
    "points",
    "success",
    "converged_starts",
}


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _parse_quarters(value, context):
    if not value.strip():
        return ()
    try:
        quarters = tuple(int(item.strip()) for item in value.split(","))
    except ValueError as exc:
        raise ValueError(f"{context} has invalid excluded_quarters") from exc
    if len(quarters) != len(set(quarters)):
        raise ValueError(f"{context} has duplicate excluded quarters")
    return tuple(sorted(quarters))


def _expected_exclusions(provenance_quarters):
    quarter_set = set(provenance_quarters)
    expected = {
        "all": (),
        "without_q3": (3,),
        "without_q6": (6,),
        "without_q10": (10,),
    }
    for season in range(4):
        expected[f"without_season{season}"] = tuple(
            quarter for quarter in provenance_quarters if quarter % 4 == season
        )
    for case in ("without_q3", "without_q6", "without_q10"):
        quarter = expected[case][0]
        if quarter not in quarter_set:
            raise ValueError(f"provenance quarters do not contain Q{quarter}")
    return expected


def _read_inputs(input_dir):
    leave_out_path = input_dir / "leave_out.csv"
    provenance_path = input_dir / "provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    if provenance.get("replay_gate_passed") is not True:
        raise ValueError("provenance replay_gate_passed must be true")

    raw_quarters = provenance.get("quarters")
    if not isinstance(raw_quarters, list) or not raw_quarters:
        raise ValueError("provenance quarters must be a non-empty list")
    if any(type(value) is not int for value in raw_quarters):
        raise ValueError("provenance quarters must contain integers")
    quarters = tuple(sorted(raw_quarters))
    if len(quarters) != len(set(quarters)):
        raise ValueError("provenance quarters contain duplicates")

    with leave_out_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = NON_DIFFERENCE_COLUMNS | {"lnlike"}
        missing_columns = sorted(required - set(reader.fieldnames or ()))
        if missing_columns:
            raise ValueError(f"leave_out.csv is missing columns: {missing_columns}")
        rows = list(reader)
    return leave_out_path, provenance_path, provenance, quarters, reader.fieldnames, rows


def validate(input_dir):
    (
        leave_out_path,
        provenance_path,
        provenance,
        quarters,
        fieldnames,
        raw_rows,
    ) = _read_inputs(Path(input_dir))
    expected_exclusions = _expected_exclusions(quarters)
    difference_fields = [name for name in fieldnames if name not in NON_DIFFERENCE_COLUMNS]
    rows_by_case = {case: {} for case in EXPECTED_CASES}

    for line_number, raw in enumerate(raw_rows, start=2):
        case = raw["case"].strip()
        if case not in rows_by_case:
            raise ValueError(f"line {line_number} has unexpected case {case!r}")
        try:
            impact = Decimal(raw["impact"])
        except InvalidOperation as exc:
            raise ValueError(f"line {line_number} has invalid impact") from exc
        if impact not in EXPECTED_IMPACTS:
            raise ValueError(f"line {line_number} has unexpected impact {impact}")
        if impact in rows_by_case[case]:
            raise ValueError(f"duplicate row for {case} at impact {impact}")

        try:
            points = int(raw["points"])
            converged_starts = int(raw["converged_starts"])
        except ValueError as exc:
            raise ValueError(f"line {line_number} has invalid integer fields") from exc
        if points <= 0:
            raise ValueError(f"line {line_number} points must be positive")
        if raw["success"].strip().lower() != "true":
            raise ValueError(f"line {line_number} success must be true")
        if converged_starts != 3:
            raise ValueError(f"line {line_number} converged_starts must equal 3")

        numeric = {}
        for field in difference_fields:
            try:
                value = float(raw[field])
            except ValueError as exc:
                raise ValueError(f"line {line_number} has non-numeric {field}") from exc
            if not math.isfinite(value):
                label = "likelihood" if field == "lnlike" else field
                raise ValueError(f"line {line_number} has non-finite {label}")
            numeric[field] = value

        exclusions = _parse_quarters(raw["excluded_quarters"], f"line {line_number}")
        rows_by_case[case][impact] = {
            "points": points,
            "excluded_quarters": exclusions,
            "numeric": numeric,
        }

    case_results = []
    difference_rows = []
    for case in EXPECTED_CASES:
        pair = rows_by_case[case]
        if set(pair) != set(EXPECTED_IMPACTS):
            missing = sorted(str(value) for value in set(EXPECTED_IMPACTS) - set(pair))
            raise ValueError(f"{case} is missing impact pair members: {missing}")
        low, high = (pair[impact] for impact in EXPECTED_IMPACTS)
        if low["points"] != high["points"]:
            raise ValueError(f"{case} has inconsistent points within its impact pair")
        if low["excluded_quarters"] != high["excluded_quarters"]:
            raise ValueError(f"{case} has inconsistent exclusions within its impact pair")
        if low["excluded_quarters"] != expected_exclusions[case]:
            raise ValueError(
                f"{case} exclusions {low['excluded_quarters']} do not match "
                f"provenance-derived {expected_exclusions[case]}"
            )

        case_results.append(
            {
                "case": case,
                "excluded_quarters": list(low["excluded_quarters"]),
                "points": low["points"],
                "lnlike_by_impact": {
                    "0": low["numeric"]["lnlike"],
                    "0.96": high["numeric"]["lnlike"],
                },
                "lnlike_high_minus_low": (
                    high["numeric"]["lnlike"] - low["numeric"]["lnlike"]
                ),
            }
        )
        difference_row = {
            "case": case,
            "excluded_quarters": ",".join(map(str, low["excluded_quarters"])),
            "points": low["points"],
        }
        for field in difference_fields:
            difference_row[f"{field}_high_minus_low"] = (
                high["numeric"][field] - low["numeric"][field]
            )
        difference_rows.append(difference_row)

    report = {
        "status": "passed",
        "scope": "Two fixed-impact conditional likelihood refits per leave-out case (impact 0 and 0.96).",
        "interpretation_limit": (
            "These are not posterior or evidence comparisons and do not prove general robustness."
        ),
        "inputs": {
            "leave_out_csv": leave_out_path.name,
            "leave_out_sha256": _sha256(leave_out_path),
            "provenance_json": provenance_path.name,
            "provenance_sha256": _sha256(provenance_path),
        },
        "validation": {
            "replay_gate_passed": provenance["replay_gate_passed"],
            "impacts": [0, 0.96],
            "required_converged_starts": 3,
            "case_count": len(case_results),
            "provenance_quarters": list(quarters),
        },
        "cases": case_results,
    }
    return report, difference_fields, difference_rows


def summarize(input_dir, output_dir=None):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir) if output_dir is not None else input_dir
    report, difference_fields, difference_rows = validate(input_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "native_leaveout_report.json"
    differences_path = output_dir / "native_leaveout_high_minus_low.csv"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    columns = ["case", "excluded_quarters", "points"] + [
        f"{field}_high_minus_low" for field in difference_fields
    ]
    with differences_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(difference_rows)
    return report_path, differences_path


def main():
    default_dir = (
        Path(__file__).resolve().parents[1]
        / "metadata"
        / "inference_assumptions_20260915"
        / "K02712_native_leaveout_20260919"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=default_dir)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    try:
        outputs = summarize(args.input_dir, args.output_dir)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        parser.error(str(exc))
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
