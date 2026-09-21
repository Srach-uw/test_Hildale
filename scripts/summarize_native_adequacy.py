"""Derive K02712 photometric chi-square per point from native quarter terms."""

import argparse
import csv
import json
import math
from pathlib import Path


REQUIRED_COLUMNS = {"component", "points", "low_lnlike", "high_lnlike"}


def summarize(input_csv, output_csv, *, provenance_json=None, convention_attested=False):
    input_csv = Path(input_csv)
    output_csv = Path(output_csv)
    if provenance_json is not None:
        provenance = json.loads(Path(provenance_json).read_text(encoding="utf-8"))
        if provenance.get("replay_gate_passed") is not True:
            raise ValueError("provenance replay_gate_passed must be true")
    elif not convention_attested:
        raise ValueError(
            "native photometric-term convention requires replay-gated provenance "
            "or explicit attestation"
        )
    with input_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = sorted(REQUIRED_COLUMNS - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"quarter contrast is missing columns: {missing}")
        raw_rows = list(reader)

    rows = []
    seen_quarters = set()
    for line_number, raw in enumerate(raw_rows, start=2):
        component = raw["component"].strip()
        if component == "ld_penalty":
            continue
        try:
            quarter = int(component)
            points = int(raw["points"])
            low_lnlike = float(raw["low_lnlike"])
            high_lnlike = float(raw["high_lnlike"])
        except ValueError as exc:
            raise ValueError(f"line {line_number} has invalid numeric data") from exc
        if not 0 <= quarter <= 17:
            raise ValueError(f"line {line_number} quarter must be in 0..17")
        if quarter in seen_quarters:
            raise ValueError(f"duplicate photometric quarter Q{quarter}")
        if points <= 0:
            raise ValueError(f"Q{quarter} points must be positive")
        if not math.isfinite(low_lnlike) or not math.isfinite(high_lnlike):
            raise ValueError(f"Q{quarter} likelihood must be finite")
        if low_lnlike > 0 or high_lnlike > 0:
            raise ValueError(f"Q{quarter} has a positive term and would imply negative chi-square")
        seen_quarters.add(quarter)
        low_chi2 = -2.0 * low_lnlike
        high_chi2 = -2.0 * high_lnlike
        rows.append(
            {
                "quarter": quarter,
                "points": points,
                "low_chi2": low_chi2,
                "low_chi2_per_point": low_chi2 / points,
                "high_chi2": high_chi2,
                "high_chi2_per_point": high_chi2 / points,
                "high_minus_low_lnlike": high_lnlike - low_lnlike,
            }
        )
    if not rows:
        raise ValueError("quarter contrast has no photometric rows")

    for prefix in ("low", "high"):
        ordered = sorted(
            rows,
            key=lambda row: (-row[f"{prefix}_chi2_per_point"], row["quarter"]),
        )
        for rank, row in enumerate(ordered, start=1):
            row[f"{prefix}_worst_rank"] = rank

    rows.sort(key=lambda row: row["quarter"])
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "quarter",
        "points",
        "low_chi2",
        "low_chi2_per_point",
        "low_worst_rank",
        "high_chi2",
        "high_chi2_per_point",
        "high_worst_rank",
        "high_minus_low_lnlike",
    ]
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main():
    directory = (
        Path(__file__).resolve().parents[1]
        / "metadata"
        / "inference_assumptions_20260915"
        / "K02712_native_profile_20260919"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=directory / "quarter_contrast.csv")
    parser.add_argument("--output", type=Path, default=directory / "photometric_adequacy.csv")
    parser.add_argument("--provenance", type=Path)
    parser.add_argument(
        "--attest-native-photometric-terms",
        action="store_true",
        help="Attest that an arbitrary input contains decomposed -0.5*r^2 photometric terms.",
    )
    args = parser.parse_args()
    provenance = args.provenance
    if provenance is None and args.input.resolve() == (directory / "quarter_contrast.csv").resolve():
        provenance = directory / "provenance.json"
    try:
        rows = summarize(
            args.input,
            args.output,
            provenance_json=provenance,
            convention_attested=args.attest_native_photometric_terms,
        )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        parser.error(str(exc))
    print(f"wrote {len(rows)} quarters to {args.output}")


if __name__ == "__main__":
    main()
