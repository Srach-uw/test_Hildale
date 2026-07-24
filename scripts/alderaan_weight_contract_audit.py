"""Audit whether ALDERAAN FITS samples are raw nested points or equalized draws."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits


def audit_file(path: Path) -> dict[str, object]:
    with fits.open(path, memmap=False) as hdul:
        if "SAMPLES" not in hdul:
            return {"file": str(path), "status": "missing_samples"}
        table = hdul["SAMPLES"]
        names = list(table.data.names or [])
        if "LN_WT" not in names:
            return {"file": str(path), "status": "missing_ln_wt"}
        ln_wt = np.asarray(table.data["LN_WT"], dtype=float)
        finite = np.isfinite(ln_wt)
        if not finite.any():
            return {"file": str(path), "status": "invalid_ln_wt"}
        weights = np.exp(ln_wt[finite] - np.max(ln_wt[finite]))
        weights /= weights.sum()
        n = int(len(ln_wt))
        ess = float(1.0 / np.sum(weights**2))
        return {
            "file": str(path),
            "status": "ok",
            "rows": n,
            "finite_ln_wt": int(finite.sum()),
            "weight_ess": ess,
            "ess_fraction": ess / n,
            "max_normalized_weight": float(weights.max()),
            "has_raw_nested_metadata": bool("NITER" in table.header and "EFF" in table.header),
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", default=None)
    parser.add_argument("--output", default="outputs/alderaan_weight_contract_audit")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    results = Path(args.results) if args.results else root / "alderaan_project" / "Results"
    files = sorted(results.rglob("*-results.fits"))
    rows = [audit_file(path) for path in files]
    good = [row for row in rows if row.get("status") == "ok"]
    summary = {
        "result_files": len(files),
        "files_with_valid_ln_wt": len(good),
        "median_ess_fraction": float(np.median([row["ess_fraction"] for row in good])) if good else None,
        "p16_ess_fraction": float(np.quantile([row["ess_fraction"] for row in good], 0.16)) if good else None,
        "p84_ess_fraction": float(np.quantile([row["ess_fraction"] for row in good], 0.84)) if good else None,
        "interpretation": "SAMPLES rows are raw nested-sampling points; faithful downstream use requires LN_WT-weighted integration or weighted-to-unweighted resampling. Equal-row treatment is diagnostic only.",
    }
    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    import pandas as pd

    pd.DataFrame(rows).to_csv(output.with_suffix(".csv"), index=False)
    output.with_suffix(".md").write_text(
        "\n".join(
            [
                "# ALDERAAN Weight Contract Audit",
                "",
                "ALDERAAN writes raw nested-sampling rows to the `SAMPLES` HDU and stores `LN_WT` separately. The cited Gilbert implementation explicitly converts weighted transit samples to unweighted samples before downstream analysis. Therefore faithful use requires weighted integration or weighted-to-unweighted resampling; treating raw rows equally is diagnostic only.",
                "",
                f"- FITS files scanned: {summary['result_files']}",
                f"- Files with valid LN_WT: {summary['files_with_valid_ln_wt']}",
                f"- Median weighted ESS / raw rows: {summary['median_ess_fraction']:.3f}",
                f"- 16th-84th percentile ESS fraction: {summary['p16_ess_fraction']:.3f} to {summary['p84_ess_fraction']:.3f}",
                "",
                "## Decision",
                "Keep both branches explicitly named: `dynesty_weighted` is the canonical weighted-integration branch; `equal_raw_rows` is a diagnostic only. Do not present the latter as a corrected posterior analysis.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
