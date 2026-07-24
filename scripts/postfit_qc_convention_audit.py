"""Compare ALDERAAN post-fit QC under weighted and raw nested rows.

Sagear's published method reports a small number of grazing and radius-precision
removals, while ALDERAAN FITS retain raw dynesty rows plus ``LN_WT``.  This
diagnostic makes that convention choice explicit and records whether it can
explain the published removal count.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from common import normalize_dynesty_weights


def result_path(target: str, source: str, archive_dir: Path, cloud_dir: Path) -> Path | None:
    if source == "original_alderaan_archive":
        path = archive_dir / f"{target}-results.fits"
        return path if path.exists() else None
    matches = sorted(cloud_dir.rglob(f"{target}-results.fits"))
    return matches[0] if matches else None


def weighted_fraction(mask: np.ndarray, weights: np.ndarray) -> float:
    total = float(np.sum(weights))
    return float(np.sum(weights[mask]) / total) if total > 0.0 else np.nan


def audit(summary: pd.DataFrame, archive_dir: Path, cloud_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (target, source), group in summary.groupby(["koi_target", "transit_fit_source"], sort=False):
        path = result_path(str(target), str(source), archive_dir, cloud_dir)
        if path is None:
            continue
        with fits.open(path, memmap=False) as hdul:
            data = hdul["SAMPLES"].data
            names = set(data.names)
            if "LN_WT" not in names:
                continue
            raw_weights = np.ones(len(data), dtype=float)
            dynesty_weights = normalize_dynesty_weights(np.asarray(data["LN_WT"], dtype=float))
            for _, planet in group.iterrows():
                index = int(planet["alderaan_planet_index"])
                ror_name = f"ROR_{index}"
                impact_name = f"IMPACT_{index}"
                if ror_name not in names or impact_name not in names:
                    continue
                ror = np.asarray(data[ror_name], dtype=float)
                impact = np.asarray(data[impact_name], dtype=float)
                valid = np.isfinite(ror) & np.isfinite(impact) & (ror > 0.0)
                if not valid.any():
                    continue
                ror = ror[valid]
                impact = impact[valid]
                weighted = dynesty_weights[valid]
                weighted /= weighted.sum()
                raw = raw_weights[valid]
                raw /= raw.sum()
                threshold = impact > (1.0 - ror)
                rows.append(
                    {
                        "kepoi_name": planet["kepoi_name"],
                        "koi_target": target,
                        "disk": planet.get("disk"),
                        "system": planet.get("system"),
                        "weighted_grazing_fraction": weighted_fraction(threshold, weighted),
                        "raw_grazing_fraction": weighted_fraction(threshold, raw),
                    }
                )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["weighted_grazing_exclude"] = result["weighted_grazing_fraction"] > 0.05
    result["raw_grazing_exclude"] = result["raw_grazing_fraction"] > 0.05
    result["convention_disagrees"] = result["weighted_grazing_exclude"] != result["raw_grazing_exclude"]
    return result


def summarize(audit_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (disk, system), group in audit_df.groupby(["disk", "system"], dropna=False):
        rows.append(
            {
                "disk": disk,
                "system": system,
                "n": len(group),
                "weighted_excluded": int(group["weighted_grazing_exclude"].sum()),
                "raw_excluded": int(group["raw_grazing_exclude"].sum()),
                "convention_disagreements": int(group["convention_disagrees"].sum()),
                "weighted_fraction_median": group["weighted_grazing_fraction"].median(),
                "raw_fraction_median": group["raw_grazing_fraction"].median(),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", default="outputs/eccentricity_posterior_summary_uniform_paired_full.csv")
    parser.add_argument("--archive-dir", default="inputs/alderaan_archive_raw_20260612/ALDERAAN_posteriors")
    parser.add_argument("--cloud-dir", default="alderaan_project/Results")
    parser.add_argument("--output", default="outputs/postfit_qc_convention_audit.csv")
    parser.add_argument("--summary-output", default="outputs/postfit_qc_convention_summary.csv")
    args = parser.parse_args()
    audit_df = audit(pd.read_csv(args.summary), Path(args.archive_dir), Path(args.cloud_dir))
    summary = summarize(audit_df)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    audit_df.to_csv(args.output, index=False)
    summary.to_csv(args.summary_output, index=False)
    print(f"audited_planets,{len(audit_df)}")
    print(f"weighted_grazing_excluded,{int(audit_df['weighted_grazing_exclude'].sum())}")
    print(f"raw_grazing_excluded,{int(audit_df['raw_grazing_exclude'].sum())}")
    print(f"convention_disagreements,{int(audit_df['convention_disagrees'].sum())}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
