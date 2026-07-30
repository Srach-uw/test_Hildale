from __future__ import annotations

"""Assemble one method-consistent posterior archive from ALDERAAN runs."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_SUMMARY_COLUMNS = {
    "kepoi_name",
    "posterior_file",
    "posterior_source",
    "impact_mode",
    "formalism",
    "include_transit_prior",
    "density_error_mode",
    "e_max",
    "qc_primary_exclude",
    "qc_reasons",
}
REQUIRED_NPZ_KEYS = {
    "e_grid",
    "omega_grid",
    "posterior",
    "include_transit_prior",
    "formalism",
    "impact_mode",
    "e_max",
}
SAMPLE_PROVENANCE_COLUMNS = [
    "kepid",
    "koi_target",
    "koi_period",
    "disk",
    "system",
    "P_thick",
    "berger_logg",
    "berger_rad",
    "berger2018_evol",
    "berger2018_bin",
    "koi_disposition",
    "koi_model_snr",
    "koi_prad",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Combine direct-importance posterior summaries while enforcing one "
            "paired-impact posterior contract."
        )
    )
    parser.add_argument("--archive", required=True, help="Summary re-extracted from the original ALDERAAN FITS archive.")
    parser.add_argument("--new", required=True, help="Summary from the missing-target ALDERAAN run.")
    parser.add_argument("--sample", required=True, help="Canonical planet sample that supplies final disk/system labels.")
    parser.add_argument("--archive-excluded", default=None)
    parser.add_argument("--new-excluded", default=None)
    parser.add_argument("--out", required=True)
    parser.add_argument("--qc-out", required=True)
    parser.add_argument("--coverage-out", required=True)
    parser.add_argument("--manifest-out", required=True)
    parser.add_argument("--overlap-out", required=True)
    parser.add_argument(
        "--skip-npz-validation",
        action="store_true",
        help="Diagnostic only. Canonical assembly validates every selected NPZ.",
    )
    return parser.parse_args()


def scalar(value: np.ndarray) -> object:
    return np.asarray(value).item()


def validate_summary(summary: pd.DataFrame, label: str) -> None:
    missing = sorted(REQUIRED_SUMMARY_COLUMNS - set(summary.columns))
    if missing:
        raise ValueError(f"{label} summary is missing required columns: {missing}")
    duplicated = summary["kepoi_name"].astype(str).duplicated(keep=False)
    if duplicated.any():
        names = sorted(summary.loc[duplicated, "kepoi_name"].astype(str).unique())
        raise ValueError(f"{label} summary contains duplicate planet IDs: {names[:10]}")
    impact_modes = sorted(summary["impact_mode"].dropna().astype(str).str.lower().unique())
    if impact_modes != ["alderaan"]:
        raise ValueError(f"{label} summary is not uniformly paired-impact: {impact_modes}")
    transit_prior = summary["include_transit_prior"].astype("boolean")
    if transit_prior.isna().any() or transit_prior.any():
        raise ValueError(f"{label} summary must have include_transit_prior=False for every row")
    method_columns = [
        "posterior_source",
        "formalism",
        "density_source",
        "density_sampling_mode",
        "density_error_mode",
        "posterior_sampling_mode",
        "period_sampling_mode",
        "nested_weight_mode",
        "e_max",
    ]
    missing_method_columns = sorted(set(method_columns) - set(summary.columns))
    if missing_method_columns:
        raise ValueError(f"{label} summary is missing method provenance: {missing_method_columns}")
    for column in method_columns:
        values = summary[column].dropna().astype(str).str.strip()
        if len(values) != len(summary) or values.eq("").any():
            raise ValueError(f"{label} summary has missing {column} provenance")
        values = values.unique()
        if len(values) != 1:
            raise ValueError(f"{label} summary mixes {column}: {sorted(values)}")
    missing_files = [p for p in summary["posterior_file"].astype(str) if not Path(p).is_file()]
    if missing_files:
        raise FileNotFoundError(f"{label} summary references {len(missing_files)} missing posterior files")


def validate_npz_contract(summary: pd.DataFrame) -> dict[str, object]:
    reference: dict[str, object] | None = None
    for path_text in summary["posterior_file"].astype(str):
        path = Path(path_text)
        with np.load(path, allow_pickle=False) as data:
            missing = sorted(REQUIRED_NPZ_KEYS - set(data.files))
            if missing:
                raise ValueError(f"Posterior product is missing {missing}: {path}")
            e_grid = np.asarray(data["e_grid"], dtype=float)
            omega_grid = np.asarray(data["omega_grid"], dtype=float)
            posterior = np.asarray(data["posterior"], dtype=float)
            current = {
                "e_grid": e_grid,
                "omega_grid": omega_grid,
                "include_transit_prior": bool(scalar(data["include_transit_prior"])),
                "formalism": str(scalar(data["formalism"])),
                "impact_mode": str(scalar(data["impact_mode"])),
                "e_max": float(scalar(data["e_max"])),
            }
        if posterior.shape != (len(e_grid), len(omega_grid)):
            raise ValueError(f"Posterior/grid shape mismatch: {path}")
        if not np.all(np.isfinite(posterior)) or np.any(posterior < 0.0) or posterior.sum() <= 0.0:
            raise ValueError(f"Posterior mass is invalid: {path}")
        if current["include_transit_prior"] or current["impact_mode"].lower() != "alderaan":
            raise ValueError(f"Posterior is not paired-impact/no-transit-prior: {path}")
        if reference is None:
            reference = current
            continue
        for key in ["formalism", "impact_mode", "e_max", "include_transit_prior"]:
            if current[key] != reference[key]:
                raise ValueError(f"Posterior metadata mismatch for {key}: {path}")
        if not np.array_equal(current["e_grid"], reference["e_grid"]):
            raise ValueError(f"Eccentricity-grid mismatch: {path}")
        if not np.array_equal(current["omega_grid"], reference["omega_grid"]):
            raise ValueError(f"Omega-grid mismatch: {path}")
    if reference is None:
        raise ValueError("No posterior products were selected")
    return {
        "formalism": reference["formalism"],
        "impact_mode": reference["impact_mode"],
        "e_max": reference["e_max"],
        "e_grid_size": len(reference["e_grid"]),
        "omega_grid_size": len(reference["omega_grid"]),
    }


def coverage_table(sample: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    available = set(summary["kepoi_name"].astype(str))
    rows = []
    for (disk, system), group in sample.groupby(["disk", "system"], dropna=False):
        covered = int(group["kepoi_name"].astype(str).isin(available).sum())
        rows.append(
            {
                "disk": disk,
                "system": system,
                "sample_planets": len(group),
                "posterior_planets": covered,
                "missing_planets": len(group) - covered,
                "coverage_fraction": covered / len(group),
            }
        )
    return pd.DataFrame(rows).sort_values(["disk", "system"]).reset_index(drop=True)


def exclusion_reason_map(paths: list[tuple[str, str | None]]) -> dict[str, str]:
    collected: dict[str, list[str]] = {}
    for source, path_text in paths:
        if not path_text or not Path(path_text).is_file():
            continue
        table = pd.read_csv(path_text)
        if not {"kepoi_name", "reason"}.issubset(table.columns):
            continue
        for _, row in table.iterrows():
            name = str(row["kepoi_name"])
            reason = f"{source}:{row['reason']}"
            collected.setdefault(name, [])
            if reason not in collected[name]:
                collected[name].append(reason)
    return {name: ";".join(reasons) for name, reasons in collected.items()}


def assemble(
    archive: pd.DataFrame,
    new: pd.DataFrame,
    sample: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    validate_summary(archive, "archive")
    validate_summary(new, "new")
    if sample["kepoi_name"].astype(str).duplicated().any():
        raise ValueError("Canonical sample contains duplicate kepoi_name values")

    archive = archive.copy()
    new = new.copy()
    archive["transit_fit_source"] = "original_alderaan_archive"
    new["transit_fit_source"] = "cloud_missing_run"
    archive["_source_rank"] = 1
    new["_source_rank"] = 0

    overlap_ids = sorted(set(archive["kepoi_name"].astype(str)) & set(new["kepoi_name"].astype(str)))
    overlap = archive[archive["kepoi_name"].astype(str).isin(overlap_ids)].merge(
        new[new["kepoi_name"].astype(str).isin(overlap_ids)],
        on="kepoi_name",
        suffixes=("_archive", "_new"),
    )

    merged = pd.concat([archive, new], ignore_index=True, sort=False)
    for column in [
        "formalism",
        "density_source",
        "density_sampling_mode",
        "density_error_mode",
        "posterior_sampling_mode",
        "period_sampling_mode",
        "nested_weight_mode",
        "e_max",
        "include_transit_prior",
    ]:
        values = merged[column].dropna().astype(str).str.strip()
        if len(values) != len(merged) or values.eq("").any() or values.str.lower().nunique() != 1:
            raise ValueError(
                "Archive and recovery summaries do not share one uniform "
                f"{column} contract: {sorted(values.unique())}"
            )
    merged = (
        merged.sort_values(["kepoi_name", "_source_rank"])
        .drop_duplicates("kepoi_name", keep="last")
        .drop(columns="_source_rank")
    )
    inherited = [column for column in SAMPLE_PROVENANCE_COLUMNS if column in merged.columns]
    merged = merged.drop(columns=inherited)
    sample_columns = ["kepoi_name", *[column for column in SAMPLE_PROVENANCE_COLUMNS if column in sample.columns]]
    labels = sample[sample_columns].copy()
    labels["kepoi_name"] = labels["kepoi_name"].astype(str)
    merged["kepoi_name"] = merged["kepoi_name"].astype(str)
    merged = labels.merge(
        merged.drop(columns=[c for c in ["kepid", "koi_target", "koi_period"] if c in merged.columns]),
        on="kepoi_name",
        how="inner",
        validate="one_to_one",
    )
    merged["posterior_source"] = "uniform_paired_direct_importance"
    merged["qc_primary_exclude"] = merged["qc_primary_exclude"].astype("boolean").fillna(True).astype(bool)
    merged["qc_reasons"] = merged["qc_reasons"].fillna("missing_qc_reason").astype(str)
    merged = merged.sort_values(["disk", "system", "koi_target", "koi_period"]).reset_index(drop=True)
    return merged, overlap


def main() -> None:
    args = parse_args()
    archive = pd.read_csv(args.archive)
    new = pd.read_csv(args.new)
    sample = pd.read_csv(args.sample)
    merged, overlap = assemble(archive, new, sample)

    contract = None if args.skip_npz_validation else validate_npz_contract(merged)
    qc = merged[~merged["qc_primary_exclude"]].reset_index(drop=True)
    coverage = coverage_table(sample, merged)

    manifest = sample.merge(
        merged[
            [
                "kepoi_name",
                "posterior_file",
                "posterior_source",
                "transit_fit_source",
                "impact_mode",
                "qc_primary_exclude",
                "qc_reasons",
            ]
        ],
        on="kepoi_name",
        how="left",
        validate="one_to_one",
    )
    manifest["posterior_status"] = np.where(
        manifest["posterior_file"].notna(), "posterior_available", "missing_after_uniform_assembly"
    )
    reasons = exclusion_reason_map(
        [
            ("original_archive", args.archive_excluded),
            ("cloud_missing", args.new_excluded),
        ]
    )
    missing = manifest["posterior_file"].isna()
    manifest.loc[missing, "qc_primary_exclude"] = True
    manifest.loc[missing, "qc_reasons"] = (
        manifest.loc[missing, "kepoi_name"].astype(str).map(reasons).fillna("missing_after_uniform_assembly")
    )

    outputs = {
        Path(args.out): merged,
        Path(args.qc_out): qc,
        Path(args.coverage_out): coverage,
        Path(args.manifest_out): manifest,
        Path(args.overlap_out): overlap,
    }
    for path, table in outputs.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(path, index=False)

    print(f"Archive rows: {len(archive)}")
    print(f"New rows: {len(new)}")
    print(f"Overlap planets (archive retained): {len(overlap)}")
    print(f"Uniform unique rows: {len(merged)}")
    print(f"Primary-QC rows: {len(qc)}")
    if contract is not None:
        print("Validated NPZ contract: " + ", ".join(f"{key}={value}" for key, value in contract.items()))
    print(coverage.to_string(index=False))


if __name__ == "__main__":
    main()
