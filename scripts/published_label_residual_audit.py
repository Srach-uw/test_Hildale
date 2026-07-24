from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


SAGEAR_RAYLEIGH = {
    "thick_singles": (0.066, 0.045, 0.096),
    "thin_singles": (0.022, 0.017, 0.029),
    "thick_multis": (0.033, 0.015, 0.065),
    "thin_multis": (0.030, 0.023, 0.031),
}

FIT_VARIANTS = {
    "strict_reconstructed_labels": (
        "rayleigh_population_fit_transit_selection_manuscript_reciprocal_"
        "UNIFORM_EQUAL_NESTED_50K_STRICT_DIAGNOSTICS.csv"
    ),
    "published_overlap_reconstructed_labels": (
        "rayleigh_population_fit_transit_selection_manuscript_reciprocal_"
        "EQUAL_NESTED_PUBLISHED_OVERLAP_RECONSTRUCTED.csv"
    ),
    "published_labels": (
        "rayleigh_population_fit_transit_selection_manuscript_reciprocal_"
        "EQUAL_NESTED_PUBLISHED_LABELS.csv"
    ),
    "published_labels_archive_only": (
        "rayleigh_population_fit_transit_selection_manuscript_reciprocal_"
        "EQUAL_NESTED_PUBLISHED_LABELS_ARCHIVE.csv"
    ),
    "published_labels_cloud_only": (
        "rayleigh_population_fit_transit_selection_manuscript_reciprocal_"
        "EQUAL_NESTED_PUBLISHED_LABELS_CLOUD.csv"
    ),
    "published_labels_logg_ge_4": (
        "rayleigh_population_fit_transit_selection_manuscript_reciprocal_"
        "EQUAL_NESTED_PUBLISHED_LABELS_LOGG4.csv"
    ),
    "published_labels_berger2018_evol_0": (
        "rayleigh_population_fit_transit_selection_manuscript_reciprocal_"
        "EQUAL_NESTED_PUBLISHED_LABELS_EVOL0.csv"
    ),
    "published_labels_outlier_floor_1e6": (
        "rayleigh_population_fit_transit_selection_manuscript_reciprocal_"
        "EQUAL_NESTED_PUBLISHED_LABELS_EPS1E6.csv"
    ),
}

REQUIRED_SUMMARY_COLUMNS = {
    "kepoi_name",
    "kepid",
    "koi_target",
    "disk",
    "system",
    "disk_reconstructed",
    "disk_published",
    "disk_label_agrees",
    "posterior_source",
    "e50",
    "zeta_median",
    "berger_logg",
    "berger_rad",
    "berger2018_evol",
    "koi_model_snr",
    "rho_true_solar",
}


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = sorted(columns - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def population_key(disk: str, system: str) -> str:
    suffix = "singles" if system == "single" else "multis"
    return f"{disk}_{suffix}"


def load_fit(path: Path, variant: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required fit is missing for {variant}: {path}")
    fit = pd.read_csv(path)
    required = {
        "population",
        "n",
        "expected_e",
        "expected_e_lo",
        "expected_e_hi",
        "boundary_flag",
    }
    require_columns(fit, required, f"fit {variant}")
    fit = fit.copy()
    fit["variant"] = variant
    refs = fit["population"].map(SAGEAR_RAYLEIGH)
    if refs.isna().any():
        unknown = sorted(fit.loc[refs.isna(), "population"].astype(str).unique())
        raise ValueError(f"Unknown populations in {path}: {unknown}")
    fit["sagear_expected_e"] = [value[0] for value in refs]
    fit["sagear_lo"] = [value[1] for value in refs]
    fit["sagear_hi"] = [value[2] for value in refs]
    fit["delta_from_sagear"] = fit["expected_e"] - fit["sagear_expected_e"]
    fit["ratio_to_sagear"] = fit["expected_e"] / fit["sagear_expected_e"]
    fit["interval_overlaps_sagear"] = (
        (fit["expected_e_hi"] >= fit["sagear_lo"])
        & (fit["expected_e_lo"] <= fit["sagear_hi"])
    )
    return fit


def numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    return pd.to_numeric(frame[column], errors="coerce")


def fraction_true(series: pd.Series) -> float:
    clean = series.dropna()
    if clean.empty:
        return np.nan
    if clean.dtype == bool:
        return float(clean.mean())
    normalized = clean.astype(str).str.strip().str.lower()
    valid = normalized.isin({"true", "false", "1", "0"})
    if not valid.any():
        return np.nan
    return float(normalized[valid].isin({"true", "1"}).mean())


def build_migration(summary: pd.DataFrame) -> pd.DataFrame:
    migration = (
        summary.groupby(
            ["disk_reconstructed", "disk_published", "system"],
            dropna=False,
        )
        .agg(planets=("kepoi_name", "size"), hosts=("kepid", "nunique"))
        .reset_index()
        .sort_values(["system", "disk_reconstructed", "disk_published"])
    )
    return migration


def enrich_qc(summary: pd.DataFrame, qc: pd.DataFrame) -> pd.DataFrame:
    require_columns(
        qc,
        {
            "kepoi_name",
            "nested_grazing_fraction",
            "planet_radius_fractional_uncertainty_approx",
            "rho_circular_to_catalog_ratio",
            "gilbert_grazing_exclude",
            "gilbert_radius_precision_exclude",
        },
        "Gilbert QC audit",
    )
    if qc["kepoi_name"].duplicated().any():
        duplicated = qc.loc[qc["kepoi_name"].duplicated(), "kepoi_name"].tolist()
        raise ValueError(f"Gilbert QC audit has duplicate planets: {duplicated[:10]}")
    return summary.merge(qc, on="kepoi_name", how="left", validate="one_to_one")


def build_composition(enriched: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_columns = ["disk", "system", "posterior_source"]
    for keys, group in enriched.groupby(group_columns, dropna=False):
        disk, system, source = keys
        e50 = numeric(group, "e50")
        zeta = numeric(group, "zeta_median")
        logg = numeric(group, "berger_logg")
        radius = numeric(group, "berger_rad")
        snr = numeric(group, "koi_model_snr")
        rho_ratio = numeric(group, "rho_circular_to_catalog_ratio")
        grazing = numeric(group, "nested_grazing_fraction")
        radius_uncertainty = numeric(
            group, "planet_radius_fractional_uncertainty_approx"
        )
        evol = numeric(group, "berger2018_evol")
        rows.append(
            {
                "population": population_key(str(disk), str(system)),
                "disk": disk,
                "system": system,
                "posterior_source": source,
                "planets": len(group),
                "hosts": group["kepid"].nunique(),
                "e50_median": e50.median(),
                "e50_p90": e50.quantile(0.90),
                "zeta_median": zeta.median(),
                "berger_logg_median": logg.median(),
                "berger_logg_p10": logg.quantile(0.10),
                "berger_radius_median": radius.median(),
                "berger_radius_p90": radius.quantile(0.90),
                "koi_model_snr_median": snr.median(),
                "rho_circular_to_catalog_median": rho_ratio.median(),
                "rho_circular_to_catalog_p90": rho_ratio.quantile(0.90),
                "nested_grazing_fraction_median": grazing.median(),
                "radius_fractional_uncertainty_median": radius_uncertainty.median(),
                "berger2018_evolved_fraction": float((evol > 0).mean()),
                "gilbert_grazing_exclude_fraction": fraction_true(
                    group["gilbert_grazing_exclude"]
                ),
                "gilbert_radius_precision_exclude_fraction": fraction_true(
                    group["gilbert_radius_precision_exclude"]
                ),
                "qc_match_fraction": float(group["raw_result_file"].notna().mean()),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["population", "posterior_source"]
    )


def build_thin_single_leverage(
    enriched: pd.DataFrame, influence: pd.DataFrame
) -> pd.DataFrame:
    require_columns(
        influence,
        {
            "kepoi_name",
            "population",
            "expected_e_full",
            "expected_e_leave_one_out",
            "fractional_shift_when_removed",
        },
        "leave-one influence",
    )
    keep = [
        "kepoi_name",
        "kepid",
        "koi_target",
        "posterior_source",
        "e16",
        "e84",
        "berger_logg",
        "berger_rad",
        "berger2018_evol",
        "koi_model_snr",
        "koi_prad",
        "rho_true_solar",
        "nested_grazing_fraction",
        "planet_radius_fractional_uncertainty_approx",
        "rho_circular_to_catalog_ratio",
        "gilbert_grazing_exclude",
        "gilbert_radius_precision_exclude",
    ]
    detail = enriched.loc[
        (enriched["disk"] == "thin") & (enriched["system"] == "single"),
        keep,
    ]
    joined = influence.loc[
        influence["population"] == "thin_singles"
    ].merge(detail, on=["kepoi_name", "kepid", "koi_target"], how="left", validate="one_to_one")
    joined["absolute_fractional_shift"] = joined[
        "fractional_shift_when_removed"
    ].abs()
    return joined.sort_values(
        ["absolute_fractional_shift", "kepoi_name"], ascending=[False, True]
    )


def format_interval(row: pd.Series) -> str:
    return (
        f"{row['expected_e']:.4f} "
        f"[{row['expected_e_lo']:.4f}, {row['expected_e_hi']:.4f}]"
    )


def make_markdown(
    fit_comparison: pd.DataFrame,
    migration: pd.DataFrame,
    composition: pd.DataFrame,
    leverage: pd.DataFrame,
) -> str:
    canonical = fit_comparison.loc[
        fit_comparison["variant"] == "published_labels"
    ].set_index("population")
    thin = canonical.loc["thin_singles"]
    archive = fit_comparison.loc[
        (fit_comparison["variant"] == "published_labels_archive_only")
        & (fit_comparison["population"] == "thin_singles")
    ].iloc[0]
    cloud = fit_comparison.loc[
        (fit_comparison["variant"] == "published_labels_cloud_only")
        & (fit_comparison["population"] == "thin_singles")
    ].iloc[0]
    reconstructed = fit_comparison.loc[
        (fit_comparison["variant"] == "published_overlap_reconstructed_labels")
        & (fit_comparison["population"] == "thin_singles")
    ].iloc[0]

    lines = [
        "# Published-label residual audit",
        "",
        "This audit isolates the remaining Sagear Table 3 discrepancy after using",
        "the exact published host classifications and the equal-raw-nested-row",
        "diagnostic that reproduces three of four Rayleigh populations.",
        "",
        "## Decisive result",
        "",
        f"- Published-label thin singles: {format_interval(thin)}.",
        "- Sagear thin singles: 0.0220 [0.0170, 0.0290].",
        f"- Published-overlap planets with reconstructed labels: {format_interval(reconstructed)}.",
        f"- Archive-only published-label thin singles: {format_interval(archive)}.",
        f"- Replacement-cloud-only published-label thin singles: {format_interval(cloud)}.",
        "",
        "The published disk labels do not solve the residual. They increase the",
        "thin-single estimate, while the replacement-cloud subset is much hotter",
        "than both the archive subset and Sagear. This localizes the unresolved",
        "difference to posterior construction and/or unpublished sample curation,",
        "not to a global Rayleigh formula or a simple label swap.",
        "",
        "## Fit variants",
        "",
        "| variant | population | N | ours (16th-84th) | Sagear | ratio | overlap |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    preferred_order = list(FIT_VARIANTS)
    ordered = fit_comparison.assign(
        _variant_order=pd.Categorical(
            fit_comparison["variant"], preferred_order, ordered=True
        )
    ).sort_values(["_variant_order", "population"])
    for _, row in ordered.iterrows():
        lines.append(
            f"| {row['variant']} | {row['population']} | {int(row['n'])} | "
            f"{format_interval(row)} | {row['sagear_expected_e']:.3f} | "
            f"{row['ratio_to_sagear']:.2f} | "
            f"{'yes' if row['interval_overlaps_sagear'] else 'no'} |"
        )

    lines.extend(
        [
            "",
            "## Disk-label migration",
            "",
            "| system | reconstructed | published | planets | hosts |",
            "|---|---|---|---:|---:|",
        ]
    )
    for _, row in migration.iterrows():
        lines.append(
            f"| {row['system']} | {row['disk_reconstructed']} | "
            f"{row['disk_published']} | {int(row['planets'])} | "
            f"{int(row['hosts'])} |"
        )

    lines.extend(
        [
            "",
            "## Source and QC composition",
            "",
            "| population | source | planets | e50 median | logg median | "
            "evolved frac | S/N median | grazing-exclude frac | radius-QC frac |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for _, row in composition.iterrows():
        lines.append(
            f"| {row['population']} | {row['posterior_source']} | "
            f"{int(row['planets'])} | {row['e50_median']:.3f} | "
            f"{row['berger_logg_median']:.3f} | "
            f"{row['berger2018_evolved_fraction']:.1%} | "
            f"{row['koi_model_snr_median']:.1f} | "
            f"{row['gilbert_grazing_exclude_fraction']:.1%} | "
            f"{row['gilbert_radius_precision_exclude_fraction']:.1%} |"
        )

    lines.extend(
        [
            "",
            "## Highest-influence thin singles",
            "",
            "| KOI | source | e50 | logg | S/N | rho_circ/rho_catalog | "
            "leave-one fractional shift |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    for _, row in leverage.head(20).iterrows():
        lines.append(
            f"| {row['kepoi_name']} | {row['posterior_source']} | "
            f"{row['e50']:.3f} | {row['berger_logg']:.3f} | "
            f"{row['koi_model_snr']:.1f} | "
            f"{row['rho_circular_to_catalog_ratio']:.2f} | "
            f"{row['fractional_shift_when_removed']:.2%} |"
        )

    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "The archive/cloud split is diagnostic, not a permissible scientific cut.",
            "Likewise, logg and Berger evolutionary cuts are sensitivities because",
            "the published methods do not state them for the planet-host sample.",
            "The remaining exact-replication inputs are the final planet inclusion",
            "table, visual rejection/refit list, exact stellar-density fields,",
            "post-model code, and confirmation of whether dynesty LN_WT was applied.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Consolidate the residual Sagear Table 3 diagnosis."
    )
    parser.add_argument("--outputs", default="outputs")
    args = parser.parse_args()
    output = Path(args.outputs).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)

    summary_path = output / "eccentricity_posterior_summary_equal_nested_published_labels.csv"
    qc_path = output / "gilbert_postfit_qc_audit_uniform_paired_full.csv"
    influence_path = output / "rayleigh_leave_one_planet_influence_EQUAL_NESTED_PUBLISHED_LABELS.csv"
    for path in (summary_path, qc_path, influence_path):
        if not path.exists():
            raise FileNotFoundError(path)

    summary = pd.read_csv(summary_path)
    require_columns(summary, REQUIRED_SUMMARY_COLUMNS, "published-label summary")
    if summary["kepoi_name"].duplicated().any():
        raise ValueError("Published-label summary contains duplicate planets")
    if not summary["published_label_available"].fillna(False).astype(bool).all():
        raise ValueError("Published-label summary contains rows without published labels")

    qc = pd.read_csv(qc_path)
    enriched = enrich_qc(summary, qc)
    influence = pd.read_csv(influence_path)

    fit_frames = [
        load_fit(output / filename, variant)
        for variant, filename in FIT_VARIANTS.items()
    ]
    fit_comparison = pd.concat(fit_frames, ignore_index=True)
    migration = build_migration(summary)
    composition = build_composition(enriched)
    leverage = build_thin_single_leverage(enriched, influence)

    fit_comparison.to_csv(
        output / "published_label_residual_fit_comparison.csv", index=False
    )
    migration.to_csv(
        output / "published_label_disk_migration.csv", index=False
    )
    composition.to_csv(
        output / "published_label_source_qc_composition.csv", index=False
    )
    leverage.to_csv(
        output / "published_label_thin_single_leverage_audit.csv", index=False
    )
    report = make_markdown(fit_comparison, migration, composition, leverage)
    (output / "published_label_residual_audit.md").write_text(
        report, encoding="utf-8"
    )

    print("Published-label residual audit complete")
    print(f"planets={len(summary)} hosts={summary['kepid'].nunique()}")
    print(
        fit_comparison.loc[
            (fit_comparison["variant"] == "published_labels")
            & (fit_comparison["population"] == "thin_singles"),
            [
                "population",
                "n",
                "expected_e",
                "expected_e_lo",
                "expected_e_hi",
                "sagear_expected_e",
            ],
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
