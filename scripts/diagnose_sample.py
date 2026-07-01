from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

from common import (
    add_audit_row,
    add_target_and_system,
    finite_columns,
    load_config,
    output_dir,
    read_angus,
    read_berger_table2,
    read_gaia_kepler,
    read_koi,
    root_path,
    write_json,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Sagear replication sample attrition audit.")
    parser.add_argument("--config", default=None)
    parser.add_argument(
        "--allow-fallback-gmm",
        action="store_true",
        help="Use a diagnostic all-Angus GMM if APOGEE chemistry is missing. This is not Sagear-equivalent.",
    )
    parser.add_argument(
        "--force-fallback-gmm",
        action="store_true",
        help="Force the diagnostic all-Angus GMM even if APOGEE chemistry exists.",
    )
    parser.add_argument(
        "--fallback-velocity",
        choices=["direct", "old_astropy"],
        default="direct",
        help=(
            "Velocity convention for the diagnostic fallback GMM. 'direct' uses Angus vy and sqrt(vx^2+vz^2); "
            "'old_astropy' reproduces the old Hilldale cylindrical projection used by its Toomre figure."
        ),
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    out_dir = output_dir()
    mode = "strict"
    missing: dict[str, str] = {}
    audit: list[dict[str, object]] = []

    koi = read_koi(cfg)
    add_audit_row(audit, "00_raw_koi_catalog", koi, "NASA Exoplanet Archive KOI table")

    df = koi[koi["koi_disposition"].isin(["CONFIRMED", "CANDIDATE"])].copy()
    add_audit_row(audit, "01_remove_false_positives", df, "koi_disposition in CONFIRMED,CANDIDATE")

    cuts = cfg["cuts"]
    df = df[
        (df["koi_period"] >= cuts["period_min_days"])
        & (df["koi_period"] <= cuts["period_max_days"])
    ].copy()
    add_audit_row(audit, "02_period_1_to_100_days", df)

    gk = read_gaia_kepler(cfg)
    keep_gaia_cols = [
        c
        for c in [
            "kepid",
            "source_id",
            "ra",
            "dec",
            "parallax",
            "pmra",
            "pmdec",
            "ruwe",
            "radial_velocity",
            "phot_g_mean_mag",
            "phot_bp_mean_mag",
            "phot_rp_mean_mag",
            "bp_rp",
            "feh",
            "teff",
            "logg",
        ]
        if c in gk.columns
    ]
    gk_keep = gk[keep_gaia_cols].drop_duplicates("kepid").rename(
        columns={"ra": "gaia_ra", "dec": "gaia_dec"}
    )
    df = df.merge(gk_keep, on="kepid", how="inner")
    add_audit_row(audit, "03_join_gaia_kepler", df)

    furlan_path = root_path(cfg, "furlan_contamination")
    if furlan_path and furlan_path.exists():
        contamination = read_contamination_table(furlan_path)
        df = apply_contamination_cut(df, contamination, cuts["contamination_max"])
        add_audit_row(audit, "04_furlan_contamination_le_5pct", df, str(furlan_path))
    else:
        missing["furlan_contamination"] = (
            "Missing Furlan et al. contamination table; contamination <= 5% cut was not applied."
        )
        add_audit_row(audit, "04_furlan_contamination_le_5pct_SKIPPED", df, "missing input")

    if "ruwe" in df.columns:
        df = df[(df["ruwe"].notna()) & (df["ruwe"] <= cuts["ruwe_max"])].copy()
    add_audit_row(audit, "05_ruwe_le_1p4", df)

    if cuts.get("parallax_min_mas") is not None:
        df = df[(df["parallax"].notna()) & (df["parallax"] > cuts["parallax_min_mas"])].copy()
        add_audit_row(audit, "06_parallax_cut", df)
    else:
        add_audit_row(audit, "06_parallax_cut_SKIPPED", df, "not in Sagear plan defaults")

    berger = read_berger_table2(cfg)
    df = df.merge(berger, on="kepid", how="inner")
    add_audit_row(audit, "07_join_berger_density", df)

    if cuts.get("stellar_teff_max") is not None:
        df = df[pd.to_numeric(df["berger_teff"], errors="coerce") < cuts["stellar_teff_max"]].copy()
        add_audit_row(audit, "07a_fgkm_teff_lt_6500", df, "FGKM upper-temperature cut from Sagear sample description")

    b18_path = root_path(cfg, "berger2018_stellar")
    if cuts.get("berger2018_bin_required") is not None:
        if b18_path and b18_path.exists():
            b18 = read_berger2018_stellar_table(b18_path)
            df = df.merge(b18[["kepid", "berger2018_bin"]], on="kepid", how="inner")
            df = df[
                pd.to_numeric(df["berger2018_bin"], errors="coerce")
                == cuts["berger2018_bin_required"]
            ].copy()
            add_audit_row(
                audit,
                "07b_berger2018_bin0",
                df,
                "Berger+2018 Bin=0 resolved-companion flag; likely hidden binary/sample cut",
            )
        else:
            missing["berger2018_stellar"] = (
                "Missing Berger+2018 stellar table; Bin=0 resolved-companion cut was not applied. "
                f"Expected path: {b18_path}"
            )
            add_audit_row(audit, "07b_berger2018_bin0_SKIPPED", df, "missing input")

    angus = read_angus(cfg)
    df = df.merge(angus, on="kepid", how="inner", suffixes=("", "_angus"))
    add_audit_row(audit, "08_join_angus_velocities", df)

    df = add_velocity_features(df)
    apogee_path = root_path(cfg, "apogee_dr17_chemical")
    if apogee_path and apogee_path.exists() and not args.force_fallback_gmm:
        df, note = classify_with_chemical_gmm(df, pd.read_csv(apogee_path), cfg)
        add_audit_row(audit, "09_sagear_apogee_calibrated_gmm", df, note)
    elif not (args.allow_fallback_gmm or args.force_fallback_gmm):
        expected = f" Expected path: {apogee_path}" if apogee_path else ""
        missing["apogee_dr17_chemical"] = (
            "Missing APOGEE DR17 chemical table; disk labels were not assigned. "
            "This is strict Sagear mode. Rerun with --allow-fallback-gmm only for diagnostics."
            + expected
        )
        df["P_thick"] = np.nan
        df["disk"] = "unclassified"
        add_audit_row(audit, "09_disk_classification_SKIPPED", df, "missing APOGEE calibration input")
    else:
        mode = "diagnostic"
        expected = f" Expected path: {apogee_path}" if apogee_path else ""
        if args.force_fallback_gmm and apogee_path and apogee_path.exists():
            missing["apogee_dr17_chemical"] = (
                "APOGEE DR17 chemical table exists, but --force-fallback-gmm requested the "
                "diagnostic all-Angus GMM. This is not Sagear-equivalent."
            )
        else:
            missing["apogee_dr17_chemical"] = (
                "Missing APOGEE DR17 chemical table; used fallback GMM on all Angus+Kepler velocities. "
                "This is diagnostic only and not Sagear-equivalent."
                + expected
            )
        if args.fallback_velocity != "direct":
            mode = f"{mode}_{args.fallback_velocity}"
        df, note = classify_with_fallback_gmm(df, cfg, args.fallback_velocity)
        add_audit_row(audit, "09_fallback_angus_allstar_gmm", df, note)

    df = add_target_and_system(df)
    add_audit_row(audit, "10_assign_single_multi_after_cuts", df)

    conv_path = root_path(cfg, "alderaan_convergence")
    if conv_path and conv_path.exists():
        conv = pd.read_csv(conv_path)
        df = apply_convergence_cut(df, conv)
        add_audit_row(audit, "11_alderaan_convergence_cut", df, str(conv_path))
    else:
        missing["alderaan_convergence"] = "Missing ALDERAAN convergence table; no convergence failures removed."
        add_audit_row(audit, "11_alderaan_convergence_cut_SKIPPED", df, "missing input")

    audit_df = pd.DataFrame(audit)
    counts = disk_counts(df)
    comparison = compare_to_sagear(counts, cfg)

    status = {
        "mode": mode,
        "sagear_equivalent": mode == "strict" and not missing,
        "missing_inputs": missing,
        "outputs": {
            "sample_audit": f"sample_audit_{mode}.csv",
            "disk_counts": f"disk_counts_{mode}.csv",
            "sagear_count_comparison": f"sagear_count_comparison_{mode}.csv",
            "canonical_sample": f"canonical_sample_{mode}.csv",
            "missing_inputs": f"missing_inputs_{mode}.json",
        },
    }

    audit_df.to_csv(out_dir / f"sample_audit_{mode}.csv", index=False)
    counts.to_csv(out_dir / f"disk_counts_{mode}.csv", index=False)
    comparison.to_csv(out_dir / f"sagear_count_comparison_{mode}.csv", index=False)
    df.to_csv(out_dir / f"canonical_sample_{mode}.csv", index=False)
    write_json(out_dir / f"missing_inputs_{mode}.json", missing)
    write_json(out_dir / f"pipeline_status_{mode}.json", status)

    print("\n=== Sample audit ===")
    print(audit_df.to_string(index=False))
    print("\n=== Disk/system counts ===")
    print(counts.to_string(index=False))
    print("\n=== Sagear count comparison ===")
    print(comparison.to_string(index=False))
    if missing:
        print("\n=== Missing inputs / non-Sagear fallbacks ===")
        for key, value in missing.items():
            print(f"- {key}: {value}")
    print(f"\nMode: {mode}")
    print(f"Sagear-equivalent run: {status['sagear_equivalent']}")
    print(f"Wrote outputs to: {out_dir}")


def read_contamination_table(path) -> pd.DataFrame:
    path = pd.io.common.stringify_path(path)
    if path.lower().endswith(".dat"):
        # Furlan+2017 VizieR J/AJ/153/71 table9.dat. The weighted
        # average primary-host planet-radius correction factor is bytes
        # 156-162; PRCF=sqrt(F_total/F_primary).
        tab = pd.read_fwf(
            path,
            colspecs=[(0, 4), (155, 162), (163, 170)],
            names=["koi_num", "furlan_prcf_avg", "furlan_prcf_avg_err"],
        )
        tab["koi_num"] = pd.to_numeric(tab["koi_num"], errors="coerce")
        tab["furlan_prcf_avg"] = pd.to_numeric(tab["furlan_prcf_avg"], errors="coerce")
        tab.loc[tab["furlan_prcf_avg"] <= 0, "furlan_prcf_avg"] = np.nan
        tab["furlan_contam"] = 1.0 - 1.0 / np.square(tab["furlan_prcf_avg"])
        return tab[["koi_num", "furlan_prcf_avg", "furlan_prcf_avg_err", "furlan_contam"]]
    return pd.read_csv(path)


def apply_contamination_cut(df: pd.DataFrame, contamination: pd.DataFrame, max_frac: float) -> pd.DataFrame:
    if {"koi_num", "furlan_contam"}.issubset(contamination.columns):
        out = df.copy()
        out["koi_num"] = (
            out["kepoi_name"].astype(str).str.extract(r"K0*(\d+)\.", expand=False).astype(float)
        )
        ctab = contamination[["koi_num", "furlan_contam", "furlan_prcf_avg"]].drop_duplicates("koi_num")
        out = out.merge(ctab, on="koi_num", how="left")
        return out[(out["furlan_contam"].isna()) | (out["furlan_contam"] <= max_frac)].copy()

    kep_col = next((c for c in contamination.columns if c.lower() in {"kepid", "kic", "kic_id"}), None)
    cont_col = next((c for c in contamination.columns if "contam" in c.lower() or "flux" in c.lower()), None)
    if kep_col is None or cont_col is None:
        raise ValueError("Furlan contamination table needs a KIC/kepid column and a contamination/flux column.")
    ctab = contamination[[kep_col, cont_col]].rename(columns={kep_col: "kepid", cont_col: "furlan_contam"})
    out = df.merge(ctab, on="kepid", how="left")
    return out[(out["furlan_contam"].isna()) | (out["furlan_contam"] <= max_frac)].copy()


def read_berger2018_stellar_table(path) -> pd.DataFrame:
    from io import StringIO
    from pathlib import Path

    lines = Path(path).read_text(errors="replace").splitlines()
    header_idx = next((i for i, line in enumerate(lines) if line.startswith("KIC\t")), None)
    if header_idx is None:
        raise ValueError(f"Could not find Berger+2018 VizieR header in {path}")
    data = "\n".join([lines[header_idx]] + lines[header_idx + 3 :])
    tab = pd.read_csv(StringIO(data), sep="\t")
    out = tab.rename(columns={"KIC": "kepid", "Teff": "berger2018_teff", "R*": "berger2018_rad", "Bin": "berger2018_bin"})
    out["kepid"] = pd.to_numeric(out["kepid"], errors="coerce")
    out = out.dropna(subset=["kepid"]).copy()
    out["kepid"] = out["kepid"].astype(int)
    if "berger2018_bin" in out:
        out["berger2018_bin"] = pd.to_numeric(out["berger2018_bin"], errors="coerce")
    return out.drop_duplicates("kepid")


def add_velocity_features(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    # Angus+2022 provides Galactocentric Cartesian velocities where vy is the
    # rotational component used as V_phi in the Sagear-style figures. Keep this
    # direct convention as the primary classifier input.
    out["V_R"] = out["vx"]
    out["V_phi"] = out["vy"]
    out["V_z"] = out["vz"]
    out["V_perp"] = np.sqrt(out["vx"] ** 2 + out["vz"] ** 2)
    out["V_phi_proxy"] = out["V_phi"]
    out["V_perp_proxy"] = out["V_perp"]

    # Also compute an explicitly position-derived cylindrical basis for
    # diagnostics. It is not the default classifier input because it moves the
    # replication counts farther from Sagear's reported sample.
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    ra_col = "gaia_ra" if "gaia_ra" in out.columns else "ra"
    dec_col = "gaia_dec" if "gaia_dec" in out.columns else "dec"
    good = finite_columns(out, [ra_col, dec_col, "parallax", "vx", "vy", "vz"]) & (out["parallax"] > 0)
    out["galcen_X_pc"] = np.nan
    out["galcen_Y_pc"] = np.nan
    out["galcen_Z_pc"] = np.nan
    out["V_R_geom"] = np.nan
    out["V_phi_geom"] = np.nan
    out["V_perp_geom"] = np.nan
    out["galcen_X_astropy_pc"] = np.nan
    out["galcen_Y_astropy_pc"] = np.nan
    out["V_R_astropy"] = np.nan
    out["V_phi_astropy"] = np.nan
    out["V_perp_astropy"] = np.nan

    if good.any():
        coord_icrs = SkyCoord(
            ra=out.loc[good, ra_col].to_numpy(float) * u.deg,
            dec=out.loc[good, dec_col].to_numpy(float) * u.deg,
            distance=(1000.0 / out.loc[good, "parallax"].to_numpy(float)) * u.pc,
            frame="icrs",
        )
        coord = coord_icrs.galactic

        l = coord.l.radian
        b = coord.b.radian
        d = coord.distance.to_value(u.pc)
        r_sun = 8340.0
        z_sun = 20.8

        x = r_sun - d * np.cos(b) * np.cos(l)
        y = -d * np.cos(b) * np.sin(l)
        z = d * np.sin(b) + z_sun
        r = np.sqrt(x * x + y * y)

        vx = out.loc[good, "vx"].to_numpy(float)
        vy = out.loc[good, "vy"].to_numpy(float)
        vz = out.loc[good, "vz"].to_numpy(float)
        v_r = (x * vx + y * vy) / r
        v_phi = (-y * vx + x * vy) / r

        out.loc[good, "galcen_X_pc"] = x
        out.loc[good, "galcen_Y_pc"] = y
        out.loc[good, "galcen_Z_pc"] = z
        out.loc[good, "V_R_geom"] = v_r
        out.loc[good, "V_phi_geom"] = v_phi
        out.loc[good, "V_perp_geom"] = np.sqrt(v_r * v_r + vz * vz)

        # Old Hilldale pipeline convention: use Astropy's default
        # Galactocentric position transform, then project the Angus Cartesian
        # velocities onto the local cylindrical basis. This gives V_phi near
        # -220 km/s for disk rotation and is the convention used by the old
        # Toomre/classification code.
        gc = coord_icrs.transform_to("galactocentric")
        x_ast = gc.x.to_value(u.pc)
        y_ast = gc.y.to_value(u.pc)
        r_ast = np.sqrt(x_ast * x_ast + y_ast * y_ast)
        v_r_ast = (x_ast * vx + y_ast * vy) / r_ast
        v_phi_ast = (x_ast * vy - y_ast * vx) / r_ast
        out.loc[good, "galcen_X_astropy_pc"] = x_ast
        out.loc[good, "galcen_Y_astropy_pc"] = y_ast
        out.loc[good, "V_R_astropy"] = v_r_ast
        out.loc[good, "V_phi_astropy"] = v_phi_ast
        out.loc[good, "V_perp_astropy"] = np.sqrt(v_r_ast * v_r_ast + vz * vz)
    return out


def classify_with_fallback_gmm(df: pd.DataFrame, cfg: dict, velocity: str = "direct") -> tuple[pd.DataFrame, str]:
    if velocity == "direct":
        x_col, y_col = "V_phi", "V_perp"
        velocity_note = "direct Angus vy and sqrt(vx^2+vz^2)"
        score_sign = -1.0
    elif velocity == "old_astropy":
        x_col, y_col = "V_phi_astropy", "V_perp_astropy"
        velocity_note = "old Hilldale/Astropy cylindrical V_phi and V_perp"
        score_sign = 1.0
    else:
        raise ValueError(f"Unknown fallback velocity convention: {velocity}")

    use = finite_columns(df, [x_col, y_col])
    x = df.loc[use, [x_col, y_col]].to_numpy()
    gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=42, n_init=20)
    gmm.fit(x)
    # Thick component: lower rotational velocity and/or higher perpendicular velocity.
    score = score_sign * gmm.means_[:, 0] + 0.25 * gmm.means_[:, 1]
    thick_i = int(np.argmax(score))
    probs = np.full((len(df), 2), np.nan)
    probs[use.to_numpy()] = gmm.predict_proba(x)
    out = df.copy()
    out["P_thick"] = probs[:, thick_i]
    out["disk"] = np.where(out["P_thick"] > cfg["cuts"]["p_thick_threshold"], "thick", "thin")
    out.loc[out["P_thick"].isna(), "disk"] = "unclassified"
    note = (
        f"Fallback GMM on all matched Angus velocities using {velocity_note}; "
        f"means={gmm.means_.round(3).tolist()}, thick_component={thick_i}"
    )
    return out, note


def classify_with_chemical_gmm(df: pd.DataFrame, apogee: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, str]:
    kep_col = next((c for c in apogee.columns if c.lower() in {"kepid", "kic", "kic_id"}), None)
    feh_col = next((c for c in apogee.columns if c.lower() in {"feh", "apogee_feh", "[fe/h]", "fe_h", "m_h"}), None)
    mg_col = next((c for c in apogee.columns if c.lower() in {"mgfe", "apogee_mgfe", "[mg/fe]", "mg_fe"}), None)
    if kep_col is None or feh_col is None or mg_col is None:
        raise ValueError("APOGEE table needs KIC/kepid, [Fe/H], and [Mg/Fe] columns.")

    cal = build_apogee_velocity_calibration(apogee, kep_col, feh_col, mg_col, cfg)
    cal = cal[finite_columns(cal, ["apogee_feh", "apogee_mgfe", "V_phi", "V_perp"])]
    high_alpha = cal["apogee_mgfe"] > (
        cfg["cuts"]["high_alpha_slope"] * cal["apogee_feh"] + cfg["cuts"]["high_alpha_intercept"]
    )
    if high_alpha.nunique() < 2:
        raise ValueError("APOGEE calibration sample does not contain both high-alpha and low-alpha stars.")

    x_cal = cal[["V_phi", "V_perp"]].to_numpy()
    gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=42, n_init=20)
    labels = gmm.fit_predict(x_cal)
    frac_high = [float(high_alpha[labels == i].mean()) for i in range(2)]
    thick_i = int(np.argmax(frac_high))

    use = finite_columns(df, ["V_phi", "V_perp"])
    probs = np.full((len(df), 2), np.nan)
    probs[use.to_numpy()] = gmm.predict_proba(df.loc[use, ["V_phi", "V_perp"]].to_numpy())

    out = df.copy()
    out["P_thick"] = probs[:, thick_i]
    out["disk"] = np.where(out["P_thick"] > cfg["cuts"]["p_thick_threshold"], "thick", "thin")
    out.loc[out["P_thick"].isna(), "disk"] = "unclassified"
    return out, f"APOGEE-calibrated GMM; calibration_n={len(cal)}, high_alpha_fraction_by_component={frac_high}"


def build_apogee_velocity_calibration(
    apogee: pd.DataFrame,
    kep_col: str,
    feh_col: str,
    mg_col: str,
    cfg: dict,
) -> pd.DataFrame:
    cal = apogee.rename(
        columns={kep_col: "kepid", feh_col: "apogee_feh", mg_col: "apogee_mgfe"}
    ).copy()
    cal = cal.drop_duplicates("kepid")

    gk = read_gaia_kepler(cfg)
    keep_gaia_cols = ["kepid"]
    for col in ["ra", "dec", "parallax"]:
        if col in gk.columns and col not in cal.columns:
            keep_gaia_cols.append(col)
    cal = cal.merge(gk[keep_gaia_cols].drop_duplicates("kepid"), on="kepid", how="inner")
    cal = cal.merge(read_angus(cfg), on="kepid", how="inner", suffixes=("", "_angus"))
    cal = add_velocity_features(cal)

    teff_col = next((c for c in ["apogee_teff", "teff", "berger_teff"] if c in cal.columns), None)
    logg_col = next((c for c in ["apogee_logg", "logg", "berger_logg"] if c in cal.columns), None)
    if teff_col:
        cal = cal[pd.to_numeric(cal[teff_col], errors="coerce") < cfg["cuts"]["apogee_teff_max"]]
    if logg_col:
        cal = cal[pd.to_numeric(cal[logg_col], errors="coerce") > cfg["cuts"]["apogee_logg_min"]]
    return cal


def apply_convergence_cut(df: pd.DataFrame, conv: pd.DataFrame) -> pd.DataFrame:
    target_col = next((c for c in conv.columns if c in {"koi_target", "target", "koi_id"}), None)
    ok_col = next((c for c in conv.columns if c.lower() in {"converged", "ok", "usable"}), None)
    if target_col is None or ok_col is None:
        raise ValueError("Convergence table needs target and converged/ok/usable columns.")
    c = conv[[target_col, ok_col]].rename(columns={target_col: "koi_target", ok_col: "alderaan_converged"})
    out = df.merge(c, on="koi_target", how="left")
    return out[(out["alderaan_converged"].isna()) | (out["alderaan_converged"].astype(bool))].copy()


def disk_counts(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for disk in ["thin", "thick", "unclassified"]:
        for system in ["single", "multi"]:
            sub = df[(df.get("disk") == disk) & (df.get("system") == system)]
            if len(sub) == 0:
                continue
            rows.append(
                {
                    "disk": disk,
                    "system": system,
                    "planets": int(len(sub)),
                    "hosts": int(sub["kepid"].nunique()),
                    "median_p_thick": float(np.nanmedian(sub["P_thick"])) if "P_thick" in sub else np.nan,
                }
            )
    return pd.DataFrame(rows)


def compare_to_sagear(counts: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    targets = cfg["sagear_targets"]

    def get(disk: str, system: str, col: str) -> int:
        row = counts[(counts["disk"] == disk) & (counts["system"] == system)]
        return int(row[col].iloc[0]) if len(row) else 0

    rows = [
        ("thin_singles_planets", get("thin", "single", "planets")),
        ("thick_singles_planets", get("thick", "single", "planets")),
        ("thin_multi_planets", get("thin", "multi", "planets")),
        ("thick_multi_planets", get("thick", "multi", "planets")),
        ("thin_multi_hosts", get("thin", "multi", "hosts")),
        ("thick_multi_hosts", get("thick", "multi", "hosts")),
    ]
    return pd.DataFrame(
        [
            {
                "metric": name,
                "ours": ours,
                "sagear": int(targets[name]),
                "delta": ours - int(targets[name]),
            }
            for name, ours in rows
        ]
    )


if __name__ == "__main__":
    main()
