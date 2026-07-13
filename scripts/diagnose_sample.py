from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from sklearn.mixture import GaussianMixture

from common import (
    add_audit_row,
    add_target_and_system,
    finite_columns,
    load_config,
    output_dir,
    read_berger2018_stellar_table,
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
    parser.add_argument(
        "--chemical-model",
        choices=["conditioned", "pooled"],
        default="conditioned",
        help=(
            "Kinematic calibration model. 'conditioned' fits one Gaussian component to each chemically "
            "defined low/high-alpha sample; 'pooled' retains the older unsupervised pooled GMM sensitivity."
        ),
    )
    parser.add_argument(
        "--chemical-velocity",
        choices=["old_astropy", "direct"],
        default="old_astropy",
        help="Velocity basis for chemical calibration; old_astropy is the cylindrical convention shown in Sagear Figure 1.",
    )
    parser.add_argument(
        "--apogee-scope",
        choices=["planet_hosts", "all_kepler"],
        default="planet_hosts",
        help="Restrict the APOGEE calibration to the selected planet-host sample, as stated in the manuscript.",
    )
    parser.add_argument(
        "--apogee-flag-policy",
        choices=["finite", "all_zero"],
        default="finite",
        help="ASPCAP quality policy. 'finite' matches the explicit manuscript text; all_zero is a conservative sensitivity.",
    )
    parser.add_argument(
        "--chemical-prior",
        choices=["equal", "empirical"],
        default="equal",
        help=(
            "Prior mixture weight for chemically conditioned components. The manuscript does not state this; "
            "equal corresponds to a relative-likelihood classification, empirical is a sensitivity."
        ),
    )
    parser.add_argument(
        "--chemical-prior-high",
        type=float,
        default=None,
        help="Diagnostic explicit high-alpha prior in (0,1); overrides --chemical-prior.",
    )
    parser.add_argument(
        "--multiplicity-basis",
        choices=["raw_koi", "filtered"],
        default="raw_koi",
        help=(
            "Define single/multi from all non-false-positive KOIs around the host (raw_koi) or only planets "
            "surviving the analysis cuts (filtered). raw_koi preserves known system architecture."
        ),
    )
    parser.add_argument(
        "--berger2018-bin0",
        action="store_true",
        help="Sensitivity only: require Berger+2018 Bin=0, a plausible but unstated extra binary/sample cut.",
    )
    parser.add_argument("--out-tag", default=None, help="Append a label to output filenames without changing method status.")
    args = parser.parse_args()

    cfg = load_config(args.config)
    if args.berger2018_bin0:
        cfg["cuts"]["berger2018_bin_required"] = 0
    out_dir = output_dir()
    mode = "strict"
    missing: dict[str, str] = {}
    audit: list[dict[str, object]] = []

    koi = read_koi(cfg)
    add_audit_row(audit, "00_raw_koi_catalog", koi, "NASA Exoplanet Archive KOI table")

    df = koi[koi["koi_disposition"].isin(["CONFIRMED", "CANDIDATE"])].copy()
    raw_host_multiplicity = df.groupby("kepid")["kepoi_name"].size()
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
            b18_cols = [
                c
                for c in ["kepid", "berger2018_bin", "berger2018_evol", "berger2018_rad"]
                if c in b18.columns
            ]
            df = df.merge(b18[b18_cols], on="kepid", how="inner")
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

    if cuts.get("berger2018_evol_required") is not None:
        if "berger2018_evol" not in df.columns and b18_path and b18_path.exists():
            b18 = read_berger2018_stellar_table(b18_path)
            b18_cols = [
                c
                for c in ["kepid", "berger2018_evol", "berger2018_rad"]
                if c in b18.columns and c not in df.columns
            ]
            if b18_cols:
                df = df.merge(b18[["kepid", *b18_cols]], on="kepid", how="inner")
        if "berger2018_evol" not in df.columns:
            missing["berger2018_evol"] = (
                "Berger+2018 Evol flag is unavailable; evolutionary-state cut was not applied. "
                "Make sure berger2018_stellar points to a table with Evol."
            )
            add_audit_row(audit, "07c_berger2018_evol_SKIPPED", df, "missing Evol column")
        else:
            df = df[
                pd.to_numeric(df["berger2018_evol"], errors="coerce")
                == cuts["berger2018_evol_required"]
            ].copy()
            add_audit_row(
                audit,
                "07c_berger2018_evol",
                df,
                "Optional Berger+2018 Evol evolutionary-state cut; diagnostic until Sagear's exact host cut is confirmed",
            )

    angus = read_angus(cfg)
    df = df.merge(angus, on="kepid", how="inner", suffixes=("", "_angus"))
    add_audit_row(audit, "08_join_angus_velocities", df)

    df = add_velocity_features(df)
    apogee_path = root_path(cfg, "apogee_dr17_chemical")
    if apogee_path and apogee_path.exists() and not args.force_fallback_gmm:
        df, note, classifier_metadata = classify_with_chemical_gmm(
            df,
            pd.read_csv(apogee_path),
            cfg,
            model=args.chemical_model,
            velocity=args.chemical_velocity,
            scope=args.apogee_scope,
            flag_policy=args.apogee_flag_policy,
            component_prior=args.chemical_prior,
            prior_high_override=args.chemical_prior_high,
            return_metadata=True,
        )
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
    if args.multiplicity_basis == "raw_koi":
        df["raw_nonfp_host_multiplicity"] = df["kepid"].map(raw_host_multiplicity)
        df["system"] = np.where(df["raw_nonfp_host_multiplicity"] > 1, "multi", "single")
        add_audit_row(audit, "10_assign_single_multi_raw_nonfp_koi", df)
    else:
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

    file_mode = f"{mode}_{args.out_tag}" if args.out_tag else mode
    if "classifier_metadata" in locals():
        classifier_metadata["sample_file_mode"] = file_mode
        write_json(out_dir / f"chemical_classifier_{file_mode}.json", classifier_metadata)
    status = {
        "mode": mode,
        "sagear_equivalent": mode == "strict" and not missing,
        "missing_inputs": missing,
        "outputs": {
            "sample_audit": f"sample_audit_{file_mode}.csv",
            "disk_counts": f"disk_counts_{file_mode}.csv",
            "sagear_count_comparison": f"sagear_count_comparison_{file_mode}.csv",
            "canonical_sample": f"canonical_sample_{file_mode}.csv",
            "missing_inputs": f"missing_inputs_{file_mode}.json",
        },
    }

    audit_df.to_csv(out_dir / f"sample_audit_{file_mode}.csv", index=False)
    counts.to_csv(out_dir / f"disk_counts_{file_mode}.csv", index=False)
    comparison.to_csv(out_dir / f"sagear_count_comparison_{file_mode}.csv", index=False)
    df.to_csv(out_dir / f"canonical_sample_{file_mode}.csv", index=False)
    write_json(out_dir / f"missing_inputs_{file_mode}.json", missing)
    write_json(out_dir / f"pipeline_status_{file_mode}.json", status)

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


def add_velocity_features(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    # Angus+2022 provides Galactocentric Cartesian velocities. Preserve those
    # aliases explicitly; vx/vy are not cylindrical V_R/V_phi away from one
    # special position. The paper-defined primary fields are filled from the
    # position-dependent cylindrical projection below.
    out["V_R_cartesian_proxy"] = out["vx"]
    out["V_phi_cartesian_proxy"] = out["vy"]
    out["V_z"] = out["vz"]
    out["V_perp_cartesian_proxy"] = np.sqrt(out["vx"] ** 2 + out["vz"] ** 2)
    out["V_phi_proxy"] = out["V_phi_cartesian_proxy"]
    out["V_perp_proxy"] = out["V_perp_cartesian_proxy"]

    # Also compute an explicitly position-derived cylindrical basis for
    # diagnostics. It is not the default classifier input because it moves the
    # replication counts farther from Sagear's reported sample.
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    ra_col = "angus_ra" if "angus_ra" in out.columns else ("gaia_ra" if "gaia_ra" in out.columns else "ra")
    dec_col = "angus_dec" if "angus_dec" in out.columns else ("gaia_dec" if "gaia_dec" in out.columns else "dec")
    if "angus_dist_pc" in out.columns:
        distance_pc = pd.to_numeric(out["angus_dist_pc"], errors="coerce")
        good = finite_columns(out, [ra_col, dec_col, "angus_dist_pc", "vx", "vy", "vz"]) & (distance_pc > 0)
    else:
        distance_pc = 1000.0 / pd.to_numeric(out["parallax"], errors="coerce")
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
    out["V_R"] = np.nan
    out["V_phi"] = np.nan
    out["V_perp"] = np.nan

    if good.any():
        coord_icrs = SkyCoord(
            ra=out.loc[good, ra_col].to_numpy(float) * u.deg,
            dec=out.loc[good, dec_col].to_numpy(float) * u.deg,
            distance=distance_pc.loc[good].to_numpy(float) * u.pc,
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
        out.loc[good, "V_R"] = v_r_ast
        out.loc[good, "V_phi"] = v_phi_ast
        out.loc[good, "V_perp"] = np.sqrt(v_r_ast * v_r_ast + vz * vz)
    return out


def classify_with_fallback_gmm(df: pd.DataFrame, cfg: dict, velocity: str = "direct") -> tuple[pd.DataFrame, str]:
    if velocity == "direct":
        x_col, y_col = "V_phi_cartesian_proxy", "V_perp_cartesian_proxy"
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


def classify_with_chemical_gmm(
    df: pd.DataFrame,
    apogee: pd.DataFrame,
    cfg: dict,
    *,
    model: str = "conditioned",
    velocity: str = "old_astropy",
    scope: str = "planet_hosts",
    flag_policy: str = "finite",
    component_prior: str = "equal",
    prior_high_override: float | None = None,
    return_metadata: bool = False,
) -> tuple[pd.DataFrame, str] | tuple[pd.DataFrame, str, dict[str, object]]:
    kep_col = next((c for c in apogee.columns if c.lower() in {"kepid", "kic", "kic_id"}), None)
    feh_col = next((c for c in apogee.columns if c.lower() in {"feh", "apogee_feh", "[fe/h]", "fe_h", "m_h"}), None)
    mg_col = next((c for c in apogee.columns if c.lower() in {"mgfe", "apogee_mgfe", "[mg/fe]", "mg_fe"}), None)
    if kep_col is None or feh_col is None or mg_col is None:
        raise ValueError("APOGEE table needs KIC/kepid, [Fe/H], and [Mg/Fe] columns.")

    cal = build_apogee_velocity_calibration(apogee, kep_col, feh_col, mg_col, cfg)
    if scope == "planet_hosts":
        cal = cal[cal["kepid"].isin(pd.to_numeric(df["kepid"], errors="coerce"))].copy()
    elif scope != "all_kepler":
        raise ValueError(f"Unknown APOGEE scope: {scope}")

    if flag_policy == "all_zero":
        flag_cols = [c for c in ["aspcapflag", "starflag", "fe_h_flag", "mg_fe_flag"] if c in cal.columns]
        for col in flag_cols:
            cal = cal[pd.to_numeric(cal[col], errors="coerce").fillna(1).eq(0)].copy()
    elif flag_policy != "finite":
        raise ValueError(f"Unknown APOGEE flag policy: {flag_policy}")

    if velocity == "old_astropy":
        x_col, y_col = "V_phi_astropy", "V_perp_astropy"
    elif velocity == "direct":
        x_col, y_col = "V_phi_cartesian_proxy", "V_perp_cartesian_proxy"
    else:
        raise ValueError(f"Unknown chemical velocity convention: {velocity}")

    cal = cal[finite_columns(cal, ["apogee_feh", "apogee_mgfe", x_col, y_col])]
    high_alpha = cal["apogee_mgfe"] > (
        cfg["cuts"]["high_alpha_slope"] * cal["apogee_feh"] + cfg["cuts"]["high_alpha_intercept"]
    )
    if high_alpha.nunique() < 2:
        raise ValueError("APOGEE calibration sample does not contain both high-alpha and low-alpha stars.")

    x_cal = cal[[x_col, y_col]].to_numpy()
    use = finite_columns(df, [x_col, y_col])
    x_all = df.loc[use, [x_col, y_col]].to_numpy()
    probs = np.full(len(df), np.nan)
    loglike_low = np.full(len(df), np.nan)
    loglike_high = np.full(len(df), np.nan)

    if model == "conditioned":
        if int((~high_alpha).sum()) < 5 or int(high_alpha.sum()) < 5:
            raise ValueError(
                "Chemically conditioned calibration requires at least five low-alpha and high-alpha stars; "
                f"got low={int((~high_alpha).sum())}, high={int(high_alpha.sum())}."
            )
        low_model = GaussianMixture(
            n_components=1, covariance_type="full", random_state=42, reg_covar=1e-6
        ).fit(x_cal[~high_alpha.to_numpy()])
        high_model = GaussianMixture(
            n_components=1, covariance_type="full", random_state=42, reg_covar=1e-6
        ).fit(x_cal[high_alpha.to_numpy()])
        if prior_high_override is not None:
            if not 0.0 < prior_high_override < 1.0:
                raise ValueError("prior_high_override must be strictly between 0 and 1")
            prior_high = float(prior_high_override)
        elif component_prior == "equal":
            prior_high = 0.5
        elif component_prior == "empirical":
            prior_high = float(high_alpha.mean())
        else:
            raise ValueError(f"Unknown chemical component prior: {component_prior}")
        low_score = low_model.score_samples(x_all)
        high_score = high_model.score_samples(x_all)
        logp = np.column_stack(
            [
                low_score + np.log1p(-prior_high),
                high_score + np.log(prior_high),
            ]
        )
        probs[use.to_numpy()] = np.exp(logp[:, 1] - logsumexp(logp, axis=1))
        loglike_low[use.to_numpy()] = low_score
        loglike_high[use.to_numpy()] = high_score
        model_note = (
            "chemically conditioned Gaussian components; "
            f"low_mean={low_model.means_[0].round(3).tolist()}, "
            f"high_mean={high_model.means_[0].round(3).tolist()}, prior_high={prior_high:.6f}"
        )
        metadata = {
            "model": "chemically_conditioned_gaussian_components",
            "velocity": velocity,
            "x_column": x_col,
            "y_column": y_col,
            "low_mean": low_model.means_[0].tolist(),
            "low_covariance": low_model.covariances_[0].tolist(),
            "high_mean": high_model.means_[0].tolist(),
            "high_covariance": high_model.covariances_[0].tolist(),
            "prior_high": prior_high,
            "calibration_n": int(len(cal)),
            "high_alpha_n": int(high_alpha.sum()),
            "scope": scope,
            "flag_policy": flag_policy,
        }
    elif model == "pooled":
        gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=42, n_init=20)
        labels = gmm.fit_predict(x_cal)
        frac_high = [float(high_alpha[labels == i].mean()) for i in range(2)]
        thick_i = int(np.argmax(frac_high))
        probs[use.to_numpy()] = gmm.predict_proba(x_all)[:, thick_i]
        model_note = (
            f"pooled unsupervised GMM; means={gmm.means_.round(3).tolist()}, "
            f"high_alpha_fraction_by_component={frac_high}, thick_component={thick_i}"
        )
        metadata = {
            "model": "pooled_unsupervised_gmm",
            "velocity": velocity,
            "x_column": x_col,
            "y_column": y_col,
            "means": gmm.means_.tolist(),
            "covariances": gmm.covariances_.tolist(),
            "weights": gmm.weights_.tolist(),
            "thick_component": thick_i,
            "calibration_n": int(len(cal)),
            "high_alpha_n": int(high_alpha.sum()),
            "scope": scope,
            "flag_policy": flag_policy,
        }
    else:
        raise ValueError(f"Unknown chemical model: {model}")

    out = df.copy()
    out["P_thick"] = probs
    out["kinematic_loglike_thin"] = loglike_low
    out["kinematic_loglike_thick"] = loglike_high
    out["disk"] = np.where(out["P_thick"] > cfg["cuts"]["p_thick_threshold"], "thick", "thin")
    out.loc[out["P_thick"].isna(), "disk"] = "unclassified"
    note = (
        f"APOGEE-calibrated model={model}; calibration_n={len(cal)}, "
        f"high_alpha_n={int(high_alpha.sum())}, velocity={velocity}, scope={scope}, "
        f"flag_policy={flag_policy}, component_prior={component_prior}, "
        f"prior_high_override={prior_high_override}; {model_note}"
    )
    return (out, note, metadata) if return_metadata else (out, note)


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
