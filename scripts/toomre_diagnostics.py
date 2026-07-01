from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from sklearn.mixture import GaussianMixture

from common import finite_columns, load_config, output_dir, read_angus, root_path, write_json
from diagnose_sample import add_velocity_features, build_apogee_velocity_calibration


def main() -> None:
    parser = argparse.ArgumentParser(description="Diagnose Toomre coordinate and disk-classifier conventions.")
    parser.add_argument("--config", default=None)
    parser.add_argument("--sample", default=None)
    parser.add_argument("--out-prefix", default="toomre_diagnostics")
    args = parser.parse_args()

    cfg = load_config(args.config)
    out_dir = output_dir()
    sample_path = Path(args.sample) if args.sample else out_dir / "canonical_sample_diagnostic.csv"
    sample = pd.read_csv(sample_path)
    sample = ensure_velocity_aliases(sample)
    if {"V_phi_astropy", "V_perp_astropy"} - set(sample.columns):
        sample = add_velocity_features(sample)

    old_angus = read_old_angus_astropy_catalog(cfg)

    apogee_path = root_path(cfg, "apogee_dr17_chemical")
    apogee = pd.read_csv(apogee_path) if apogee_path and apogee_path.exists() else None
    calibration = None
    high_alpha = None
    if apogee is not None:
        calibration = build_apogee_velocity_calibration(apogee, "kepid", "feh", "mgfe", cfg)
        calibration = ensure_velocity_aliases(calibration)
        calibration = calibration[finite_columns(calibration, ["apogee_feh", "apogee_mgfe", "V_phi", "V_perp"])]
        high_alpha = calibration["apogee_mgfe"] > (
            cfg["cuts"]["high_alpha_slope"] * calibration["apogee_feh"] + cfg["cuts"]["high_alpha_intercept"]
        )

    variants: list[dict[str, object]] = []
    variants.append(run_unsupervised_variant(sample, sample, "planet_host_gmm_direct", "V_phi", "V_perp", cfg))
    variants.append(run_unsupervised_variant(sample, sample, "planet_host_gmm_geom", "V_phi_geom", "V_perp_geom", cfg))
    variants.append(
        run_unsupervised_variant(sample, sample, "planet_host_gmm_old_astropy", "V_phi_astropy", "V_perp_astropy", cfg)
    )
    if old_angus is not None:
        old_train = old_angus[
            finite_columns(old_angus, ["V_phi_astropy", "V_perp_astropy"])
            & (np.abs(old_angus["V_phi_astropy"] + 220.0) < 200.0)
            & (old_angus["V_perp_astropy"] < 300.0)
        ].copy()
        variants.append(
            run_unsupervised_variant(
                sample,
                old_train,
                "old_pipeline_kicwide_gmm_astropy",
                "V_phi_astropy",
                "V_perp_astropy",
                cfg,
            )
        )
    if calibration is not None and high_alpha is not None:
        variants.append(
            run_calibration_unsupervised_variant(
                sample, calibration, high_alpha, "apogee_calibration_gmm_direct", "V_phi", "V_perp", cfg
            )
        )
        variants.append(
            run_calibration_unsupervised_variant(
                sample, calibration, high_alpha, "apogee_calibration_gmm_geom", "V_phi_geom", "V_perp_geom", cfg
            )
        )
        variants.append(
            run_calibration_unsupervised_variant(
                sample,
                calibration,
                high_alpha,
                "apogee_calibration_gmm_old_astropy",
                "V_phi_astropy",
                "V_perp_astropy",
                cfg,
            )
        )
        variants.extend(run_supervised_prior_sweep(sample, calibration, high_alpha, "V_phi", "V_perp", cfg))
        variants.extend(
            run_supervised_component_grid(sample, calibration, high_alpha, "V_phi", "V_perp", cfg)
        )

    variants_df = pd.DataFrame(variants)
    variants_path = out_dir / f"{args.out_prefix}_classifier_variants.csv"
    variants_df.to_csv(variants_path, index=False)

    coord_summary = coordinate_summary(sample, calibration, high_alpha)
    coord_summary_path = out_dir / f"{args.out_prefix}_coordinate_summary.csv"
    coord_summary.to_csv(coord_summary_path, index=False)

    plot_path = out_dir / f"{args.out_prefix}_plot.png"
    make_toomre_plot(sample, cfg, plot_path, old_angus)

    metadata = {
        "sample": str(sample_path),
        "classifier_variants": str(variants_path),
        "coordinate_summary": str(coord_summary_path),
        "plot": str(plot_path),
        "interpretation_hint": (
            "Sagear-style Toomre display uses x=-vy from Angus+2022 so the rotating disk appears at negative V_phi. "
            "Classifier counts are very sensitive to direct Angus vs position-derived geometric cylindrical velocities."
        ),
    }
    write_json(out_dir / f"{args.out_prefix}_metadata.json", metadata)

    print(variants_df.to_string(index=False))
    print(f"\nWrote: {variants_path}")
    print(f"Wrote: {coord_summary_path}")
    print(f"Wrote: {plot_path}")


def ensure_velocity_aliases(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "V_phi" not in out and "vy" in out:
        out["V_phi"] = out["vy"]
    if "V_perp" not in out and {"vx", "vz"}.issubset(out.columns):
        out["V_perp"] = np.sqrt(out["vx"] ** 2 + out["vz"] ** 2)
    if "V_phi_geom" not in out and "V_phi" in out:
        out["V_phi_geom"] = np.nan
    if "V_perp_geom" not in out and "V_perp" in out:
        out["V_perp_geom"] = np.nan
    if "V_phi_astropy" not in out and "V_phi" in out:
        out["V_phi_astropy"] = np.nan
    if "V_perp_astropy" not in out and "V_perp" in out:
        out["V_perp_astropy"] = np.nan
    return out


def read_old_angus_astropy_catalog(cfg: dict) -> pd.DataFrame | None:
    """Parse Angus+2022 raw table and reproduce the old pipeline coordinates."""
    raw_path = Path(cfg["_root"]) / "angus_velocities.dat.gz"
    if not raw_path.exists():
        return None

    rows = []
    with gzip.open(raw_path, "rt") as f:
        for line in f:
            try:
                kepid = int(line[7:15])
                ra = float(line[36:43])
                dec = float(line[50:56])
                plx = float(line[63:69])
                dist_text = line[76:83].strip()
                dist_pc = float(dist_text) if dist_text else np.nan
                vx_c = line[175:182].strip()
                vy_c = line[195:202].strip()
                vz_c = line[214:221].strip()
                vx = float(vx_c) if vx_c else float(line[183:189])
                vy = float(vy_c) if vy_c else float(line[203:208])
                vz = float(vz_c) if vz_c else float(line[222:228])
                rows.append((kepid, ra, dec, plx, dist_pc, vx, vy, vz, bool(vx_c)))
            except Exception:
                continue
    ang = pd.DataFrame(rows, columns=["kepid", "ra", "dec", "plx", "dist_pc", "vx", "vy", "vz", "has_rv"])
    if ang.empty:
        return None
    bad_dist = ~np.isfinite(ang["dist_pc"])
    ok_plx = bad_dist & np.isfinite(ang["plx"]) & (ang["plx"] > 0)
    ang.loc[ok_plx, "dist_pc"] = 1000.0 / ang.loc[ok_plx, "plx"]
    good = finite_columns(ang, ["ra", "dec", "dist_pc", "vx", "vy", "vz"]) & (ang["dist_pc"] > 0)
    ang["V_R_astropy"] = np.nan
    ang["V_phi_astropy"] = np.nan
    ang["V_perp_astropy"] = np.nan

    if good.any():
        from astropy.coordinates import Galactocentric, SkyCoord
        import astropy.units as u

        sc = SkyCoord(
            ra=ang.loc[good, "ra"].to_numpy(float) * u.deg,
            dec=ang.loc[good, "dec"].to_numpy(float) * u.deg,
            distance=ang.loc[good, "dist_pc"].to_numpy(float) * u.pc,
        )
        gc = sc.transform_to(Galactocentric())
        x = gc.x.to_value(u.pc)
        y = gc.y.to_value(u.pc)
        r = np.sqrt(x * x + y * y)
        vx = ang.loc[good, "vx"].to_numpy(float)
        vy = ang.loc[good, "vy"].to_numpy(float)
        vz = ang.loc[good, "vz"].to_numpy(float)
        v_r = (x * vx + y * vy) / r
        v_phi = (x * vy - y * vx) / r
        ang.loc[good, "V_R_astropy"] = v_r
        ang.loc[good, "V_phi_astropy"] = v_phi
        ang.loc[good, "V_perp_astropy"] = np.sqrt(v_r * v_r + vz * vz)
    return ang


def run_unsupervised_variant(
    sample: pd.DataFrame,
    train: pd.DataFrame,
    name: str,
    x_col: str,
    y_col: str,
    cfg: dict,
) -> dict[str, object]:
    use_train = finite_columns(train, [x_col, y_col])
    use_sample = finite_columns(sample, [x_col, y_col])
    if use_train.sum() < 20 or use_sample.sum() < 20:
        return variant_row(name, "skipped_insufficient_finite", np.nan, sample, np.full(len(sample), np.nan), cfg)

    gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=42, n_init=20)
    x_train = train.loc[use_train, [x_col, y_col]].to_numpy()
    gmm.fit(x_train)
    score = -gmm.means_[:, 0] + 0.25 * gmm.means_[:, 1]
    thick_i = int(np.argmax(score))

    probs = np.full(len(sample), np.nan)
    probs[use_sample.to_numpy()] = gmm.predict_proba(sample.loc[use_sample, [x_col, y_col]].to_numpy())[:, thick_i]
    note = f"{x_col},{y_col}; means={gmm.means_.round(3).tolist()}; thick_component={thick_i}"
    return variant_row(name, note, thick_i, sample, probs, cfg)


def run_calibration_unsupervised_variant(
    sample: pd.DataFrame,
    calibration: pd.DataFrame,
    high_alpha: pd.Series,
    name: str,
    x_col: str,
    y_col: str,
    cfg: dict,
) -> dict[str, object]:
    use_cal = finite_columns(calibration, [x_col, y_col])
    use_sample = finite_columns(sample, [x_col, y_col])
    if use_cal.sum() < 20 or use_sample.sum() < 20:
        return variant_row(name, "skipped_insufficient_finite", np.nan, sample, np.full(len(sample), np.nan), cfg)

    gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=42, n_init=20)
    x_cal = calibration.loc[use_cal, [x_col, y_col]].to_numpy()
    labels = gmm.fit_predict(x_cal)
    high = high_alpha.loc[use_cal].reset_index(drop=True)
    frac_high = [float(high[labels == i].mean()) if np.any(labels == i) else np.nan for i in range(2)]
    thick_i = int(np.nanargmax(frac_high))

    probs = np.full(len(sample), np.nan)
    probs[use_sample.to_numpy()] = gmm.predict_proba(sample.loc[use_sample, [x_col, y_col]].to_numpy())[:, thick_i]
    note = (
        f"{x_col},{y_col}; cal_n={int(use_cal.sum())}; means={gmm.means_.round(3).tolist()}; "
        f"frac_high={np.round(frac_high, 3).tolist()}; thick_component={thick_i}"
    )
    return variant_row(name, note, thick_i, sample, probs, cfg)


def run_supervised_prior_sweep(
    sample: pd.DataFrame,
    calibration: pd.DataFrame,
    high_alpha: pd.Series,
    x_col: str,
    y_col: str,
    cfg: dict,
) -> list[dict[str, object]]:
    use_cal = finite_columns(calibration, [x_col, y_col])
    use_sample = finite_columns(sample, [x_col, y_col])
    high = high_alpha.loc[use_cal]
    cal = calibration.loc[use_cal]
    if high.sum() < 10 or (~high).sum() < 10 or use_sample.sum() < 20:
        return []

    thin_g = GaussianMixture(n_components=1, covariance_type="full", random_state=42)
    thick_g = GaussianMixture(n_components=1, covariance_type="full", random_state=42)
    thin_g.fit(cal.loc[~high, [x_col, y_col]].to_numpy())
    thick_g.fit(cal.loc[high, [x_col, y_col]].to_numpy())

    x_sample = sample.loc[use_sample, [x_col, y_col]].to_numpy()
    log_l_thin = thin_g.score_samples(x_sample)
    log_l_thick = thick_g.score_samples(x_sample)
    rows = []
    for p_thick_prior in [0.25, 0.35, 0.40, 0.425, 0.45, 0.50]:
        probs = np.full(len(sample), np.nan)
        log_num = log_l_thick + np.log(p_thick_prior)
        log_den = logsumexp(
            np.vstack([log_l_thin + np.log(1.0 - p_thick_prior), log_num]),
            axis=0,
        )
        probs[use_sample.to_numpy()] = np.exp(log_num - log_den)
        note = (
            f"single Gaussian per chemical class; p_thick_prior={p_thick_prior}; "
            f"thin_mean={thin_g.means_.round(3).tolist()}; thick_mean={thick_g.means_.round(3).tolist()}"
        )
        rows.append(variant_row(f"supervised_gaussian_prior_{p_thick_prior:.3f}", note, np.nan, sample, probs, cfg))
    return rows


def run_supervised_component_grid(
    sample: pd.DataFrame,
    calibration: pd.DataFrame,
    high_alpha: pd.Series,
    x_col: str,
    y_col: str,
    cfg: dict,
) -> list[dict[str, object]]:
    """Fit separate mixtures to chemical classes and report the best prior.

    This tests the interpretation that Sagear's "GMM calibrated using chemical
    abundances" means a supervised kinematic likelihood for the low-alpha and
    high-alpha calibration samples, not a fully unsupervised GMM on all APOGEE
    dwarfs.
    """
    use_cal = finite_columns(calibration, [x_col, y_col])
    use_sample = finite_columns(sample, [x_col, y_col])
    high = high_alpha.loc[use_cal]
    cal = calibration.loc[use_cal]
    if high.sum() < 10 or (~high).sum() < 10 or use_sample.sum() < 20:
        return []

    rows = []
    for n_thin, n_thick in [(1, 1), (2, 1), (1, 2), (2, 2)]:
        thin_g = GaussianMixture(n_components=n_thin, covariance_type="full", random_state=42, n_init=10)
        thick_g = GaussianMixture(n_components=n_thick, covariance_type="full", random_state=42, n_init=10)
        thin_g.fit(cal.loc[~high, [x_col, y_col]].to_numpy())
        thick_g.fit(cal.loc[high, [x_col, y_col]].to_numpy())

        x_sample = sample.loc[use_sample, [x_col, y_col]].to_numpy()
        log_l_thin = thin_g.score_samples(x_sample)
        log_l_thick = thick_g.score_samples(x_sample)
        best: tuple[float, float, np.ndarray] | None = None
        for p_thick_prior in np.linspace(0.25, 0.55, 61):
            probs = np.full(len(sample), np.nan)
            log_num = log_l_thick + np.log(p_thick_prior)
            log_den = logsumexp(
                np.vstack([log_l_thin + np.log(1.0 - p_thick_prior), log_num]),
                axis=0,
            )
            probs[use_sample.to_numpy()] = np.exp(log_num - log_den)
            row = variant_row(
                f"supervised_mix_{n_thin}thin_{n_thick}thick_prior_{p_thick_prior:.3f}",
                "",
                np.nan,
                sample,
                probs,
                cfg,
            )
            score = (
                abs(row["delta_thin_singles"])
                + abs(row["delta_thick_singles"])
                + abs(row["delta_thin_multi_planets"])
                + abs(row["delta_thick_multi_planets"])
            )
            if best is None or score < best[0]:
                best = (float(score), float(p_thick_prior), probs)
        assert best is not None
        _, best_prior, best_probs = best
        note = (
            f"separate chemical-class mixtures; n_thin={n_thin}; n_thick={n_thick}; "
            f"best_prior={best_prior:.3f}; thin_means={thin_g.means_.round(3).tolist()}; "
            f"thick_means={thick_g.means_.round(3).tolist()}"
        )
        rows.append(
            variant_row(
                f"supervised_mix_{n_thin}thin_{n_thick}thick_best_prior",
                note,
                np.nan,
                sample,
                best_probs,
                cfg,
            )
        )
    return rows


def variant_row(
    name: str,
    note: str,
    thick_component: float,
    sample: pd.DataFrame,
    p_thick: np.ndarray,
    cfg: dict,
) -> dict[str, object]:
    disk = np.where(p_thick > cfg["cuts"]["p_thick_threshold"], "thick", "thin").astype(object)
    disk[~np.isfinite(p_thick)] = "unclassified"

    def planets(d: str, s: str) -> int:
        return int(((disk == d) & (sample["system"] == s)).sum())

    def hosts(d: str, s: str) -> int:
        return int(sample.loc[(disk == d) & (sample["system"] == s), "kepid"].nunique())

    targets = cfg["sagear_targets"]
    return {
        "variant": name,
        "thin_singles": planets("thin", "single"),
        "thick_singles": planets("thick", "single"),
        "thin_multi_planets": planets("thin", "multi"),
        "thick_multi_planets": planets("thick", "multi"),
        "thin_multi_hosts": hosts("thin", "multi"),
        "thick_multi_hosts": hosts("thick", "multi"),
        "delta_thin_singles": planets("thin", "single") - int(targets["thin_singles_planets"]),
        "delta_thick_singles": planets("thick", "single") - int(targets["thick_singles_planets"]),
        "delta_thin_multi_planets": planets("thin", "multi") - int(targets["thin_multi_planets"]),
        "delta_thick_multi_planets": planets("thick", "multi") - int(targets["thick_multi_planets"]),
        "median_p_thick": float(np.nanmedian(p_thick)),
        "thick_component": thick_component,
        "note": note,
    }


def coordinate_summary(
    sample: pd.DataFrame,
    calibration: pd.DataFrame | None,
    high_alpha: pd.Series | None,
) -> pd.DataFrame:
    rows = []
    rows.extend(summary_rows(sample, "planet_sample", sample.get("disk"), sample.get("system")))
    if calibration is not None and high_alpha is not None:
        chem_label = np.where(high_alpha, "chemical_high_alpha", "chemical_low_alpha")
        rows.extend(summary_rows(calibration, "apogee_calibration", chem_label, None))
    return pd.DataFrame(rows)


def summary_rows(df: pd.DataFrame, source: str, group_a, group_b) -> list[dict[str, object]]:
    rows = []
    labels_a = pd.Series(group_a if group_a is not None else "all", index=df.index, dtype="object")
    labels_b = pd.Series(group_b if group_b is not None else "all", index=df.index, dtype="object")
    for (a, b), sub_idx in pd.DataFrame({"a": labels_a, "b": labels_b}).groupby(["a", "b"]).groups.items():
        sub = df.loc[sub_idx]
        for x_col, y_col, convention in [
            ("V_phi", "V_perp", "direct_angus_positive_rotation"),
            ("V_phi_geom", "V_perp_geom", "position_derived_geometric"),
            ("V_phi_astropy", "V_perp_astropy", "old_pipeline_astropy_galactocentric"),
        ]:
            if x_col not in sub or y_col not in sub:
                continue
            finite = finite_columns(sub, [x_col, y_col])
            if finite.sum() == 0:
                continue
            rows.append(
                {
                    "source": source,
                    "group_a": a,
                    "group_b": b,
                    "convention": convention,
                    "n": int(finite.sum()),
                    "median_vphi_positive": float(np.nanmedian(sub.loc[finite, x_col])),
                    "median_vphi_sagear_plot": float(np.nanmedian(-sub.loc[finite, x_col])),
                    "median_vperp": float(np.nanmedian(sub.loc[finite, y_col])),
                    "p90_vperp": float(np.nanpercentile(sub.loc[finite, y_col], 90)),
                }
            )
    return rows


def make_toomre_plot(sample: pd.DataFrame, cfg: dict, path: Path, old_angus: pd.DataFrame | None) -> None:
    angus = read_angus(cfg)
    angus = angus[np.isfinite(angus[["vx", "vy", "vz"]]).all(axis=1)].copy()
    rng = np.random.default_rng(42)
    if len(angus) > 35000:
        angus = angus.iloc[rng.choice(np.arange(len(angus)), size=35000, replace=False)].copy()
    angus["x_sagear"] = -angus["vy"]
    angus["y_perp"] = np.sqrt(angus["vx"] ** 2 + angus["vz"] ** 2)

    sample = sample.copy()
    sample["x_sagear"] = -sample["V_phi"]
    sample["y_perp"] = sample["V_perp"]

    ncols = 3 if old_angus is not None else 2
    fig, axes = plt.subplots(1, ncols, figsize=(7 * ncols, 6), constrained_layout=True)
    if ncols == 2:
        axes = np.asarray(axes)
    panels = [
        (axes[0], "Sagear-style display: x = -Angus vy", "x_sagear"),
        (axes[1], "Internal direct-Angus classifier: x = +vy", "V_phi"),
    ]
    if old_angus is not None:
        panels.append((axes[2], "Old pipeline Astropy cylindrical", "V_phi_astropy"))

    for ax, title, x_col in panels:
        if x_col == "V_phi_astropy":
            bg = old_angus[finite_columns(old_angus, ["V_phi_astropy", "V_perp_astropy"])]
            if len(bg) > 35000:
                bg = bg.iloc[rng.choice(np.arange(len(bg)), size=35000, replace=False)].copy()
            bg_x = bg["V_phi_astropy"]
            bg_y = bg["V_perp_astropy"]
            host_y = sample["V_perp_astropy"]
        else:
            bg_x = angus["x_sagear"] if x_col == "x_sagear" else angus["vy"]
            bg_y = angus["y_perp"]
            host_y = sample["y_perp"]
        ax.scatter(bg_x, bg_y, s=2, c="0.72", alpha=0.18, linewidths=0, label="Angus Kepler stars")
        for disk, color, marker in [("thin", "#008b8b", "+"), ("thick", "#9b0000", "+")]:
            sub = sample[sample["disk"] == disk]
            ax.scatter(
                sub[x_col],
                host_y.loc[sub.index],
                s=32,
                c=color,
                marker=marker,
                alpha=0.82,
                linewidths=1.2,
                label=f"{disk} planet hosts",
            )
        ax.set_title(title)
        ax.set_xlabel(r"$V_\phi$ [km s$^{-1}$]")
        ax.set_ylabel(r"$\sqrt{V_R^2 + V_Z^2}$ [km s$^{-1}$]")
        ax.set_ylim(0, 200)
        if x_col in {"x_sagear", "V_phi_astropy"}:
            ax.set_xlim(-330, -100)
        else:
            ax.set_xlim(100, 330)
        ax.grid(alpha=0.18)
        ax.legend(loc="upper left", frameon=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
