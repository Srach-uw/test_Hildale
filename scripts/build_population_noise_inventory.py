"""Build the fixed-cohort source inventory for the population noise audit."""

import argparse
import hashlib
import json
import os
import os
from pathlib import Path

import pandas as pd
from astropy.io import fits


SUMMARY_SHA256 = "020008a34eec7993affb8bb67ff553e32b9e4c5de413800e79b152d7273b9bae"
EXPECTED_EXCLUSIONS = {"K00972.02", "K04388.02", "K02706.01", "K03165.01"}
SOURCE_CLASSES = {
    "inputs/alderaan_archive_raw_20260612/ALDERAAN_posteriors": "source_original_archive",
    "alderaan_project/Results/sagear_missing": "source_recovery_local",
    "tmp/published_inventory_recovery_postprocess_20260726_final/extracted/Results/sagear_published_inventory_missing": "source_recovery_published",
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _bool_series(values, label):
    normalized = values.astype(str).str.strip().str.lower()
    if not normalized.isin(["true", "false"]).all():
        raise ValueError(f"invalid {label} values")
    return normalized.eq("true")


def _fits_identity(path):
    header = fits.getheader(path, 0)
    return str(header.get("TARGET", "")).strip(), int(header.get("NPL", -1))


def _deduplicate_existing(paths):
    unique = []
    for path in paths:
        if not path.is_file():
            continue
        if any(os.path.samefile(path, prior) for prior in unique):
            continue
        unique.append(path.resolve())
    return unique


def _companion_paths(target, source_file, staging_file, companion_roots):
    locations = [Path(source_file).parent, Path(staging_file).parent]
    for root in companion_roots:
        root = Path(root)
        locations.extend([root, root / target])
    locations = list(dict.fromkeys(path.resolve() for path in locations))

    def find(name):
        return _deduplicate_existing([location / name for location in locations])

    return find, locations


def _presence(paths):
    if not paths:
        return "missing", ""
    if len(paths) == 1:
        return "present_unverified_source", str(paths[0])
    return "ambiguous_present_unverified_source", "|".join(map(str, paths))


def _parse_alternatives(row):
    files = str(row.get("candidate_files", "")).split("|")
    counts = str(row.get("candidate_planet_counts", "")).split("|")
    count_by_file = {}
    for item in counts:
        if ":" in item:
            path, count = item.rsplit(":", 1)
            count_by_file[str(Path(path).resolve())] = int(count)
    selected = str(Path(row["selected_file"]).resolve())
    return [
        (Path(path), count_by_file.get(str(Path(path).resolve())))
        for path in files
        if path and str(Path(path).resolve()) != selected
    ]


def _write_conflicts(output_dir, conflicts, scope):
    output_dir.mkdir(parents=True, exist_ok=True)
    columns = ["koi_target", "conflict_type", "detail"]
    pd.DataFrame(conflicts, columns=columns).to_csv(output_dir / "conflicts.csv", index=False)
    (output_dir / "scope.json").write_text(json.dumps(scope, indent=2) + "\n", encoding="utf-8")


def build_inventory(
    summary_path,
    source_manifest_path,
    leverage_dir,
    output_dir,
    *,
    companion_roots=(),
    expected_summary_sha256=SUMMARY_SHA256,
    expected_planets=2119,
    expected_targets=1583,
    expected_exclusions=EXPECTED_EXCLUSIONS,
):
    summary_path = Path(summary_path)
    source_manifest_path = Path(source_manifest_path)
    leverage_dir = Path(leverage_dir)
    output_dir = Path(output_dir)
    actual_summary_hash = sha256(summary_path)
    if actual_summary_hash != expected_summary_sha256:
        raise ValueError("fixed summary SHA-256 mismatch")

    summary = pd.read_csv(summary_path)
    if summary["kepoi_name"].duplicated().any():
        raise ValueError("fixed summary has duplicate planet IDs")
    if "alderaan_planet_index" not in summary:
        raise ValueError("fixed summary is missing alderaan_planet_index")
    excluded = _bool_series(summary["qc_primary_exclude"], "qc_primary_exclude")
    excluded_ids = set(summary.loc[excluded, "kepoi_name"].astype(str))
    if excluded_ids != set(expected_exclusions):
        raise ValueError(f"unexpected primary exclusions: {sorted(excluded_ids)}")
    selected = summary.loc[~excluded].copy()
    if len(selected) != expected_planets or selected["koi_target"].nunique() != expected_targets:
        raise ValueError("fixed selected planet/target counts do not match the gate")

    manifest = pd.read_csv(source_manifest_path)
    if manifest["koi_target"].duplicated().any():
        raise ValueError("source manifest has duplicate targets")
    target_context = selected.groupby("koi_target", sort=True).agg(
        qc_selected_planet_count=("kepoi_name", "size"),
        selected_planet_ids=("kepoi_name", lambda values: "|".join(sorted(map(str, values)))),
        selected_planet_indices=("alderaan_planet_index", lambda values: "|".join(map(str, sorted(map(int, values))))),
        disk=("disk", lambda values: "|".join(sorted(set(map(str, values))))),
        system=("system", lambda values: "|".join(sorted(set(map(str, values))))),
    ).reset_index()
    all_ids = summary.groupby("koi_target")["kepoi_name"].agg(
        lambda values: "|".join(sorted(map(str, values)))
    )
    target_context["summary_planet_ids_including_primary_exclusions"] = target_context[
        "koi_target"
    ].map(all_ids)
    joined = target_context.merge(manifest, on="koi_target", how="left", validate="one_to_one")
    if joined["selected_file"].isna().any():
        missing = joined.loc[joined["selected_file"].isna(), "koi_target"].tolist()
        raise ValueError(f"source manifest is missing fixed targets: {missing[:10]}")

    conflicts = []
    inventory_rows = []
    for raw in joined.to_dict("records"):
        target = str(raw["koi_target"])
        selected_root = str(raw["selected_root"])
        source_class = SOURCE_CLASSES.get(selected_root, "unverified_source_root")
        if source_class == "unverified_source_root":
            conflicts.append({"koi_target": target, "conflict_type": "unknown_source_root", "detail": selected_root})

        source_file = Path(raw["selected_file"])
        staging_file = Path(raw["staged_file"])
        expected_sha = str(raw["selected_sha256"]).lower()
        required_npl = int(raw["required_npl"])
        source_state = "missing"
        source_hash = ""
        fits_target = ""
        fits_npl = -1
        if source_file.is_file():
            source_hash = sha256(source_file)
            fits_target, fits_npl = _fits_identity(source_file)
            source_state = "verified"
            if source_hash != expected_sha:
                source_state = "hash_mismatch"
                conflicts.append({"koi_target": target, "conflict_type": "source_hash_mismatch", "detail": str(source_file)})
            if fits_target != target or fits_npl < required_npl:
                source_state = "identity_mismatch"
                conflicts.append({"koi_target": target, "conflict_type": "source_fits_identity", "detail": f"TARGET={fits_target};NPL={fits_npl};expected_NPL={required_npl}"})
        else:
            conflicts.append({"koi_target": target, "conflict_type": "source_result_missing", "detail": str(source_file)})

        staging_state = "missing"
        if staging_file.is_file() and source_file.is_file():
            staging_hash = source_hash if os.path.samefile(source_file, staging_file) else sha256(staging_file)
            staging_state = "verified_same_canonical_sha" if staging_hash == expected_sha else "hash_mismatch"
            if staging_hash != expected_sha:
                conflicts.append({"koi_target": target, "conflict_type": "staging_hash_mismatch", "detail": str(staging_file)})

        alternative_paths = []
        alternative_states = []
        for alternative, recorded_npl in _parse_alternatives(raw):
            alternative_paths.append(str(alternative))
            if not alternative.is_file():
                alternative_states.append("missing")
                continue
            alt_target, alt_npl = _fits_identity(alternative)
            alt_hash = sha256(alternative)
            if alt_target != target:
                state = "target_identity_mismatch"
            elif alt_npl != fits_npl or (recorded_npl is not None and recorded_npl < required_npl):
                state = "alternative_different_system_size"
            elif alt_hash == expected_sha:
                state = "alternative_same_content"
            else:
                state = "conflicting_full_system_fit"
                conflicts.append({"koi_target": target, "conflict_type": state, "detail": str(alternative)})
            alternative_states.append(state)

        find, _ = _companion_paths(target, source_file, staging_file, companion_roots)
        lc_paths = find(f"{target}_lc_filtered.fits")
        sc_paths = find(f"{target}_sc_filtered.fits")
        parameter_paths = find(f"{target}_transit_parameters.csv")
        timing_paths = []
        missing_timing = []
        ambiguous_timing = []
        full_joint_npl = fits_npl if fits_npl >= required_npl else required_npl
        selected_indices = [int(value) for value in str(raw["selected_planet_indices"]).split("|")]
        if (
            len(selected_indices) != len(set(selected_indices))
            or any(index < 0 or index >= full_joint_npl for index in selected_indices)
        ):
            conflicts.append({
                "koi_target": target,
                "conflict_type": "selected_planet_index_incompatible",
                "detail": f"indices={selected_indices};FITS_NPL={fits_npl}",
            })
        for planet_index in range(full_joint_npl):
            paths = find(f"{target}_{planet_index:02d}_quick.ttvs")
            if not paths:
                missing_timing.append(planet_index)
            elif len(paths) > 1:
                ambiguous_timing.append(planet_index)
            timing_paths.extend(paths)
        for label, paths in (("lc", lc_paths), ("sc", sc_paths), ("parameters", parameter_paths)):
            if len(paths) > 1:
                conflicts.append({"koi_target": target, "conflict_type": "ambiguous_companion_location", "detail": f"{label}:{'|'.join(map(str, paths))}"})
        if ambiguous_timing:
            conflicts.append({"koi_target": target, "conflict_type": "ambiguous_timing_location", "detail": str(ambiguous_timing)})

        lc_state, lc_value = _presence(lc_paths)
        sc_state, sc_value = _presence(sc_paths)
        parameters_state, parameters_value = _presence(parameter_paths)
        all_timing_present = not missing_timing and not ambiguous_timing
        if not lc_paths or not parameter_paths or missing_timing:
            readiness = "blocked_missing_companions"
        else:
            readiness = "blocked_unverified_companion_provenance"
        inventory_rows.append({
            "koi_target": target,
            "disk": raw["disk"],
            "system": raw["system"],
            "selected_planet_count_qc": int(raw["qc_selected_planet_count"]),
            "manifest_required_planet_count": required_npl,
            "full_joint_planet_count": full_joint_npl,
            "joint_fit_membership_state": (
                "fits_has_additional_joint_planets" if full_joint_npl > required_npl else "matches_manifest_required_count"
            ),
            "selected_planet_ids": raw["selected_planet_ids"],
            "selected_planet_indices": raw["selected_planet_indices"],
            "summary_planet_ids_including_primary_exclusions": raw["summary_planet_ids_including_primary_exclusions"],
            "source_class": source_class,
            "selected_root": selected_root,
            "source_result_file": str(source_file),
            "source_result_state": source_state,
            "expected_result_sha256": expected_sha,
            "actual_result_sha256": source_hash,
            "fits_target": fits_target,
            "fits_npl": fits_npl,
            "staging_file": str(staging_file),
            "staging_state": staging_state,
            "alternative_fit_count": len(alternative_paths),
            "alternative_fit_paths": "|".join(alternative_paths),
            "alternative_fit_states": "|".join(alternative_states),
            "lc_filtered_state": lc_state,
            "lc_filtered_paths": lc_value,
            "sc_filtered_state": sc_state,
            "sc_filtered_paths": sc_value,
            "parameters_state": parameters_state,
            "parameters_paths": parameters_value,
            "timing_files_present": len(timing_paths),
            "timing_files_required": full_joint_npl,
            "all_full_system_timing_present": all_timing_present,
            "missing_timing_indices": "|".join(map(str, missing_timing)),
            "timing_paths": "|".join(map(str, timing_paths)),
            "companion_source_match": "unverified_no_companion_hash_manifest",
            "replay_readiness": readiness,
        })

    scope = {
        "status": "failed_conflicts" if conflicts else "passed_inventory",
        "fixed_summary": str(summary_path),
        "fixed_summary_sha256": actual_summary_hash,
        "source_manifest": str(source_manifest_path),
        "source_manifest_sha256": sha256(source_manifest_path),
        "selected_planets": len(selected),
        "selected_targets": selected["koi_target"].nunique(),
        "companion_roots": [str(Path(path)) for path in companion_roots],
        "limitations": "Presence is not source provenance. No target is replay-ready without mapped companion hashes and a verified likelihood convention. Pilot rows are candidates, not authorized fits.",
        "npl_interpretation": "FITS NPL is the full jointly fit system size and may exceed manifest required_npl or QC-selected planet count; timing inventory uses FITS NPL.",
        "superseded_false_conflict_note": "An earlier local run incorrectly treated FITS NPL greater than manifest required_npl as a source conflict. Those were audit-code false positives, not source-data failures; refreshed outputs use the corrected rule.",
    }
    _write_conflicts(output_dir, conflicts, scope)
    if conflicts:
        raise RuntimeError(f"inventory stopped with {len(conflicts)} conflicts; see conflicts.csv")

    inventory = pd.DataFrame(inventory_rows).sort_values("koi_target")
    inventory.to_csv(output_dir / "target_inventory.csv", index=False)
    source_summary = inventory.groupby("source_class", sort=True).agg(
        targets=("koi_target", "size"),
        planets=("selected_planet_count_qc", "sum"),
        verified_results=("source_result_state", lambda values: int((values == "verified").sum())),
        targets_with_lc=("lc_filtered_state", lambda values: int(values.str.startswith("present").sum())),
        targets_with_sc=("sc_filtered_state", lambda values: int(values.str.startswith("present").sum())),
        targets_with_parameters=("parameters_state", lambda values: int(values.str.startswith("present").sum())),
        targets_with_all_timings=("all_full_system_timing_present", "sum"),
        replay_ready=("replay_readiness", lambda values: int((values == "ready").sum())),
    ).reset_index()
    source_summary.to_csv(output_dir / "source_summary.csv", index=False)

    leverage_files = sorted(leverage_dir.glob("*_leverage.csv"))
    if len(leverage_files) != 4:
        raise ValueError("expected exactly four bounded leverage ledgers")
    leverage = pd.concat([pd.read_csv(path) for path in leverage_files], ignore_index=True)
    if set(leverage["kepoi_name"].astype(str)) != set(selected["kepoi_name"].astype(str)):
        raise ValueError("leverage ledgers do not match fixed selected membership")
    leverage = leverage.merge(inventory[["koi_target", "source_class", "replay_readiness"]], on="koi_target", validate="many_to_one")
    pilot_rows = []
    for (disk, system, source_class), group in leverage.groupby(["disk", "system", "source_class"], sort=True):
        ranked = group.sort_values(["log_contrast_fit_over_paper", "kepoi_name"], ascending=[False, True])
        high = ranked.iloc[0]
        controls = group[group["koi_target"] != high["koi_target"]].copy()
        controls["absolute_contrast"] = controls["log_contrast_fit_over_paper"].abs()
        control = controls.sort_values(["absolute_contrast", "kepoi_name"]).iloc[0] if len(controls) else None
        for role, row in (("canonical_leverage", high), ("low_contrast_control", control)):
            if row is None:
                continue
            pilot_rows.append({
                "disk": disk,
                "system": system,
                "source_class": source_class,
                "candidate_role": role,
                "koi_target": row["koi_target"],
                "priority_planet": row["kepoi_name"],
                "log_contrast_fit_over_paper": row["log_contrast_fit_over_paper"],
                "replay_readiness": row["replay_readiness"],
                "authorization_state": "candidate_not_authorized",
            })
    pd.DataFrame(pilot_rows).to_csv(output_dir / "pilot_candidates.csv", index=False)
    return inventory, source_summary, pd.DataFrame(pilot_rows)


def main():
    repo = Path(__file__).resolve().parents[1]
    research_root = os.environ.get("SAGEAR_RESEARCH_ROOT")
    research = Path(research_root).expanduser() if research_root else None
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=research is None,
                        default=None if research is None else research / "outputs/eccentricity_posterior_summary_SOURCE_FAITHFUL_RMS_PERIASTRON_GILBERT_QC_20260807.csv")
    parser.add_argument("--source-manifest", type=Path, required=research is None,
                        default=None if research is None else research / "tmp/exact_inventory_berger2020_fixed_dynesty_pairedperiod_rms_periastron_20260807/raw_result_manifest.csv")
    parser.add_argument("--leverage-dir", type=Path, default=repo / "metadata/inference_assumptions_20260915/full_priorities_20260916")
    parser.add_argument("--output-dir", type=Path, default=repo / "metadata/population_noise_audit_20260920")
    parser.add_argument("--companion-root", type=Path, action="append", default=[])
    args = parser.parse_args()
    roots = args.companion_root or ([] if research is None else [
        research / "alderaan_project/Data",
        research / "alderaan_project/Results/sagear_missing",
        research / "inputs/alderaan_archive_raw_20260612/ALDERAAN_posteriors",
        research / "tmp/published_inventory_recovery_postprocess_20260726_final/extracted/Results/sagear_published_inventory_missing",
    ])
    inventory, sources, pilots = build_inventory(
        args.summary,
        args.source_manifest,
        args.leverage_dir,
        args.output_dir,
        companion_roots=roots,
    )
    print(f"wrote {len(inventory)} targets, {len(sources)} source strata, {len(pilots)} pilot candidates")


if __name__ == "__main__":
    main()
