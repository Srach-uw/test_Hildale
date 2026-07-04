# Gap Hunt Findings — 2026-07-04

Audit of config-vs-code gaps, short-cadence scoping, and disk-classification methodology vs. the paper's prose (arXiv:2509.23973). Read-only investigation; **no script fixes applied** — proposals below are for a future retry/second pass. The live GCP run was not touched (verified the staged bundle's `run_one_target.sh` is byte-identical to `cloud/run_one_target.sh`, i.e. the fixed wget version).

## Findings table

| # | Finding | Severity / scope | Fixable? | Proposed fix |
|---|---------|------------------|----------|--------------|
| F1 | **Short cadence never used.** `config.json` `alderaan.include_short_cadence_when_available: true` is consumed by *nothing* (grep across both repos). All download paths hardcode `-c long` and `--use_sc` is never passed to `detrend_and_estimate_ttvs.py`. Paper §3.1 explicitly integrates models "over the Kepler long- and short-cadence exposure times." **Measured: 11/50 sampled missing-batch targets (22%) have SC files on MAST** (95% CI ≈ 12–36% → ~70–210 of the 592; median 15 SC files among SC targets). | Medium-high: ~130 of 592 retry targets; fidelity gap vs paper for the whole sample; plausible rescue for "Over 50% of transits flagged low quality" failures | Yes (retry pass) | See "F1 patch" below |
| F2 | **Bundle generator still emits the curl redirect bug.** `scripts/prepare_gcp_missing_alderaan_bundle.py:194` and `scripts/cloud_prepare.py:173` write `--cmdtype curl` into generated `run_one_target.sh`; `scripts/alderaan_batch.py:211` same for local downloads. `cloud/run_one_target.sh:35-39` documents why curl silently saves an HTML redirect page instead of FITS (http→https 301, no `-L`). The staged bundle was hand-fixed to wget, so the live run is safe — but **regenerating the bundle reintroduces corrupt downloads silently**. | High if bundle is ever regenerated; zero for live run | Yes, trivial | Change `--cmdtype curl` → `--cmdtype wget` at `prepare_gcp_missing_alderaan_bundle.py:194`, `cloud_prepare.py:173`, `alderaan_batch.py:211` |
| F3 | **Canonical classifier ≠ paper's described classifier — and the paper's literal method demonstrably cannot reproduce its own counts.** Canonical sample (`canonical_sample_old_astropy_rawcc.csv`) uses the *fallback* sample-trained unsupervised GMM (old-astropy cylindrical velocities), not the APOGEE-chemical-calibrated GMM the paper describes (§2.1). New tests today (see "F3 evidence"): every literal implementation of the paper's prose lands L1 = 633–943 vs paper category counts, while the canonical fallback sits at L1 = 69. | High conceptually, but **blocked**: paper prose under-determines their implementation | No (without author code) | Document as known divergence (this file). Do not tune further — consistent with `sagear_total_mismatch_audit.md` closed conclusion |
| F4 | **`classify_with_chemical_gmm` uses a Cartesian velocity proxy, paper specifies cylindrical.** `diagnose_sample.py:431-433` classifies on `V_phi`(=vy), `V_perp`(=√(vx²+vz²)); paper converts to Galactocentric cylindrical. `toomre_diagnostics_coordinate_summary.csv`: thick-star median V_perp is 78–81 km/s (Cartesian) vs 53–57 km/s (cylindrical) — material. All `classifier_threshold_diagnostics.py` chem variants (lines 94-105) also used the Cartesian columns. | Medium (only matters if chemical GMM path is ever promoted to canonical) | Yes | Parameterize `classify_with_chemical_gmm` with `x_col/y_col` like `classify_with_fallback_gmm` (diagnose_sample.py:378) and pass `V_phi_astropy/V_perp_astropy`. Note per F3 this still won't match paper counts |
| F5 | **Two sample cuts applied that the paper never states.** (a) Berger+2018 `Bin=0` resolved-companion cut (`diagnose_sample.py:126-139`, config `berger2018_bin_required: 0`) removes **200 planets** (audit stage 07b: 2970→2770). (b) Sample-level `Teff < 6500` on Berger Teff (stage 07a) removes 169; the paper states 6500 K only for the *APOGEE calibration* sample (its main sample is just "FGKM dwarfs", so (b) is a defensible inference). Also `diagnose_sample.py:108` drops RUWE=NaN rows; paper only says "filter out RUWE > 1.4". | Medium: 200 planets for (a); (b)/(c) minor | Decision needed, not a bug | Review whether 07b was added purely to chase counts. If yes, consider dropping it in a clearly-labeled sensitivity rerun; do not silently change the canonical sample mid-run |
| F6 | **Limb-darkening prior means differ in provenance (minor).** Vendored ALDERAAN applies N(µ, 0.1) on (u1,u2) (`alderaan/dynesty_helpers.py:229-232`, `sig_ld_sq=0.01`) — mechanics exactly match paper Table 1. But µ comes from KOI DR25 catalog coeffs (`alderaan_batch.py:181-186`), while the paper derives µ from Gaia + PHOENIX atmosphere models (Husser+2013, Parviainen & Aigrain 2015). Config fallback `default_limbdark_1/2` (0.3/0.2) triggered on **0 of 767** live catalog rows, so the config knob is effectively dead but harmless. | Low: both µ sources are model-derived for the same star; σ=0.1 dominates | Optional | If pursuing fidelity later: compute µ via `ldtk` with Berger+2020 Teff/logg/[Fe/H] in `build_alderaan_catalog`. Not worth blocking a retry pass |
| F7 | **Config internal contradiction:** `alderaan.cadence: "long"` (honored by `build_alderaan_needed_manifest.py:478`) coexists with `include_short_cadence_when_available: true` (honored by nothing). | Low (documentation hazard) | Yes | When applying F1, make the SC flag real; otherwise delete the key so config matches reality |

Verified-clean items (no gap): period 1–100 d, CONFIRMED+CANDIDATE FP removal, Furlan >5% contamination (keep-if-unreported semantics match paper wording), RUWE ≤ 1.4, high-α line [Mg/Fe] = −0.08·[Fe/H] + 0.14 (config exact match), APOGEE cal cuts Teff<6500/logg>4.0, `p_thick_threshold` 0.5 (paper states 0.5 explicitly; the 0.61 in `classifier_threshold_diagnostics_best.csv` is post-hoc tuning, not the paper's method), `posterior_resample_size`/e/ω grid sizes wired in `extract_eccentricity_posteriors.py:80-158`, null cuts (`impact_max`, `parallax_min_mas`) properly skipped with audit notes. **Also verified: the "use direct Gaia velocities where RV available" requirement (paper §2.1) is already satisfied** — `angus_velocities_parsed.csv` carries `has_rv` and its vx/vy/vz come from the Angus table's *measured* columns when RV exists (spot-checked 50/50 rows against the raw `angus_velocities.dat.gz` fixed-width measured fields).

## F1 evidence and patch (proposed — NOT applied)

Sampling: every 12th row of `targets_missing_launchable.csv` (n=50), MAST directory listing `https://archive.stsci.edu/missions/kepler/lightcurves/<kic4>/<kic9>/`, counting `_slc.fits`. Result: 11/50 with SC (22%); SC-rich targets are mostly low-numbered KOIs (K00108, K00142, K00179, K00281, K00370, K00708, K01175, K01860, K01929, K02414, K04637).

**Recommendation: worth fixing for the retry pass, but selectively.** SC is ~30× the data volume per quarter → materially slower detrend + fit. Best ROI: after the live run finishes, take the *failed* targets, check SC availability for exactly those (same directory-listing probe, minutes to run), and re-run only failed-AND-SC-available targets with the patch below. Do not blanket-enable SC for targets that already succeeded on long cadence unless a full-fidelity rerun is planned.

Patch to `cloud/run_one_target.sh` (then copy to bundle per established pattern):

```bash
# after line 40 (`bash "get_${KEPID}_lc.sh" || true`), inside the same pushd:
python "$ALDERAAN_REPO/bin/get_kepler_data.py" "$KEPID" -c short -t lightcurve -o "get_${KEPID}_sc.sh" --cmdtype wget
bash "get_${KEPID}_sc.sh" || true
```

and line 44 becomes:

```bash
python bin/detrend_and_estimate_ttvs.py --mission "$MISSION" --target "$TARGET" --run_id "$RUN_ID" --project_dir "$PROJECT_DIR" --data_dir "$DATA_DIR/" --catalog "$CATALOG_NAME" --use_sc True
```

Safety notes, verified in vendored source: with `--use_sc True` and zero SC files present, `detrend_and_estimate_ttvs.py:242-248` gets an empty list from `io.read_mast_files` and proceeds long-cadence-only, so the flag is safe unconditionally. **argparse gotcha:** `--use_sc` is `type=bool` (line 123-127), so `--use_sc False` evaluates *True* (`bool("False") is True`) — never pass the flag with a "False" value; omit it instead. Downstream stages need no changes: `fit_transit_shape_simultaneous_nested.py:198-204` auto-loads `<TARGET>_sc_filtered.fits` when it exists.

## F3 evidence (new tests, 2026-07-04)

Paper §2.1 prose: crossmatch with APOGEE DR17, split at the high/low-α line, plot both subsamples in (Vϕ, √(Vr²+Vz²)) **cylindrical** space, "fit a two-component GMM to the sample", apply to all hosts, thick if Pthick > 0.5. Implemented literally against the current 2474-planet sample (system labels held fixed at canonical rawcc definition):

| Variant (all at threshold 0.5) | thin_s | thick_s | thin_m | thick_m | L1 vs paper |
|---|---|---|---|---|---|
| Pooled unsupervised GMM, cylindrical (astropy) | 1326 | 55 | 1035 | 58 | 747 |
| Pooled unsupervised GMM, cylindrical (geometric) | 1325 | 56 | 1035 | 58 | 745 |
| Per-subsample Gaussians, empirical w_thick=0.06, cylindrical | 1365 | 16 | 1085 | 8 | 925 |
| Per-subsample Gaussians, w=0.5, cylindrical | 905 | 476 | 766 | 327 | 633 |
| **Canonical: sample-trained fallback GMM, old-astropy** | **1109** | **272** | **860** | **233** | **69** |
| Paper | 1121 | 275 | 862 | 207 | — |

The pooled GMM in cylindrical velocities latches onto a halo-like component (mean V_perp ≈ 132 km/s, only 26.7% high-α purity) instead of a thick-disk component; per-subsample Gaussians with the calibration's true 6% high-α weight classify nearly everything thin. No literal reading of the prose comes close. This extends the closed `sagear_total_mismatch_audit.md` conclusion with direct evidence: the classification difference is not a threshold, velocity-frame, or system-definition knob — the paper's actual implementation differs from its prose in some undocumented way (component weights, initialization, KIC-wide training set, or frame details). Exact reconciliation remains blocked; the canonical fallback-GMM sample is the defensible surrogate and its divergence is now quantified and documented.

## Synthesis of pre-existing diagnostics (Task 2, already on disk)

- `system_definition_counts.csv`: 12 classifier × single/multi-definition combos tested. Canonical (old_astropy + raw CONFIRMED/CANDIDATE count>1) is the best at L1=69; "after-all-cuts" definition moves 45 planets and is worse overall despite hitting thin_singles=1121 exactly under the direct classifier. Paper is silent on the definition — genuine ambiguity, already fully mapped.
- `classifier_threshold_diagnostics_best.csv`: tuned thresholds (0.61 direct; chem variants with priors 0.285–0.425) reach L1≈51–55 — but these are post-hoc fits, not the paper's stated 0.5, and were rightly not adopted.
- `classifier_disagreement_diagnostics_pairwise_summary.csv`: 47–640 planets flip between classifier variants; classifier choice dominates the +26 thick-multi excess.
- `toomre_classifier_grid_counts.csv`: KIC-background-trained GMMs (closest to the paper's "apply to the entire Kepler sample" language) are much worse (L1 273–419) than sample-trained.
- `sagear_total_mismatch_audit.md`: paper's own macros are internally inconsistent (2465 vs 2474); exact planet list unrecoverable. Closed; nothing here reopens it.

## Recommended action order for a future session

1. **F2** (3-line fix, prevents silent data corruption on any regeneration) — do first.
2. **F1** patch + failed-target SC probe → targeted SC retry batch after the live run completes.
3. **F7** config cleanup alongside F1.
4. **F5(a)** decide on the Bin=0 cut explicitly; record the decision either way.
5. F4/F6 only if a fidelity-focused rerun is ever commissioned; F3 is documentation-only.
