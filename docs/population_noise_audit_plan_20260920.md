# Population residual-noise audit plan

## Direction and boundary

Stop K02712 leave-out variants. Its saved-sample replay establishes fidelity to
one archived native likelihood, not adequacy of the photometric noise model.
This plan specifies the next local audit; it does not authorize a broad scanner,
cloud work, canonical edits, error rescaling, or population refits.

## Fixed membership gate

Use only the strict source-faithful posterior summary,
`eccentricity_posterior_summary_SOURCE_FAITHFUL_RMS_PERIASTRON_GILBERT_QC_20260807.csv`,
at SHA-256 `020008a34eec7993affb8bb67ff553e32b9e4c5de413800e79b152d7273b9bae`.
Supply its location explicitly with `--summary`, or set
`SAGEAR_RESEARCH_ROOT` for the local audit defaults.
It has 2,123 unique planets; exclude only the four rows already marked
`qc_primary_exclude=True` for `importance_ess_below_threshold`, leaving 2,119
planets in 1,583 targets. Do not convert target counts into planet counts or
treat planets in one joint target fit as independent photometric datasets.

The strict summary already follows the documented post-fit cuts: the 342
inventory omissions are those flagged by the grazing rule (nested weight above
`b > 1-Rp/Rstar` exceeds 5%) and/or approximate planet-radius fractional
uncertainty above 20%. Preserve unique `kepoi_name`, target/planet indexing,
the four population labels, `impact_mode=alderaan`, dynesty `LN_WT`, homogeneous
method provenance, and `include_transit_prior=False`. Abort on a membership,
hash, QC, period-index, or source-path discrepancy.

## Available provenance and unknowns

The paired-period extraction's `raw_result_manifest.csv` maps each staged
result FITS to a selected SHA-256 and source. Among the fixed 2,119 planets:

| Source root | Planets | Targets |
| --- | ---: | ---: |
| `inputs/alderaan_archive_raw_20260612/ALDERAAN_posteriors` | 1,590 | 1,164 |
| `alderaan_project/Results/sagear_missing` | 418 | 330 |
| `tmp/published_inventory_recovery_postprocess_20260726_final/extracted/Results/sagear_published_inventory_missing` | 111 | 89 |

The staged result FITS and posterior NPZ paths are available. Population-wide
availability and hashes of corresponding filtered light curves, errors,
quarters, transit catalogs and timing files remain **unknown**. Result FITS alone
cannot reconstruct native masks or residuals. Missing original light curves
therefore prevents a population-wide residual audit, source-specific likelihood
replay, and checks of quoted-error calibration or time correlation. Do not infer
light-curve availability from result-FITS availability.

## Residual specification

After parent review, first build a read-only target manifest, not the scanner.
For each of 1,583 targets record member planets and fit indices; original or
recovery source; result, light-curve, catalog and timing paths; SHA-256; cadence;
quarter coverage; and an explicit likelihood convention.

For targets with complete attested inputs, replay representative stored
high-weight samples and require saved-likelihood agreement before interpreting
likelihood-derived quantities. Verify separately for each source class whether
and where `LN_LIKE = -0.5 sum(r^2)` plus penalties or constants. Never treat
`-2*LN_LIKE` as chi-square merely because the column exists.

For each target and quarter, compute from normalized photometric residuals:

- point count, `sum(r^2)/N`, mean, median, RMS and robust MAD scale;
- fractions with `|r| > 3` and `|r| > 5`;
- lag-1 and short-lag correlations, with cadence and gap handling fixed;
- in-transit versus local out-of-transit summaries using the native mask.

These are target-level measurements. Join them to hosted planets only for
linkage, retaining a target key and clustered summaries. Stratify by archive
versus each recovery source, cadence, single/multi status and quarter coverage.
Link planet-level results to strict-summary density mismatch
(`rho_circular_to_catalog_ratio` and uncertainty fields) and the four
`full_priorities_20260916/*_leverage.csv` ledgers. Report target-clustered
associations; do not duplicate one residual series as independent evidence for
each planet.

## Decision gates

1. **Identity:** exact summary hash, 2,119 planets, 1,583 targets, four named
   exclusions, unique planet IDs and unchanged population membership.
2. **Source:** complete hashed photometric inputs and verified likelihood
   convention for every reported source stratum; otherwise report coverage and
   stop short of a population-wide claim.
3. **Pilot:** run a small predeclared source/population-stratified target sample
   before any broad scan. Freeze metrics and thresholds without reference to
   Table 3.
4. **Interpretation:** separate localized target/source failures from widespread
   residual inadequacy. Association with density mismatch or leverage is
   diagnostic, not causation or a new exclusion rule.
5. **Escalation:** if inadequacy is widespread, require a small end-to-end
   photometric injection/recovery experiment with matched cadence, masks,
   selection and known true geometry before population reinterpretation. The
   existing controls add Gaussian noise in log density and do not validate
   residual-noise handling. Do not rescale errors to chase Table 3.

Operational next step after parent approval: produce the 1,583-target source
availability manifest and a named pilot list only. No residual scanner or fit
launch belongs in that step.
