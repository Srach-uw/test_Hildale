# Hilldale / Sagear Replication Handoff Report

Generated: 2026-06-30

Latest operational update: 2026-07-01. The newest cloud-run status and decision log are appended at the bottom and supersede older intermediate count snapshots where they differ.

## Short Answer

We can make the pipeline much more faithful before running any new ALDERAAN jobs, but I would not call it "perfect" yet.

The pre-ALDERAAN sample construction is now close to Sagear. The largest old sample-size mismatch was real and was mostly fixed by adding:

- `Berger Teff < 6500 K`, matching the stated FGKM sample.
- Berger+2018 `Bin = 0`, likely excluding resolved-companion/binary systems.

After those cuts, the current canonical diagnostic sample is:

| quantity | current | Sagear target / macro | note |
|---|---:|---:|---|
| planets | 2474 | 2465 or 2474 | Sagear manuscript macros are internally inconsistent by 9 planets |
| hosts | 1857 | 1888 | still 31 low |
| thin singles | 1121 | 1121 | exact |
| thick singles | 305 | 275 | +30 |
| thin multi planets | 822 | 862 | -40 |
| thick multi planets | 226 | 207 | +19 |
| thick multi hosts | 95 | 98 | -3 |

So: sample construction is now strong; disk classification is close but not exact; eccentricity/population replication is not done until ALDERAAN posterior processing is done.

Important correction: ALDERAAN posterior FITS files already exist locally in:

`<alderaan-results-dir>`

This folder contains:

- 1692 `*-results.fits` files
- about 4.39 GB of outputs
- coverage of 1768 / 2474 planets in the corrected current sample
- 706 planets still missing

Coverage by current disk/system group:

| group | covered by existing FITS | missing |
|---|---:|---:|
| thin singles | 779 | 342 |
| thick singles | 163 | 142 |
| thin multis | 671 | 151 |
| thick multis | 155 | 71 |

This means we should process the existing ALDERAAN posteriors first. We only need new ALDERAAN runs for the missing systems or for a controlled validation/refit if the existing files are not Sagear-equivalent.

## What Changed Since the Recent File Review

### 1. Corrected Sagear Target Counts

The previous working assumption that Sagear had `378` thick singles was wrong. The paper/manuscript macros and the saved `thin_thick_rayleigh.png` indicate:

- thin singles: 1121
- thick singles: 275
- thin multi planets: 862
- thick multi planets: 207
- thin multi hosts: 394
- thick multi hosts: 98

The `378` number refers to a thick-host total macro/paragraph value, not thick-single planets.

The Sagear manuscript itself appears internally inconsistent:

| check | value A | value B | delta |
|---|---:|---:|---:|
| `allthinplanets + allthickplanets` vs `allplanets` | 2474 | 2465 | +9 |
| subgroup planet macros vs `allplanets` | 2465 | 2465 | 0 |
| thick subgroup planets vs `allthickplanets` | 482 | 491 | -9 |
| `allthinstars + allthickstars` vs `allstars` | 1893 | 1888 | +5 |
| subgroup host-like macros vs `allstars` | 1888 | 1888 | 0 |

Therefore a current sample of 2474 planets is not obviously wrong; it matches one set of Sagear's own macros.

### 2. Added Missing Sample Cuts

The current pipeline now applies:

- false-positive removal
- period `1-100 d`
- Furlan contamination `<= 5%`
- RUWE `<= 1.4`
- Berger stellar-density join
- `Berger Teff < 6500 K`
- Berger+2018 `Bin = 0`
- Angus velocity availability

The key attrition stages are:

| stage | planets | hosts |
|---|---:|---:|
| after Berger density join | 3139 | 2339 |
| after `Teff < 6500` | 2970 | 2208 |
| after Berger+2018 `Bin = 0` | 2770 | 2071 |
| after Angus velocities | 2474 | 1857 |

This is the biggest genuine improvement over the older 2305-planet sample.

### 3. Investigated the Toomre / Classifier Problem

Elena's comment that the Toomre diagram itself looked wrong was important.

The old zip contained a classification pipeline using an Astropy Galactocentric cylindrical velocity convention where disk rotation appears near:

`V_phi ~ -220 km/s`

The current diagnostic pipeline originally had a direct Angus-velocity convention where `vy` is near `+220 km/s`, and the plot needed a sign flip to resemble Sagear.

I ported the old coordinate convention into the current pipeline as:

- `V_phi_astropy`
- `V_perp_astropy`

But the old classifier is not the answer by itself. The old saved sample is much too small:

| sample | planets | hosts |
|---|---:|---:|
| old zip `sample.csv` | 1716 | 1278 |
| current corrected sample | 2474 | 1857 |

Old zip counts:

| group | old count |
|---|---:|
| thin singles | 875 |
| thick singles | 100 |
| thin multis | 671 |
| thick multis | 70 |

So the old zip is useful as a coordinate/sign reference, not as a sample or final classifier to adopt.

Classifier comparison on the corrected current sample:

| classifier | thin singles | thick singles | thin multi planets | thick multi planets | thick multi hosts |
|---|---:|---:|---:|---:|---:|
| direct Angus planet-host GMM | 1121 | 305 | 822 | 226 | 95 |
| geometric cylindrical planet-host GMM | 1148 | 278 | 819 | 229 | 96 |
| old-Astropy planet-host GMM | 1146 | 280 | 823 | 225 | 95 |
| old-pipeline KIC-wide Astropy GMM | 1031 | 395 | 742 | 306 | 127 |
| Sagear target | 1121 | 275 | 862 | 207 | 98 |

Interpretation:

- The old coordinate convention helps explain the Toomre visual mismatch.
- The old KIC-wide classifier overclassifies thick-disk planets.
- The closest current single-count match is old-Astropy/geometric planet-host GMM.
- The closest exact thin-single match is direct Angus planet-host GMM.
- No tested classifier reproduces all groups exactly.

### 4. Threshold Sweep Is Diagnostic, Not Proof

A threshold sweep showed that `P_thick > 0.61` can force the direct-GMM thick-single count to 275. But it also moves other groups:

| classifier | threshold | thin singles | thick singles | thin multi planets | thick multi planets | L1 planet delta |
|---|---:|---:|---:|---:|---:|---:|
| direct GMM | 0.610 | 1151 | 275 | 839 | 209 | 55 |
| geometric GMM | 0.525 | 1154 | 272 | 828 | 220 | 83 |
| chem 1-thin/1-thick | 0.490 | 1121 | 305 | 858 | 190 | 51 |
| chem 2-thin/2-thick | 0.460 | 1121 | 305 | 846 | 202 | 51 |

This is useful for diagnosis, but it is not a real replication unless Sagear's actual probability construction and threshold justify it.

### 5. Eccentricity Formula / Posterior Logic

The current `sagear_reproduction` code now follows the correct conceptual split:

- individual `e, omega` posteriors should not include the transit geometric prior by default;
- the hierarchical population likelihood should apply the reciprocal transit-probability correction:

`(1 - e^2) / (1 + e sin omega)`

This direction is mathematically correct. However, the earlier review is right that this does not by itself prove the final science result. The formula sanity checks are synthetic checks, not a reproduction of Sagear Table 2.

### 6. Current Eccentricity Diagnostics Are Not Final

The existing diagnostic eccentricity summary based on old `e_photo` point estimates is:

| population | n with e | mean e | median e | q95 e | Sagear Rayleigh `<e>` |
|---|---:|---:|---:|---:|---:|
| thin singles | 958 | 0.0323 | 0.0233 | 0.0855 | 0.022 |
| thick singles | 262 | 0.0351 | 0.0231 | 0.0891 | 0.066 |
| thin multis | 742 | 0.0357 | 0.0241 | 0.1001 | 0.030 |
| thick multis | 205 | 0.0362 | 0.0252 | 0.0961 | 0.033 |

These are not final replication values. They are triage only. The `e_photo` point-estimate workflow cannot replace Sagear's posterior-based hierarchical inference.

## Answer To: "Can We Do Everything Perfectly Until ALDERAAN Is Needed?"

Mostly, but with one important caveat.

### Yes, before ALDERAAN we can do these well

- Build the planet sample with clear attrition accounting.
- Apply the likely hidden FGKM and Berger+2018 binary cuts.
- Match the total planet count to within the manuscript's own ambiguity.
- Recreate diagnostic Toomre plots with correct sign/convention variants.
- Quantify exactly how many planets move under different classifiers.
- Prepare ALDERAAN target lists and identify missing systems.
- Process existing ALDERAAN FITS outputs into `e, omega` posterior grids.

### Not perfectly resolved before ALDERAAN

- Exact disk classifier remains ambiguous.
- Sagear's manuscript counts are internally inconsistent.
- The strict APOGEE chemical GMM currently underproduces thick-disk planets badly.
- The best classifier depends on whether the target is exact thin singles, thick singles, multis, or Toomre visual similarity.
- Final eccentricity distributions cannot be claimed until posterior processing and hierarchical fits are rerun.

So the correct statement is:

> We can make the pre-ALDERAAN sample and diagnostics very strong, but not mathematically identical to Sagear without knowing Sagear's exact disk-probability construction and convergence-removal table.

## Answer To: "Are ALDERAAN Results Already There?"

Yes. Existing ALDERAAN FITS results are already downloaded:

`<alderaan-results-dir>`

There are `1692` system FITS files. They cover `1768` of the `2474` planets in the corrected current sample.

Missing from existing FITS:

- 342 thin singles
- 142 thick singles
- 151 thin multi planets
- 71 thick multi planets

This is enough to start a serious posterior-based reproduction immediately. It is not enough for complete Sagear coverage.

## What The Earlier Review Got Right

The earlier review was right on these points:

- The two sample cuts are the most important practical discovery so far.
- The threshold sweep should not be merged as a final classifier.
- The formula sanity check is not a final validation.
- No new ALDERAAN fits have been run in the current pipeline.
- The thin-single discrepancy is not solved merely by sample cuts.

The earlier claim that there were no ALDERAAN posteriors available became outdated after checking the local machine. The files existed locally but were not wired cleanly into the canonical `sagear_reproduction` workflow.

## What Needed Revision In The Earlier Plan

The earlier plan suggested merging only the two cuts into the older `HILDALE_RESEARCH_SUMMER` pipeline and rerunning the old forward model.

I agree with the spirit, but I would add two steps:

1. Process the existing 1692 ALDERAAN FITS against the corrected sample before deciding to run new cloud jobs.
2. Do not trust the old classifier blindly. The old KIC-wide Astropy classifier overclassifies thick-disk systems relative to Sagear.

## Recommended Next Actions

### Immediate

1. Add a flat-results-dir option to `extract_eccentricity_posteriors.py`, or mirror the existing FITS files into the expected project layout.
2. Reprocess the existing 1692 ALDERAAN `*-results.fits` files against `canonical_sample_diagnostic.csv`.
3. Use period matching inside FITS systems, not just sorted planet index, because multi-planet order mismatches can silently corrupt posterior assignment.
4. Run `hierarchical_rayleigh.py` on the resulting joint posterior grids.
5. Compare Table 2-style Rayleigh values for covered-only populations.

### Then

6. Identify which missing 706 planets matter most for the thin-single discrepancy.
7. Prioritize missing/high-impact/high-e thin singles for new ALDERAAN runs.
8. Only after that, scale cloud ALDERAAN to the rest of the missing sample.

### For Merging

Merge into the main research pipeline:

- `Teff < 6500`
- Berger+2018 `Bin = 0`
- count-target correction: thick singles = 275, not 378
- Toomre sign/convention diagnostics
- existing ALDERAAN coverage audit

Do not merge as final science:

- `P_thick > 0.61` as a claimed Sagear classifier
- any final eccentricity claim from `e_photo` point estimates
- synthetic formula sanity checks as evidence of Table 2 replication

## Bottom Line

The current work did not fully solve the Sagear replication. It did materially improve the upstream sample and clarified that the Toomre/classifier issue is real.

The new best diagnosis is:

1. Sample construction is now close.
2. Disk classification is close but not exact.
3. Existing ALDERAAN posteriors are available for about 71.5% of the corrected sample.
4. The next scientific step is not to run all ALDERAAN from scratch; it is to process the existing FITS correctly, rerun hierarchical fits, and then decide what missing systems need new ALDERAAN.

## Update After Processing Existing ALDERAAN FITS

I implemented and ran the existing-posterior path after writing the first version of this report.

Code changes:

- `extract_eccentricity_posteriors.py` now accepts `--results-dir` for a flat folder of `*-results.fits`.
- Planet matching is now by ALDERAAN transit-time period versus KOI period, not by sorted row index.
- Berger stellar-density errors are now treated as log10 absolute errors in solar-density units, not as upper/lower log-density bounds.
- The default impact-parameter treatment is now `--impact-mode geometric`, marginalizing over `b ~ U(0, 1+Rp/R*)` because the ALDERAAN impact posterior is often prior-dominated.
- `hierarchical_rayleigh.py` was optimized by collapsing each 2-D `e,omega` posterior over omega once before the sigma sweep.
- Rayleigh fits now write separate outputs for transit-selection and no-selection diagnostics.

Processed files:

- ALDERAAN source: `<alderaan-results-dir>`
- Systems checked with FITS: 1278
- Systems missing FITS: 579
- Planet posterior grids extracted: 1729

Coverage after real extraction:

| group | sample planets | posterior extracted | no FITS | FITS but no matched/usable posterior |
|---|---:|---:|---:|---:|
| thin singles | 1121 | 779 | 342 | 0 |
| thick singles | 305 | 163 | 142 | 0 |
| thin multis | 822 | 637 | 151 | 34 |
| thick multis | 226 | 150 | 71 | 5 |

Coverage detail file:

`<repository-root>\outputs\alderaan_existing_coverage_detail.csv`

Rayleigh fit from existing extracted posteriors, with reciprocal transit-selection correction:

| population | N | current `<e>` | Sagear `<e>` |
|---|---:|---:|---:|
| thick singles | 163 | 0.109 | 0.066 |
| thin singles | 779 | 0.139 | 0.022 |
| thick multis | 150 | 0.006 | 0.033 |
| thin multis | 637 | 0.006 | 0.030 |

No-transit-selection diagnostic:

| population | N | `<e>` |
|---|---:|---:|
| thick singles | 163 | 0.135 |
| thin singles | 779 | 0.172 |
| thick multis | 150 | 0.006 |
| thin multis | 637 | 0.006 |

Interpretation:

- The existing ALDERAAN/Gilbert posteriors do not reproduce Sagear Table 2 under the current extraction.
- Singles come out much too eccentric, especially thin singles.
- Multis collapse to the low-e boundary in the Rayleigh fit.
- This strengthens the idea that Sagear's result depends on either her exact ALDERAAN reductions, exact posterior formalism, exact sample/classifier/convergence cuts, or some combination of these.
- It also means the next step is not blindly running 300 new systems. First, validate the extraction against a small set of known systems and compare `T14`, `Rp/Rs`, `b`, and duration ratio against Sagear/Gilbert expectations.

Immediate next diagnostic:

1. Pick 10 systems: low-e multis, high-e thin singles, thick singles, and systems with large zeta tails.
2. For each, compare:
   - KOI period and ALDERAAN period;
   - KOI duration and ALDERAAN `DUR14`;
   - catalog impact, ALDERAAN impact, and geometric-impact marginalized result;
   - Berger density prior;
   - resulting `zeta = T_obs / T_circ`.
3. Reproduce the old project's forward-model result on the corrected sample before deciding whether the hierarchical posterior extraction is the right final estimator.

## Update After Shape And Sensitivity Diagnostics

I added and ran:

- `alderaan_shape_diagnostics.py`
- `rayleigh_sensitivity_diagnostics.py`

Shape diagnostics output:

- `<repository-root>\outputs\alderaan_shape_diagnostics.csv`
- `<repository-root>\outputs\alderaan_shape_diagnostics_flagged.csv`
- `<repository-root>\outputs\alderaan_shape_diagnostics.png`
- `<repository-root>\outputs\alderaan_shape_diagnostics.md`

Main shape result:

| group | N | median zeta | q16 zeta | q84 zeta | median e50 | median ALDERAAN/KOI duration |
|---|---:|---:|---:|---:|---:|---:|
| thick multis | 150 | 1.068 | 0.761 | 1.172 | 0.218 | 0.994 |
| thick singles | 163 | 0.969 | 0.617 | 1.183 | 0.281 | 1.001 |
| thin multis | 637 | 1.050 | 0.755 | 1.164 | 0.216 | 0.995 |
| thin singles | 779 | 0.962 | 0.691 | 1.159 | 0.277 | 0.999 |

Interpretation:

- KOI and ALDERAAN durations agree at the group median, so the mismatch is not an obvious duration-unit or planet-matching bug.
- Stellar-density fractional errors are modest at the median, roughly 9-11%.
- High `e50` is strongly tied to `zeta` tails, not to broad density priors.
- Thin and thick singles have more short-zeta systems than multis.

Correlation of `e50` with `zeta_median` is strongly negative in all groups:

| group | corr(e50, zeta_median) |
|---|---:|
| thick multis | -0.902 |
| thick singles | -0.854 |
| thin multis | -0.871 |
| thin singles | -0.692 |

Fraction with `e50 > 0.5`:

| group | all planets | among zeta < 0.7 |
|---|---:|---:|
| thick multis | 0.113 | 0.895 |
| thick singles | 0.196 | 0.865 |
| thin multis | 0.110 | 0.908 |
| thin singles | 0.154 | 0.892 |

Rayleigh sensitivity:

| filter | thick singles | thin singles | thick multis | thin multis |
|---|---:|---:|---:|---:|
| all | 0.109 | 0.139 | 0.006 | 0.006 |
| zeta >= 0.7 | 0.103 | 0.134 | 0.006 | 0.006 |
| 0.7 <= zeta <= 1.3 | 0.006 | 0.006 | 0.006 | 0.006 |
| no duration shift >25% | 0.107 | 0.137 | 0.006 | 0.006 |
| ALDERAAN b < 0.85 | 0.103 | 0.139 | 0.006 | 0.006 |
| SNR >= 20 | 0.113 | 0.122 | 0.006 | 0.006 |
| strict clean | 0.006 | 0.006 | 0.006 | 0.006 |

This says the high-e signal in the existing posteriors is driven by duration-ratio tails. Removing only short-zeta systems does not erase the singles signal because long-zeta tails can also require eccentricity. Removing both tails forces all populations to the circular boundary.

Best current diagnosis:

> The remaining discrepancy is not basic sample construction, not a gross ALDERAAN duration mismatch, and not simply broad Berger density errors. It is specifically about how the duration-ratio tails enter the posterior/population inference, and whether Sagear's own ALDERAAN reductions or convergence/quality cuts suppress those tails relative to the Gilbert FITS files on disk.

Concrete next validation set:

Use `alderaan_shape_diagnostics_flagged.csv` and inspect/refit the top flagged systems, especially:

- high-e thin singles: `K03232.01`, `K02519.01`, `K06878.01`, `K01793.01`, `K03104.01`;
- high-e thin multis: `K00352.02`, `K00780.02`, `K03384.01`, `K01955.04`;
- high-e thick singles: `K03358.01`, `K07638.01`, `K02617.01`, `K02426.01`;
- high-e thick multis: `K00282.03`, `K03437.02`, `K01240.02`.

For each one, compare the light curve fit visually if possible. The question is no longer "do we have ALDERAAN files?" It is "are these specific duration-ratio tails present in Sagear's ALDERAAN run, and did she cut them?"

## Update After Refit Manifest And Dossiers

I added and ran:

- `alderaan_refit_plan.py`
- `alderaan_target_dossiers.py`

These turn the vague next step into concrete run lists and visual dossiers.

New outputs:

- `<repository-root>\outputs\alderaan_refit_full_manifest.csv`
- `<repository-root>\outputs\alderaan_refit_validation_targets.csv`
- `<repository-root>\outputs\alderaan_missing_targets_all.csv`
- `<repository-root>\outputs\alderaan_refit_plan.md`
- `<repository-root>\outputs\alderaan_target_dossiers\index.md`
- `<repository-root>\outputs\alderaan_target_dossiers\*.png`
- `<repository-root>\outputs\alderaan_shape_pathology_rates.csv`

ALDERAAN project files updated:

- `<repository-root>\alderaan_project\sagear_refit_validation_targets.csv`
- `<repository-root>\alderaan_project\Catalogs\sagear_refit_validation_catalog.csv`
- `<repository-root>\alderaan_project\Catalogs\sagear_validation_catalog.csv`
- `<repository-root>\alderaan_project\Scripts\download_validation_lightcurves.ps1`
- `<repository-root>\alderaan_project\Scripts\run_validation_alderaan.ps1`

The validation target list contains:

| reason | systems |
|---|---:|
| suspicious existing ALDERAAN duration-ratio tail refits | 40 |
| missing ALDERAAN high-SNR validation systems | 32 |
| clean existing ALDERAAN controls | 20 |

Full missing ALDERAAN list:

| group | missing systems | missing planets |
|---|---:|---:|
| thick multis | 30 | 71 |
| thick singles | 142 | 142 |
| thin multis | 65 | 151 |
| thin singles | 342 | 342 |

One error spotted and fixed:

- The refit plan initially wrote `sagear_refit_validation_catalog.csv`, but the existing ALDERAAN run-script generator expected `sagear_validation_catalog.csv`.
- I now write both names so the generated `run_validation_alderaan.ps1` is usable without manual editing.

Another important environment finding:

- The ALDERAAN source repository is represented here as `<alderaan-source-dir>`.
- The active Anaconda environment does **not** have the required ALDERAAN dependencies:
  - missing `pymc3`
  - missing `exoplanet`
  - missing `dynesty`
  - missing `batman`
  - missing `celerite2`
  - missing `ldtk`
  - missing `arviz`

So the true remaining blocker is now very concrete:

> Create/activate the ALDERAAN conda environment, then run the generated validation/refit scripts.

Additional pathology-rate finding:

| group | ROR shift >25% | duration shift >25% | ALDERAAN b > 1 | e50 > 0.5 |
|---|---:|---:|---:|---:|
| thick multis | 0.060 | 0.040 | 0.000 | 0.113 |
| thick singles | 0.043 | 0.049 | 0.006 | 0.196 |
| thin multis | 0.047 | 0.044 | 0.003 | 0.110 |
| thin singles | 0.042 | 0.030 | 0.006 | 0.154 |

The top dossier, `K03232.01`, is especially instructive:

- ALDERAAN duration is shifted relative to KOI.
- ALDERAAN `Rp/Rs` is wildly larger than KOI.
- ALDERAAN impact posterior is partly beyond 1.
- The derived eccentricity posterior peaks high.

This is not a tiny modeling nuance. It looks like the kind of case where Sagear's convergence/visual-quality cuts could matter a lot.

Current pre-ALDERAAN status:

- Sample audit: done and close to Sagear.
- Disk-count diagnosis: done; exact classifier still ambiguous but quantified.
- Existing ALDERAAN coverage audit: done.
- Existing posterior extraction: done.
- Population fit on existing posteriors: done, not Sagear.
- Shape-tail diagnostics: done.
- Sensitivity cuts: done.
- Missing/refit/control manifests: done.
- Visual dossiers for top suspicious systems: done.
- Generated ALDERAAN validation catalog and scripts: done.
- Remaining work: ALDERAAN environment + actual ALDERAAN reruns/refits.

## Pre-ALDERAAN Status Table

Asked: are we good on everything before ALDERAAN, and how many planets/eccentricities are in each category?

Short answer:

- sample construction is close and usable;
- disk classification is close but not exact;
- existing ALDERAAN/Gilbert posteriors have been processed, but they do not reproduce Sagear;
- final remaining blocker is ALDERAAN validation/refits and Sagear-equivalent convergence/quality cuts.

| population | Sagear planets | current corrected sample | current hosts | existing ALDERAAN posteriors | missing/no posterior | old point-estimate median e | old point-estimate mean e | current ALDERAAN Rayleigh `<e>` | Sagear Rayleigh `<e>` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| thin singles | 1121 | 1121 | 1121 | 779 | 342 | 0.023 | 0.032 | 0.139 | 0.022 |
| thick singles | 275 | 305 | 305 | 163 | 142 | 0.023 | 0.035 | 0.109 | 0.066 |
| thin multis | 862 | 822 | 336 | 637 | 185 | 0.024 | 0.036 | 0.006 | 0.030 |
| thick multis | 207 | 226 | 95 | 150 | 76 | 0.025 | 0.036 | 0.006 | 0.033 |

Interpretation:

- Thin-single count matches Sagear exactly.
- Thick-single and multi counts are still offset by tens of planets, mostly from the unresolved disk-classifier/convergence details.
- The old `e_photo` medians look deceptively close to Sagear, but those are not the posterior-based measurement.
- The currently processed existing ALDERAAN posteriors produce high singles and near-circular multis, so the existing FITS plus current extraction are not a Sagear reproduction.
- This is why the prepared ALDERAAN refit/validation batch matters.

## Toomre/Classifer Recheck, June 30

I regenerated Sagear-style Toomre diagnostics instead of relying on the older generic scatter plots.

New outputs:

| file | purpose |
|---|---|
| `<repository-root>\outputs\toomre_sagear_style_direct.png` | direct Angus velocity Toomre plot with `P_thick` background |
| `<repository-root>\outputs\toomre_sagear_style_old_astropy.png` | old Hilldale/Astropy cylindrical convention plot |
| `<repository-root>\outputs\toomre_sagear_reference_comparison.png` | side-by-side with Sagear reference |
| `<repository-root>\outputs\toomre_classifier_grid_with_reference.png` | reference plus 4 classifier convention panels |
| `<repository-root>\outputs\toomre_classifier_grid_counts.csv` | count comparison for those classifier variants |

Important error caught:

- My first classifier-grid run filtered direct positive Angus `vy` as if it were already negative plotted `V_phi`.
- That trained the KIC-wide direct panel on the wrong velocity tail.
- I fixed the filter so it applies after the display-sign convention.

Classifier-count result:

| variant | thin singles | thick singles | thin multis | thick multis | L1 planet-count error |
|---|---:|---:|---:|---:|---:|
| direct Angus, planet-sample GMM | 1121 | 305 | 822 | 226 | 89 |
| old Astropy, planet-sample GMM | 1146 | 280 | 823 | 225 | 87 |
| direct Angus, KIC-wide GMM | 1074 | 352 | 777 | 271 | 273 |
| old Astropy, KIC-wide GMM | 1034 | 392 | 744 | 304 | 419 |

Interpretation:

- The professor's Toomre concern is valid: the direct-Angus and old-Astropy conventions produce visibly different probability fields.
- The old-Astropy planet-sample GMM is marginally closer in total count error and much closer for thick singles.
- But it does not solve the replication: thin/multi counts remain off, and no classifier variant reproduces all four Sagear groups.
- The old KIC-wide classifier from the old project makes the Toomre figure look plausible but badly overclassifies thick-disk planets in this corrected sample.

I made this classifier choice explicit in `diagnose_sample.py`:

```powershell
python diagnose_sample.py --force-fallback-gmm --fallback-velocity direct
python diagnose_sample.py --force-fallback-gmm --fallback-velocity old_astropy
```

The old-Astropy diagnostic sample is now:

`<repository-root>\outputs\canonical_sample_diagnostic_old_astropy.csv`

## Existing-ALDERAAN Sensitivity To Classifier Choice

I relabeled the already-extracted existing ALDERAAN posterior summary with the old-Astropy disk labels and reran the Rayleigh population fit.

Output:

`<repository-root>\outputs\pre_alderaan_classifier_eccentricity_comparison.csv`

| population | Sagear planets | direct planets | direct ALDERAAN posteriors | direct Rayleigh `<e>` | old-Astropy planets | old-Astropy ALDERAAN posteriors | old-Astropy Rayleigh `<e>` | Sagear Rayleigh `<e>` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| thin singles | 1121 | 1121 | 779 | 0.139 | 1146 | 788 | 0.139 | 0.022 |
| thick singles | 275 | 305 | 163 | 0.109 | 280 | 154 | 0.112 | 0.066 |
| thin multis | 862 | 822 | 637 | 0.006 | 823 | 632 | 0.006 | 0.030 |
| thick multis | 207 | 226 | 150 | 0.006 | 225 | 155 | 0.006 | 0.033 |

Interpretation:

- Changing to the Toomre-better old-Astropy labels barely changes the existing-ALDERAAN eccentricity result.
- Therefore the current mismatch is not fixed by disk-label convention alone.
- The strongest remaining non-ALDERAAN issue is exact Sagear convergence/quality-cut knowledge, especially for the multi sample and the 9-planet inconsistency in Sagear's reported totals.
- The strongest ALDERAAN issue remains posterior coverage and suspicious shape-tail cases. Existing local FITS are not enough to reproduce Sagear.

Updated bottom line:

- Pre-ALDERAAN sample/cut work is now close and quantified.
- Toomre/classifier ambiguity is real but bounded.
- Classifier changes do not rescue the eccentricity mismatch.
- The next scientifically clean step is still the validation/refit ALDERAAN batch, plus a convergence/visual-quality cut table.

## Multiplicity/System-Definition Diagnosis, June 30

I found a major upstream issue after the Toomre work:

> The single/multi label should probably not be assigned only after all final planet cuts.

The previous pipeline labeled a system as "single" if only one planet survived all cuts. But Sagear's counts are better matched if a system is labeled "multi" when the raw KOI confirmed/candidate system has more than one planet.

New diagnostic script:

`<repository-root>\system_definition_diagnostics.py`

New outputs:

| file | purpose |
|---|---|
| `<repository-root>\outputs\system_definition_counts.csv` | count comparison for different single/multi definitions |
| `<repository-root>\outputs\system_definition_moved_planets.csv` | planets moved between single/multi labels |
| `<repository-root>\outputs\system_definition_moved_planet_eccentricity_summary.csv` | eccentricity summary of moved planets |

Best count variant so far:

| classifier + system definition | thin singles | thick singles | thin multis | thick multis | thin multi hosts | thick multi hosts | L1 error |
|---|---:|---:|---:|---:|---:|---:|---:|
| old-Astropy disk + raw confirmed/candidate multiplicity | 1109 | 272 | 860 | 233 | 373 | 103 | 69 |
| direct disk + raw confirmed/candidate multiplicity | 1087 | 294 | 856 | 237 | 370 | 106 | 121 |
| old-Astropy disk + after-all-cuts multiplicity | 1146 | 280 | 823 | 225 | 336 | 95 | 148 |
| direct disk + after-all-cuts multiplicity | 1121 | 305 | 822 | 226 | 336 | 95 | 150 |

This is the strongest pre-ALDERAAN improvement so far.

Why it matters:

- In the best-label variant, thin multis move from 823 to 860 planets, close to Sagear's 862.
- Thick singles move from 280 to 272 planets, close to Sagear's 275.
- Thick multi hosts move from 95 to 103, close to Sagear's 98.
- Thin multi hosts remain low: 373 vs Sagear's 394.
- Thick multi planets remain high: 233 vs Sagear's 207.

The moved planets are not random:

| variant | disk | moved from | moved to | planets | with posteriors | median e50 |
|---|---|---|---|---:|---:|---:|
| direct rawcc | thin | single | multi | 34 | 30 | 0.232 |
| direct rawcc | thick | single | multi | 11 | 7 | 0.213 |
| old-Astropy rawcc | thin | single | multi | 37 | 32 | 0.245 |
| old-Astropy rawcc | thick | single | multi | 8 | 5 | 0.171 |

Interpretation:

- The after-all-cuts multiplicity definition was hiding a high-e tail inside the single-planet sample.
- Moving those systems back to "multi" makes thin multis agree much better with Sagear.
- This does not fix thin singles; they remain too eccentric.

## Best Current Pre-ALDERAAN Checkpoint

I created a checkpoint script:

`<repository-root>\pre_alderaan_checkpoint.py`

It uses the current best pre-ALDERAAN labels by default:

- old-Astropy disk classifier;
- raw confirmed/candidate KOI multiplicity for single/multi labels.

New outputs:

| file | purpose |
|---|---|
| `<repository-root>\outputs\canonical_sample_old_astropy_rawcc.csv` | best current sample-label table |
| `<repository-root>\outputs\eccentricity_posterior_summary_old_astropy_rawcc.csv` | existing posteriors relabeled with best labels |
| `<repository-root>\outputs\pre_alderaan_best_population_table.csv` | best current population/eccentricity table |
| `<repository-root>\outputs\pre_alderaan_best_missing_posteriors.csv` | missing posterior planets under best labels |
| `<repository-root>\outputs\pre_alderaan_best_top_outliers.csv` | highest-e / most suspicious existing posterior rows |
| `<repository-root>\outputs\pre_alderaan_best_eccentricity_diagnostics.png` | e50 and zeta diagnostic plot |
| `<repository-root>\outputs\pre_alderaan_best_summary.md` | compact checkpoint report |

Best current table:

| population | sample planets | sample hosts | posterior planets | missing posteriors | coverage | median e50 | Rayleigh `<e>` | Sagear `<e>` | delta |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| thin singles | 1109 | 1109 | 756 | 353 | 0.682 | 0.281 | 0.140 | 0.022 | +0.118 |
| thick singles | 272 | 272 | 149 | 123 | 0.548 | 0.274 | 0.117 | 0.066 | +0.051 |
| thin multis | 860 | 373 | 664 | 196 | 0.772 | 0.219 | 0.041 | 0.030 | +0.011 |
| thick multis | 233 | 103 | 160 | 73 | 0.687 | 0.213 | 0.006 | 0.033 | -0.027 |

This is the most honest current status:

- Thin multis are now close to Sagear.
- Thick singles have good counts but eccentricity remains high.
- Thin singles are still badly high.
- Thick multis have decent counts but the Rayleigh fit collapses to the lower bound.

## Quality/Convergence Proxy Diagnosis

I reran the shape/outlier sensitivity under the best labels.

Output:

`<repository-root>\outputs\rayleigh_sensitivity_diagnostics_old_astropy_rawcc.csv`

Important results:

| filter | thick singles `<e>` | thin singles `<e>` | thick multis `<e>` | thin multis `<e>` |
|---|---:|---:|---:|---:|
| all | 0.116 | 0.140 | 0.006 | 0.041 |
| zeta >= 0.7 | 0.111 | 0.134 | 0.006 | 0.041 |
| no duration shift >25% | 0.111 | 0.140 | 0.006 | 0.039 |
| ALDERAAN b < 0.85 | 0.114 | 0.140 | 0.006 | 0.041 |
| SNR >= 20 | 0.122 | 0.122 | 0.006 | 0.039 |
| strict clean | 0.006 | 0.006 | 0.006 | 0.006 |

Interpretation:

- Simple quality cuts do not reproduce Sagear.
- Removing all zeta tails makes every population circular, which is not a valid fix.
- Thin singles remain high even after duration-shift, impact, and zeta >= 0.7 cuts.
- Thick multis remain at the lower bound under every simple proxy.
- Therefore the remaining issue is not one obvious catalog mismatch or one obvious filter.

I also checked no-transit-selection under the best labels:

| population | no-selection Rayleigh `<e>` |
|---|---:|
| thick singles | 0.142 |
| thin singles | 0.175 |
| thick multis | 0.006 |
| thin multis | 0.045 |

So the transit-selection correction is not causing the remaining mismatch; removing it makes singles worse.

Updated diagnosis:

1. **Fixed/strongly improved:** single/multi definition.
2. **Mostly bounded:** disk classifier convention.
3. **Still unresolved:** exact Sagear convergence/quality cuts.
4. **Still unresolved:** existing ALDERAAN posterior coverage and posterior-shape behavior.
5. **Most likely next necessary step:** run the prepared ALDERAAN validation/refit batch, then compare fitted transit-shape posteriors directly to the existing local FITS.

## Total Planet Count Mismatch Audit, July 1

Question asked:

> Why is our total 2474 planets while Sagear says 2465? Can we identify the planets?

I wrote a targeted audit:

`<repository-root>\sagear_total_mismatch_audit.py`

New outputs:

| file | purpose |
|---|---|
| `<repository-root>\outputs\sagear_total_mismatch_macro_audit.csv` | compares our total to Sagear total macros |
| `<repository-root>\outputs\sagear_total_mismatch_category_audit.csv` | compares every subgroup count |
| `<repository-root>\outputs\sagear_total_mismatch_implied_removed.csv` | calculates the implied hidden bucket |
| `<repository-root>\outputs\sagear_total_mismatch_candidate_hosts.csv` | current thick hosts ranked for inspection |
| `<repository-root>\outputs\sagear_total_mismatch_candidate_planets.csv` | current thick planets ranked for inspection |
| `<repository-root>\outputs\sagear_total_mismatch_audit.md` | compact explanation |

Key result:

| quantity | value | delta vs our 2474 |
|---|---:|---:|
| our current best sample rows | 2474 | 0 |
| Sagear `allplanets` macro | 2465 | +9 |
| Sagear `allthinplanets + allthickplanets` | 2474 | 0 |
| Sagear subgroup total: thin/thick singles + multis | 2465 | +9 |

This means our total of 2474 is not arbitrary. It exactly equals one set of Sagear macros:

```text
allthinplanets + allthickplanets = 1983 + 491 = 2474
```

But another set of Sagear macros says:

```text
allplanets = 2465
allthinsingles + allthicksingles + allthinmultiplanets + allthickmultiplanets
= 1121 + 275 + 862 + 207 = 2465
```

So the paper/source has an internal count inconsistency.

The inconsistency is specifically in the thick-disk aggregate:

| implied bucket | count | calculation |
|---|---:|---|
| thick planets in disk total but not final subgroups | 9 | `allthickplanets - allthicksingles - allthickmultiplanets = 491 - 275 - 207` |
| thick hosts in disk total but not final subgroups | 5 | `allthickstars - allthicksingles - allthickmultistars = 378 - 275 - 98` |

Interpretation:

- Sagear likely had an intermediate sample with 2474 planets.
- Then nine thick-disk planets from five thick hosts were removed from the final subgroup analysis.
- The most plausible reason is the later ALDERAAN visual/convergence removal mentioned in the methods.
- But the LaTeX source archive contains no final planet list and no convergence/removal table.
- Therefore the exact nine planet IDs cannot be recovered from the local Sagear source files alone.

Important: deleting nine planets from our current sample cannot make the subgroup counts match by itself.

Current best subgroup deltas:

| metric | ours | Sagear | delta |
|---|---:|---:|---:|
| thin singles | 1109 | 1121 | -12 |
| thick singles | 272 | 275 | -3 |
| thin multis | 860 | 862 | -2 |
| thick multis | 233 | 207 | +26 |
| thin total planets | 1969 | 1983 | -14 |
| thick total planets | 505 | 491 | +14 |
| thin multi hosts | 373 | 394 | -21 |
| thick multi hosts | 103 | 98 | +5 |

So there are two separate things:

1. **Total mismatch:** likely a Sagear macro/update inconsistency around nine removed thick planets.
2. **Composition mismatch:** our disk labels still put too many planets into thick multis and too few into thin categories.

Candidate-host inspection:

The audit ranks current thick-disk hosts by missing posterior / low SNR / shape-pathology proxies. These are not confirmed Sagear removals, only candidates to inspect if we are looking for hidden convergence removals.

Top current candidate hosts include:

| kepid | candidate planets | system | planet names | missing posterior planets |
|---:|---:|---|---|---:|
| 9458613 | 5 | multi | K00707.01,K00707.02,K00707.03,K00707.04,K00707.05 | 5 |
| 5511081 | 4 | multi | K01930.01,K01930.02,K01930.03,K01930.04 | 4 |
| 6359320 | 3 | multi | K01127.02,K01127.03,K01127.01 | 3 |
| 11450414 | 3 | multi | K01992.02,K01992.03,K01992.01 | 3 |
| 5965819 | 3 | multi | K03319.03,K03319.01,K03319.02 | 3 |

Bottom line:

- I cannot honestly name "the nine Sagear planets" from the files we have, because the source archive does not contain the final table.
- I can say the nine are very likely a paper-internal omitted/removed thick-disk bucket.
- To recover them exactly, we need either Sagear's final sample/convergence table or to rerun ALDERAAN and create our own convergence/visual-quality removal table.

## Final Pre-ALDERAAN Hardening Pass, July 1

Goal:

Create a defensible planet-level list of everything that still needs ALDERAAN, plus a launch-ready KOI-system list that will not silently drop or crash on bad catalog inputs.

New scripts/edits:

| file | purpose |
|---|---|
| `<repository-root>\build_alderaan_needed_manifest.py` | builds full needed/runnable/unseeded planet manifests, runnable target manifest, ALDERAAN catalog, and runner scripts |
| `<repository-root>\validate_alderaan_needed_manifest.py` | validates duplicates, target coverage, required catalog fields, and catalog/target consistency |
| `<repository-root>\alderaan_batch.py` | relaxed the catalog seed requirement so missing impact no longer blocks a run; period, epoch, depth, and duration still must be valid |

Important correction:

I found that 27 planets in the "needs ALDERAAN" list have no usable KOI depth. ALDERAAN reads `catalog.depth` directly in `detrend_and_estimate_ttvs.py`, `analyze_autocorrelated_noise.py`, and `fit_transit_shape_simultaneous_nested.py`, so those planets cannot be honestly launched without recovering or manually estimating a transit-depth seed.

I checked local fallback catalogs in the ALDERAAN repo and older project zip. A few historical rows have conflicting old TCE fit values, but not enough to recover the set cleanly, and some durations differ strongly from the current cumulative catalog. I therefore did not silently fill them. They are preserved in a separate blocked list.

Final counts:

| quantity | count |
|---|---:|
| current best sample planets | 2474 |
| existing extracted ALDERAAN posteriors | 1729 |
| planets that scientifically need ALDERAAN/refit | 1065 |
| launch-ready needed planets | 1038 |
| blocked/unseeded needed planets | 27 |
| all KOI systems with at least one needed planet | 913 |
| runnable KOI systems for batch execution | 889 |
| ALDERAAN catalog rows after system expansion | 1225 |

Needed planets by current best population:

| population | all needed | runnable now | blocked missing seed |
|---|---:|---:|---:|
| thick multis | 100 | 98 | 2 |
| thick singles | 156 | 153 | 3 |
| thin multis | 302 | 299 | 3 |
| thin singles | 507 | 488 | 19 |

Needed planets by action:

| action | all needed | runnable now |
|---|---:|---:|
| missing ALDERAAN results file | 706 | 682 |
| bad existing results-file match / bad zeta extraction | 39 | 36 |
| hard shape pathology refit | 109 | 109 |
| high-e / extreme-zeta tail refit | 204 | 204 |
| low-SNR tail stress refit | 7 | 7 |

Runnable target systems by priority:

| priority | target type | runnable targets |
|---:|---|---:|
| 1 | required missing/bad extraction | 593 |
| 2 | required hard-shape refit | 99 |
| 3 | recommended tail refit | 191 |
| 4 | recommended low-priority tail refit | 6 |

Validation result:

`validate_alderaan_needed_manifest.py` passes all checks:

- no duplicate planet IDs in the sample/status/needed tables;
- every runnable needed planet maps to a runnable KOI target;
- every runnable needed planet appears in the ALDERAAN catalog;
- every target in the execution CSV is justified by at least one runnable needed planet;
- required ALDERAAN catalog fields are finite/positive;
- catalog `npl` matches rows per target.

Runner improvements:

The generated project scripts are now compact and parameterized rather than thousands of hard-coded lines:

| script | use |
|---|---|
| `<repository-root>\alderaan_project\Scripts\download_needed_lightcurves.ps1` | Windows/PowerShell download helper |
| `<repository-root>\alderaan_project\Scripts\run_needed_alderaan.ps1` | Windows/PowerShell serial runner with `MaxPriority`, `Start`, `Limit`, and phase selection |
| `<repository-root>\alderaan_project\Scripts\download_needed_lightcurves.sh` | Linux/GCP download helper |
| `<repository-root>\alderaan_project\Scripts\run_needed_alderaan_parallel.sh` | Linux/GCP parallel runner using `JOBS`, `MAX_PRIORITY`, `START`, `LIMIT`, and `PHASE` |

I also checked the ALDERAAN nested-sampling script. It contains `ncores = cpu_count() - 2`, but `USE_MULTIPRO = False`, so the fit is currently serial per target. I still added `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, and `NUMEXPR_NUM_THREADS=1` to the generated runners to prevent hidden BLAS/Aesara oversubscription when running many targets in parallel.

Deliverables copied into the project output folder:

| file | meaning |
|---|---|
| `<user-home>\Documents\project\2026-06-29\8\outputs\alderaan_needed_planets_best.csv` | full planet-level scientific need list |
| `<user-home>\Documents\project\2026-06-29\8\outputs\alderaan_runnable_needed_planets_best.csv` | launch-ready planet-level subset |
| `<user-home>\Documents\project\2026-06-29\8\outputs\alderaan_unseeded_needed_planets_best.csv` | needed planets blocked by missing transit seeds |
| `<user-home>\Documents\project\2026-06-29\8\outputs\alderaan_needed_targets_best.csv` | runnable KOI-system execution list |
| `<user-home>\Documents\project\2026-06-29\8\outputs\alderaan_needed_targets_all_best.csv` | all target systems, including blocked missing-seed targets |
| `<user-home>\Documents\project\2026-06-29\8\outputs\alderaan_needed_catalog_best.csv` | ALDERAAN input catalog for runnable targets |
| `<user-home>\Documents\project\2026-06-29\8\outputs\alderaan_all_planets_status_best.csv` | full 2474-row status/flags table |
| `<user-home>\Documents\project\2026-06-29\8\outputs\alderaan_needed_manifest.md` | readable manifest summary |
| `<user-home>\Documents\project\2026-06-29\8\outputs\alderaan_needed_validation.md` | validation report |

Current best next step:

Run ALDERAAN first on the priority-1 runnable target set, not the optional tail-refit set. That is 593 KOI systems covering 718 launch-ready tier-1 planet rows. The 27 missing-depth planets should be handled separately through manual transit-seed recovery or excluded with an explicit documented cut.

### Historical Transit-Seed Recovery Audit

I added one more audit:

`<repository-root>\audit_unseeded_transit_recovery.py`

Outputs:

| file | purpose |
|---|---|
| `<user-home>\Documents\project\2026-06-29\8\outputs\alderaan_unseeded_historical_seed_audit.csv` | row-level check of local current/historical KOI/TCE catalogs for the 27 blocked planets |
| `<user-home>\Documents\project\2026-06-29\8\outputs\alderaan_unseeded_historical_seed_audit.md` | readable summary |

Result:

| catalog | matched rows | valid historical depth rows | valid depth + duration-consistent rows |
|---|---:|---:|---:|
| current cumulative 2026 | 27 | 0 | 0 |
| ALDERAAN cumulative 2024 | 27 | 0 | 0 |
| DR22 Mullally Q1-Q16 | 16 | 1 | 0 |
| DR24 Coughlin Q1-Q17 | 11 | 3 | 0 |
| DR25 Thompson Q1-Q17 | 4 | 0 | 0 |
| merged planets full | 27 | 0 | 0 |

Interpretation:

Some old TCE rows contain depth values for a few blocked planets, but none pass the simple automatic-fill rule requiring both a positive historical depth and duration agreement within 25 percent of the current row. Therefore I did not autofill them into the ALDERAAN catalog. They should remain a manual seed-recovery/exclusion subset.

## Existing Posterior Inventory Clarification, July 1

Important correction/clarification:

We are **not** starting from zero on ALDERAAN. The local posterior archive already covers most of the current best sample.

| population | total planets | extracted posteriors | coverage | missing posteriors |
|---|---:|---:|---:|---:|
| thick multis | 233 | 160 | 68.7% | 73 |
| thick singles | 272 | 149 | 54.8% | 123 |
| thin multis | 860 | 664 | 77.2% | 196 |
| thin singles | 1109 | 756 | 68.2% | 353 |
| **ALL** | **2474** | **1729** | **69.9%** | **745** |

The earlier "1065 need ALDERAAN" number combines:

| bucket | planets |
|---|---:|
| missing extracted posterior | 745 |
| existing posterior but flagged for review/refit/stress-test | 320 |
| total on broad ALDERAAN-needed manifest | 1065 |

Operational queues now written:

| queue | planets | targets | thin singles | thick singles | thin multis | thick multis |
|---|---:|---:|---:|---:|---:|---:|
| baseline usable existing posteriors | 1409 | 1069 | 602 | 116 | 558 | 133 |
| existing posterior flagged review/refit | 320 | 302 | 154 | 33 | 106 | 27 |
| missing extracted posterior total | 745 | 617 | 353 | 123 | 196 | 73 |
| missing launchable ALDERAAN now | 718 | 592 | 334 | 120 | 193 | 71 |
| missing unseeded manual seed needed | 27 | 27 | 19 | 3 | 3 | 2 |

New files:

| file | meaning |
|---|---|
| `<user-home>\Documents\project\2026-06-29\8\outputs\existing_eccentricity_posteriors_best.csv` | all 1729 extracted posterior summaries with FITS file references |
| `<user-home>\Documents\project\2026-06-29\8\outputs\posterior_coverage_status_best.csv` | coverage by population |
| `<user-home>\Documents\project\2026-06-29\8\outputs\posterior_status_counts_best.csv` | counts by missing / usable / flagged |
| `<user-home>\Documents\project\2026-06-29\8\outputs\posterior_operational_queues_best.csv` | queue counts for planning |
| `<user-home>\Documents\project\2026-06-29\8\outputs\posterior_inventory_and_queues_best.md` | readable version |
| `<user-home>\Documents\project\2026-06-29\8\outputs\posterior_queue_baseline_usable_best.csv` | 1409 existing posterior rows not currently flagged |
| `<user-home>\Documents\project\2026-06-29\8\outputs\posterior_queue_existing_flagged_for_review_best.csv` | 320 existing posterior rows worth inspecting/refitting |
| `<user-home>\Documents\project\2026-06-29\8\outputs\posterior_queue_missing_launchable_alderaan_best.csv` | 718 missing-posterior rows launchable now |
| `<user-home>\Documents\project\2026-06-29\8\outputs\posterior_queue_missing_unseeded_best.csv` | 27 missing-posterior rows blocked by missing transit depth |

Existing flagged-posterior reasons:

| flag | flagged existing posterior planets |
|---|---:|
| `flag_zeta_lt_0p7` | 242 |
| `flag_e50_gt_0p5` | 239 |
| `flag_low_snr_lt_20` | 178 |
| `flag_ror_shift_25pct` | 79 |
| `flag_duration_shift_25pct` | 65 |
| `flag_zeta_gt_1p3` | 23 |
| `flag_impact_gt_1` | 8 |
| `flag_results_file_bad_match` | 0 |

Interpretation:

- The posterior archive is highly useful and covers about 70% of the sample.
- The most defensible immediate baseline is the 1409 existing posteriors not flagged by automated checks.
- The 320 flagged existing posteriors are not "missing"; they are suspicious or scientifically high-leverage and should be inspected/refit selectively.
- The true missing-posterior ALDERAAN queue is 745 planets, of which 718 are launch-ready.
- Individual posterior medians (`e50`) are often high, especially in flagged tails, so the archive should be used through hierarchical posterior-sample inference rather than by over-interpreting per-planet medians.

## Cloud-Ready ALDERAAN Preparation, July 1

Question:

> Do we have to run ALDERAAN on Google Cloud, and can everything else be done first?

Answer:

For the full missing-posterior run, yes, cloud compute is the practical route. But the local preparation is now complete enough that the cloud stage should be mostly execution, not analysis design.

New/updated local code:

| file | purpose |
|---|---|
| `<repository-root>\prepare_gcp_missing_alderaan_bundle.py` | creates the cloud bundle from the true missing-posterior launch-ready queue |
| `<repository-root>\extract_eccentricity_posteriors.py` | patched to accept `--run-id`, cloud-specific summary paths, and posterior subdirectories |
| `<repository-root>\merge_posterior_summaries.py` | merges new cloud posterior summaries into the existing 1729-row posterior archive |
| `<repository-root>\postprocess_missing_cloud_results.ps1` | local one-command postprocessing after the cloud tarball is retrieved |

Generated cloud bundle:

`<repository-root>\cloud_missing_batch`

Bundle validation:

`validate_bundle.py` reports:

```text
VALIDATION OK: 592 targets, 767 catalog rows
```

Cloud bundle counts:

| quantity | count |
|---|---:|
| true missing-posterior launch-ready planets | 718 |
| KOI systems to run | 592 |
| ALDERAAN catalog rows after whole-system expansion | 767 |
| thick multi missing planets | 71 |
| thick single missing planets | 120 |
| thin multi missing planets | 193 |
| thin single missing planets | 334 |

Shard layout:

| shard | targets | missing planets |
|---:|---:|---:|
| 0 | 148 | 185 |
| 1 | 148 | 237 |
| 2 | 148 | 148 |
| 3 | 148 | 148 |

Copied cloud planning files:

| file | purpose |
|---|---|
| `<user-home>\Documents\project\2026-06-29\8\outputs\pre_cloud_alderaan_checklist.md` | concise execution checklist |
| `<user-home>\Documents\project\2026-06-29\8\outputs\cloud_missing_manifest.md` | bundle manifest |
| `<user-home>\Documents\project\2026-06-29\8\outputs\README_GCP_MISSING.md` | GCP command guide |
| `<user-home>\Documents\project\2026-06-29\8\outputs\targets_missing_launchable.csv` | 592 target execution list |
| `<user-home>\Documents\project\2026-06-29\8\outputs\shard_summary_missing_launchable.csv` | shard sizes |

Recommended cloud run order:

1. Create VM and copy `cloud_missing_batch`.
2. Run `python validate_bundle.py` on the VM.
3. Run one shard first:

```bash
TARGET_CSV=shards/targets_shard_000.csv JOBS=8 bash run_batch.sh
```

4. Retrieve that shard and run:

```powershell
.\postprocess_missing_cloud_results.ps1 -TarPath "C:\path\to\alderaan_results_sagear_missing_*.tar.gz"
```

5. If extraction/merge works, run the full queue:

```bash
JOBS=30 bash run_batch.sh
```

Current recommendation:

Do **not** spend cloud time on the 320 existing-but-flagged posteriors yet. First fill the 718 true missing launch-ready posteriors. Then merge, rerun hierarchical fits, and only then decide whether the 320 suspicious existing posteriors need refits.

## 2026-07-01: Current Operational Decision

The best next step is not to rerun everything. It is to run only the true missing-posterior ALDERAAN queue first, in a staged cloud execution.

Why:

- We already have useful extracted eccentricity posterior summaries for 1729 / 2474 planets.
- The scientifically missing set is 745 planets, of which 718 are launch-ready in ALDERAAN right now.
- The remaining 27 missing rows are blocked by missing or unreliable catalog transit seeds; local historical catalogs did not provide a safe automatic recovery.
- There are 320 existing posterior rows flagged for review or possible refit, but spending cloud time on those before filling the true missing set would be premature.

The current cloud bundle is:

`<repository-root>\cloud_missing_batch`

The ready-to-upload zip is:

`<user-home>\Documents\project\2026-06-29\8\outputs\cloud_missing_batch_ready_for_gcp.zip`

Validation was rerun on 2026-07-01:

```text
VALIDATION OK: 592 targets, 767 catalog rows
```

Important counts:

| quantity | count |
|---|---:|
| missing launch-ready planet rows | 718 |
| target systems to run | 592 |
| ALDERAAN catalog rows after whole-system expansion | 767 |
| existing posterior rows already extracted | 1729 |
| existing posterior rows flagged for later review/refit | 320 |
| missing-depth blocked rows | 27 |

Cost-control / risk-control decision:

- Use a Spot `e2-standard-32` VM first, because the workload is embarrassingly parallel and resumable by KOI target.
- Use a 150 GB SSD boot disk even though the expected data volume is smaller, because ALDERAAN intermediate products and environment setup are more expensive to lose than the extra one-day disk cost.
- Run one shard first with `JOBS=8`, retrieve it, and verify local extraction/merge before launching the full queue.
- Only after the missing posteriors merge cleanly should we consider rerunning the 320 existing-but-flagged posterior cases.

Local environment check:

- `gcloud` is not currently available on the Windows PATH in this project session.
- Therefore the recommended launch path is Google Cloud Shell, not local Windows Cloud SDK setup.
- A Cloud Shell quickstart has been written to:

`<user-home>\Documents\project\2026-06-29\8\outputs\cloud_shell_quickstart.md`

Stop-cost requirement:

After the tarball is packed and copied out, delete the VM:

```bash
gcloud compute instances delete alderaan-missing-e2-32 --zone $ZONE
```

2026-07-01 billing safety update:

- The VM creation script now sets `MAX_RUN_DURATION=20h` by default and passes `--max-run-duration "$MAX_RUN_DURATION"` to `gcloud compute instances create`.
- This stops compute automatically after 20 hours, but it does not eliminate all possible costs because stopped VM disks can still incur storage charges until deletion.
- A no-surprise-cost checklist has been written to:

`<user-home>\Documents\project\2026-06-29\8\outputs\gcp_no_charge_safety_checklist.md`

- Do not create the VM unless the Google Cloud Billing page confirms an active Free Trial / credit-covered billing account, a university/grant billing account, or you explicitly accept possible real charges.

What success looks like after the first shard:

- ALDERAAN finishes for a nontrivial number of targets without environment errors.
- The returned `alderaan_results_sagear_missing_*.tar.gz` can be unpacked locally.
- `postprocess_missing_cloud_results.ps1` produces a new `eccentricity_posterior_summary_sagear_missing.csv`.
- The new summary merges into the existing archive without duplicated planet IDs or broken category labels.
- The thin-single posterior coverage increases before any refit of already-existing posteriors.

## 2026-07-13 factorial-validation partial snapshot

The repaired factorial runner was resumed after the VM's 24-hour automatic
stop. Existing nonempty FITS are detected and skipped, so the restart did not
repeat the 47 completed fits. The cloud matrix continues independently.

A status-gated snapshot was taken while the interrupted `K01127` reference-LD
fit was rerunning:

| snapshot item | value |
|---|---:|
| completed original-LD long-cadence systems | 24 |
| completed reference-LD long-cadence systems | 23 |
| completed FITS total | 47 / 82 |
| archive size | 118,322,978 bytes |
| SHA-256 | `e8d9a772b30fe6a7e23a54989ba4dc71e100beefc92d1dc552fe35ed95784bc8` |

The local and VM checksums matched exactly. The archive is stored at:

`<validation-artifact-dir>\alderaan_factorial_partial_20260713T102557Z.tar.gz`

The factorial comparator previously required all 82 expected FITS, which made
the planned concurrent partial analysis impossible. It now has an explicit
`--allow-incomplete` mode. Strict full-matrix behavior remains the default;
partial mode extracts only discovered FITS, compares only complete within-planet
pairs, preserves the full discovery/missing audit, and marks practical-effect
claims unassessable until the repeatability arm exists. The focused comparator
suite passes 8/8 tests in both the Git checkout and live pipeline.

The 23 matched systems contain 31 matched planets. All 31 passed direct
importance-extraction QC. Reference minus original limb-darkening results are:

| quantity | median signed shift | 95% system-bootstrap interval | median absolute shift | 95th percentile absolute shift |
|---|---:|---:|---:|---:|
| eccentricity posterior median | +0.00168 | [-0.00103, +0.00523] | 0.01568 | 0.05548 |
| zeta posterior median | -0.00151 | [-0.00464, +0.00837] | 0.00899 | 0.10507 |
| impact parameter median | +0.00944 | [-0.00089, +0.02143] | 0.02763 | 0.14118 |
| duration T14 (hours) | -0.00529 | [-0.01458, +0.01934] | 0.02058 | 0.30409 |
| Rp/Rstar | +0.000155 | [+0.000053, +0.000249] | 0.000188 | 0.01061 |

For the 13 thin-single planets in this targeted validation subset, the median
eccentricity shift is +0.00168 and the median absolute shift is 0.00604. This
does not support a coherent limb-darkening shift large enough, by itself, to
explain the thin-single population discrepancy. Individual systems remain
sensitive: the largest absolute eccentricity shifts currently include
`K00283.03` (0.080), `K02533.03` (0.072), `K01001.02` (0.039), and
`K02109.01` (0.039). `K02533.03` also changes duration by 2.22 hours and must
be inspected as a high-leverage case rather than allowed to define a global
limb-darkening conclusion.

This is provisional, not the final factorial result. There is no repeat-seed
arm in the snapshot, so no observed shift can yet be distinguished from
ALDERAAN nested-sampling variability. The cadence, reference-LD plus
short-cadence, repeatability, and printed-prior arms remain necessary and are
still running in GCP.

Partial analysis products are under:

`<repository-root>\outputs\factorial_validation_partial_20260713T102557Z`

## 2026-07-13 multiplicity contract and current four-bin result

The professor's notebook diagnosis is correct for the historical workflow:
single/multi architecture must be determined before planet-level quality cuts.
The current canonical sample already does this from all KOIs with archive
disposition `CONFIRMED` or `CANDIDATE`. A direct data audit established:

| multiplicity audit | count |
|---|---:|
| canonical planet rows | 2474 |
| rows missing a raw non-FP host count | 0 |
| current labels disagreeing with the pre-cut raw count | 0 |
| planets that an after-cut recount would mislabel | 45 |
| affected thin-multi planets/hosts | 41 / 41 |
| affected thick-multi planets/hosts | 4 / 4 |

`common.add_target_and_system` now accepts an explicit host-indexed pre-cut
multiplicity series and fails closed if any KIC is missing. The canonical
sample builder passes this series directly; the filtered recount remains only
an explicit sensitivity option. Three dedicated regression tests were added,
and the full live test suite passes 23/23.

The current method-consistent population result is the 710-planet uniform
direct-posterior subset under the count-calibrated q=0.535 classifier context.
It is diagnostic because coverage is incomplete and nonrandom. Under Sagear's
printed reciprocal transit-selection convention, the Rayleigh expected-value
comparison is:

| population | current fit N | current mean e (16th-84th) | Sagear full N | Sagear mean e (16th-84th) | current/Sagear |
|---|---:|---:|---:|---:|---:|
| thin singles | 304 | 0.335 (0.319-0.353) | 1121 | 0.022 (0.017-0.029) | 15.2x |
| thick singles | 108 | 0.288 (0.263-0.315) | 275 | 0.066 (0.045-0.096) | 4.4x |
| thin multis | 222 | 0.122 (0.108-0.136) | 862 | 0.030 (0.023-0.031) | 4.1x |
| thick multis | 69 | 0.117 (0.086-0.154) | 207 | 0.033 (0.015-0.065) | 3.5x |

The corresponding q=0.535 sample inventory has 1096 thin singles, 285 thick
singles, 886 thin-multi planets, and 207 thick-multi planets, but only the fit-N
subsets above have uniform paired-impact direct posteriors. These results must
not be represented as a final full-sample replication.

The same after-cut recount pattern was also found in the published-label
overlap diagnostic. It did not feed the hierarchical fit, but it biased the
reported overlap-bin counts. After preserving the canonical architecture, the
broad no-binary-cut overlap contains 1052 thin singles, 803 thin-multi planets,
269 thick singles, and 205 thick-multi planets; 39 rows would have been labeled
differently by the overlap-only recount. This supersedes the historical
1082/773/278/196 overlap counts.

Two related fail-open defaults were removed during the same audit. The
system-definition sensitivity no longer invents a `single` label when a KIC is
absent from its multiplicity map, and the legacy posterior merger no longer
treats unknown system completeness as known complete. Missing multiplicity now
raises an error; unknown completeness is explicitly recorded and excluded by
primary QC. These changes harden diagnostics and legacy merges but do not alter
the current 710-row uniform direct-posterior result.

An aggressive combined quality sensitivity (confirmed, high-S/N, dwarf,
small-planet filters) gives thin singles 0.022 on only 49 objects, but the fit
is boundary-sensitive and the multis remain high. It diagnoses where the
tension is concentrated; it is not a defensible replacement for Sagear's
stated sample.

At 2026-07-13 10:46 UTC the factorial VM was running with 47/82 FITS complete.
It restarted at 2026-07-13 10:01:17 UTC with a 24-hour maximum runtime, so its
next automatic stop is approximately 2026-07-14 10:01:17 UTC (15:31:17 IST).
Changing the maximum runtime requires stopping the VM. The cost-controlled
default is therefore to let it run, then restart and use the repaired
FITS-aware runner, which skips every completed result and grants a fresh
24-hour window.

## 2026-07-13 scientific interim hardening

Two additional diagnostics were completed locally while the factorial VM continued.
They do not consume cloud compute and do not alter the running matrix.

First, the 31 matched original/reference limb-darkening planets were analyzed with
target-system clustered bootstrap uncertainty, leave-system-out leverage, impact and
posterior-width sensitivity cuts, and exploratory parameter correlations. The global
median change in `e50` remains +0.00168. Removing the highest-leverage system changes
that median by only -0.00018, requiring impact below 0.8 changes it by -0.00060, and
1,000 random leave-10%-of-systems trials have a 95th percentile absolute median shift
of 0.00083. The current thick-multi subgroup remains underpowered because it contains
only one matched system until K01127 finishes.

Second, the current 710-planet direct posterior population input was tested by treating
each host system as the resampling unit. Under the printed reciprocal selection mode:

| population | planets / hosts | host-bootstrap e16-e84 | fraction of leave-10%-host trials within 5% | largest one-host shift |
|---|---:|---:|---:|---:|
| thin singles | 304 / 304 | 0.313-0.357 | 97.5% | 2.5%, K04921 |
| thick singles | 108 / 108 | 0.252-0.320 | 83.6% | 8.4%, K04943 |
| thin multis | 222 / 101 | 0.099-0.139 | 67.6% | 8.8%, K02857 |
| thick multis | 69 / 28 | 0.072-0.181 | 31.9% | 34.3%, K01992 |

This sharpens the diagnosis. The high thin-single result is not caused by one host, but
the current multi results, especially thick multis, are not robust to system-level
dependence. A forward-normalized selection sensitivity produces the same qualitative
conclusion. The complete synthesis is in
`<repository-root>\docs\scientific_interim_assessment.md`.
