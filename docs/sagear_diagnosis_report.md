# Sagear Replication Diagnosis Snapshot

> **2026-07-10 superseding audit:** The historical sections below are retained
> as an audit trail, but several earlier claims are no longer current. In
> particular, the 1680-row archive is not a uniform paired-impact product, the
> old Toomre script refit an unrelated unsupervised GMM, and a 20-fit one-arm
> validation was insufficient. The following section is the current decision
> record.

## Deep Recovery Audit (2026-07-10)

### Bottom line

This is not a dead end, but the current mixed 2395-planet result is not a
publishable replication. The 547 recovered ALDERAAN FITS files support 710
planet posteriors under an exact sample-level importance extractor. The other
1680 archived grids were generated with synthetic geometric impact draws and
cannot be repaired without their original paired transit samples or a rerun.

The remaining disagreement is not primarily the thin/thick label. On the 710
uniform direct posteriors, changing from a method-first classifier to a
count-calibrated classifier moves memberships substantially but leaves singles
very warm. The recovered transit fits or their inputs differ materially from
the fits used by Sagear.

### Newly established facts

- The four published planet bins are 1121 thin singles, 275 thick singles, 862
  thin multi planets, and 207 thick multi planets. They sum exactly to 2465.
- `378` is an inconsistent total-thick-host macro used incorrectly in prose.
  The coherent host arithmetic is 275 single hosts + 98 multi hosts = 373
  thick hosts, and 1515 thin + 373 thick = 1888 total hosts.
- The final published paper appeared as AJ 172:42 on 2026-06-23. Its conclusion
  explicitly identifies Berger et al. 2020 density priors, resolving the
  submitted-source Berger-2018/Berger-2020 conflict in favor of Berger 2020 for
  the primary reconstruction. Berger 2018 flags remain labeled sensitivities.
- Figure 1 and the methods require Galactocentric cylindrical `V_phi` and
  `sqrt(V_R^2+V_Z^2)`, with the GMM calibrated using high/low-alpha chemistry.
- The raw Angus table contains the EDR3 positions and Bailer-Jones distances
  paired with its Cartesian velocities. The former classifier mixed those
  velocities with another Gaia table and inverse-parallax distances; this is
  now corrected.
- The previous Sagear-style Toomre script ignored the chemical classifier and
  fit a fresh unsupervised GMM to planet hosts. It now applies the serialized
  chemically calibrated means, covariances, and class prior to the Angus field.

### Sample and classifier reconstruction

The explicit-paper sample without an extra Berger binary flag has 2642 planets
around 1970 hosts before transit-fit convergence removal. Adding the plausible
but unstated Berger 2018 `Bin=0` sensitivity gives 2474 planets around 1857
hosts. The latter is numerically close to 2465 but is not promoted solely for
that reason.

Multiplicity is now defined by default from all non-false-positive KOIs around
the host before the period and stellar cuts. This preserves known system
architecture and matches the meaning of single/multi-transit systems better
than relabeling a multi as a single when one companion falls outside 1-100 d.

The chemical classifier is highly prior-sensitive. For the `Bin=0`, raw-KOI
multiplicity hypothesis:

| high-alpha prior | thin single | thick single | thin multi | thick multi |
|---:|---:|---:|---:|---:|
| 0.500 | 1183 | 198 | 950 | 143 |
| 0.535, count-calibrated diagnostic | 1096 | 285 | 886 | 207 |
| Sagear | 1121 | 275 | 862 | 207 |

The 0.535 row is not an inferred physical truth. It demonstrates that an
undocumented 0.035 change in component prior moves hundreds of stars. Exact
classifier reproduction therefore requires Sagear's GMM weights or final
planet table. The new classifier sweep is in
`classifier_reconstruction_methods.csv` and
`classifier_reconstruction_prior_sweep.csv`.

The corrected comparison plot is
`toomre_sagear_reference_comparison.png`; its middle panel now uses the actual
chemical model. The old diagnostic panels are retained but must not be cited as
Sagear-equivalent.

### Exact eccentricity extraction

`extract_eccentricity_posteriors_direct.py` now implements the commented
Sagear/MacDougall equation directly:

- paired `T14`, `Rp/R*`, and `b` rows from ALDERAAN;
- normalized dynesty nested-sampling weights;
- uniform `e` and `omega` proposals;
- Berger 2020 density likelihood;
- no zeta KDE and no transit prior in the individual posterior;
- deterministic RNG, one-to-one period matching, ESS, and explicit exclusions;
- common `e,omega` grids with `omega endpoint=False`.

Six focused tests pass, including circular injections at `b=0,0.5,0.8`. The
split-error sensitivity now includes the Gaussian scale normalization. The
extractor intentionally refuses to consume the old diagnostic sample by
default; an explicit sample is required until classification is resolved.

The 547 successful FITS systems yielded 710 matched planet posteriors. Seven
rows have importance ESS below 100 and are flagged, not silently deleted.
There are 2195 explicit sample rows without a corresponding recovered FITS in
the broad mapping table.

For the 651 planets overlapping the earlier paired-zeta extraction, direct and
old median eccentricities correlate at 0.982. Therefore the extreme individual
posteriors are not caused by the direct extractor replacing the zeta KDE.

### Uniform-subset population results

These are diagnostic because they cover only the 547 recovered systems, not a
representative complete Sagear sample. They use Sagear's printed reciprocal
transit-probability likelihood and a uniform prior on Rayleigh sigma.

| classifier context | population | n after ESS QC | inferred mean e |
|---|---|---:|---:|
| explicit-paper/equal prior | thin singles | 334 | 0.331 |
| explicit-paper/equal prior | thick singles | 78 | 0.291 |
| count-calibrated q=0.535 | thin singles | 304 | 0.335 |
| count-calibrated q=0.535 | thick singles | 108 | 0.288 |
| count-calibrated q=0.535 | thin multis | 222 | 0.122 |
| count-calibrated q=0.535 | thick multis | 69 | 0.117 |

The classifier changes the membership but not the conclusion: recovered
singles are much warmer than Sagear. Omitting transit selection or using the
scientifically generative forward-normalized selection also leaves them warm,
so the discrepancy is upstream of that convention.

The recovered ALDERAAN shapes are not obviously nonsensical in aggregate:
median fitted/catalog duration ratio is 1.004 and median fitted impact is 0.51.
However, the thin-single zeta distribution is shifted and enough individual
systems favor short durations relative to Berger density to warm the
hierarchical result. Top-k removal is not robust, and the full mixed product
must remain rejected.

### Quality-cut diagnosis, not a replacement sample

Labeled sensitivities show where the tension is concentrated:

| diagnostic slice | thin-single mean e | thick-single mean e |
|---|---:|---:|
| all direct subset | 0.335 | 0.288 |
| Berger logg >= 4 | 0.336 | 0.237 |
| confirmed KOIs only | 0.262 | 0.196 |
| KOI S/N >= 10 | 0.316 | 0.249 |
| planet radius < 3.5 Earth radii | 0.245 | 0.172 |
| confirmed + S/N + logg + small | 0.022 | 0.039 |

The final row is boundary-sensitive and retains only 49 thin singles and 26
thick singles. It is evidence that weak candidates, large planets, and
non-dwarf hosts concentrate the discrepancy; it is not permission to claim a
replication by discarding most of Sagear's stated sample.

### Hierarchical inference repairs

`hierarchical_rayleigh.py` now:

- defaults to the reciprocal convention printed in Sagear's commented HBM;
- keeps a separate forward-normalized mode, which passes a synthetic recovery
  test and is scientifically generative but is not silently substituted for
  the manuscript method;
- rejects mixed posterior sources, non-paired impact modes, missing QC, corrupt
  grids, and mixed grids in canonical runs;
- reports posterior medians and 16/84 intervals under uniform sigma, plus map
  and boundary flags;
- writes per-planet leverage, random leave-10%-out, and top-k diagnostics.

A misplaced diagnostic-output block that would crash correctly tagged direct
summaries was found and fixed. Contract regression tests now cover this path.

### Factorial ALDERAAN validation before any full rerun

The old 20-target one-arm validation has been replaced by an 82-fit paired
factorial design built around 24 systems / 34 seedable planets:

- 20 targeted high-leverage/stable-control systems;
- four deterministic random controls, one per population;
- original versus reference limb-darkening centers on all 24;
- long-only versus long+short cadence on the nine audited SC systems;
- eight same-configuration repeat fits with a different seed;
- eight Sagear-Table-1 prior fits in a separate patched ALDERAAN clone.

All arms use the same pinned public commit
`7443dff16b7f9092e14a6f0cc1f8948d457c9e0b`, deterministic target seeds,
full known KOI systems (including seedable companions outside 1-100 d), and
per-target provenance. The public-code and manuscript-prior hypotheses remain
separate. Bash syntax, Python compilation, full-system `npl`, and bundle
validation pass.

`compare_factorial_validation.py` performs strict FITS discovery, exact direct
extraction, paired shape/eccentricity comparisons, deterministic system-cluster
bootstrap summaries, and calibration of effects against observed repeat-run
variation. Sixteen related tests pass. It will fail clearly on incomplete arms
instead of manufacturing a comparison.

### What must happen next

1. Run the 82-fit factorial validation, not the old 20-fit script.
2. Analyze it locally with `compare_factorial_validation.py`.
3. Select the full-run configuration only if limb darkening, cadence, or Table
   1 priors move the high-leverage tail by more than repeatability variation.
4. Then rerun the approximately 1680 geometric-impact archive systems with the
   winning configuration. The old NPZ grids cannot be repaired in place.
5. Re-extract every successful FITS uniformly with the direct extractor, attach
   one declared classifier context, and require leave-10%-out stability before
   comparing Table 2.
6. Request the artifacts listed in `author_clarification_request.md`; receiving
   the final Sagear table or paired posteriors could eliminate much of the
   otherwise necessary rerun.

There is no honest path from the present mixed archive to a final replicated
Table 2 without either broader ALDERAAN reruns or author products. There is,
however, a concrete and tested path forward; the project is constrained, not
dead.

This is the current state after correcting the Sagear thick-single target count, adding the likely hidden FGKM/Berger-binary sample cuts, rebuilding the APOGEE crossmatch, and regenerating the Toomre/eccentricity/ALDERAAN diagnostics.

## Corrections Made

The previous local target count for thick singles was wrong. The manuscript macros and `thin_thick_rayleigh.png` indicate:

- thin singles: 1121
- thick singles: 275
- thin multi planets: 862
- thick multi planets: 207
- thin multi hosts: 394
- thick multi hosts: 98

The earlier `378` value is the total thick-host macro/paragraph value, not the thick-single planet count.

The manuscript macros are internally inconsistent:

| check | left | right | delta |
|---|---:|---:|---:|
| `allthinplanets + allthickplanets` vs `allplanets` | 2474 | 2465 | +9 |
| subgroup planet macros vs `allplanets` | 2465 | 2465 | 0 |
| thick subgroup planet macros vs `allthickplanets` | 482 | 491 | -9 |
| `allthinstars + allthickstars` vs `allstars` | 1893 | 1888 | +5 |
| subgroup host-like macros vs `allstars` | 1888 | 1888 | 0 |

Therefore the current corrected sample total of 2474 planets matches one set of Sagear disk-total macros. The 9-planet gap to `allplanets=2465` should not be overfit without external clarification or the original analysis table.

## Sample Audit

Two additional sample cuts were promoted into the canonical audit:

- `berger_teff < 6500`, matching the stated FGKM sample description;
- Berger+2018 `Bin = 0`, a likely resolved-companion/binary exclusion inherited from the Berger 2018 stellar catalog.

These cuts reduce the current sample from 2795 to 2474 planets before ALDERAAN convergence/visual-fit removal. Sagear reports 2465 planets, so the sample-size mismatch is now only 9 planets.

Current diagnostic sample:

| stage | planets | hosts |
|---|---:|---:|
| after Berger density join | 3139 | 2339 |
| after `Teff < 6500` | 2970 | 2208 |
| after Berger+2018 `Bin=0` | 2770 | 2071 |
| after Angus velocities | 2474 | 1857 |

Sagear reports 2465 planets / 1888 stars. Planet count is now very close; host count is still 31 low, suggesting catalog-version, duplicate/multiplicity, or velocity-availability differences remain.

## Disk Counts

The strict pooled APOGEE GMM is still not Sagear-equivalent:

| group | ours | Sagear | delta |
|---|---:|---:|---:|
| thin singles | 1332 | 1121 | +211 |
| thick singles | 94 | 275 | -181 |
| thin multi planets | 987 | 862 | +125 |
| thick multi planets | 61 | 207 | -146 |
| thin multi hosts | 406 | 394 | +12 |
| thick multi hosts | 25 | 98 | -73 |

The diagnostic all-Angus planet-host GMM is now very close:

| group | ours | Sagear | delta |
|---|---:|---:|---:|
| thin singles | 1121 | 1121 | 0 |
| thick singles | 305 | 275 | +30 |
| thin multi planets | 822 | 862 | -40 |
| thick multi planets | 226 | 207 | +19 |
| thin multi hosts | 336 | 394 | -58 |
| thick multi hosts | 95 | 98 | -3 |

This means the largest previous problem was upstream sample construction, not eccentricity.

## Toomre/Classifier Diagnosis

The professor's Toomre concern was real. The plot must be compared with Sagear's sign/display convention for the rotational velocity. After the sample correction, the best diagnostic classifier variants are:

| variant | thin singles | thick singles | thin multi planets | thick multi planets | thick multi hosts |
|---|---:|---:|---:|---:|---:|
| planet-host GMM, direct Angus velocities | 1121 | 305 | 822 | 226 | 95 |
| planet-host GMM, geometric cylindrical velocities | 1148 | 278 | 819 | 229 | 96 |
| supervised chemical Gaussian, prior 0.425 | 1150 | 276 | 875 | 173 | 73 |
| supervised mix 2 thin / 2 thick, best prior | 1146 | 280 | 861 | 187 | 78 |

No current chemically calibrated model reproduces all groups exactly. The exact Sagear GMM convention and prior/mixing treatment remain the main unresolved methodological detail.

### Old Hilldale Zip Check

The old project zip was extracted read-only into `external/old_hildale_project_zip/`. Its classification code was useful, but its saved sample is not a target to revert to:

- old saved `sample.csv`: 1716 planets / 1278 hosts;
- current corrected diagnostic sample: 2474 planets / 1857 hosts;
- overlap: 1545 planets;
- old disk counts: 875 thin singles, 100 thick singles, 671 thin multi planets, 70 thick multi planets.

I ported the old coordinate convention into the current pipeline as `V_phi_astropy` and `V_perp_astropy`. This reproduces the old Astropy Galactocentric cylindrical projection where rotating disk stars sit near `V_phi = -220 km/s`.

Expanded classifier variants:

| variant | thin singles | thick singles | thin multi planets | thick multi planets | thick multi hosts |
|---|---:|---:|---:|---:|---:|
| planet-host direct Angus GMM | 1121 | 305 | 822 | 226 | 95 |
| planet-host geometric cylindrical GMM | 1148 | 278 | 819 | 229 | 96 |
| planet-host old-Astropy GMM | 1146 | 280 | 823 | 225 | 95 |
| old-pipeline KIC-wide Astropy GMM | 1031 | 395 | 742 | 306 | 127 |

Conclusion: the old zip confirms a coordinate/sign convention issue for the Toomre plot, but the old KIC-wide GMM does not fix the Sagear counts. It overclassifies thick-disk planets relative to Sagear by +120 thick singles and +99 thick multi planets. The closest old-code-inspired option is the planet-host old-Astropy GMM, which nearly matches thick singles but still leaves multis low/high by tens of planets.

## Threshold Sweep

A threshold sweep was added to test whether the remaining mismatch is only a `P_thick` boundary issue.

Best rows:

| classifier | threshold | thin singles | thick singles | thin multi planets | thick multi planets | L1 planet delta |
|---|---:|---:|---:|---:|---:|---:|
| planet-host direct GMM | 0.610 | 1151 | 275 | 839 | 209 | 55 |
| planet-host geometric GMM | 0.525 | 1154 | 272 | 828 | 220 | 83 |
| chem 1-thin/1-thick prior 0.425 | 0.490 | 1121 | 305 | 858 | 190 | 51 |
| chem 2-thin/2-thick prior 0.350 | 0.460 | 1121 | 305 | 846 | 202 | 51 |

Interpretation: threshold choice matters, but it does not fully reproduce Sagear. The probability field and/or remaining sample membership still differ.

Pairwise classifier disagreement:

| comparison | changed planets | changed systems |
|---|---:|---:|
| direct GMM `P>0.50` vs direct GMM `P>0.61` | 47 | 38 |
| direct GMM `P>0.61` vs geometric GMM `P>0.525` | 230 | 183 |
| direct GMM `P>0.61` vs chemical 2-thin/2-thick | 247 | 185 |

Interpretation: changing only the threshold moves a small boundary set. Changing coordinate/probability construction moves hundreds of planets. The key unresolved issue is the exact kinematic probability construction, not just the threshold.

## Eccentricity Triage

No local ALDERAAN posterior products are present yet. These eccentricity diagnostics use the older `e_photo` point estimates only.

With the corrected diagnostic sample:

| population | n with e | mean e | median e | q95 e | Sagear Rayleigh `<e>` |
|---|---:|---:|---:|---:|---:|
| thin singles | 958 | 0.0323 | 0.0233 | 0.0855 | 0.022 |
| thick singles | 262 | 0.0351 | 0.0231 | 0.0891 | 0.066 |
| thin multis | 742 | 0.0357 | 0.0241 | 0.1001 | 0.030 |
| thick multis | 205 | 0.0362 | 0.0252 | 0.0961 | 0.033 |

The thin-single median remains close to Sagear, while the mean is inflated by a high-e tail. The tail is concentrated in short-period and often high-impact systems, so ALDERAAN validation should focus there.

## Formula Audit

The eccentricity posterior/hierarchical scripts were corrected before ALDERAAN scaling:

- individual `e, omega` posteriors are now extracted without multiplying by the geometric transit prior by default;
- the Rayleigh population fit now reads the full joint posterior grid instead of only marginal `e_pdf`;
- the hierarchical likelihood now applies Sagear's reciprocal geometric transit-probability factor `(1-e^2)/(1+e sin omega)` over the joint grid.

Formula sanity checks passed:

- transit-selection correction is finite and positive;
- periastron correction is smaller than apastron correction at fixed `e=0.5`;
- `zeta=1` prefers lower eccentricities than a short-duration `zeta=0.6` case;
- synthetic low-e posteriors recover a reasonable Rayleigh expected eccentricity.

## Berger 2018

After the new cuts, the diagnostic sample has complete Berger+2018 stellar matches:

- 2474 / 2474 planet rows matched;
- 1857 / 1857 hosts matched.

Berger 2018 availability is therefore not the remaining mismatch.

## ALDERAAN Status

The validation project and cloud batch were regenerated from the corrected sample. Existing downloaded ALDERAAN/Gilbert FITS were also found and processed.

Validation batch:

- 20 systems total;
- 8 high-e thin-single stress-test systems;
- 3 normal representatives each for thin singles, thick singles, thin multis, and thick multis.

Cloud batch:

- 210 systems;
- 300 planet rows;
- balanced to 75 planet rows per disk/multiplicity bin.

Existing local ALDERAAN results:

- path: `<alderaan-results-dir>`;
- 1692 FITS result files;
- 1278 corrected-sample systems with result files;
- 579 corrected-sample systems without result files;
- 1729 planet posterior grids extracted;
- 706 planets have no result file;
- 39 planets have a result file but no period-matched/usable posterior.

Coverage after extraction:

| group | sample planets | posterior extracted | missing/no usable posterior |
|---|---:|---:|---:|
| thin singles | 1121 | 779 | 342 |
| thick singles | 305 | 163 | 142 |
| thin multis | 822 | 637 | 185 |
| thick multis | 226 | 150 | 76 |

The extractor now matches ALDERAAN planets to KOIs by period from `TTIMES_XX`, treats Berger density uncertainties as log10 absolute errors, and marginalizes impact parameter geometrically by default.

Rayleigh fit to the currently extracted existing posteriors:

| population | N | current `<e>` with transit correction | Sagear `<e>` |
|---|---:|---:|---:|
| thick singles | 163 | 0.109 | 0.066 |
| thin singles | 779 | 0.139 | 0.022 |
| thick multis | 150 | 0.006 | 0.033 |
| thin multis | 637 | 0.006 | 0.030 |

These are not a Sagear reproduction. They are a new diagnostic showing that the existing ALDERAAN/Gilbert posteriors plus current posterior extraction still disagree strongly with Sagear, especially for thin singles and multis.

ALDERAAN shape diagnostics:

| group | N | median zeta | q16 zeta | q84 zeta | median e50 | median ALDERAAN/KOI duration |
|---|---:|---:|---:|---:|---:|---:|
| thick multis | 150 | 1.068 | 0.761 | 1.172 | 0.218 | 0.994 |
| thick singles | 163 | 0.969 | 0.617 | 1.183 | 0.281 | 1.001 |
| thin multis | 637 | 1.050 | 0.755 | 1.164 | 0.216 | 0.995 |
| thin singles | 779 | 0.962 | 0.691 | 1.159 | 0.277 | 0.999 |

Key interpretation:

- KOI and ALDERAAN durations agree at the group median, so the disagreement is not an obvious duration-unit or period-matching bug.
- Median Berger density fractional errors are modest, roughly 9-11%.
- High `e50` is strongly tied to `zeta` tails.
- The systems to inspect/refit are listed in `alderaan_shape_diagnostics_flagged.csv`.

Rayleigh sensitivity to shape cuts:

| filter | thick singles | thin singles | thick multis | thin multis |
|---|---:|---:|---:|---:|
| all | 0.109 | 0.139 | 0.006 | 0.006 |
| zeta >= 0.7 | 0.103 | 0.134 | 0.006 | 0.006 |
| 0.7 <= zeta <= 1.3 | 0.006 | 0.006 | 0.006 | 0.006 |
| no duration shift >25% | 0.107 | 0.137 | 0.006 | 0.006 |
| ALDERAAN b < 0.85 | 0.103 | 0.139 | 0.006 | 0.006 |
| SNR >= 20 | 0.113 | 0.122 | 0.006 | 0.006 |
| strict clean | 0.006 | 0.006 | 0.006 | 0.006 |

Removing only short-zeta systems does not erase the singles signal. Removing both short and long zeta tails collapses every population to the circular boundary. This makes the remaining problem very specific: how Sagear's exact ALDERAAN reductions, convergence/quality cuts, and posterior formalism treated the duration-ratio tails.

ALDERAAN refit planning is now concrete:

- `alderaan_refit_full_manifest.csv`
- `alderaan_refit_validation_targets.csv`
- `alderaan_missing_targets_all.csv`
- `alderaan_refit_plan.md`
- `alderaan_target_dossiers/index.md`

The validation batch contains:

| reason | systems |
|---|---:|
| suspicious existing ALDERAAN duration-ratio tail refits | 40 |
| missing ALDERAAN high-SNR validation systems | 32 |
| clean existing ALDERAAN controls | 20 |

The full missing ALDERAAN list contains 579 systems / 706 planets.

Visual dossiers were generated for the 40 most suspicious existing systems. The top case, `K03232.01`, shows a strongly pathological existing fit: ALDERAAN `Rp/Rs` is wildly shifted relative to KOI, impact extends beyond 1, and the derived eccentricity posterior peaks high.

The ALDERAAN repo exists locally, but the active conda environment is missing required packages (`pymc3`, `exoplanet`, `dynesty`, `batman`, `celerite2`, `ldtk`, `arviz`). The generated validation catalog and PowerShell scripts are ready; the remaining execution blocker is the ALDERAAN environment and actual reruns/refits.

## Current Best Diagnosis

The replication is now much healthier:

1. Sample construction is close to Sagear after adding FGKM and Berger+2018 binary cuts.
2. The Toomre/disk-classifier convention remains the biggest non-ALDERAAN mismatch.
3. Existing ALDERAAN/Gilbert FITS cover most, but not all, of the corrected sample.
4. Processing those FITS does not reproduce Sagear; the next high-value action is a system-level validation of `T14`, `Rp/Rs`, impact treatment, stellar density, and `zeta` for representative outliers before running new ALDERAAN at scale.

## 2026-07-07 Recovery Plan Implementation

I implemented the recovery-plan fixes in the live canonical folder:

`<repository-root>`

The important code changes are:

- `extract_eccentricity_posteriors.py`
  - Default extraction now uses paired ALDERAAN impact samples (`--impact-mode alderaan`) rather than synthetic geometric impact draws.
  - Stores both fitted ALDERAAN impact summaries and actually used impact summaries.
  - Uses stable deterministic RNG for dynesty resampling.
  - Uses `omega_grid endpoint=False`.
  - Makes the zeta KDE/support range consistent with the actual `e <= 0.95` eccentricity grid.
- `common.py`
  - Replaced the approximate circular-duration calculation with the exact circular transit-duration formula.
  - Added a NumPy-compatible `trapezoid` helper.
- `merge_posterior_summaries.py`
  - Writes a full QC manifest for every sample planet.
  - Recomputes zeta support flags for old and new rows.
  - Writes a primary-QC usable merged table.
  - Explicitly records replaced rows where new ALDERAAN results supersede older archive rows.
- `hierarchical_rayleigh.py`
  - Replaced the earlier point/grid-max fit with a posterior-summary Rayleigh fit.
  - Uses forward transit probability with population normalization.
  - Rejects mixed `include_transit_prior=True` posteriors in canonical transit-corrected runs.
  - Reports posterior median and 16/84 intervals for `<e>`, not just the grid maximum.
  - Adds per-planet leverage, top-leverage removal, leave-10%-out stability, and zeta/e50 ECDF diagnostics.
- `formula_sanity_checks.py`
  - Adds tests for transit probability, transit normalization, exact duration, zeta posterior behavior, and low-e synthetic Rayleigh recovery.
- `impact_mode_comparison.py`
  - New diagnostic comparing old geometric-impact extraction against paired ALDERAAN-impact extraction.

### Validation Status

The lightweight validation checks pass:

- `py_compile` passed for the patched scripts.
- `formula_sanity_checks.py` passed.
- The merge/QC manifest rerun completed with no pandas FutureWarnings.

### Paired-Impact Extraction Coverage

The corrected paired-impact extraction used the existing 547 ALDERAAN result FITS. It produced 705 new posterior summaries and excluded 7 planets whose zeta values were outside the `e <= 0.95` model support.

After merging these 705 rows into the existing 1729-row archive, the merged posterior table contains 2393 unique planets:

| disk | system | sample planets | posterior planets | missing planets | coverage |
|---|---|---:|---:|---:|---:|
| thick | multi | 233 | 225 | 8 | 0.966 |
| thick | single | 272 | 255 | 17 | 0.938 |
| thin | multi | 860 | 849 | 11 | 0.987 |
| thin | single | 1109 | 1064 | 45 | 0.959 |

Total sample remains 2474 planets. This is still slightly above Sagear's 2465, but the 9-planet difference is not large enough to explain the eccentricity mismatch.

### Primary Rayleigh Results After Fixes

These are the canonical merged paired-impact results before QC exclusion:

| population | N | `<e>` | 16th | 84th | Sagear `<e>` |
|---|---:|---:|---:|---:|---:|
| thick singles | 255 | 0.290 | 0.275 | 0.306 | 0.066 |
| thin singles | 1064 | 0.287 | 0.280 | 0.295 | 0.022 |
| thick multis | 225 | 0.071 | 0.061 | 0.082 | 0.033 |
| thin multis | 849 | 0.092 | 0.085 | 0.098 | 0.030 |

Primary-QC filtered results:

| population | N | `<e>` | 16th | 84th | Sagear `<e>` |
|---|---:|---:|---:|---:|---:|
| thick singles | 248 | 0.230 | 0.210 | 0.251 | 0.066 |
| thin singles | 1027 | 0.256 | 0.248 | 0.264 | 0.022 |
| thick multis | 221 | 0.071 | 0.061 | 0.082 | 0.033 |
| thin multis | 835 | 0.091 | 0.085 | 0.098 | 0.030 |

This is not publishable as a Sagear replication. The mechanical pipeline is healthier, but the inferred eccentricities are still much too high, especially thin singles.

### Source Split Diagnostic

Running the same Rayleigh fit separately on old archive rows and newly re-extracted ALDERAAN rows shows the problem is not only in the old archive:

| source | thick singles | thin singles | thick multis | thin multis |
|---|---:|---:|---:|---:|
| old archive only | 0.156 | 0.189 | 0.0004 boundary | 0.046 |
| new paired-impact extraction only | 0.358 | 0.390 | 0.143 | 0.155 |

The newly extracted ALDERAAN systems are even higher-e than the archive. That means the paired-impact correction alone is not the solution.

### Impact-Mode Comparison

Paired ALDERAAN impact samples reduce some extreme thin-single posterior medians, but not enough to recover Sagear:

| disk | system | N compared | geometric median e50 | paired ALDERAAN median e50 | high-e count geometric | high-e count paired |
|---|---|---:|---:|---:|---:|---:|
| thick | multi | 69 | 0.232 | 0.242 | 0 | 0 |
| thick | single | 106 | 0.316 | 0.298 | 13 | 12 |
| thin | multi | 222 | 0.250 | 0.256 | 5 | 4 |
| thin | single | 308 | 0.384 | 0.368 | 51 | 34 |

This confirms the earlier extractor default was wrong, but also confirms it was not the dominant source of the Sagear mismatch.

### QC Manifest And Leverage

The QC manifest has one row for each of the 2474 sample planets:

| disk | system | rows | primary-QC excluded |
|---|---|---:|---:|
| thick | multi | 233 | 12 |
| thick | single | 272 | 24 |
| thin | multi | 860 | 25 |
| thin | single | 1109 | 82 |

The most common review flag is `zeta_tail_review` (499 rows). There are 74 missing-after-merge rows.

Top-leverage removal shows the fits are still fragile:

- Thick singles: removing the top 5 leverage planets shifts `<e>` by more than 5%; removing the top 10 collapses the fit near the circular boundary.
- Thin singles: top 20 removal is still within 5%, but top 50 collapses near the circular boundary.
- Thick multis: removing the top 3 already collapses near the circular boundary.
- Thin multis: top 3 removal fails the 5% stability threshold.

This is a major scientific warning. The posterior mass being fed into the population model is controlled by zeta/eccentricity tails, not by a stable Sagear-like population signal.

### Current Best Diagnosis After Recovery

The recovery plan fixed real implementation errors, but it did not recover Sagear Table 2.

The current evidence points to a deeper mismatch in how the eccentricity posteriors are being constructed from ALDERAAN transit-shape results, or how the ALDERAAN shape parameters in these FITS should be interpreted relative to Sagear's actual post-model. The next useful work is not a blind full rerun. It is a targeted astrophysical/posterior audit:

1. Pick a small set of high-leverage planets from `rayleigh_per_planet_leverage_PAIRED_EXACT_qcprimary.csv`.
2. For each one, compare ALDERAAN `T14`, `Rp/Rs`, `b`, period, and stellar-density prior directly against the KOI catalog and Sagear assumptions.
3. Verify whether the ALDERAAN result FITS variables used here are the same transit-shape variables Sagear used.
4. Re-derive one planet by hand from the MacDougall/Sagear duration-ratio formalism to confirm the exact zeta-to-`e,omega` transformation and transit-selection correction.
5. Only after that, decide whether the existing 547 FITS can be salvaged or whether selected systems need to be refit with a stricter ALDERAAN configuration.

### Key Outputs From This Recovery Run

- `outputs/eccentricity_posterior_summary_sagear_missing_paired_exact.csv`
- `outputs/eccentricity_posterior_summary_merged_paired_exact.csv`
- `outputs/eccentricity_posterior_summary_merged_paired_exact_qcprimary.csv`
- `outputs/eccentricity_posterior_qc_manifest_paired_exact.csv`
- `outputs/eccentricity_posterior_coverage_merged_paired_exact.csv`
- `outputs/rayleigh_population_fit_transit_selection_PAIRED_EXACT_primary.csv`
- `outputs/rayleigh_population_fit_transit_selection_PAIRED_EXACT_qcprimary.csv`
- `outputs/rayleigh_per_planet_leverage_PAIRED_EXACT_qcprimary.csv`
- `outputs/rayleigh_top_leverage_out_PAIRED_EXACT_qcprimary.csv`
- `outputs/impact_mode_comparison_paired_exact_summary.csv`
- `outputs/formula_sanity_checks/formula_sanity_checks.json`

## 2026-07-07 High-Leverage And Sample-Quality Audit

After the paired-impact recovery run still failed to reproduce Sagear, I audited the planets with the highest leverage in the Rayleigh likelihood. This was meant to answer a very specific question:

Are the high eccentricities caused by our shortcut posterior extraction, or are they already implied by the ALDERAAN transit-shape posteriors and stellar-density priors?

### Direct Post-Model Importance Sampling Check

I added:

- `post_model_formalism_audit.py`

This script takes the top-leverage planets and recomputes eccentricity constraints using a direct MacDougall/Sagear-style post-model importance sampler. Instead of using the current zeta-KDE shortcut, it draws `(e, omega)` from the prior, combines them with paired ALDERAAN `(T14, Rp/Rs, b)` samples, computes the circular/eccentric stellar density implied by the transit shape, and weights each draw by the Berger stellar-density prior.

The audit output files are:

- `outputs/post_model_formalism_audit.csv`
- `outputs/post_model_formalism_audit_summary.csv`
- `outputs/post_model_formalism_audit_current_vs_direct.png`

Summary for the top 15 leverage planets per population:

| population | audited rows | direct-FITS rows | median current e50 | median direct-importance e50 | current minus direct |
|---|---:|---:|---:|---:|---:|
| thick multis | 15 | 8 | 0.278 | 0.335 | -0.013 |
| thick singles | 15 | 8 | 0.426 | 0.461 | -0.003 |
| thin multis | 15 | 9 | 0.410 | 0.433 | -0.010 |
| thin singles | 15 | 13 | 0.410 | 0.442 | -0.005 |

Interpretation:

- The direct post-model importance sampler does not reduce the high-e tail.
- In the audited high-leverage new-ALDERAAN planets, direct importance generally gives eccentricities as high as or higher than the zeta-KDE extractor.
- Therefore the main high-e problem is not the paired-impact fix or the zeta-KDE shortcut.
- The high-e tail is already implied by the ALDERAAN transit shapes when compared to the stellar-density prior.

### Circular Transit Density Mismatch

The direct audit also computes the stellar density implied by the ALDERAAN transit shape at `e=0`, then compares it to the Berger density prior.

Key result:

- 37/37 direct-high-e audited new-ALDERAAN rows have a large circular-transit-density mismatch against the Berger density prior.
- 34/60 audited rows are zeta-tail cases.
- 12/60 are grazing or super-grazing ALDERAAN fits.
- 4/60 have large `Rp/Rs` shifts.
- 3/60 have period-mismatch flags.

This is now the best physical diagnosis of the eccentricity problem:

The population fit is high because many individual planets have ALDERAAN transit shapes whose circular-orbit stellar density is inconsistent with the Berger stellar density. The eccentricity posterior is doing what the photoeccentric formalism tells it to do: it moves to high `e` to reconcile transit shape and stellar density.

That does not mean the high eccentricities are astrophysically real. It means the remaining discrepancy is probably in sample/stellar-density choice, ALDERAAN fit quality, or Sagear's visual/convergence vetting.

### Evolved-Star Sensitivity

I added:

- `dwarf_cut_diagnostics.py`

This surfaced an important upstream problem: the pipeline was using the Berger+2018 `Bin=0` flag but was not carrying or optionally filtering on Berger+2018 `Evol`, even though the downloaded Berger+2018 table includes it.

I patched:

- `diagnose_sample.py`
- `config.json`

The pipeline now carries `berger2018_evol` and `berger2018_rad` through the sample builder and supports an optional config key:

```json
"berger2018_evol_required": null
```

I left it as `null`, not `0`, because applying `Evol=0` changes the sample too aggressively unless we confirm Sagear's exact evolved-star rule.

Observed evolved-star issue in the current sample:

| diagnostic cut | sample planets | hosts | removed planets |
|---|---:|---:|---:|
| current no dwarf cut | 2474 | 1857 | 0 |
| Berger logg >= 4.0 | 2261 | 1667 | 213 |
| Berger radius <= 2 | 2333 | 1731 | 141 |
| rho linear >= 0.05 solar | 2432 | 1818 | 42 |
| Berger logg >= 4.0 and radius <= 2 | 2261 | 1667 | 213 |
| Berger+2018 Evol == 0 | 2046 | 1514 | 428 |

The current sample includes obvious evolved/giant cases that can force high photoeccentricities. Examples include `K04943.01`, `K05532.01`, `K01241.02` / Kepler-56 b, and other low-density stars.

However:

- A strict `Evol=0` cut drops the sample to 2046 planets, far below Sagear's 2465.
- A simpler `berger_logg >= 4.0 and berger_rad <= 2` cut drops the sample to 2261 planets, also far below Sagear.
- Therefore these cuts are scientifically important diagnostics, but I am not calling them canonical without confirming Sagear's exact host-star selection.

### Rayleigh Sensitivity To Sample-Quality Filters

I generated filtered posterior summaries and ran QC-primary Rayleigh fits for several sensitivity cuts. Outputs:

- `outputs/dwarf_cut_diagnostics.csv`
- `outputs/population_filter_sensitivity_counts.csv`
- `outputs/population_filter_sensitivity_rayleigh_results.csv`
- `outputs/population_filter_sensitivity_rayleigh_pivot.csv`

Rayleigh `<e>` sensitivity summary:

| filter | thick singles | thin singles | thick multis | thin multis |
|---|---:|---:|---:|---:|
| paired exact QC-primary | 0.230 | 0.256 | 0.071 | 0.091 |
| dwarf logg/radius QC-primary | 0.194 | 0.251 | 0.064 | 0.072 |
| small planets only | 0.155 | 0.207 | 0.068 | 0.072 |
| dwarf + small planets | 0.152 | 0.203 | 0.065 | 0.072 |
| confirmed only | 0.136 | 0.200 | 0.071 | 0.092 |
| confirmed or score >= 0.5 | 0.161 | 0.235 | 0.071 | 0.091 |
| dwarf + confirmed or score >= 0.5 | 0.158 | 0.228 | 0.064 | 0.071 |
| dwarf + small + confirmed/score >= 0.5 | 0.135 | 0.195 | 0.065 | 0.071 |

Interpretation:

- Evolved-star, giant-planet, and low-reliability candidate filters all move the results in the right direction.
- None of them gets close to Sagear's thin-single value of 0.022.
- The thin-single discrepancy remains the dominant failure.
- The current evidence says the remaining problem is not one isolated bug. It is a combination of density-prior/source mismatch, sample-quality/fit-quality vetting, and perhaps Sagear's exact ALDERAAN postprocessing choices.

### New Best Diagnosis

The replication is now blocked by one of the following, in priority order:

1. **Stellar-density prior mismatch.** The manuscript says the eccentricity inference uses Berger+2018 densities, while this pipeline is still using Berger+2020 `table2.dat.gz` density values. Berger+2018 is present locally for radii/evolution flags, but not as a complete density-prior table in the current pipeline.
2. **Unreproduced Sagear visual/convergence vetting.** Sagear says less than 2 percent of the sample was removed after visual inspection and convergence checks. Our high-leverage planets include many grazing, large-radius, or low-score candidate-like cases that could plausibly be removed by that vetting.
3. **Evolved-star handling is ambiguous.** The current sample contains evolved stars, but a strict Berger+2018 `Evol=0` cut makes the count too small. We need Sagear's actual host-star file or exact rule.
4. **ALDERAAN shape interpretation still needs target-level validation.** Direct importance sampling confirms the high-e tails are implied by ALDERAAN shapes plus density priors, so the next audit should inspect the actual light-curve fits for the top leverage systems, not just their posterior CSVs.

### Immediate Next Actions

The most valuable next tasks are:

1. Find or reconstruct Sagear's exact stellar-density prior table. Do not assume Berger+2020 is equivalent to Berger+2018.
2. Get the exact list of Sagear planets removed for ALDERAAN non-convergence/visual-quality reasons, or reproduce that vetting with a documented rule.
3. For the top thin-single leverage systems, inspect the ALDERAAN fit plots/light curves directly, especially systems with `koi_score=0`, very large `Rp`, grazing `b`, or density ratios far from 1.
4. Only after those checks decide whether a new cloud ALDERAAN run is scientifically worth it.

## 2026-07-08 Go/No-Go Assessment

### Direct Answer

This is not a dead end yet, but it is a dead end for blind local tinkering. The remaining gap is no longer explained by obvious code bugs, incomplete ALDERAAN coverage, or simple sample-count mismatch.

The current project can still become scientifically useful in two ways:

1. A faithful replication, if we obtain or reconstruct Sagear's exact stellar-density priors and ALDERAAN fit-quality rejection decisions.
2. A rigorous independent replication attempt, if those hidden choices remain unavailable, explicitly showing that the published result is not reproducible from the public method description plus available catalogs alone.

### Public And Local Source Checks

I checked the public arXiv manuscript, ALDERAAN repository, MacDougall et al. photoeccentric formalism reference, and the local source tarballs.

Findings:

- The public paper reports a sample of 2465 planets / candidates around 1888 stars, with 1121 thin singles, 378 thick singles, 862 thin multi planets, and 207 thick multi planets.
- The paper says eccentricities are constrained from ALDERAAN transit-shape posteriors plus stellar-density priors.
- The methods section says the stellar-density sample/prior comes from Berger et al. 2018.
- Elsewhere in the manuscript/conclusion and comments in the source mention Berger et al. 2020 density values.
- The paper states that each transit fit was visually inspected and that less than 2 percent of the sample was removed for non-convergent posteriors.
- The local Sagear source tarball contains manuscript source and figures, not the analysis code, posterior samples, density-prior table, or removed-target list.

That combination matters because the failure mode in our current run is exactly density-prior sensitive: many high-leverage ALDERAAN posteriors imply circular transit densities that strongly disagree with the stellar density adopted by our extractor. Changing sample filters helps, but not enough.

### What We Should Stop Doing

- Do not spend more cloud time rerunning the full ALDERAAN batch until the density-prior and fit-vetting issues are resolved.
- Do not treat the current Rayleigh numbers as publishable Sagear replication results.
- Do not silently remove high-e planets just to force agreement. Any veto needs a documented, reproducible rule or a match to Sagear's own rejection list.
- Do not assume Berger+2020 densities are interchangeable with the Berger+2018 density priors implied by the methods text.

### What Is Still Worth Doing

1. Obtain the missing artifacts directly from the authors or repository owners:
   - exact stellar-density prior table used for every KOI,
   - exact ALDERAAN results/posterior files if shareable,
   - target list removed by visual inspection/non-convergence,
   - population-fit code or likelihood normalization used for Table 2.

2. Reconstruct the density-prior branch locally:
   - build density priors from Berger+2018 radii plus available masses where possible,
   - compare Berger+2018-like densities against Berger+2020 densities for the top leverage systems,
   - rerun only the post-model eccentricity extraction from existing ALDERAAN FITS, not the photometry.

3. Audit the top leverage planets by hand:
   - inspect ALDERAAN transit fit plots/light curves,
   - record grazing, poor TTV, large-radius, low-score, and density-mismatch flags,
   - test whether a documented QC rule removes the high-e tail without destroying the published sample counts.

4. Only after those checks:
   - rerun ALDERAAN for the 45 failed targets if they are needed,
   - rerun the hierarchical fits across Rayleigh, beta, modified beta, and half-Gaussian models.

### Blunt Status

We are not out of ideas. We are out of high-confidence fixes that can be made using only the current local files.

The strongest current diagnosis is:

- The pipeline mechanics are now mostly sound.
- The broad sample counts are close enough that they are not the main cause.
- The disk classification is no longer the obvious dominant failure.
- The eccentricity posteriors are still not Sagear-equivalent.
- The mismatch appears dominated by the stellar-density prior / fit-quality-vetting layer, not by the population Rayleigh code.

So the honest decision is: keep going, but change tactics. The next move is not more brute-force ALDERAAN. The next move is to recover or reconstruct the exact density-prior and ALDERAAN rejection information.

## 2026-07-08 Additional Recovery Work

### Added Shared Berger+2018 Parser

I promoted the Berger+2018 VizieR parser into `common.py` as `read_berger2018_stellar_table()`, and updated `diagnose_sample.py` to use the shared parser instead of a local duplicate.

Reason: the project now uses Berger+2018 radii/evolution flags in several diagnostics, so the parser should not silently diverge between scripts.

### Density-Prior Reconstruction Test

New script:

- `density_prior_diagnostics.py`

New outputs:

- `outputs/density_prior_reconstruction_diagnostics_QC_PRIMARY.csv`
- `outputs/density_prior_reconstruction_summary_QC_PRIMARY.csv`
- `outputs/density_prior_top_leverage_QC_PRIMARY.csv`
- `outputs/density_prior_required_vs_b18like_QC_PRIMARY.png`

Test performed:

- Constructed a conservative "Berger+2018-like" density proxy using Berger+2018 radius and Berger+2020 mass:
  - `rho_b18like = M_Berger2020 / R_Berger2018^3`
- Compared that against the current Berger+2020 density table.
- For every posterior, estimated the density shift required to make the median transit shape circular:
  - `rho_required / rho_current = zeta_median^-3`

Key result:

| population | n | median required rho/current | median Berger2018-radius rho/current | fraction needing rho factor > 2 |
|---|---:|---:|---:|---:|
| thick singles | 248 | 1.095 | 0.996 | 0.302 |
| thin singles | 1027 | 1.185 | 1.011 | 0.294 |
| thick multis | 221 | 0.815 | 0.975 | 0.172 |
| thin multis | 835 | 0.867 | 0.993 | 0.184 |

For the top thin-single leverage planets, the required density shifts are much larger:

| KOI | e50 | zeta | rho needed/current | Berger2018-radius rho/current | notes |
|---|---:|---:|---:|---:|---|
| K00125.01 | 0.306 | 0.827 | 1.77 | 1.08 | confirmed, giant-radius planet |
| K01465.01 | 0.364 | 0.782 | 2.09 | 1.01 | candidate, score 0, grazing |
| K00856.01 | 0.378 | 0.773 | 2.16 | 0.99 | evolved flag 1, grazing, giant-radius |
| K00791.01 | 0.409 | 0.748 | 2.39 | 1.01 | confirmed, evolved flag 1 |
| K01455.01 | 0.426 | 0.736 | 2.51 | 0.87 | candidate, score 0, grazing |
| K03320.01 | 0.444 | 0.708 | 2.82 | 1.03 | candidate, score 0, giant-radius |

Conclusion:

- The Berger+2018-vs-Berger+2020 radius difference is not big enough to explain the high-e tail.
- Density-prior ambiguity is still important, but only if Sagear used a materially different stellar-density prior than either Berger+2020 table density or this Berger+2018-radius/Berger+2020-mass reconstruction.
- The problem is therefore not solved by simply swapping in Berger+2018 radii.

### Posterior Source-Mix Audit

New script:

- `posterior_source_diagnostics.py`

New outputs:

- `outputs/posterior_source_population_summary_QC_PRIMARY.csv`
- `outputs/posterior_source_population_counts_QC_PRIMARY.csv`
- `outputs/alderaan_needed_to_replace_archive_targets_QC_PRIMARY.csv`
- `outputs/posterior_source_high_leverage_mix_QC_PRIMARY.csv`

Key count:

| population | existing archive / geometric | new ALDERAAN / paired impact |
|---|---:|---:|
| thick singles | 149 | 99 |
| thin singles | 752 | 275 |
| thick multis | 156 | 65 |
| thin multis | 623 | 212 |

Interpretation:

- The merged QC-primary table is still a mixed-source table:
  - 1680 old archive rows use `impact_mode=geometric`.
  - 651 new rows use `impact_mode=alderaan`.
- The old archive `.npz` files do not contain raw paired `T14, b, Rp/Rs` ALDERAAN samples. They contain only final `e, omega` grids and summary fields.
- Therefore the old archive rows cannot be repaired into Sagear-equivalent paired-impact posteriors without raw ALDERAAN FITS or rerunning ALDERAAN.
- However, the top leverage thin-single failure is dominated by new ALDERAAN rows, not archive rows. Archive replacement alone will not fix the main discrepancy.

Source-only Rayleigh fits:

| source subset | thick singles | thin singles | thick multis | thin multis |
|---|---:|---:|---:|---:|
| new ALDERAAN only | 0.358 | 0.390 | 0.143 | 0.155 |
| old archive only | 0.156 | 0.189 | 0.0004 boundary | 0.046 |
| merged QC-primary | 0.230 | 0.256 | 0.071 | 0.091 |

Conclusion:

- The new ALDERAAN subset makes the eccentricity discrepancy worse, not better.
- The archive mixture is still a method-equivalence problem, but it is not the dominant explanation for why thin singles are too eccentric.

### QC / Visual-Vetting Sensitivity Test

New script:

- `qc_vetting_sensitivity.py`

New outputs:

- `outputs/qc_vetting_sensitivity_counts_QC_VETTING.csv`
- `outputs/qc_vetting_sensitivity_rayleigh_QC_VETTING.csv`
- `outputs/qc_vetting_sensitivity_rayleigh_pivot_QC_VETTING.csv`
- filtered posterior summaries named `outputs/eccentricity_posterior_summary_filter_*.csv`

Purpose:

- Test whether plausible transparent vetting rules could recover Sagear-like Rayleigh values.
- In particular, test the manuscript statement that visual/non-convergence removals were `<2%`.

Key results:

| rule | thick singles | thin singles | thick multis | thin multis | comment |
|---|---:|---:|---:|---:|---|
| base QC-primary | 0.230 | 0.256 | 0.071 | 0.091 | current best merged run |
| drop top 2% leverage per population | 0.209 | 0.248 | 0.0003 boundary | 0.0002 boundary | Sagear-scale removal does not fix thin singles |
| no e84 > 0.9 review | 0.193 | 0.243 | 0.070 | 0.087 | small improvement only |
| no grazing b >= 0.9 | 0.223 | 0.226 | 0.0017 | 0.085 | not enough for thin singles |
| confirmed or score >= 0.5 | 0.161 | 0.235 | 0.071 | 0.091 | not enough |
| planet radius < 4 R_Earth | 0.158 | 0.206 | 0.067 | 0.074 | not enough |
| no zeta-tail review | 0.122 | 0.116 | 0.0004 boundary | 0.0011 | removes 18-29% of planets, far beyond Sagear's stated <2% |
| transparent strict fit QC | 0.111 | 0.114 | 0.0004 boundary | 0.0011 | removes 21-36% of planets |
| transparent strict astro QC | 0.0004 boundary | 0.084 | 0.0004 boundary | 0.0011 | removes 23-42% of planets |

Conclusion:

- Sagear-scale `<2%` visual/convergence removal cannot explain the thin-single discrepancy.
- Transparent aggressive QC can reduce thin singles, but only by removing roughly one-third of the sample, and it often collapses other populations to the sigma-grid floor.
- Therefore "Sagear visually removed a few bad ALDERAAN fits" is not a sufficient explanation.
- The high-e tail is either:
  1. a deeper difference in the eccentricity/post-model formalism,
  2. a deeper difference in the transit-shape inputs/posteriors,
  3. a hidden sample/planet-selection difference much larger than described,
  4. or a problem in the way our ALDERAAN-derived `zeta` is converted into population likelihoods.

### Updated Diagnosis After These Tests

The project is not at a dead end, but the failure is now narrower:

- Berger+2018 radius substitution does not fix the density issue.
- Archive-vs-new posterior mixing is real but does not explain the high-e tail.
- Sagear-scale visual vetting does not fix thin singles.
- Aggressive QC can lower eccentricities, but only by removing too much data and creating boundary artifacts.

The next scientifically meaningful branch is to audit the **post-model eccentricity formalism itself** against the MacDougall/Sagear equations and ALDERAAN's actual fitted parameter definitions, especially:

1. whether `DUR14` is exactly first-to-fourth contact duration in days for every ALDERAAN output,
2. whether the current `zeta = T14 / T14_circ` convention is inverted relative to Sagear/MacDougall,
3. whether the population fit is using posterior samples/evidence in a way that double-counts or omits the single-planet posterior prior,
4. whether Sagear's hierarchical likelihood used the full `e, omega` posterior grids or a different summary/integration convention.

### ALDERAAN Parameter Definition Check

I checked the local ALDERAAN source under `external/alderaan`.

Relevant findings:

- `combine_umbrella_samples.py` labels `DUR14` as:
  - unit: `days`
  - description: `Transit duration: 1st-4th contact`
- `alderaan/astro.py` separately describes the total transit duration as contact I-IV.
- The Sagear manuscript also describes the fitted shape parameter as full first-to-fourth contact duration `T14`.

Conclusion:

- There is no evidence that our use of `DUR14` as first-to-fourth contact duration in days is wrong.
- The high-e tail is not explained by an hour/day unit mistake or by using ingress/egress duration instead of total duration.

### Transit-Selection Formalism Audit

New script:

- `selection_correction_audit.py`

Patched script:

- `hierarchical_rayleigh.py`

New outputs:

- `outputs/selection_correction_rayleigh_comparison_QC_PRIMARY.csv`
- `outputs/selection_correction_rayleigh_pivot_QC_PRIMARY.csv`
- `outputs/rayleigh_population_fit_transit_selection_manuscript_reciprocal_PAIRED_EXACT_qcprimary_MANUSCRIPT_RECIPROCAL.csv`

Reason:

- The manuscript text says the hierarchical likelihood accounts for non-uniform transit probability.
- The commented equation in `main.tex` explicitly applies the reciprocal transit-probability factor:
  - `(1-e^2)/(1 + e sin omega)`
- The older pipeline used a forward transit-probability weighting plus a population normalizer:
  - `(1 + e sin omega)/(1-e^2)`, followed by a normalization.

I added explicit selection modes:

- `legacy_forward_norm`: old pipeline behavior.
- `manuscript_reciprocal`: direct implementation of the commented Sagear HBM equation.
- `manuscript_reciprocal_with_norm`: diagnostic reciprocal mode with an additional population normalizer.
- `none`: no transit-selection correction.

Results:

| selection mode | thick singles | thin singles | thick multis | thin multis |
|---|---:|---:|---:|---:|
| legacy forward + normalization | 0.230 | 0.256 | 0.071 | 0.091 |
| manuscript reciprocal | 0.178 | 0.217 | 0.068 | 0.083 |
| manuscript reciprocal + normalization | 0.197 | 0.240 | 0.070 | 0.087 |
| no selection correction | 0.225 | 0.264 | 0.071 | 0.090 |

Conclusion:

- This was a real method mismatch in the prior pipeline.
- The manuscript reciprocal correction moves the results in the right direction.
- It does **not** solve the replication discrepancy:
  - thin singles remain at `0.217`, still far above Sagear's `0.022`;
  - thick singles remain at `0.178`, still far above Sagear's `0.066`.
- Future Sagear-comparison Rayleigh runs should report the manuscript reciprocal mode, not only the older legacy-forward mode.

### Updated Diagnosis After Selection Audit

We have now eliminated or weakened several explanations:

- Not a `DUR14` unit/definition error.
- Not fixed by Berger+2018 radius substitution.
- Not fixed by replacing old archive rows alone.
- Not fixed by Sagear-scale `<2%` visual vetting.
- Not fixed by switching to the manuscript reciprocal transit-selection correction, although that is a real improvement.

The remaining plausible root causes are now narrower and more uncomfortable:

1. Our individual eccentricity posteriors are fundamentally broader/higher-e than Sagear's, even for new ALDERAAN rows.
2. The ALDERAAN run configuration/catalog inputs may differ from Sagear's in a way that changes `T14`, `b`, or `Rp/Rs`.
3. Sagear's unpublished posterior extraction may include additional priors, clipping, convergence checks, or target-specific fit vetos not captured in the public manuscript.
4. The old archive rows are not method-equivalent and still need raw ALDERAAN FITS or reruns, but they are not the whole problem.

The next useful check is therefore not another population-fit tweak. It is a target-level ALDERAAN-shape audit:

- compare ALDERAAN `DUR14`, `b`, and `Rp/Rs` against KOI catalog values for the highest-leverage new ALDERAAN systems;
- identify whether the high-e tail comes from real ALDERAAN shape shifts, grazing solutions, giant-planet candidates, or fit pathologies;
- inspect whether the ALDERAAN project used the same input catalog and limb-darkening assumptions Sagear likely used.

### ALDERAAN Shape Audit Results

Updated script:

- `alderaan_shape_diagnostics.py`

New outputs:

- `outputs/alderaan_shape_diagnostics_all_current_matched.csv`
- `outputs/alderaan_shape_diagnostics_flagged_all_current_matched.csv`
- `outputs/alderaan_shape_diagnostics_top120_leverage_current.csv`
- `outputs/alderaan_shape_diagnostics_flagged_top120_leverage_current.csv`
- `outputs/alderaan_shape_diagnostics_top120_leverage_current.md`
- `outputs/alderaan_shape_diagnostics_top120_leverage_current.png`

Reason:

- The population eccentricities are high because many individual posteriors prefer high eccentricity.
- We needed to distinguish whether this is caused by ALDERAAN transit-shape drift relative to the KOI catalog, impact/grazing pathologies, or a coherent photoeccentric density mismatch.

All matched new-ALDERAAN rows:

| population | n | median zeta | zeta 16-84 | median e50 | median ALDERAAN/KOI duration |
|---|---:|---:|---:|---:|---:|
| thick multis | 65 | 1.089 | 0.836-1.289 | 0.239 | 1.002 |
| thick singles | 99 | 0.923 | 0.590-1.166 | 0.293 | 1.005 |
| thin multis | 212 | 1.055 | 0.742-1.279 | 0.252 | 1.009 |
| thin singles | 275 | 0.920 | 0.512-1.271 | 0.338 | 1.000 |

Top-leverage subset:

| population | n | median zeta | median e50 | median ALDERAAN/KOI duration |
|---|---:|---:|---:|---:|
| thick multis | 4 | 2.300 | 0.515 | 1.085 |
| thick singles | 8 | 0.708 | 0.456 | 1.005 |
| thin multis | 13 | 1.253 | 0.421 | 1.015 |
| thin singles | 55 | 0.638 | 0.526 | 0.998 |

Flag rates among flagged top-leverage rows:

| flag | count/fraction |
|---|---:|
| short zeta | 39/69 = 56.5% |
| long zeta | 20/69 = 29.0% |
| duration shift >25% | 6/69 = 8.7% |
| radius-ratio shift >25% | 17/69 = 24.6% |
| high impact | 28/69 = 40.6% |
| large stellar-density uncertainty | 1/69 = 1.4% |
| high e50 | 37/69 = 53.6% |

Conclusion:

- The high-e tail is not mostly caused by ALDERAAN duration medians drifting away from KOI catalog durations. The ALDERAAN/KOI duration ratio is near unity for the full current matched sample and also for high-leverage thin singles.
- The problematic thin-single tail is instead a coherent `zeta` tail: the fitted transit shape, when interpreted under the current stellar-density prior, implies a circular density that is too different from the catalog stellar density.
- This points back to either:
  1. the stellar-density prior,
  2. the ALDERAAN shape fit configuration/catalog priors,
  3. or an unpublished target-level vetting/selection difference.

### Limb-Darkening Catalog Audit

New script:

- `limb_darkening_audit.py`

New outputs:

- `outputs/limb_darkening_catalog_offsets_QC_PRIMARY.csv`
- `outputs/limb_darkening_catalog_offset_summary_QC_PRIMARY.csv`
- `outputs/limb_darkening_population_context_QC_PRIMARY.csv`
- `outputs/limb_darkening_population_summary_QC_PRIMARY.csv`
- `outputs/limb_darkening_top_leverage_QC_PRIMARY.csv`

Reason:

- Sagear states that ALDERAAN used system-level quadratic limb-darkening priors `(u1, u2) ~ Normal(mu, 0.1)`, with `mu` derived from Gaia/Berger stellar parameters and stellar atmosphere models.
- The local ALDERAAN source confirms the code applies a Gaussian penalty with `sigma=0.1` around the catalog `limbdark_1/limbdark_2` values.
- Our batch catalogs used KOI catalog `koi_ldm_coeff1/2` values. I compared them to ALDERAAN's bundled `kepler_dr25_gaia_dr2_crossmatch.csv`, which appears to be the closer Gaia/Berger-style reference catalog.

Catalog-level offsets against the ALDERAAN bundled reference:

| catalog | matched rows | median du1 | median du2 | frac max(|du|)>0.05 | frac max(|du|)>0.10 |
|---|---:|---:|---:|---:|---:|
| initial cloud batch | 277 | -0.103 | +0.120 | 0.906 | 0.650 |
| missing fixed batch | 578 | -0.128 | +0.142 | 0.964 | 0.858 |
| missing original batch | 578 | -0.128 | +0.141 | 0.958 | 0.853 |
| local needed catalog | 1036 | -0.107 | +0.125 | 0.931 | 0.744 |

For new ALDERAAN rows in the fixed missing batch:

| population | n | median e50 | median zeta | median max(|du|) | frac max(|du|)>0.10 |
|---|---:|---:|---:|---:|---:|
| thick multis | 60 | 0.230 | 1.086 | 0.142 | 0.800 |
| thick singles | 78 | 0.265 | 1.008 | 0.142 | 0.872 |
| thin multis | 176 | 0.238 | 1.051 | 0.137 | 0.813 |
| thin singles | 203 | 0.306 | 0.948 | 0.149 | 0.911 |

Top-100 leverage subset:

| population | n | median max(|du|) | frac max(|du|)>0.10 | median e50 | median zeta |
|---|---:|---:|---:|---:|---:|
| thick multis | 4 | 0.087 | 0.500 | 0.409 | 1.558 |
| thick singles | 8 | 0.143 | 0.750 | 0.445 | 0.708 |
| thin multis | 18 | 0.127 | 0.611 | 0.387 | 1.274 |
| thin singles | 70 | 0.146 | 0.800 | 0.593 | 0.587 |

Conclusion:

- This is the strongest remaining non-population-fit lead.
- It is not proof that limb darkening caused the discrepancy, because limb darkening is partially fit and high-e is not perfectly correlated with `du`.
- But it is a real Sagear-equivalence problem: our ALDERAAN run almost certainly did not use the same limb-darkening prior centers described in the manuscript.
- The next efficient experiment is not a full rerun. It is a validation rerun of roughly 10-20 high-leverage systems using Gaia/Berger/LDTk-style limb-darkening centers, then comparing `T14`, `b`, `Rp/Rs`, `zeta`, and `e50` against the existing runs.

### Blunt Current Status

No, this is not where the project gives up.

But yes, the cheap local fixes are mostly exhausted. We have already tested and/or corrected:

- sample counts and disk classification,
- Berger+2018 binary/radius/evolution diagnostics,
- paired ALDERAAN impact samples,
- exact duration formula,
- zeta support and QC manifests,
- transit-selection weighting,
- source mixing between old archive and new ALDERAAN rows,
- aggressive QC and outlier sensitivity,
- ALDERAAN `DUR14` units/definition,
- ALDERAAN-vs-KOI duration shifts.

The remaining credible blockers are now upstream of the population fit:

1. **ALDERAAN input-prior mismatch**, especially limb darkening.
2. **Exact stellar-density prior table**, because the manuscript is inconsistent/ambiguous between Berger+2018 and Berger+2020 language.
3. **Unpublished Sagear target-level vetting**, because removing only a few bad fits cannot reproduce the published thin-single value, but a hidden target-selection difference still could.
4. **Access to Sagear's actual per-planet posterior table**, if available from the authors, because that would immediately show whether the discrepancy is in our ALDERAAN transit fits or our population layer.

Decision:

- We are not at a scientific dead end.
- We are at the point where further progress requires either targeted reruns with corrected ALDERAAN priors, exact external data from Sagear/authors, or both.
- Continuing to tune the hierarchical population fit without changing/checking individual posteriors would be misleading.

### Cloud Bundle Reproducibility Fix

Patched files:

- `cloud_missing_batch/run_batch.sh`
- `cloud_missing_batch/validate_bundle.py`

Problem:

- The fixed missing-target catalog existed as `sagear_missing_catalog_FIXED.csv`, but `run_batch.sh` still copied and validated the old `sagear_missing_catalog.csv` by default.
- The old catalog has 31 systems with inconsistent system-level `limbdark_1` and 31 systems with inconsistent system-level `limbdark_2`, which ALDERAAN can reject or handle differently depending on live patches.

Fix:

- `run_batch.sh` now prefers `sagear_missing_catalog_FIXED.csv` when present.
- It still copies the chosen source into the ALDERAAN project as `sagear_missing_catalog.csv`, so downstream scripts do not need to change.
- `validate_bundle.py` now fails fast if any target has more than one unique `limbdark_1` or `limbdark_2` value across planets in the same system.

Validation:

- Fixed catalog:
  - `VALIDATION OK: 592 targets, 767 catalog rows`
- Original catalog:
  - fails with `31 target(s) have inconsistent system-level limbdark_1`
  - fails with `31 target(s) have inconsistent system-level limbdark_2`

This does not repair already-finished FITS files, but it prevents the old catalog bug from contaminating future reruns.

### Limb-Darkening Validation Rerun Bundle

New scripts:

- `build_limb_darkening_validation_batch.py`
- `compare_limb_darkening_validation.py`

New bundle folder:

- `cloud_ld_validation_batch/`

Portable bundle zips:

- `<downloads-dir>\Hildale_LD_Validation_Batch_20260708.zip`
- `<validation-artifact-dir>\Hildale_LD_Validation_Batch_20260708.zip`

New bundle files:

- `targets_ld_reference_validation.csv`
- `sagear_ld_reference_catalog.csv`
- `ld_reference_validation_selection.csv`
- `run_ld_validation.sh`
- `pack_ld_validation_results.sh`
- `README_LD_VALIDATION.md`

Reason:

- The limb-darkening audit found that the existing ALDERAAN input prior centers differ strongly from the ALDERAAN bundled Kepler-Gaia reference values.
- Rather than launching another large cloud job, I built a small A/B validation batch.
- The validation reruns the same systems under a new run id, `sagear_ld_reference_validation`, using reference limb-darkening centers and preserving the existing `sagear_missing` FITS for direct comparison.

Important selection guard:

- The ALDERAAN bundled reference catalog contains 88 rows, about 2.3%, with exact `u1 = 0.1`, `u2 = 0.1`.
- Those look like fallback/sentinel values, not robust stellar-atmosphere priors.
- The validation builder now excludes those sentinel rows by default.
- The final validation catalog has 20 targets, 25 planet rows, and zero `0.1/0.1` sentinel rows.

Selected targets:

| priority | target | driver | population | reason | e50 | zeta | max LD offset |
|---:|---|---|---|---|---:|---:|---:|
| 1 | K00856 | K00856.01 | thin singles | high-leverage short-zeta | 0.378 | 0.773 | 0.170 |
| 2 | K00791 | K00791.01 | thin singles | high-leverage short-zeta | 0.409 | 0.748 | 0.146 |
| 3 | K01553 | K01553.01 | thin singles | high-leverage short-zeta | 0.473 | 0.682 | 0.161 |
| 4 | K00716 | K00716.01 | thin singles | high-leverage short-zeta | 0.526 | 0.638 | 0.161 |
| 5 | K01299 | K01299.01 | thin singles | high-leverage short-zeta | 0.558 | 0.607 | 0.101 |
| 6 | K00815 | K00815.01 | thin singles | high-leverage short-zeta | 0.598 | 0.577 | 0.113 |
| 7 | K01787 | K01787.01 | thin singles | high-leverage short-zeta | 0.617 | 0.558 | 0.152 |
| 8 | K00846 | K00846.01 | thin singles | high-leverage short-zeta | 0.635 | 0.540 | 0.160 |
| 9 | K00890 | K00890.01 | thick singles | high-leverage short-zeta | 0.438 | 0.716 | 0.167 |
| 10 | K00064 | K00064.01 | thick singles | high-leverage short-zeta | 0.452 | 0.701 | 0.107 |
| 11 | K02109 | K02109.01 | thick singles | high-leverage short-zeta | 0.700 | 0.462 | 0.185 |
| 12 | K01001 | K01001.01 | thin multis | high-leverage multi | 0.353 | 1.280 | 0.182 |
| 13 | K02714 | K02714.01 | thin multis | high-leverage multi | 0.420 | 1.436 | 0.141 |
| 14 | K00283 | K00283.01 | thin multis | high-leverage multi | 0.493 | 0.651 | 0.133 |
| 15 | K02533 | K02533.01 | thick multis | high-leverage multi | 0.457 | 3.037 | 0.100 |
| 16 | K00680 | K00680.01 | thin singles | long-zeta control | 0.380 | 1.323 | 0.164 |
| 17 | K02712 | K02712.01 | thin singles | long-zeta control | 0.859 | 3.351 | 0.174 |
| 18 | K00428 | K00428.01 | thin singles | stable LD-offset control | 0.154 | 0.947 | 0.227 |
| 19 | K04382 | K04382.01 | thin multis | stable LD-offset control | 0.191 | 1.068 | 0.193 |
| 20 | K00319 | K00319.01 | thin singles | stable LD-offset control | 0.107 | 0.977 | 0.187 |

Validation performed locally:

- Bundle validator:
  - `VALIDATION OK: 20 targets, 25 catalog rows`
- Sentinel check:
  - final validation catalog has `0` rows with exact `u1=u2=0.1`
- Python compile checks:
  - `build_limb_darkening_validation_batch.py` passed
  - `compare_limb_darkening_validation.py` passed
- Bash syntax check:
  - not run locally because this Windows environment does not have `bash` installed.

Cloud command:

```bash
cd ~/sagear_ld_validation_batch
JOBS=4 nohup bash run_ld_validation.sh > ld_validation_run.log 2>&1 &
```

After cloud results are copied back, run:

```bash
python extract_eccentricity_posteriors.py --sample outputs/canonical_sample_old_astropy_rawcc.csv --run-id sagear_ld_reference_validation --impact-mode alderaan --posterior-subdir eccentricity_posteriors_ld_reference_validation --summary-out outputs/eccentricity_posterior_summary_ld_reference_validation.csv --coverage-out outputs/eccentricity_posterior_coverage_ld_reference_validation.csv --excluded-out outputs/eccentricity_posterior_excluded_ld_reference_validation.csv
python compare_limb_darkening_validation.py
```

Interpretation gate:

- If the high-leverage systems move toward `zeta ~ 1` and lower `e50`, then the previous ALDERAAN limb-darkening prior centers were a material non-Sagear mismatch.
- If they do not move, limb darkening is probably not the main culprit, and the next branch is exact stellar-density priors or Sagear author data/posteriors.

## 2026-07-12: Factorial Batch Reliability Repair

The cloud factorial bundle was subsequently expanded to the current
24-target/34-catalog-row design (82 target fits across six arms). This
supersedes the older 20-target/25-row description above.

The first cloud attempt produced valid FITS files but exposed a runner defect:

- `summarize_progress.sh` counted only status files already present, so
  never-started rows were omitted and an arm could print `Batch complete` with
  fewer results than requested.
- `run_one_target.sh` marked a target `running` but had no signal/exit trap;
  a celerite2 numerical exception or VM SIGHUP could therefore leave an
  indistinguishable stale status.
- GNU Parallel received SIGHUP when the VM's original 24-hour maximum run
  duration was reached. That interruption is operational, not a scientific
  result.

The repaired bundle now:

- derives the expected target set from the CSV and writes a complete
  per-target manifest;
- accepts an arm as complete only when every requested target has a nonempty
  `*-results.fits` file;
- records `never_started`, `stale_running`, `interrupted_*`, `failed_exit_*`,
  missing-light-curve, and missing-result outcomes separately;
- uses target-specific append-only stdout/stderr logs;
- resumes only targets without a valid FITS result, preserving successful
  output; and
- continues the remaining factorial arms while returning a nonzero matrix
  status if any arm is incomplete.

The repair has passed local Bash syntax checks and an offline accounting test
covering completed, stale-running, and never-started targets. The repaired
bundle must be deployed before interpreting the factorial comparison.
# 2026-07-13: Published AJ Record Audit

The final journal article (AJ 172, 42; DOI 10.3847/1538-3881/ae71bf) and its machine-readable tables became available after the original replication work. This materially improves what can be reproduced without inference or reverse engineering.

## Published classification ground truth

The journal machine-readable Table 1 contains exactly 1,888 unique KIC hosts with Gaia DR3 identifiers, inferred and measured Galactocentric cylindrical velocities, `Pthick`, and the final disk assignment. Its internally consistent counts are 1,515 thin hosts and 373 thick hosts; 585 hosts have measured velocities and 1,303 use inferred velocities. Every label is exactly reproduced by `Pthick > 0.5`.

The prose statement of 378 thick hosts is a publication typo: 1,515 + 378 is inconsistent with 1,888, whereas 1,515 + 373 = 1,888. The subgroup counts also require 373 thick hosts: 275 thick singles + 98 thick multi hosts.

The live pipeline now parses this table with a strict count/threshold contract in `common.read_sagear2026_kinematic_hosts`. `published_sagear_audit.py` writes an authoritative host CSV, host reconciliation tables, publication-relabeled planet tables, a count summary, and a Toomre diagram made directly from Sagear's published velocities and labels.

## Classification discrepancy found

For `canonical_sample_strict_raw_no_bin.csv`, 1,764 hosts overlap the published table. Only 84.58% of overlapping labels agree; 272 hosts disagree. There are also 206 hosts only in our sample and 124 only in the published host table. This proves that the reconstructed classifier was not Sagear-equivalent and was a major upstream source of the incorrect Toomre diagram.

After applying published labels to the overlapping broad sample while preserving
the canonical pre-cut KOI architecture, the retained counts are 1,052 thin
singles, 803 thin multi planets, 269 thick singles, and 205 thick multi planets.
Sagear reports 1,121, 862, 275, and 207, respectively. The previous overlap
audit incorrectly recounted planets after intersection and reported
1,082/773/278/196; 39 planet rows disagreed with their pre-cut architecture.
With that diagnostic fixed, both thick populations are close, while the thin
deficits remain a planet-sample/crossmatch problem rather than solely a disk
classifier problem.

## Eccentricity-method audit against final publication

The final paper cites ALDERAAN v0.1.0 (Zenodo DOI 10.5281/zenodo.19208536). The cited tag commit differs from our pinned ALDERAAN commit only by a license file, so the completed cloud fits do not need to be discarded on version grounds.

The published Sagear et al. 2026 M-dwarf methods paper confirms the core direct extractor choices: paired posterior samples of `T14`, `Rp/Rs`, and `b`; uniform proposals over `e in [0,0.95]` and `omega in [-pi/2,3pi/2]`; the exact MacDougall density equation; stellar-density likelihood importance weights; and reciprocal geometric transit-probability correction in the hierarchical likelihood. These match the current direct extractor and `manuscript_reciprocal` population mode.

Remaining non-equivalences or unsafe defaults:

- The final disk paper explicitly says Berger et al. (2018) stellar densities, but Berger 2018 publishes radii rather than the homogeneous log-density column consumed by our code. Our current `table2.dat.gz` is Berger et al. (2020). This is unresolved publication ambiguity and must be tested or clarified with the authors, not silently guessed.
- Missing density errors currently fall back to an undocumented 13% uncertainty. Canonical replication should exclude these rows or label them as sensitivity-only.
- Low importance-ESS rows are flagged but included unless `--exclude-qc-primary` is supplied. Canonical runs must exclude them.
- The direct extractor and hierarchical fitter have different default summary filenames, creating a risk of fitting stale products.
- The deterministic sigma-grid posterior correctly implements a uniform `sigma_R in (0,1]` prior but is not the paper's stated NumPyro two-chain/1,000-step implementation with `Rhat < 1.05`. It is an excellent cross-check, not yet an exact computational reproduction.
- The reported untruncated Rayleigh mean is inconsistent with the grid-truncated model for broad/boundary solutions. This is negligible near Sagear's small fitted sigmas but must be corrected for diagnostics.

## New canonical artifacts

- `outputs/sagear2026_published_kinematic_hosts.csv`
- `outputs/sagear2026_publication_audit_summary.csv`
- `outputs/sagear2026_host_reconciliation_*.csv`
- `outputs/sagear2026_published_relabel_*.csv`
- `outputs/toomre_sagear2026_published_truth.png`

These published-label products supersede the reconstructed GMM labels for the primary replication. The reconstructed classifier remains useful only as a reproducibility and sensitivity diagnostic.
