# Deep Audit Findings  -  2026-07-04

Expert-level audit of the Sagear et al. (arXiv:2509.23973) replication pipeline. Read-only
investigation. **No script fixes applied**; every fix below is a specific proposal for a human to
review and apply on a future retry/full run. The live GCP run (`alderaan-full-e2-32`) was not
touched. Builds on `gap_hunt_findings_2026-07-04.md` and the two bugs (6, 9) fixed live today.

---

## Task A  -  Bug-6-shaped (missing required argparse args) and Bug-9-shaped (per-row system property) issues

### A.1 Required-argument matrix (verified by reading each `bin/*.py` argparse block)

| Script | required=True args |
|---|---|
| `detrend_and_estimate_ttvs.py` | mission, target, project_dir, **data_dir**, **catalog**, run_id |
| `analyze_autocorrelated_noise.py` | mission, target, run_id, project_dir, **data_dir**, **catalog** |
| `fit_transit_shape_simultaneous_nested.py` | mission, target, run_id, project_dir |
| `get_kepler_data.py` | positional `kepids` only (all flags have defaults) |

`analyze_autocorrelated_noise.py` genuinely uses these: `args.data_dir`/`args.catalog` (lines
118–119), `io.parse_catalog(...)` (line 187). Its top-level `except SystemExit: warnings.warn(...)`
(lines 124–125) swallows the argparse error  -  the exact bug-6 failure mode (silent, then a
confusing downstream `NameError` because MISSION/TARGET never get assigned).

### A.2 Findings

| # | Script / function | Issue | Severity / scope | Fixable? | Proposed fix (line-level) |
|---|---|---|---|---|---|
| A1 | **`cloud/run_one_target.sh`** (live) | All 3 invocations pass every required arg. **CLEAN**  -  the two prior fixes closed it. | none |  -  | none |
| A2 | **`scripts/alderaan_batch.py::write_run_script` line 259** | `analyze_autocorrelated_noise.py` invoked WITHOUT `--data_dir` and `--catalog` (both required=True). Bug-6 recurrence. | High if this generator is ever used (currently drives the validation-target run script only); it will fail exactly like the K00179 smoke-test crash | Yes, trivial | Line 259: append `--data_dir "$ProjectDir\Data\" --catalog sagear_validation_catalog.csv` (mirror the detrend line 258 args) |
| A3 | **`scripts/prepare_gcp_missing_alderaan_bundle.py` line 200** | Generated `run_one_target.sh` writes `analyze_autocorrelated_noise.py` WITHOUT `--data_dir`/`--catalog`. | High if bundle is regenerated from scratch (the live bundle was hand-fixed, so live run is safe); zero for current run | Yes, trivial | Line 200: append `--data_dir "$DATA_DIR/" --catalog "$CATALOG_NAME"` |
| A4 | **`scripts/cloud_prepare.py` line 179** | Same omission in its generated `run_one_target.sh`. | High if regenerated; zero for live run | Yes, trivial | Line 179: append `--data_dir "$DATA_DIR/" --catalog "$CATALOG_NAME"` |
| A5 | **`scripts/build_alderaan_needed_manifest.py` lines 569 and 694** | Both the PowerShell (`*>> $log`) and bash generated runners omit `--data_dir`/`--catalog` on `analyze_autocorrelated_noise.py`. | High if regenerated | Yes, trivial | Line 569: add `--data_dir "$ProjectDir\Data\" --catalog $CatalogName` before `*>> $log`. Line 694: add `--data_dir "$PROJECT_DIR/Data/" --catalog "$CATALOG_NAME"` before `>> "$log"` |

**Root-cause note:** this omission is *systemic*  -  it appears in every generator EXCEPT the one that
was manually patched (`cloud/run_one_target.sh`). Anyone who regenerates a bundle reintroduces
bug 6. Same class as the curl→wget regression (F2). Recommend fixing all four generators together.

`fit_transit_shape_simultaneous_nested.py` (only mission/target/run_id/project_dir required) is
correctly invoked everywhere. `get_kepler_data.py` needs only the positional kepid  -  clean everywhere.

### A.3 Bug-9-shaped audit (system-wide catalog columns)

`io.parse_catalog()` enforces cross-planet consistency on exactly four columns: **kic_id, npl,
limbdark_1, limbdark_2** (io.py lines 46–60), plus a per-row NaN check on period/epoch/depth/
duration/impact (lines 63–66). Cross-check vs `build_alderaan_catalog()` (alderaan_batch.py 152–201):

| Column | System-wide? | Computed once/system? | Verdict |
|---|---|---|---|
| kic_id | yes | yes (from `p["kepid"]`, constant within a `koi_target` group) | OK |
| npl | yes | yes (`len(grp)`, line 161) | OK |
| limbdark_1/2 | yes (stellar) | **yes, post bug-9 fix** (lines 171–182 compute once per group) | OK  -  fix is correct |
| period/epoch/depth/duration | per-planet | per-row | correct (not consistency-checked) |
| impact | per-planet | per-row, defaulted to 0.5 if NaN (line 196) | correct; also pre-empts the parse_catalog NaN-on-impact crash |

**No further bug-9 columns.** All four ALDERAAN-checked columns are correctly system-scoped.

**Defense-in-depth gap (A6, low severity):** the bundle validator
(`prepare_gcp_missing_alderaan_bundle.py` lines 345–358) checks column presence, NaNs, positivity,
and **npl** consistency, but does **NOT** validate `limbdark_1`/`limbdark_2` consistency across a
system  -  the very thing that crashed the live run. If bug 9 ever regressed via a different catalog
path, the validator would pass it through to a mid-run ValueError. Proposed fix: after line 358, add
a per-`koi_id` check that `limbdark_1`/`limbdark_2` each have `nunique() == 1`.

---

## Task B  -  the two new live-run errors: input-data pathology vs our bug

Both errors originate in **`detrend_and_estimate_ttvs.py`** (the FIRST pipeline stage; it is the
only stage using PyMC3  -  the transit fit uses dynesty, not PyMC3). So all 3 failing targets died in
detrending, before the transit-shape fit.

### B.1 `AttributeError: ... Normal has no finite default value ... Pass testval` (1 occurrence)

**Mechanism.** PyMC3 `Normal`s are built with data-derived `mu` that can go non-finite:
- line 775 `log_dur = pm.Normal("log_dur", mu=np.log(durs), ...)`
- line 809 / 1197 `log_jit = pm.Normal(..., mu=np.log(vbq), ...)` (vbq = per-quarter variance)

If `durs` or `vbq` contains a 0 (or negative), `np.log(...) → -inf/nan → "no finite default"`.

**Verdict: INPUT-DATA / pathological target, NOT our bug.** Our catalog cannot supply a zero
duration: `build_alderaan_catalog` filters `koi_duration > 0` and `koi_depth > 0` (line 158). The
non-finite value arises *inside* ALDERAAN: at lines 856–857 `durs` is **recomputed from the MAP
transit-shape fit** (`durs = shape_map["dur"]`), and the second PyMC3 (TTV) model at line 1189/1197
consumes that refit value. For a degenerate target (no real transit signal, or a quarter with a
single/near-constant point → variance 0), the MAP collapses `dur→0` or `vbq→0`, and `np.log`
poisons the next model. This is ALDERAAN correctly failing on data it cannot fit  -  legitimate
attrition. Not fixable from our side without modifying vendored ALDERAAN (which we must not do).

### B.2 `TypeError: bad operand type for unary ~: 'NoneType'` (2 occurrences)

**Mechanism.** `refit_mask_sc` and `refit_mask_lc` are set to `None` when `sc`/`lc` (the short- and
long-cadence light-curve objects) are `None` (lines 1109, 1122). They then flow into mask arithmetic
(lines 1133, 1146) and downstream boolean operations. When a target has effectively **no usable
long-cadence data surviving the quality cuts**, `lc` is `None`, `refit_mask_lc` stays `None`, and a
`~mask` / arithmetic-on-None raises this TypeError.

**Verdict: INPUT-DATA / pathological target, NOT our bug  -  with one caveat.**
- Reasoning it is not our bug: SC is never downloaded in our pipeline (findings F1), so `sc` and
  therefore `refit_mask_sc = None` for **every** target. If an *unguarded* `~refit_mask_sc` were the
  cause, ALL ~592 targets would fail, not 2. They don't. So these 2 are targets where even **long
  cadence** produced no usable mask (`lc is None`)  -  a data-quality pathology specific to those KOIs
  (e.g. "Over 50% of transits flagged low quality"-type targets). Legitimate attrition.
- Caveat / honest limit: I could not read the actual failing targets' KIDs or their downloaded
  FITS (no run-log access), so I cannot 100% exclude a per-target download corruption (e.g. a
  wget/redirect edge case leaving a truncated FITS that parses to an empty `lc`). Given the live
  bundle uses the wget fix, corruption is unlikely, but **not fully verifiable from here.** If a
  retry is run, capture the 3 failing KIDs and confirm each has a non-empty `*_lc.fits` with real
  flux; if the FITS is fine, it is confirmed data-quality attrition.

**Bottom line for Task B:** all 3 errors are best explained as legitimate ALDERAAN attrition on
pathological targets, not a catalog/invocation bug of the bug-9 type. No fix warranted; recommend
only that the retry pass log the failing KIDs so B.2's caveat can be closed.

---

## Task C  -  20-year-veteran methodology review (5 items chosen)

Chose the 5 highest-leverage checks that (a) are primary axes of the paper's analysis and (b) are
easy to get subtly wrong in re-implementation. Each was run against
`outputs/canonical_sample_old_astropy_rawcc.csv` (2474 planets, 1857 hosts) and the current KOI
table `cumulative_2026.02.11_22.33.58.csv`.

### C1  -  False-positive exclusion uses the CURRENT, authoritative disposition (checked; **correct**)
- Check: mapped every sampled `kepoi_name` to `koi_disposition` in the current cumulative table.
- Found: **1703 CONFIRMED + 771 CANDIDATE, 0 FALSE POSITIVE.** Filtering is on `koi_disposition`
  (the archive's current disposition), not the stale `koi_pdisposition` (diagnose_sample.py line 57).
- Subtlety confirmed handled: **3 sampled KOIs have `koi_pdisposition == FALSE POSITIVE` but
  `koi_disposition == CANDIDATE`** (later re-promoted by the archive). Using `koi_disposition` keeps
  them correctly. Had the pipeline used pdisposition, it would have wrongly dropped 3 real planets.
  **Non-issue  -  the right column was chosen.** KOI table date (2026.02.11) is recent.

### C2  -  Duplicate / renamed KOI handling (checked; **clean**)
- Check: `kepid`+`period`(3 dp) combos with >1 distinct `kepoi_name`; and duplicate `kepoi_name`.
- Found: **0 collisions, 0 duplicate KOIs; 2474 rows all unique `kepoi_name`.** No evidence of a
  reused/split/merged KOI leaking a duplicate planet into the sample. **Non-issue.**

### C3  -  Whole-system completeness / single-vs-multi labeling (checked; **fragile code, but non-issue in practice**)
- The label is assigned at **step 10, AFTER all sample cuts** (`add_target_and_system`,
  common.py:116, counts `kepoi_name` per `kepid` on the *already-cut* df; audit stage literally
  named `10_assign_single_multi_after_cuts`). In principle a system that entered with ≥2
  CONFIRMED/CANDIDATE KOIs but lost a sibling to a cut (contamination, RUWE, Teff, period 1–100 d,
  Bin=0) would be **mislabeled `single`**  -  a direct corruption of the paper's primary single/multi
  axis (eccentricity differs by category).
- Measured: of the **1381 single hosts** in the canonical sample, how many currently have >1
  CONFIRMED/CANDIDATE KOI in the cumulative table? **0.** So in the canonical sample the "singles"
  are genuine singles; the post-cut labeling did **not** manufacture false singles here.
- **Verdict: non-issue for the current canonical sample, but the code pattern is fragile.** Proposed
  hardening (only if a fidelity rerun is commissioned): compute `system` from the RAW
  CONFIRMED+CANDIDATE KOI roster per `kepid` (before cuts) and carry it through, so single/multi is
  invariant to which siblings survive photometric/kinematic cuts. Sagear's classification is
  system-level, so the raw roster is the more defensible definition. **Not** a silent mid-run change.

### C4  -  Grazing / high-impact transits (checked; **fidelity flag, not a code bug**)
- Photoeccentric duration→eccentricity inversion degrades badly for grazing geometries (b→1, the
  duration/impact/eccentricity degeneracy blows up).
- Found in sample: median b = 0.387; **11.2% have b > 0.9; 1.5% (≈37 planets) have b > 1.0**
  (formally grazing); **27 planets have NaN `koi_impact`** → defaulted to b = 0.5 seed in the catalog
  (alderaan_batch.py:196). The seed only initializes ALDERAAN's b prior (b is refit), so the default
  is defensible, but for the ~37 grazing planets the per-planet eccentricity posterior will be
  weakly constrained/biased.
- **Verdict: not a code bug** (nothing in our prose says to exclude grazers, and the paper does not
  clearly state a b cut either), but a **fidelity item**: confirm whether the paper's §3 excludes or
  down-weights high-b systems; if so, add a `koi_impact < ~0.9` (or fitted-b) filter or flag column
  in the population fit. Recommend flagging b>1 planets in the eccentricity table for a sensitivity
  test rather than silently dropping. **Worth a human decision.**

### C5  -  Stellar density prior construction (the photoeccentric core) (checked; **correct**)
- The eccentricity measurement hinges on the stellar-density prior. `extract_eccentricity_posteriors.py`
  draws `rho` from Berger+2020 log-density with asymmetric errors (`draw_rho_log`, lines 311–333).
- Ambiguity worth checking: `absolute_density_error` (line 336) returns `10**value` for the stored
  `rho_log_upper/lower`. That is only correct if Berger+2020 stores the error as *log10 of the
  absolute error in solar-density units* (not as the ± on log-rho).
- Verified against the data (`alderaan_needed_catalog_rows_best.csv`): `rho_log` median −0.27 (rho ≈
  0.54 solar); `rho_log_upper/lower` are **always negative** (median −1.18), so `10**rho_log_upper`
  ≈ 0.066 ≈ **12% of rho**  -  physically sensible, and consistent with the `0.13*rho` fallback used
  when the error is missing (line 317). The stored-value interpretation is therefore **correct**;
  `10**value` recovers a sane absolute density error. **Non-issue  -  good defensive design.**

### Task C summary
| Item | Result |
|---|---|
| C1 FP exclusion / current disposition | Correct (right column; 3 pdisp≠disp cases handled correctly) |
| C2 Duplicate/renamed KOI | Clean (0 duplicates) |
| C3 Single/multi completeness | Fragile code, 0 actual mislabels in canonical sample; harden only on rerun |
| C4 Grazing/high-b | ~37 grazing (b>1) + 27 NaN-b planets; fidelity flag, needs human decision on b cut |
| C5 Stellar-density prior | Correct (Berger+2020 log-error interpretation verified against data) |

**Not checked** (out of scope this pass): exact transit-probability geometric correction weighting,
TTV-system special handling in the population fit, and the half-Gaussian/monotonic-Beta
hyperprior forms.

---

## Consolidated recommended fixes (for human review  -  none applied)
1. **A2–A5:** add `--data_dir`/`--catalog` to the `analyze_autocorrelated_noise.py` line in all four
   generators (fix alongside the curl→wget F2 fix  -  same "regeneration reintroduces a live bug" class).
2. **A6:** add a limbdark consistency check to the bundle validator (defense in depth for bug 9).
3. **B:** no fix; log the 3 failing KIDs on retry to confirm B.2's data-quality verdict.
4. **C3:** compute single/multi from the pre-cut raw KOI roster (only in a labeled fidelity rerun).
5. **C4:** decide on a grazing-transit (b) cut/flag; add a sensitivity test rather than silent drop.
