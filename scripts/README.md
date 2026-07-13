# Sagear Replication Pipeline

Canonical pipeline for diagnosing and reproducing the Sagear Kepler thin/thick disk eccentricity result.

This folder is deliberately separate from the exploratory notebooks and older scripts. The first goal is not to force a match; it is to make every mismatch visible.

Run commands from the repository root in an environment created from
`requirements.txt`. The publication-only audit is self-contained because the official
machine-readable tables are bundled under `reference/data/`:

```powershell
python scripts/published_sagear_audit.py
python -m pytest -q scripts/test_published_sagear_audit.py
```

Other analyses still require the large catalogs and ALDERAAN products excluded from
Git. Override `paths.research_root` in a local config instead of editing committed
defaults. The final published disk labels are canonical; reconstructed GMM labels are
diagnostic only.

## Files

- `config.json`: paths, thresholds, Sagear target counts, and missing optional inputs.
- `common.py`: shared catalog parsing, fixed-width Berger reader, disk helpers.
- `diagnose_sample.py`: attrition audit and disk-count comparison.
- `prepare_external_inputs.py`: downloads Furlan and prepares the APOGEE DR17 Kepler crossmatch.
- `alderaan_batch.py`: ALDERAAN project setup, validation target selection, and command generation.
- `extract_eccentricity_posteriors.py`: converts ALDERAAN transit-shape results into `e, omega` posterior grids.
- `hierarchical_rayleigh.py`: first-pass Rayleigh population fit from posterior grids.
- `toomre_diagnostics.py`: compares Toomre coordinate/sign conventions and disk-classifier variants.
- `classifier_threshold_diagnostics.py`: sweeps `P_thick` thresholds for selected classifier probability fields.
- `target_consistency_diagnostics.py`: parses Sagear manuscript macros and reports count inconsistencies.
- `classifier_disagreement_diagnostics.py`: lists systems whose labels change across near-match classifier choices.
- `eccentricity_diagnostics.py`: checks eccentricity distributions, high-e outliers, and quality-cut sensitivity.
- `catalog_diagnostics.py`: compares the current stellar catalog joins against Berger+2018.
- `formula_sanity_checks.py`: deterministic checks for photoeccentric/posterior and Rayleigh-selection formulas.

## First Run

For a full local-data audit, create a local config pointing at the research data root:

```powershell
python scripts/diagnose_sample.py --config path\to\local_config.json
```

That strict command is the Sagear-equivalence gate. If Furlan, APOGEE, or ALDERAAN convergence inputs are missing, it marks the sample as unclassified/missing instead of inventing labels.

Strict outputs:

- `sagear_reproduction/outputs/sample_audit_strict.csv`
- `sagear_reproduction/outputs/disk_counts_strict.csv`
- `sagear_reproduction/outputs/canonical_sample_strict.csv`
- `sagear_reproduction/outputs/pipeline_status_strict.json`

## External Inputs

Furlan+2017 table 9 is now configured at `data/furlan2017_table9.dat`. If it ever needs to be refreshed:

```powershell
& "C:\Users\shres\anaconda3\python.exe" sagear_reproduction\prepare_external_inputs.py --furlan
```

The APOGEE DR17 allStarLite file is about 1.7 GB. To download it and build the Kepler/APOGEE chemical crossmatch used by the strict GMM classifier:

```powershell
& "C:\Users\shres\anaconda3\python.exe" sagear_reproduction\prepare_external_inputs.py --apogee
```

That writes `data/apogee_dr17_kepler_crossmatch.csv`, which `diagnose_sample.py` already expects.

Berger+2018 `J/ApJ/866/99/table1` is now used for an explicit `Bin=0` resolved-companion/binary cut. The local file is `data/berger2018_table1_min.tsv`.

For ALDERAAN validation target selection with the old all-Angus diagnostic classifier, run the explicit fallback:

```powershell
& "C:\Users\shres\anaconda3\python.exe" sagear_reproduction\diagnose_sample.py --force-fallback-gmm
```

Diagnostic outputs:

- `sagear_reproduction/outputs/sample_audit_diagnostic.csv`
- `sagear_reproduction/outputs/disk_counts_diagnostic.csv`
- `sagear_reproduction/outputs/canonical_sample_diagnostic.csv`
- `sagear_reproduction/outputs/pipeline_status_diagnostic.json`

## ALDERAAN Stage

ALDERAAN requires its own conda environment. The current base environment is missing the core ALDERAAN stack (`pymc3`, `exoplanet`, `dynesty`, `batman`, `celerite2`, `ldtk`).

```powershell
git clone https://github.com/gjgilbert/alderaan external\alderaan
conda env create -n alderaan -f external\alderaan\environment.yml
```

In this workspace, a direct `conda env create` attempt timed out after 10 minutes and the partial environment was removed. Use `sagear_reproduction/alderaan_project/Scripts/setup_alderaan_env.ps1` from a normal terminal so the solve/install can run to completion.

Then build a validation-batch project. By default this now includes representatives from each disk/multiplicity bin
plus high-e thin-single stress-test systems from `eccentricity_diagnostics.py`:

```powershell
& "C:\Users\shres\anaconda3\python.exe" sagear_reproduction\alderaan_batch.py prepare --n-per-bin 3 --n-high-e 8
```

This writes catalogs and command scripts under `sagear_reproduction/alderaan_project/`.
By default this uses `canonical_sample_diagnostic.csv` so validation batches can still be chosen while APOGEE/Furlan are missing. Pass `--sample` to use a fully strict sample once those inputs exist.

## GCP Batch Stage

To prepare a Linux/GCP bundle for around 300 planet rows:

```powershell
& "C:\Users\shres\anaconda3\python.exe" sagear_reproduction\cloud_prepare.py --max-planets 300 --jobs 30
```

This writes `sagear_reproduction/cloud_batch/`, including target CSVs, an ALDERAAN catalog, a VM setup script, and a parallel batch runner.
By default this also uses `canonical_sample_diagnostic.csv`; pass `--sample` when the strict Sagear sample is available.

## Diagnostics

To diagnose the Toomre/classifier mismatch:

```powershell
& "C:\Users\shres\anaconda3\python.exe" sagear_reproduction\toomre_diagnostics.py
```

Key outputs:

- `sagear_reproduction/outputs/toomre_diagnostics_plot.png`
- `sagear_reproduction/outputs/toomre_diagnostics_classifier_variants.csv`

To test whether the remaining disk-count mismatch can be explained by a classifier threshold choice:

```powershell
& "C:\Users\shres\anaconda3\python.exe" sagear_reproduction\classifier_threshold_diagnostics.py
```

Key outputs:

- `sagear_reproduction/outputs/classifier_threshold_diagnostics_best.csv`
- `sagear_reproduction/outputs/classifier_threshold_diagnostics_grid.csv`

To audit manuscript target-count consistency:

```powershell
& "C:\Users\shres\anaconda3\python.exe" sagear_reproduction\target_consistency_diagnostics.py
```

To identify the systems/planets driving classifier disagreements:

```powershell
& "C:\Users\shres\anaconda3\python.exe" sagear_reproduction\classifier_disagreement_diagnostics.py
```

Key outputs:

- `sagear_reproduction/outputs/target_consistency_diagnostics_checks.csv`
- `sagear_reproduction/outputs/classifier_disagreement_diagnostics_pairwise_summary.csv`
- `sagear_reproduction/outputs/classifier_disagreement_diagnostics_changed_planets.csv`

To diagnose the thin-single eccentricity distribution and outliers:

```powershell
& "C:\Users\shres\anaconda3\python.exe" sagear_reproduction\eccentricity_diagnostics.py
```

If ALDERAAN posterior summaries are absent, this uses the old `e_photo` point estimates and marks the results as triage-only.

To run formula-level sanity checks for the photoeccentric posterior and Rayleigh hierarchical selection correction:

```powershell
& "C:\Users\shres\anaconda3\python.exe" sagear_reproduction\formula_sanity_checks.py
```

Key output:

- `sagear_reproduction/outputs/formula_sanity_checks/formula_sanity_checks.json`

To compare Berger+2018 catalog availability/radii against the current joined stellar catalog:

```powershell
& "C:\Users\shres\anaconda3\python.exe" sagear_reproduction\catalog_diagnostics.py
```

Key outputs:

- `sagear_reproduction/outputs/catalog_diagnostics_summary.csv`
- `sagear_reproduction/outputs/catalog_diagnostics_missing_by_population.csv`
- `sagear_reproduction/outputs/catalog_diagnostics_radius_comparison.png`
- `sagear_reproduction/outputs/sagear_diagnosis_report.md`

## Important Current Gaps

The local folder now has Furlan+2017 and an APOGEE DR17 Kepler crossmatch. The current strict audit still does not reproduce Sagear's disk counts:

- strict APOGEE-calibrated pooled GMM currently gives `94` thick singles and `61` thick multi planets, far below Sagear's `275` and `207`.
- forced diagnostic all-Angus GMM gives `1121` thin singles, `305` thick singles, `822` thin multi planets, and `226` thick multi planets. This is much closer, but still not Sagear-equivalent.

The sample audit now includes two important cuts that were missing from the first pass:

- `berger_teff < 6500`, matching the paper's FGKM sample description.
- Berger+2018 `Bin=0`, a likely resolved-companion/binary exclusion inherited from the stellar catalog.

Together these reduce the pre-ALDERAAN sample to `2474` planets, close to Sagear's `2465`.

The remaining blockers for a faithful Sagear sample/classification match are:

- the exact Sagear GMM implementation/convention used to turn APOGEE high/low-alpha calibration into `P_thick`;
- ALDERAAN convergence/results outputs for the final posterior-quality cut.

Additional diagnostics now indicate:

- Berger+2018 availability is not the remaining sample-count culprit: after the new cuts, all `2474` diagnostic planet rows match the downloaded Berger+2018 table.
- Current-vs-Berger+2018 stellar radii are typically within a few percent, so catalog radius differences are unlikely to explain the large thin-single eccentricity discrepancy by themselves.
- The old `e_photo` thin-single median is close to Sagear, but its mean is elevated by a short-period/high-impact tail; this is exactly what the ALDERAAN validation batch should test.
- `extract_eccentricity_posteriors.py` now leaves individual `e,omega` posteriors without a geometric transit prior by default, and `hierarchical_rayleigh.py` applies the reciprocal transit-probability correction over the joint posterior grid.
- Sagear's manuscript macros are internally inconsistent: disk-total planet macros sum to `2474`, while subgroup/all-planets macros sum to `2465`. The current corrected sample lands at `2474`, so this 9-planet discrepancy should not be overfit without external clarification.

Until the classifier convention and ALDERAAN outputs are resolved, the audit remains non-Sagear-equivalent.
