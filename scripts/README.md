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
python scripts/prepare_external_inputs.py --furlan
```

The APOGEE DR17 allStarLite file is about 1.7 GB. To download it and build the Kepler/APOGEE chemical crossmatch used by the strict GMM classifier:

```powershell
python scripts/prepare_external_inputs.py --apogee
```

That writes `data/apogee_dr17_kepler_crossmatch.csv`, which `diagnose_sample.py` already expects.

Berger+2018 `J/ApJ/866/99/table1` is now used for an explicit `Bin=0` resolved-companion/binary cut. The local file is `data/berger2018_table1_min.tsv`.

For ALDERAAN validation target selection with the old all-Angus diagnostic classifier, run the explicit fallback:

```powershell
python scripts/diagnose_sample.py --force-fallback-gmm
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
python scripts/alderaan_batch.py prepare --n-per-bin 3 --n-high-e 8
```

This writes catalogs and command scripts under `sagear_reproduction/alderaan_project/`.
By default this uses `canonical_sample_diagnostic.csv` so validation batches can still be chosen while APOGEE/Furlan are missing. Pass `--sample` to use a fully strict sample once those inputs exist.

## GCP Batch Stage

To prepare a Linux/GCP bundle for around 300 planet rows:

```powershell
python scripts/cloud_prepare.py --max-planets 300 --jobs 30
```

This writes `sagear_reproduction/cloud_batch/`, including target CSVs, an ALDERAAN catalog, a VM setup script, and a parallel batch runner.
By default this also uses `canonical_sample_diagnostic.csv`; pass `--sample` when the strict Sagear sample is available.

## Diagnostics

To diagnose the Toomre/classifier mismatch:

```powershell
python scripts/toomre_diagnostics.py
```

Key outputs:

- `sagear_reproduction/outputs/toomre_diagnostics_plot.png`
- `sagear_reproduction/outputs/toomre_diagnostics_classifier_variants.csv`

To test whether the remaining disk-count mismatch can be explained by a classifier threshold choice:

```powershell
python scripts/classifier_threshold_diagnostics.py
```

Key outputs:

- `sagear_reproduction/outputs/classifier_threshold_diagnostics_best.csv`
- `sagear_reproduction/outputs/classifier_threshold_diagnostics_grid.csv`

To audit manuscript target-count consistency:

```powershell
python scripts/target_consistency_diagnostics.py
```

To identify the systems/planets driving classifier disagreements:

```powershell
python scripts/classifier_disagreement_diagnostics.py
```

Key outputs:

- `sagear_reproduction/outputs/target_consistency_diagnostics_checks.csv`
- `sagear_reproduction/outputs/classifier_disagreement_diagnostics_pairwise_summary.csv`
- `sagear_reproduction/outputs/classifier_disagreement_diagnostics_changed_planets.csv`

To diagnose the thin-single eccentricity distribution and outliers:

```powershell
python scripts/eccentricity_diagnostics.py
```

If ALDERAAN posterior summaries are absent, this uses the old `e_photo` point estimates and marks the results as triage-only.

To run formula-level sanity checks for the photoeccentric posterior and Rayleigh hierarchical selection correction:

```powershell
python scripts/formula_sanity_checks.py
```

Key output:

- `sagear_reproduction/outputs/formula_sanity_checks/formula_sanity_checks.json`

To compare Berger+2018 catalog availability/radii against the current joined stellar catalog:

```powershell
python scripts/catalog_diagnostics.py
```

Key outputs:

- `sagear_reproduction/outputs/catalog_diagnostics_summary.csv`
- `sagear_reproduction/outputs/catalog_diagnostics_missing_by_population.csv`
- `sagear_reproduction/outputs/catalog_diagnostics_radius_comparison.png`
- `docs/replication_status.md`

## Current Status

The count and classifier discussion that originally followed this heading
described the superseded 2,474-planet reconstruction. It is retained in the
dated July audit documents, not as current guidance.

The public reconstruction now contains exactly 2,465 planets and uses the
published host labels directly:

- 1,109 thin singles;
- 269 thick singles;
- 878 thin multis;
- 209 thick multis.

Multiplicity is frozen from the full eligible system before planet-level fit
cuts. The direct extractor preserves paired ALDERAAN transit-shape and period
samples and applies dynesty weights. The canonical hierarchy uses the
forward transit probability with population normalization; manuscript
reciprocal modes remain explicit sensitivities rather than defaults.

The current blocker is not the old 2,474-versus-2,465 count ambiguity. The
full 2,465-planet circular-density audit localizes the discrepancy upstream of
the population fit, while the remaining missing provenance is the exact
stellar-density input, final visual-QC exclusions, planet-level posterior
export, and Table 3 population implementation. See
`docs/replication_status.md` and
`metadata/public_reconstruction_20260727/README.md` for current numbers.
