# ALDERAAN Factorial Validation Results

This directory is the curated release of the completed 82-fit ALDERAAN factorial
validation matrix. The binary result FITS are stored with Git LFS. They are exact copies
of the `*-results.fits` files selected from the hash-verified archive
`alderaan_factorial_FULL_82_20260715T073853Z.tar.gz`.

Archive SHA-256:

```text
dcb1a64b75bb1fd67f7afc65f63a008f2c2b153277343ae23e342c22ef9b0859
```

## Contents

| arm | FITS files | purpose |
|---|---:|---|
| `original_lc` | 24 | Baseline long-cadence run with the original catalog limb-darkening inputs. |
| `reference_lc` | 24 | Long-cadence run with the Sagear-reference limb-darkening inputs. |
| `original_lcsc` | 9 | Original limb darkening with long and available short cadence. |
| `reference_lcsc` | 9 | Reference limb darkening with long and available short cadence. |
| `original_lc_repeat` | 8 | Repeated baseline fits for nested-sampling variability. |
| `paper_priors_original_lc` | 8 | Baseline long-cadence fits using the printed paper-prior sensitivity. |

Each FITS file corresponds to one KOI target system, not necessarily one planet. Multi-
planet systems were fitted simultaneously by ALDERAAN.

## Layout

- `results/<arm>/`: the 82 binary ALDERAAN result FITS files.
- `provenance/input_catalogs/`: the catalog inputs supplied to the original and reference
  limb-darkening arms, plus the full system inventory.
- `provenance/target_sets/`: target selections for the limb-darkening, cadence, and repeat
  validation tests.
- `provenance/status_manifests/`: completion records for all six arms.
- `provenance/run_spec/`: the matrix runner used to schedule the six arms.
- `provenance/selection/`: the target-selection ledger. Its `posterior_file` values are
  posterior basenames only; an irrelevant local machine path prefix was removed before
  publication. No numerical, target, or quality-control field was changed.
- `provenance/SHA256SUMS.txt`: checksums for every file in this release other than itself.

## Verification

After obtaining LFS objects, verify the data by computing SHA-256 values against
`provenance/SHA256SUMS.txt`. On Linux or macOS:

```bash
cd data/alderaan_factorial_validation_20260715
sha256sum -c provenance/SHA256SUMS.txt
```

The release intentionally excludes Kepler light curves, detrended and filtered light
curves, dynesty checkpoints, Gaussian-process intermediates, and VM logs. Those files are
regenerable or duplicate the information in the result FITS and provenance manifests. The
source archive remains locally preserved and hash-verified.

## Scientific Use

These files establish the controlled ALDERAAN comparison. They do not themselves resolve
the Sagear replication discrepancy. Analysis should compare matched arms while keeping
target-system correlations, quality-control flags, and the distinction between literal
published-method replication and a generative population model explicit.
