# Sagear replication walkthrough

Updated: 2026-09-20

This supporting guide follows Sagear et al.'s analysis in order and points to
the corresponding code and evidence in this repository. It distinguishes a
paper-stated method, a public-data reconstruction, and an input that the
public record does not expose.

## 1. Build the Kepler host and planet sample

The paper begins with Kepler planet hosts and applies its
stated false-positive, period, contamination, and Gaia-quality rules.

The reconstruction uses [`scripts/diagnose_sample.py`](../scripts/diagnose_sample.py).
It
reads the required local catalogs and writes attrition after every join and
cut. [`scripts/published_sagear_audit.py`](../scripts/published_sagear_audit.py)
independently audits the machine-readable article tables. The current compact
inventory in [`metadata/public_reconstruction_20260727/`](../metadata/public_reconstruction_20260727/)
has 2,465 planets, matching the paper's headline total.

A matching total does not prove matching membership. The paper
does not provide an identifier-level ledger of final accepted and rejected
planets.

## 2. Assign Galactic disk membership

Sagear et al. use Galactic velocities and an APOGEE chemical
calibration to assign a thick-disk probability. A host is thick when that
probability exceeds 0.5. The Toomre diagram is a velocity-classification check,
not an eccentricity measurement.

The public-catalog reconstruction is in
[`scripts/diagnose_sample.py`](../scripts/diagnose_sample.py). Its velocity and
classifier checks are in [`scripts/toomre_diagnostics.py`](../scripts/toomre_diagnostics.py),
[`scripts/toomre_sagear_style.py`](../scripts/toomre_sagear_style.py), and
[`scripts/classifier_reconstruction_audit.py`](../scripts/classifier_reconstruction_audit.py).

The catalog reconstruction did not reliably reproduce the published labels.
The primary population comparison therefore uses the machine-readable host
table directly through [`scripts/published_sagear_audit.py`](../scripts/published_sagear_audit.py).
It contains 1,888 hosts: 1,515 thin and 373 thick. The manuscript text reports
378 thick hosts, which conflicts with the table total and subgroup arithmetic.

Published labels are matched in the primary comparison. The
independent classifier remains useful as a diagnostic of the Toomre mismatch.

## 3. Freeze single and multi architecture before fit cuts

A single or multi label describes the known system
architecture. A multi should remain a multi if one companion later fails a
quality rule.

[`scripts/diagnose_sample.py`](../scripts/diagnose_sample.py)
uses the `raw_koi` multiplicity basis and freezes the system count before
individual planet cuts. [`scripts/published_sagear_audit.py`](../scripts/published_sagear_audit.py)
preserves that label while applying published disk labels. The regression
evidence is in [`scripts/test_multiplicity_contract.py`](../scripts/test_multiplicity_contract.py)
and [`metadata/recovery_preflight_20260724/`](../metadata/recovery_preflight_20260724/).

This corrected an exploratory bookkeeping error that could turn surviving
members of a multi-planet system into apparent singles.

## 4. Fit transit shapes with ALDERAAN

ALDERAAN fits Kepler light curves at the system level. Its
samples of duration, radius ratio, impact parameter, and period are correlated.
The later eccentricity calculation needs each sampled tuple, not independent
median values.

The runner and recovery helpers are under
[`cloud/`](../cloud/). The release includes 82 returned validation FITS files
and their provenance under
[`data/alderaan_factorial_validation_20260715/`](../data/alderaan_factorial_validation_20260715/).
Large light curves and full posterior archives are deliberately excluded; see
[`data_availability.md`](data_availability.md).

The 82-fit matrix gives valid tests of selected limb-darkening, radius-ratio
prior, and sampler-seed choices. Its intended short-cadence comparison is not
valid because the archived runners did not pass `--use_sc True` and their logs
show no processed short-cadence data. The limitation is documented in
[`corrected_cadence_validation_plan.md`](corrected_cadence_validation_plan.md).

## 5. Derive planet-level eccentricity posteriors

Photoeccentric inference compares the density implied by a
transit duration with a stellar-density prior, while marginalizing over
eccentricity and orbital orientation.

[`scripts/extract_eccentricity_posteriors_direct.py`](../scripts/extract_eccentricity_posteriors_direct.py)
implements the finite-duration importance calculation. It keeps duration,
radius ratio, impact, and period paired row by row and uses Dynesty `LN_WT`
weights. [`scripts/test_direct_importance_extractor.py`](../scripts/test_direct_importance_extractor.py)
checks those contracts. Earlier point estimates and geometric impact draws are
diagnostic-only.

The exact stellar-density ledger is not publicly reproducible.
The paper cites Berger et al. (2018) in its Methods and Berger et al. (2020)
elsewhere. The public 2018 table lacks the density, mass, and covariance
information needed to reconstruct the original per-planet prior. See
[`uncertainty_and_density_findings_20260810.md`](uncertainty_and_density_findings_20260810.md).

## 6. Combine planets into four populations

Sagear et al. combine eccentricity information for thin
singles, thick singles, thin multis, and thick multis with a Rayleigh model and
a transit-selection correction.

[`scripts/hierarchical_rayleigh.py`](../scripts/hierarchical_rayleigh.py)
implements a Rayleigh grid model with explicit selection modes. It rejects
posterior files that would double count a transit prior. Its core contracts are
tested in [`scripts/test_hierarchical_contract.py`](../scripts/test_hierarchical_contract.py)
and [`scripts/test_hierarchical_outlier_floor.py`](../scripts/test_hierarchical_outlier_floor.py).

The paper describes a NumPyro implementation, but its exact input order and
configuration are not public. This is therefore a tested reconstruction of the
stated model, not a claim of identical software.

## 7. Current comparison with Table 3

The current public QC values are versioned in
[`metadata/public_reconstruction_20260727/population_comparison.csv`](../metadata/public_reconstruction_20260727/population_comparison.csv).

| Population | Current fit N | Current mean | Paper N | Paper mean |
| --- | ---: | ---: | ---: | ---: |
| Thin singles | 1,104 | 0.232 (0.225-0.239) | 1,121 | 0.022 (0.017-0.029) |
| Thick singles | 268 | 0.208 (0.194-0.223) | 275 | 0.066 (0.045-0.096) |
| Thin multis | 873 | 0.094 (0.088-0.103) | 862 | 0.030 (0.023-0.031) |
| Thick multis | 209 | 0.133 (0.121-0.147) | 207 | 0.033 (0.015-0.065) |

The mismatch is present in every group. The small remaining count differences
cannot plausibly account for it alone. The current audit places the difference
before the population fit: available transit-density posteriors and adopted
stellar-density inputs do not produce the narrow low-eccentricity population
signal reported in the paper.

## 8. Additional diagnostic work

| Check | Why it matters | Code or evidence |
| --- | --- | --- |
| Published-label audit | Separates disk classification from eccentricity inference. | [`published_sagear_audit.py`](../scripts/published_sagear_audit.py) |
| Toomre checks | Tests velocity and chemistry-calibration assumptions. | [`toomre_diagnostics.py`](../scripts/toomre_diagnostics.py) |
| Dynesty-weight audit | Prevents raw nested rows from being treated as equal posterior draws. | [`nested_weighting_forensic_audit.py`](../scripts/nested_weighting_forensic_audit.py) |
| Gilbert control | Tests the extraction and hierarchy against an independent low-eccentricity result. | [`gilbert_real_alderaan_control.py`](../scripts/gilbert_real_alderaan_control.py) |
| Factorial comparisons | Measures selected prior, limb-darkening, and seed sensitivity within systems. | [`compare_factorial_validation.py`](../scripts/compare_factorial_validation.py) |
| Robustness checks | Measures leverage and leave-out behavior. | [`population_system_robustness.py`](../scripts/population_system_robustness.py) |
| Injection preflight | Freezes a small auditable photometric recovery test. | [`pilot_injection_preflight.md`](pilot_injection_preflight.md) |

## 9. What would close the remaining gap

The repository has established the primary implementation contracts: published
host labels, pre-cut multiplicity, paired transit samples, Dynesty weights, and
explicit selection handling. It has also excluded invalid ways of making the
numbers appear closer, including equal-row nested weighting and the mislabeled
cadence arm.

The smallest high-value author products are:

1. the final accepted planet ledger with KOI IDs and rejection reasons;
2. the stellar-density rows or posterior samples supplied to the
   photoeccentric calculation.

These products would identify whether the remaining difference enters through
membership, stellar priors, or downstream eccentricity inference. For the
current boundary and all supporting evidence, read
[`replication_match_ledger.md`](replication_match_ledger.md),
[`current_inference_audit.md`](current_inference_audit.md), and
[`final_public_data_boundary.md`](final_public_data_boundary.md).
