# Reproducibility map

This map links the repository's main scientific claims to their code, tests,
and compact evidence. It is a review aid, not a substitute for the methods.

| Claim or contract | Implementation | Verification or evidence |
| --- | --- | --- |
| Published host labels are canonical | `scripts/published_sagear_audit.py` | `reference/data/`; `metadata/public_reconstruction_20260727/` |
| Multiplicity precedes planet cuts | `scripts/diagnose_sample.py` | `scripts/test_recovery_preflight.py`; `metadata/recovery_preflight_20260724/` |
| Transit samples stay paired | `scripts/extract_eccentricity_posteriors_direct.py` | `scripts/test_direct_importance_extractor.py` |
| Dynesty points use `LN_WT` | `scripts/common.py`; direct extractor | `scripts/test_nested_weighting_forensic_audit.py`; `metadata/final_forensic_20260808/` |
| Selection modes remain explicit | `scripts/hierarchical_rayleigh.py` | hierarchy and outlier-floor tests |
| Factorial arms compare within planet | `scripts/compare_factorial_validation.py` | `scripts/test_compare_factorial_validation.py`; `metadata/factorial_validation_20260715/` |
| Gilbert provides a real-data control | `scripts/gilbert_real_alderaan_control.py` | `metadata/final_forensic_20260808/` |
| Density inputs drive the mismatch | density and residual audits | `metadata/uncertainty_calibration_20260810/`; `docs/uncertainty_and_density_findings_20260810.md` |
| Public inputs do not fix the Table 3 path | final boundary audits | `metadata/final_public_boundary_20260813/`; `docs/final_public_data_boundary.md` |

## Review order

1. Read `README.md` and `docs/replication_status.md`.
2. Check `docs/final_public_data_boundary.md` for closed and open leads.
3. Run the test suite and `scripts/check_professor_release.py`.
4. Inspect the compact metadata directories before requesting excluded data.
