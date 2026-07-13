# Repo Contents

This layout was independently reviewed and verified against the real local data
(2026-07-02) before publishing; see "Independent Verification" in `README.md`.

## Include

- `scripts/`: reproducible pipeline scripts and diagnostics.
- `cloud/`: ALDERAAN/GCP helper scripts for the missing-posterior run.
- `docs/`: concise reports and operational checklists.
- `metadata/`: small CSV summaries needed to understand current state.
- `reference/`: submitted paper/figures, original photoeccentric script, and the final
  AJ machine-readable Tables 1-3 under `reference/data/`.
- `legacy/`: superseded earlier-attempt results, kept for provenance.

## Exclude

- Raw Kepler light curves.
- ALDERAAN posterior FITS files.
- Downloaded `.tar.gz` source/data archives.
- Full generated posterior directories.
- Cloud result tarballs.
- Large raw catalog products unless explicitly needed and small enough for GitHub.

## Review Status (completed 2026-07-02)

1. Scripts compile from this layout (`py_compile` / `bash -n` pass). ✔
2. `.gitignore` excludes heavy/generated artifacts (spot-checked with `git check-ignore`). ✔
3. README states current status with verification results and known caveats. ✔
4. Cloud scripts reviewed: SPOT + `--max-run-duration` cap, no hidden billable operations;
   `pack_results.sh` hardened to tar only existing paths. ✔
5. Regeneration documented: external inputs via `scripts/prepare_external_inputs.py`,
   bundle via `scripts/prepare_gcp_missing_alderaan_bundle.py`, validation via
   `cloud/validate_bundle.py`. ✔
