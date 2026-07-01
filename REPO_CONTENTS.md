# New Repo Candidate Contents

This folder is the clean candidate layout Claude should review and then use for a new GitHub repo.

## Include

- `scripts/`: reproducible pipeline scripts and diagnostics.
- `cloud/`: ALDERAAN/GCP helper scripts for the missing-posterior run.
- `docs/`: concise reports and operational checklists.
- `metadata/`: small CSV summaries needed to understand current state.

## Exclude

- Raw Kepler light curves.
- ALDERAAN posterior FITS files.
- Downloaded `.tar.gz` source/data archives.
- Full generated posterior directories.
- Cloud result tarballs.
- Large raw catalog products unless explicitly needed and small enough for GitHub.

## Review Before Publishing

Claude should verify:

1. The scripts still run from this new layout or update path handling.
2. The `.gitignore` excludes heavy/generated artifacts.
3. The README states current status without overclaiming.
4. The cloud scripts are safe enough and do not hide billable operations.
5. The repository can be recreated from documented local inputs.
