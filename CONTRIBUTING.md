# Contributing

Changes should preserve the distinction between reproduction, sensitivity
analysis, and exploratory diagnosis.

## Before opening a change

```bash
python -m pytest -q
python scripts/check_professor_release.py
git diff --check
```

## Scientific expectations

- State the source and version of every external catalog.
- Assign multiplicity before planet-level exclusions.
- Preserve paired posterior samples and their nested-sampling weights.
- Keep source-faithful and alternative model conventions in separate outputs.
- Record exclusions with stable identifiers and explicit reasons.
- Add a focused test when changing a formula, join, prior, or quality rule.
- Do not tune thresholds after viewing the published answer.

## Repository hygiene

Commit compact evidence and manifests. Keep light curves, full posterior
archives, credentials, personal paths, and chronological worklogs outside the
repository. Use Git LFS only for the curated validation FITS already covered by
the repository's data-availability policy.

Third-party article and data products retain their original terms. New project
code is covered by the MIT license.
