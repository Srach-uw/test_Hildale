# Single-target circular-injection runner

This runner is for the first gated circular recovery only. It accepts an
isolated fixture created by `build_diagnostic_alderaan_fixture.py`, copies it
into a fresh project directory, and runs the three public ALDERAAN stages.
Unlike the production runner, it never calls `get_kepler_data.py` or any
download command.

Run it only in the pinned Linux ALDERAAN environment after the fixture,
catalog, and input manifests have been checked. Use a new `PROJECT_DIR` for
every run. The output is a diagnostic result, not a canonical posterior and
not evidence that the population replication has succeeded.

```bash
FIXTURE_DIR=/path/to/fixture \
PROJECT_DIR=/path/to/fresh_project \
ALDERAAN_REPO=$HOME/alderaan_sagear_pinned \
RUN_ID=diagnostic_k00367_circular \
bash cloud/diagnostic_injection/run_one_fixture.sh K00367 4815520
```

Before interpreting a result, verify the fixture and result hashes in
`provenance_sha256.txt`, confirm the result FITS is nonempty, inspect the fit,
and compare its recovered circular density with the injection specification.
Before launching the runner, execute
`scripts/verify_diagnostic_fixture_reader.py` in the same pinned environment.
Its audit must retain every injection-eligible cadence and no masked real
transit cadence. It reports any extra ALDERAAN-reader cadences explicitly.
