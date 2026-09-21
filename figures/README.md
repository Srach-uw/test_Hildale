# Figures

`population_comparison.png` compares the published Table 3 Rayleigh means with
the deterministic-QC public reconstruction. It is generated from the compact,
versioned table at `metadata/public_reconstruction_20260727/population_comparison.csv`:

```bash
python scripts/plot_population_comparison.py
```

The image describes a replication discrepancy. It is not a new population
measurement and should be read with `docs/replication_status.md` and
`docs/current_inference_audit.md`.

Published figure copies are retained under `reference/` for source comparison.
