# Final Forensic Checks

These compact tables record the final public-data checks completed on
2026-08-08. The full per-planet ledgers are regenerable and are intentionally
not tracked here.

The independent Gilbert control uses real ALDERAAN nested samples, the released
Gilbert stellar catalog, paired transit-shape samples, dynesty `LN_WT` weights,
and the published small-planet quality cuts. Its Beta mean eccentricity is
0.0487 for the full small-planet sample and 0.0684 for observed singles,
consistent with the corresponding published values near 0.05 and 0.073.

The Sagear residual audit shows that the reconstructed 2,465 planets have the
right total but not the paper's final population membership. The public
reconstruction contains 1,109 thin singles, 269 thick singles, 878 thin
multis, and 209 thick multis, compared with 1,121, 275, 862, and 207 in the
paper. The identities of the paper's final visual or posterior-quality
rejections have not been released.

The nested-weighting audit compares two explicit diagnostics. Proper dynesty
weighting disagrees with the Sagear means and favors Beta-like shapes. Treating
raw nested points equally moves some means toward the paper but favors a
half-Gaussian shape, contrary to the paper's reported Rayleigh model ordering.
Equal-row results are therefore a forensic clue, not a valid replacement for
nested-sampling weights.

Files in this directory are produced by the corresponding scripts in
`scripts/`. See `docs/replication_status.md` for interpretation and the
remaining author-dependent inputs.
