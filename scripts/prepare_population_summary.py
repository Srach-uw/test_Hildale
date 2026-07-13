from __future__ import annotations

"""Attach one audited sample/classifier context to uniform posterior products."""

import argparse
from pathlib import Path

import pandas as pd


def attach_context(posteriors: pd.DataFrame, sample: pd.DataFrame) -> pd.DataFrame:
    required_post = {"kepoi_name", "posterior_file", "posterior_source", "impact_mode"}
    required_sample = {"kepoi_name", "koi_target", "kepid", "disk", "system"}
    missing_post = sorted(required_post - set(posteriors.columns))
    missing_sample = sorted(required_sample - set(sample.columns))
    if missing_post:
        raise ValueError(f"Posterior summary lacks required provenance fields: {missing_post}")
    if missing_sample:
        raise ValueError(f"Sample context lacks required fields: {missing_sample}")
    if posteriors["kepoi_name"].duplicated().any():
        raise ValueError("Posterior summary contains duplicate kepoi_name rows")
    context_columns = ["kepoi_name", "koi_target", "kepid", "disk", "system"]
    context_columns += [
        column
        for column in (
            "berger_logg",
            "berger_rad",
            "berger2018_evol",
            "berger2018_bin",
            "rho_log",
            "koi_disposition",
            "koi_model_snr",
            "koi_prad",
        )
        if column in sample.columns
    ]
    context = sample[context_columns].copy()
    if context["kepoi_name"].duplicated().any():
        raise ValueError("Sample context contains duplicate kepoi_name rows")
    drop = [column for column in context_columns if column != "kepoi_name" and column in posteriors]
    out = posteriors.drop(columns=drop).merge(context, on="kepoi_name", how="inner", validate="one_to_one")
    out["sample_context_rows_dropped"] = len(posteriors) - len(out)
    return out.sort_values(["disk", "system", "koi_target", "kepoi_name"]).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Attach audited disk/system labels to a uniform posterior summary.")
    parser.add_argument("--posteriors", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = attach_context(pd.read_csv(args.posteriors), pd.read_csv(args.sample))
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    print(f"Wrote {len(out)} uniformly constructed posterior rows with audited context: {path}")
    print(out.groupby(["disk", "system"]).size().to_string())


if __name__ == "__main__":
    main()
