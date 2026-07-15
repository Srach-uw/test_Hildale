"""Optional sensitivity patch matching the priors printed in Sagear Table 1.

Public ALDERAAN uses C0,C1 ~ N(0,0.1) and log-uniform Rp/R*. Sagear Table 1
prints N(0,1) and uniform Rp/R*. Because it is unknown whether the manuscript
or the private analysis code is authoritative, this patch is for a labeled
validation arm only and must not replace the pinned public-code arm.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("repo")
    args = parser.parse_args()
    path = Path(args.repo).resolve() / "alderaan" / "dynesty_helpers.py"
    text = path.read_text(encoding="utf-8")
    replacements = [
        ("norm_ppf(u_[0 + npl * 5], 0.0, 0.1)", "norm_ppf(u_[0 + npl * 5], 0.0, 1.0)"),
        ("norm_ppf(u_[1 + npl * 5], 0.0, 0.1)", "norm_ppf(u_[1 + npl * 5], 0.0, 1.0)"),
        ("loguniform_ppf(u_[2 + npl * 5], 1e-5, 0.99)", "uniform_ppf(u_[2 + npl * 5], 1e-5, 0.99)"),
    ]
    for old, new in replacements:
        if new in text:
            continue
        if old not in text:
            raise RuntimeError(f"Expected prior-transform context not found: {old}")
        text = text.replace(old, new, 1)
    path.write_text(text, encoding="utf-8")
    print(f"Applied Sagear Table 1 prior sensitivity patch to {path}")


if __name__ == "__main__":
    main()
