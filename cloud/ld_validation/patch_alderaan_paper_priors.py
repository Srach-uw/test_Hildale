"""Optional sensitivity patch matching the priors printed in Sagear Table 1.

Public ALDERAAN uses C0,C1 ~ N(0,0.1) and log-uniform Rp/R*. Sagear Table 1
prints N(0,1) and uniform Rp/R*. Because it is unknown whether the manuscript
or the private analysis code is authoritative, this patch is for a labeled
validation arm only and must not replace the pinned public-code arm.
"""

from __future__ import annotations

import argparse
import re
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
        old_pattern = r"(?<![\w])" + re.escape(old)
        new_pattern = r"(?<![\w])" + re.escape(new)
        old_count = len(re.findall(old_pattern, text))
        new_count = len(re.findall(new_pattern, text))
        if old_count == 0 and new_count == 1:
            continue
        if old_count != 1 or new_count != 0:
            raise RuntimeError(f"Expected prior-transform context not found: {old}")
        text = re.sub(old_pattern, lambda _: new, text, count=1)
    path.write_text(text, encoding="utf-8")
    print(f"Applied Sagear Table 1 prior sensitivity patch to {path}")


if __name__ == "__main__":
    main()
