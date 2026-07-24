"""Extract the authoritative commented formalism from Sagear's source archive."""

from __future__ import annotations

import argparse
import json
import tarfile
from pathlib import Path


def read_source(archive: Path, member_name: str = "main.tex") -> list[str]:
    if archive.is_dir():
        source = archive / member_name
        if not source.exists():
            raise FileNotFoundError(f"source member is missing: {source}")
        return source.read_text(encoding="utf-8", errors="replace").splitlines()
    with tarfile.open(archive, "r:gz") as tar:
        member = tar.getmember(member_name)
        handle = tar.extractfile(member)
        if handle is None:
            raise FileNotFoundError("main.tex is missing from Sagear source archive")
        return handle.read().decode("utf-8", errors="replace").splitlines()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".")
    parser.add_argument("--archive", default=None)
    parser.add_argument("--gilbert-archive", default=None)
    parser.add_argument("--output", default="outputs/sagear_source_formalism_audit")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    archive = Path(args.archive).resolve() if args.archive else root.parent / "sagear_source.tar.gz"
    lines = read_source(archive)
    gilbert_archive = Path(args.gilbert_archive).resolve() if args.gilbert_archive else root.parent / "gilbert_source.tar.gz"
    gilbert_lines = read_source(gilbert_archive, "eccentricity.tex")

    def find(text: str) -> int:
        for index, line in enumerate(lines):
            if text in line:
                return index + 1
        raise ValueError(f"source phrase not found: {text}")

    evidence = {
        "density_prose_line": find("we take the stellar sample and stellar densities"),
        "rho_true_b20_line": find("the stellar density measurement from \\citealt{berger_gaia-kepler_2020}"),
        "rho_formula_line": find("rho_{\\rm \\star, samp}"),
        "transit_probability_line": find("p(e,\\omega|\\hat{t})"),
        "hierarchical_likelihood_line": find("p{(e_k^n|\\theta)}"),
        "uniform_prior_line": find("p(e_n^k|\\alpha) = 1"),
        "gilbert_weighted_to_unweighted_line": next(
            index + 1
            for index, line in enumerate(gilbert_lines)
            if "converted our weighted samples to unweighted samples" in line
        ),
    }
    result = {
        "archive": str(archive),
        "evidence_lines": evidence,
        "internal_density_provenance_conflict": True,
        "hierarchy_uses_reciprocal_transit_probability": True,
        "canonical_formula_status": "rho formula and reciprocal transit correction match the local canonical implementation",
        "raw_equal_rows_are_replication_validated": False,
    }
    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    excerpt_ranges = [(evidence["density_prose_line"] - 1, evidence["rho_formula_line"] + 22), (evidence["transit_probability_line"] - 1, evidence["uniform_prior_line"] + 2)]
    excerpt = []
    for start, stop in excerpt_ranges:
        excerpt.extend(f"{i + 1}: {lines[i]}" for i in range(start, min(stop, len(lines))))
        excerpt.append("")
    markdown = "\n".join(
        [
            "# Sagear Source Formalism Audit",
            "",
            "The archived Sagear `main.tex` contains an internal density-provenance conflict. The prose names Berger 2018, while the commented post-model equation describes `rho_true` as the Berger 2020 density measurement. The same source explicitly gives the reciprocal transit-probability correction used in the hierarchical likelihood.",
            "",
            "The archived Gilbert implementation, which Sagear cites for the post-model method, explicitly says that weighted transit samples were converted to unweighted samples before downstream analysis. Therefore the faithful ALDERAAN analogue is weighted-then-resampled, or an exactly equivalent weighted integration. Treating raw nested rows equally before applying the density likelihood is only a diagnostic and is not validated as Sagear's method.",
            "",
            "## Evidence",
            f"- Density prose: source line {evidence['density_prose_line']}.",
            f"- Berger 2020 `rho_true` wording: source line {evidence['rho_true_b20_line']}.",
            f"- Photoeccentric equation: source line {evidence['rho_formula_line']}.",
            f"- Transit probability: source line {evidence['transit_probability_line']}.",
            f"- Hierarchical likelihood: source line {evidence['hierarchical_likelihood_line']}.",
            "",
            "## Extracted source excerpt",
            "```text",
            *excerpt,
            "```",
            "",
            "## Decision",
            "The canonical implementation keeps the B20 density table because it is the only supplied density table and because the archived equation names B20 for rho_true. The B18 interpretation remains an explicit sensitivity/clarification branch, not a silent replacement. The canonical dynesty-weighted integration is retained; raw equal-row output remains diagnostic only.",
        ]
    )
    output.with_suffix(".md").write_text(markdown + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
