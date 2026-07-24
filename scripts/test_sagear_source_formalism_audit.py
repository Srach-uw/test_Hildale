from pathlib import Path

from sagear_source_formalism_audit import read_source


ROOT = Path(__file__).resolve().parents[1]


def test_archived_source_contains_density_conflict_and_hbm_equation():
    lines = "\n".join(read_source(ROOT / "reference" / "sources" / "sagear_source.tar.gz"))
    assert "berger_revised_2018" in lines
    assert "berger_gaia-kepler_2020" in lines
    assert "p(e,\\omega|\\hat{t})" in lines
    assert "p{(e_k^n|\\theta)}" in lines


def test_gilbert_source_requires_weighted_then_unweighted_samples():
    lines = "\n".join(read_source(ROOT / "reference" / "sources" / "gilbert_source.tar.gz", "eccentricity.tex"))
    assert "converted our weighted samples to unweighted samples" in lines
