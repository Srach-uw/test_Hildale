import numpy as np
import pandas as pd

from gilbert_postfit_qc_audit import (
    build_result_index,
    catalog_density_solar,
    grazing_exceeds_limit,
    radius_fractional_uncertainty,
    result_file,
    weighted_fraction,
)


def test_gilbert_analysis_notebook_grazing_limit_is_five_percent() -> None:
    assert grazing_exceeds_limit(0.051)
    assert not grazing_exceeds_limit(0.05)
    assert not grazing_exceeds_limit(np.nan)


def test_weighted_fraction_uses_nested_weights() -> None:
    mask = np.array([False, True, True])
    weights = np.array([0.8, 0.1, 0.1])
    assert np.isclose(weighted_fraction(mask, weights), 0.2)


def test_radius_fractional_uncertainty_combines_independent_terms() -> None:
    ror = np.array([0.09, 0.10, 0.11])
    weights = np.ones(3) / 3
    _, median, _, frac = radius_fractional_uncertainty(
        ror,
        weights,
        stellar_radius=1.0,
        radius_err_hi=0.1,
        radius_err_lo=0.1,
    )
    assert np.isclose(median, 0.10)
    assert 0.10 < frac < 0.15


def test_catalog_density_prefers_direct_solar_value() -> None:
    row = pd.Series({"rho_true_solar": 0.7, "rho_log": -2.0})
    assert np.isclose(catalog_density_solar(row), 0.7)


def test_catalog_density_falls_back_to_log_value() -> None:
    row = pd.Series({"rho_true_solar": np.nan, "rho_log": -0.3})
    assert np.isclose(catalog_density_solar(row), 10**-0.3)


def test_result_index_falls_back_across_roots(tmp_path) -> None:
    archive = tmp_path / "archive"
    project = tmp_path / "project" / "nested"
    archive.mkdir()
    project.mkdir(parents=True)
    result = project / "K00001-results.fits"
    result.write_bytes(b"result")

    index = build_result_index(archive, tmp_path / "project")
    found, conflict, count = result_file("K00001", index)

    assert found == result.resolve()
    assert conflict is False
    assert count == 1


def test_result_file_prefers_largest_duplicate(tmp_path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()
    small = first / "K00001-results.fits"
    large = second / "K00001-results.fits"
    small.write_bytes(b"a")
    large.write_bytes(b"larger")

    found, conflict, count = result_file(
        "K00001", build_result_index(first, second)
    )

    assert found == large.resolve()
    assert conflict is True
    assert count == 2
