from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest


MODULE = Path(__file__).parents[2] / "scripts" / "direct_chord_eccentricity_diagnostic.py"
SPEC = importlib.util.spec_from_file_location("direct_chord", MODULE)
direct_chord = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(direct_chord)


def test_duration_is_circular_at_zero_eccentricity() -> None:
    circular = direct_chord.finite_duration_days(31.58, 0.042, 0.2, 0.0, 0.0, 1.02)
    assert 0 < circular < 1


def test_higher_photoeccentric_velocity_shortens_duration() -> None:
    circular = direct_chord.finite_duration_days(31.58, 0.042, 0.2, 0.0, 0.0, 1.02)
    fast = direct_chord.finite_duration_days(31.58, 0.042, 0.2, 0.4, np.pi / 2, 1.02)
    assert fast < circular


def test_invalid_chord_is_rejected() -> None:
    with pytest.raises(direct_chord.DirectChordError):
        direct_chord.finite_duration_days(10.0, 0.1, 1.2, 0.0, 0.0, 1.0)


def test_native_mapping_preserves_layout() -> None:
    direct = direct_chord.direct_prior_transform(np.full(8, 0.5))
    native = direct_chord.native_alderaan_vector(direct, 20.0, 1.0)
    assert native.shape == (7,)
    assert np.isclose(native[0], direct[0])
    assert np.isclose(native[1], direct[1])
    assert np.isclose(native[2], direct[2])
    assert np.isclose(native[3], direct[3])
    assert np.isclose(native[-2], direct[-2])
    assert np.isclose(native[-1], direct[-1])


def test_direct_transform_uses_normal_ephemeris_support() -> None:
    center = direct_chord.direct_prior_transform(np.full(8, 0.5))
    assert np.isclose(center[0], 0.0)
    assert np.isclose(center[1], 0.0)
    assert 0.0 < center[4] < direct_chord.ECCENTRICITY_MAX


def test_preflight_hashes_a_complete_one_planet_fixture(tmp_path: Path) -> None:
    result = tmp_path / "Results" / "run" / "K00001"
    result.mkdir(parents=True)
    (result / "K00001_transit_parameters.csv").write_text(
        ",koi_id,kic_id,npl\n0,K00001,1,1\n", encoding="utf-8"
    )
    (result / "K00001_lc_filtered.fits").write_bytes(b"fixture")
    (result / "K00001_00_quick.ttvs").write_text("0 1 1\n1 2 2\n", encoding="utf-8")
    manifest = direct_chord.fixture_preflight(tmp_path, "run", "K00001", 1.0)
    assert manifest["catalog_planets"] == 1
    assert len(manifest["required_inputs"]) == 3
