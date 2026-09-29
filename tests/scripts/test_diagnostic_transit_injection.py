from __future__ import annotations

import hashlib
import importlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import diagnostic_transit_injection as injection  # noqa: E402


def valid_spec(**changes) -> injection.InjectionSpec:
    values = dict(
        period_days=10.0,
        rho_star_solar=1.15,
        radius_ratio=0.03,
        impact_parameter=0.4,
        t0_days=100.0,
        limb_darkening_u1=0.3,
        limb_darkening_u2=0.2,
        exposure_minutes=29.4,
        supersample_factor=15,
        seed=17,
    )
    values.update(changes)
    return injection.InjectionSpec(**values)


def write_input(path: Path, *, mask_index: np.ndarray | None = None, **optional: np.ndarray) -> Path:
    time = np.linspace(99.7, 100.3, 11, dtype=np.float64)
    flux = np.linspace(0.98, 1.02, 11, dtype=np.float64)
    np.savez(
        path,
        time=time,
        flux=flux,
        mask_index=np.array([2, 4, 7, 9], dtype=np.int64) if mask_index is None else mask_index,
        **optional,
    )
    return path


def test_circular_geometry_round_trip_is_consistent():
    a_over_rstar = injection.circular_a_over_rstar(12.5, 0.87)
    assert a_over_rstar > 1
    assert injection.circular_rho_star_solar(12.5, a_over_rstar) == pytest.approx(0.87, rel=1e-12)


def test_circular_duration_matches_the_fixture_geometry():
    period_days = 12.5
    rho_star_solar = 0.87
    radius_ratio = 0.04
    impact_parameter = 0.5
    a_over_rstar = injection.circular_a_over_rstar(period_days, rho_star_solar)
    sin_inclination = np.sqrt(1.0 - (impact_parameter / a_over_rstar) ** 2)
    argument = np.sqrt((1.0 + radius_ratio) ** 2 - impact_parameter**2) / (a_over_rstar * sin_inclination)
    expected = period_days / np.pi * np.arcsin(argument)
    assert injection.circular_first_to_fourth_duration_days(
        period_days, rho_star_solar, radius_ratio, impact_parameter
    ) == pytest.approx(expected, rel=1e-14)


def test_renderer_supports_the_public_orbital_batman_contract(monkeypatch):
    captured = {}

    class Params:
        t0 = per = rp = a = inc = ecc = w = u = limb_dark = None

    class Model:
        def __init__(self, params, time, **kwargs):
            captured["params"] = params
            captured["time"] = time
            captured["kwargs"] = kwargs

        def light_curve(self, _params):
            return np.ones(len(captured["time"]))

    class Batman:
        __version__ = "public-fake"
        TransitParams = Params
        TransitModel = Model

    spec = valid_spec()
    monkeypatch.setattr(injection, "_batman_module", lambda: Batman)
    model, version = injection.exposure_integrated_circular_model(spec, np.array([99.9, 100.0, 100.1]))
    params = captured["params"]
    a_over_rstar = injection.circular_a_over_rstar(spec.period_days, spec.rho_star_solar)
    assert version == "public-fake"
    assert np.array_equal(model, np.ones(3))
    assert params.a == pytest.approx(a_over_rstar)
    assert params.inc == pytest.approx(np.degrees(np.arccos(spec.impact_parameter / a_over_rstar)))
    assert params.ecc == 0.0
    assert params.w == 90.0
    assert captured["kwargs"] == {
        "supersample_factor": spec.supersample_factor,
        "exp_time": spec.exposure_minutes / 1440.0,
    }


def test_renderer_supports_the_pinned_alderaan_batman_contract(monkeypatch):
    captured = {}

    class Params:
        __slots__ = ("t0", "per", "rp", "b", "T14", "u", "limb_dark")

    class Model:
        def __init__(self, params, time, **kwargs):
            captured["params"] = params
            captured["time"] = time
            captured["kwargs"] = kwargs

        def light_curve(self, _params):
            return np.ones(len(captured["time"]))

    class Batman:
        __version__ = "alderaan-fake"
        TransitParams = Params
        TransitModel = Model

    spec = valid_spec()
    monkeypatch.setattr(injection, "_batman_module", lambda: Batman)
    model, version = injection.exposure_integrated_circular_model(spec, np.array([99.9, 100.0, 100.1]))
    params = captured["params"]
    assert version == "alderaan-fake"
    assert np.array_equal(model, np.ones(3))
    assert params.b == spec.impact_parameter
    assert params.T14 == pytest.approx(
        injection.circular_first_to_fourth_duration_days(
            spec.period_days, spec.rho_star_solar, spec.radius_ratio, spec.impact_parameter
        )
    )
    assert not hasattr(params, "a")
    assert captured["kwargs"] == {
        "supersample_factor": spec.supersample_factor,
        "exp_time": spec.exposure_minutes / 1440.0,
    }


@pytest.mark.parametrize(
    "changes, message",
    [
        ({"period_days": 0.0}, "period_days"),
        ({"rho_star_solar": np.nan}, "rho_star_solar"),
        ({"impact_parameter": 1.04, "radius_ratio": 0.03}, "transiting limit"),
        ({"radius_ratio": 1.0}, "radius_ratio"),
        ({"limb_darkening_u1": 0.8, "limb_darkening_u2": 0.3}, "limb-darkening"),
        ({"supersample_factor": 0}, "supersample_factor"),
    ],
)
def test_invalid_specifications_fail_clearly(changes, message):
    with pytest.raises(injection.InjectionInputError, match=message):
        valid_spec(**changes).validate()


def test_dry_run_writes_deterministic_immutable_manifest(tmp_path):
    source = write_input(tmp_path / "input.npz")
    first = tmp_path / "first"
    second = tmp_path / "second"
    injection.run_injection(valid_spec(), source, first, dry_run=True)
    injection.run_injection(valid_spec(), source, second, dry_run=True)

    one = (first / injection.MANIFEST_NAME).read_bytes()
    two = (second / injection.MANIFEST_NAME).read_bytes()
    assert one == two
    payload = json.loads(one)
    assert payload["mode"] == "preflight"
    assert payload["software"]["batman"] == "not_loaded_preflight"
    assert payload["input"]["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    with pytest.raises(injection.InjectionInputError, match="immutable manifest"):
        injection._write_json_once(first / injection.MANIFEST_NAME, payload)


def test_contract_rejects_nonfinite_or_nonpreserving_masks(tmp_path):
    source = write_input(tmp_path / "bad.npz", mask_index=np.array([2, 2], dtype=np.int64))
    with pytest.raises(injection.InjectionInputError, match="strictly increasing"):
        injection.load_input_contract(source)
    np.savez(source, time=np.array([1.0, np.nan]), flux=np.ones(2), mask_index=np.array([0], dtype=np.int64))
    with pytest.raises(injection.InjectionInputError, match="finite"):
        injection.load_input_contract(source)
    np.savez(source, time=np.array([1.0, 0.9]), flux=np.ones(2), mask_index=np.array([0], dtype=np.int64))
    with pytest.raises(injection.InjectionInputError, match="strictly increasing"):
        injection.load_input_contract(source)


def test_render_preserves_time_and_mask_identity_with_fake_batman(tmp_path, monkeypatch):
    source = write_input(tmp_path / "input.npz")
    selected_times = np.linspace(99.7, 100.3, 11, dtype=np.float64)[[2, 4, 7, 9]]

    class Params:
        t0 = per = rp = a = inc = ecc = w = u = limb_dark = None

    class Model:
        def __init__(self, _params, time, **kwargs):
            self.time = time
            self.kwargs = kwargs

        def light_curve(self, _params):
            return np.where(np.isin(self.time, selected_times), 0.99, 1.0)

    class Batman:
        __version__ = "fake-1"
        TransitParams = Params
        TransitModel = Model

    monkeypatch.setattr(injection, "_batman_module", lambda: Batman)
    output = tmp_path / "rendered"
    injection.run_injection(valid_spec(), source, output)
    with np.load(source) as original, np.load(output / injection.OUTPUT_NAME) as rendered:
        assert np.array_equal(rendered["time"], original["time"])
        assert np.array_equal(rendered["mask_index"], original["mask_index"])
        outside = np.setdiff1d(np.arange(len(original["flux"])), original["mask_index"])
        assert np.array_equal(rendered["flux"][outside], original["flux"][outside])
        assert np.allclose(rendered["flux"][original["mask_index"]], original["flux"][original["mask_index"]] * 0.99)
    manifest = json.loads((output / injection.MANIFEST_NAME).read_text())
    assert manifest["software"]["batman"] == "fake-1"
    assert manifest["output"]["sha256"] == hashlib.sha256((output / injection.OUTPUT_NAME).read_bytes()).hexdigest()
    assert manifest["injected_transit_coverage"] == {
        "all_in_transit_cadences_eligible": True,
        "eligible_grid_count": 4,
        "eligible_in_transit_count": 4,
        "full_grid_count": 11,
        "in_transit_full_grid_count": 4,
        "in_transit_tolerance": injection.IN_TRANSIT_TOLERANCE,
        "masked_in_transit_count": 0,
    }


def test_render_rejects_an_injected_transit_partially_removed_by_mask(tmp_path, monkeypatch):
    source = write_input(tmp_path / "input.npz", mask_index=np.array([2, 4, 7], dtype=np.int64))

    class Params:
        t0 = per = rp = a = inc = ecc = w = u = limb_dark = None

    class Model:
        def __init__(self, _params, time, **_kwargs):
            self.time = time

        def light_curve(self, _params):
            return np.where(np.isin(self.time, [99.94, 100.12, 100.24]), 0.99, 1.0)

    class Batman:
        __version__ = "fake-1"
        TransitParams = Params
        TransitModel = Model

    monkeypatch.setattr(injection, "_batman_module", lambda: Batman)
    with pytest.raises(injection.InjectionInputError, match="would remove 1 of 3 modeled in-transit cadences"):
        injection.run_injection(valid_spec(), source, tmp_path / "rendered")


def test_render_requires_batman_but_dry_run_does_not(tmp_path, monkeypatch):
    source = write_input(tmp_path / "input.npz")
    real_import = importlib.import_module

    def missing(name: str):
        if name == "batman":
            raise ModuleNotFoundError(name)
        return real_import(name)

    monkeypatch.setattr(injection.importlib, "import_module", missing)
    with pytest.raises(injection.BatmanDependencyError, match="batman-package"):
        injection.run_injection(valid_spec(), source, tmp_path / "rendered")
    injection.run_injection(valid_spec(), source, tmp_path / "preflight", dry_run=True)


def test_render_preserves_optional_provenance_arrays_and_reports_hashes(tmp_path, monkeypatch):
    point_count = 11
    selected_times = np.linspace(99.7, 100.3, point_count, dtype=np.float64)[[2, 4, 7, 9]]
    optional = {
        "flux_err": np.linspace(0.001, 0.002, point_count, dtype=np.float32),
        "cadence": np.arange(100, 100 + point_count, dtype=np.int32),
        "quality": np.arange(point_count, dtype=np.uint32),
        "quarter": np.full(point_count, 7, dtype=np.int16),
        "source_file_index": np.zeros(point_count, dtype=np.int16),
        "source_row_index": np.arange(point_count, dtype=np.int64),
        "known_transit_mask": np.array([False, True] * 5 + [False]),
        "source_files": np.array(["kplr000000001-llc.fits"], dtype="<U24"),
    }
    source = write_input(tmp_path / "input.npz", **optional)

    class Params:
        t0 = per = rp = a = inc = ecc = w = u = limb_dark = None

    class Model:
        def __init__(self, _params, time, **_kwargs):
            self.time = time

        def light_curve(self, _params):
            return np.where(np.isin(self.time, selected_times), 0.99, 1.0)

    class Batman:
        __version__ = "fake-1"
        TransitParams = Params
        TransitModel = Model

    monkeypatch.setattr(injection, "_batman_module", lambda: Batman)
    output = tmp_path / "rendered"
    injection.run_injection(valid_spec(), source, output)
    with np.load(source, allow_pickle=False) as original, np.load(output / injection.OUTPUT_NAME, allow_pickle=False) as rendered:
        for name, values in optional.items():
            assert rendered[name].dtype == original[name].dtype
            assert rendered[name].shape == original[name].shape
            assert rendered[name].tobytes() == original[name].tobytes()
            assert np.array_equal(rendered[name], values)
    manifest = json.loads((output / injection.MANIFEST_NAME).read_text())
    for name, values in optional.items():
        report = manifest["output"]["optional_arrays"][name]
        assert report["dtype"] == values.dtype.str
        assert report["shape"] == list(values.shape)
        assert report["sha256"] == hashlib.sha256(np.ascontiguousarray(values).tobytes()).hexdigest()


@pytest.mark.parametrize("name", injection.POINT_ALIGNED_OPTIONAL_ARRAYS)
def test_contract_rejects_misaligned_optional_point_arrays(tmp_path, name):
    source = write_input(tmp_path / "bad_optional.npz", **{name: np.arange(10, dtype=np.int64)})
    with pytest.raises(injection.InjectionInputError, match=rf"optional array {name} must have length matching time"):
        injection.load_input_contract(source)


def test_contract_allows_source_mapping_with_non_point_length(tmp_path):
    source = write_input(
        tmp_path / "source_mapping.npz",
        source_files=np.array(["quarter1.fits", "quarter2.fits"], dtype="<U16"),
    )
    contract = injection.load_input_contract(source)
    assert np.array_equal(contract["source_files"], np.array(["quarter1.fits", "quarter2.fits"], dtype="<U16"))
