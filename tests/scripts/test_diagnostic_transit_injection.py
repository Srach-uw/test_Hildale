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


def write_input(path: Path, *, mask_index: np.ndarray | None = None) -> Path:
    time = np.linspace(99.7, 100.3, 11, dtype=np.float64)
    flux = np.linspace(0.98, 1.02, 11, dtype=np.float64)
    np.savez(path, time=time, flux=flux, mask_index=np.array([2, 4, 7, 9], dtype=np.int64) if mask_index is None else mask_index)
    return path


def test_circular_geometry_round_trip_is_consistent():
    a_over_rstar = injection.circular_a_over_rstar(12.5, 0.87)
    assert a_over_rstar > 1
    assert injection.circular_rho_star_solar(12.5, a_over_rstar) == pytest.approx(0.87, rel=1e-12)


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

    class Params:
        pass

    class Model:
        def __init__(self, _params, time, **kwargs):
            self.time = time
            self.kwargs = kwargs

        def light_curve(self, _params):
            return np.full(len(self.time), 0.99)

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
