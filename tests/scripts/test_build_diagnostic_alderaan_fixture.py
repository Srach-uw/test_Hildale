import ast
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from build_diagnostic_alderaan_fixture import FixtureInputError, build_fixture


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_fixture_builder_avoids_python39_incompatible_path_write_text_newline():
    source = (SCRIPTS / "build_diagnostic_alderaan_fixture.py").read_text(encoding="utf-8")
    calls = [
        node
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "write_text"
        and any(keyword.arg == "newline" for keyword in node.keywords)
    ]
    assert not calls


def write_source(path):
    columns = [
        fits.Column(name="TIME", format="D", array=np.array([1.0, 2.0, 3.0, 4.0])),
        fits.Column(name="PDCSAP_FLUX", format="D", array=np.array([100.0, 101.0, 102.0, 103.0])),
        fits.Column(name="PDCSAP_FLUX_ERR", format="D", array=np.ones(4)),
        fits.Column(name="CADENCENO", format="K", array=np.arange(4)),
        fits.Column(name="SAP_QUALITY", format="J", array=np.zeros(4, dtype=np.int32)),
    ]
    primary = fits.PrimaryHDU()
    primary.header["OBSMODE"] = "long cadence"
    primary.header["KEPLERID"] = 42
    lightcurve = fits.BinTableHDU.from_columns(columns, name="LIGHTCURVE")
    fits.HDUList([primary, lightcurve]).writeto(path)


def write_inputs(tmp_path):
    source = tmp_path / "kplr000000042-test_llc.fits"
    write_source(source)
    prepared = tmp_path / "prepared.npz"
    injection = tmp_path / "injection.npz"
    np.savez_compressed(prepared, time=np.array([1., 2., 3., 4.]), flux=np.array([100., 101., 102., 103.]), source_file_index=np.zeros(4, dtype=np.int32), source_row_index=np.arange(4, dtype=np.int64), known_transit_mask=np.array([True, False, False, True]), source_files=np.array([str(source.resolve())]))
    np.savez_compressed(injection, time=np.array([1., 2., 3., 4.]), flux=np.array([90., 80., 70., 60.]))
    spec = {"period_days": 10.0, "rho_star_solar": 1.0, "radius_ratio": 0.1, "impact_parameter": 0.2, "t0_days": 2.0, "limb_darkening_u1": 0.3, "limb_darkening_u2": 0.2, "exposure_minutes": 29.4, "supersample_factor": 15, "seed": 1}
    prepared_manifest = tmp_path / "prepared.json"
    injection_manifest = tmp_path / "injection.json"
    prepared_manifest.write_text(json.dumps({"target": "K00001", "kic": 42, "output": {"sha256": sha256(prepared)}, "cadence": {"exposure_minutes": 29.4, "exposure_match_tolerance_minutes": 0.05}, "source_files": [{"path": source.name, "sha256": sha256(source)}]}))
    injection_manifest.write_text(json.dumps({"input": {"sha256": sha256(prepared)}, "output": {"sha256": sha256(injection)}, "specification": spec}))
    catalog = tmp_path / "catalog.csv"
    pd.DataFrame([{ "koi_id": "K00001", "kic_id": 42, "npl": 1, "period": 9., "epoch": 1., "depth": 1., "duration": 2., "impact": 0.1, "limbdark_1": 0.1, "limbdark_2": 0.1 }]).to_csv(catalog)
    return prepared, prepared_manifest, injection, injection_manifest, catalog


def test_builds_isolated_fixture_and_masks_known_transits(tmp_path):
    prepared, prepared_manifest, injection, injection_manifest, catalog = write_inputs(tmp_path)
    result = build_fixture(prepared_path=prepared, prepared_manifest=prepared_manifest, injection_path=injection, injection_manifest=injection_manifest, catalog_path=catalog, target="K00001", output_dir=tmp_path / "fixture")
    assert result.is_file()
    with fits.open(tmp_path / "fixture" / "Data" / "kplr000000042-test_llc.fits") as hdul:
        flux = hdul[1].data["PDCSAP_FLUX"]
        assert np.isnan(flux[0]) and np.isnan(flux[3])
        assert flux[1] == 80.0 and flux[2] == 70.0
        assert hdul[0].header["HILDDIAG"]
    row = pd.read_csv(tmp_path / "fixture" / "diagnostic_catalog.csv", index_col=0).iloc[0]
    assert row.period == 10.0 and row.epoch == 2.0 and row.impact == 0.2
    assert row.depth == pytest.approx(10000.0)
    manifest = json.loads(result.read_text())
    assert manifest["known_transits_are_nan_masked"] is True
    assert manifest["results_written"] is False


def test_rejects_non_single_planet_catalog(tmp_path):
    prepared, prepared_manifest, injection, injection_manifest, catalog = write_inputs(tmp_path)
    table = pd.read_csv(catalog, index_col=0)
    table.loc[table.index[0], "npl"] = 2
    table.to_csv(catalog)
    with pytest.raises(FixtureInputError, match="npl=1"):
        build_fixture(prepared_path=prepared, prepared_manifest=prepared_manifest, injection_path=injection, injection_manifest=injection_manifest, catalog_path=catalog, target="K00001", output_dir=tmp_path / "fixture")


def test_rejects_prepared_flux_that_does_not_match_source(tmp_path):
    prepared, prepared_manifest, injection, injection_manifest, catalog = write_inputs(tmp_path)
    with np.load(prepared) as source:
        payload = {name: source[name] for name in source.files}
    payload["flux"] = payload["flux"].copy()
    payload["flux"][1] += 1
    np.savez_compressed(prepared, **payload)
    source_file = Path(np.load(prepared)["source_files"][0])
    prepared_manifest.write_text(json.dumps({"target": "K00001", "kic": 42, "output": {"sha256": sha256(prepared)}, "cadence": {"exposure_minutes": 29.4, "exposure_match_tolerance_minutes": 0.05}, "source_files": [{"path": source_file.name, "sha256": sha256(source_file)}]}))
    injection_record = json.loads(injection_manifest.read_text())
    injection_record["input"]["sha256"] = sha256(prepared)
    injection_manifest.write_text(json.dumps(injection_record))
    with pytest.raises(FixtureInputError, match="prepared flux"):
        build_fixture(prepared_path=prepared, prepared_manifest=prepared_manifest, injection_path=injection, injection_manifest=injection_manifest, catalog_path=catalog, target="K00001", output_dir=tmp_path / "fixture")


def test_rejects_injection_from_a_different_prepared_input(tmp_path):
    prepared, prepared_manifest, injection, injection_manifest, catalog = write_inputs(tmp_path)
    manifest = json.loads(injection_manifest.read_text())
    manifest["input"]["sha256"] = "0" * 64
    injection_manifest.write_text(json.dumps(manifest))
    with pytest.raises(FixtureInputError, match="does not identify"):
        build_fixture(prepared_path=prepared, prepared_manifest=prepared_manifest, injection_path=injection, injection_manifest=injection_manifest, catalog_path=catalog, target="K00001", output_dir=tmp_path / "fixture")


def test_rejects_inconsistent_declared_exposure(tmp_path):
    prepared, prepared_manifest, injection, injection_manifest, catalog = write_inputs(tmp_path)
    manifest = json.loads(injection_manifest.read_text())
    manifest["specification"]["exposure_minutes"] = 1.0
    injection_manifest.write_text(json.dumps(manifest))
    with pytest.raises(FixtureInputError, match="exposure"):
        build_fixture(prepared_path=prepared, prepared_manifest=prepared_manifest, injection_path=injection, injection_manifest=injection_manifest, catalog_path=catalog, target="K00001", output_dir=tmp_path / "fixture")
