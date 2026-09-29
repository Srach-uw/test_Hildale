from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest
from astropy.io import fits


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import diagnostic_transit_injection as injection  # noqa: E402
import prepare_diagnostic_injection_input as preparation  # noqa: E402


def install_fake_lightkurve(monkeypatch):
    calls = []

    class Flags:
        @staticmethod
        def create_quality_mask(quality, *, bitmask):
            calls.append((np.asarray(quality).copy(), bitmask))
            return np.asarray(quality) == 0

    class FakeLightkurve:
        __version__ = "2.1.1"
        KeplerQualityFlags = Flags

    monkeypatch.setattr(preparation, "_load_pinned_lightkurve", lambda: FakeLightkurve)
    return calls


def write_pdcsap(
    path: Path,
    *,
    kic: int = 123456,
    quarter: int = 4,
    quarter_location: str = "extension",
    obsmode: str = "long cadence",
    obsmode_location: str = "primary",
    timedel_days: float = 0.02043359821692,
    timedel_location: str = "extension",
    time=None,
    quality=None,
):
    time = np.array([100.0, 100.05, 100.10, 100.15, 100.20]) if time is None else np.asarray(time)
    quality = (
        np.array([0, 0, 128, 0, 0], dtype=np.int32)
        if quality is None and len(time) == 5
        else np.zeros(len(time), dtype=np.int32)
        if quality is None
        else np.asarray(quality)
    )
    columns = [
        fits.Column(name="TIME", format="D", array=time),
        fits.Column(name="PDCSAP_FLUX", format="D", array=np.ones(len(time))),
        fits.Column(name="PDCSAP_FLUX_ERR", format="D", array=np.full(len(time), 0.001)),
        fits.Column(name="CADENCENO", format="K", array=np.arange(1000, 1000 + len(time))),
        fits.Column(name="SAP_QUALITY", format="J", array=quality),
    ]
    primary = fits.PrimaryHDU()
    primary.header["KEPLERID"] = kic
    primary.header["BJDREFI"] = 2454833
    primary.header["BJDREFF"] = 0.0
    primary.header["TIMEUNIT"] = "d"
    primary.header["TIMESYS"] = "TDB"
    if obsmode_location in {"primary", "both"}:
        primary.header["OBSMODE"] = obsmode
    if timedel_location in {"primary", "both"}:
        primary.header["TIMEDEL"] = timedel_days
    if quarter_location in {"primary", "both"}:
        primary.header["QUARTER"] = quarter
    lightcurve = fits.BinTableHDU.from_columns(columns, name="LIGHTCURVE")
    if quarter_location in {"extension", "both"}:
        lightcurve.header["QUARTER"] = quarter
    lightcurve.header["BJDREFI"] = 2454833
    lightcurve.header["BJDREFF"] = 0.0
    lightcurve.header["TIMEUNIT"] = "d"
    lightcurve.header["TIMESYS"] = "TDB"
    if obsmode_location in {"extension", "both"}:
        lightcurve.header["OBSMODE"] = obsmode
    if timedel_location in {"extension", "both"}:
        lightcurve.header["TIMEDEL"] = timedel_days
    fits.HDUList([primary, lightcurve]).writeto(path)
    return path


def timing_hdu(index: int, model: list[float], *, values=None):
    model_values = np.asarray(model, dtype=float)
    if values is None:
        values = {
            "INDEX": np.arange(len(model_values)),
            "TTIME": model_values + 0.001,
            "MODEL": model_values,
            "OUT_PROB": np.zeros(len(model_values)),
            "OUT_FLAG": np.zeros(len(model_values)),
        }
    return fits.BinTableHDU.from_columns(
        [
            fits.Column(name="INDEX", format="D", array=np.asarray(values["INDEX"], dtype=float)),
            fits.Column(name="TTIME", format="D", array=np.asarray(values["TTIME"], dtype=float)),
            fits.Column(name="MODEL", format="D", array=np.asarray(values["MODEL"], dtype=float)),
            fits.Column(name="OUT_PROB", format="D", array=np.asarray(values["OUT_PROB"], dtype=float)),
            fits.Column(name="OUT_FLAG", format="D", array=np.asarray(values["OUT_FLAG"], dtype=float)),
        ],
        name=f"TTIMES_{index:02d}",
    )


def write_result(path: Path, *, target: str = "K00367", models=((100.05, 101.05), (100.15, 101.15))):
    primary = fits.PrimaryHDU()
    primary.header["TARGET"] = target
    primary.header["NPL"] = len(models)
    fits.HDUList([primary, *[timing_hdu(index, list(model)) for index, model in enumerate(models)]]).writeto(path)
    return path


def write_catalog(path: Path, *, target: str = "K00367", kic: int = 123456, npl: int = 2):
    rows = [
        f"{target},{kic},{npl},{index + 1}.0,{100.0 + index},6.0"
        for index in range(npl)
    ]
    path.write_text(
        "koi_id,kic_id,npl,period,epoch,duration\n" + "\n".join(rows) + "\n",
        encoding="utf-8",
    )
    return path


def prepare_arguments(pdcsap, result, catalog, **changes):
    values = dict(
        pdcsap_paths=pdcsap,
        result_path=result,
        system_catalog_path=catalog,
        target="K00367",
        kic=123456,
        duration_days=0.25,
        duration_multiple=1.0,
        timing_uncertainty_margin_days=0.01,
        cadence_mode="long",
        exposure_minutes=29.4,
        timing_time_reference="BKJD",
        timing_timeunit="d",
        timing_timesys="TDB",
    )
    values.update(changes)
    return preparation.prepare_input(**values)


def test_preparation_preserves_provenance_and_satisfies_injection_contract(tmp_path, monkeypatch):
    calls = install_fake_lightkurve(monkeypatch)
    first = write_pdcsap(tmp_path / "q4.fits", quarter=4)
    second = write_pdcsap(tmp_path / "q5.fits", quarter=5, time=np.array([102.25, 102.30]))
    result = write_result(tmp_path / "result.fits")
    catalog = write_catalog(tmp_path / "system.csv")

    arrays, manifest = prepare_arguments([first, second], result, catalog)
    output, manifest_path = preparation.write_prepared_input(arrays, manifest, tmp_path / "prepared")

    assert len(calls) == 2
    assert all(bitmask == "default" for _, bitmask in calls)
    assert arrays["time"].tolist() == [100.0, 100.05, 100.15, 100.2, 102.25, 102.3]
    assert arrays["known_transit_mask"].tolist() == [True, True, True, True, False, False]
    assert arrays["mask_index"].tolist() == [4, 5]
    assert arrays["source_file_index"].tolist() == [0, 0, 0, 0, 1, 1]
    assert arrays["source_row_index"].tolist() == [0, 1, 3, 4, 0, 1]
    assert np.array_equal(arrays["flux"], np.ones(6))
    assert injection.load_input_contract(output)["mask_index"].tolist() == [4, 5]

    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert payload["cadence_counts"] == {
        "raw_total": 7,
        "quality_pass_total": 6,
        "finite_total": 7,
        "retained_total": 6,
        "known_transit_excluded_total": 4,
        "injection_eligible_total": 2,
    }
    assert payload["known_transit_exclusion"]["result_fits"]["npl"] == 2
    timing_tables = payload["known_transit_exclusion"]["result_fits"]["timing_tables"]
    assert len(timing_tables) == 2
    assert set(timing_tables[0]) == {
        "extension", "index", "model", "out_flag", "out_prob", "planet_index", "ttime"
    }
    assert payload["output"]["sha256"] == hashlib.sha256(output.read_bytes()).hexdigest()
    assert payload["system_catalog"]["npl"] == 2
    assert payload["known_transit_exclusion"]["half_width_days"] == pytest.approx(0.26)
    assert payload["time_convention"] == {
        "time_reference": "BKJD", "timeunit": "d", "timesys": "TDB", "bjdref": 2454833.0, "timezero": 0.0
    }
    assert payload["cadence"] == {
        "mode": "long",
        "exposure_minutes": 29.4,
        "observed_mode": "long",
        "timedel_days": pytest.approx(0.02043359821692),
        "observed_exposure_minutes": pytest.approx(29.4243814323648),
        "exposure_match_tolerance_minutes": 0.05,
    }


def test_primary_only_quarter_zero_is_retained_and_contradictions_fail(tmp_path, monkeypatch):
    install_fake_lightkurve(monkeypatch)
    primary_only = write_pdcsap(
        tmp_path / "q0-primary.fits", quarter=0, quarter_location="primary"
    )
    record = preparation.read_pdcsap_file(primary_only, 123456, 0)
    assert record["quarter"] == 0
    assert np.array_equal(record["quarter_array"], np.zeros(record["retained_count"], dtype=np.int16))

    contradictory = write_pdcsap(
        tmp_path / "contradictory-quarter.fits", quarter=0, quarter_location="both"
    )
    with fits.open(contradictory, mode="update") as hdul:
        hdul[1].header["QUARTER"] = 1
    with pytest.raises(preparation.PreparationInputError, match="contradictory primary and LIGHTCURVE QUARTER"):
        preparation.read_pdcsap_file(contradictory, 123456, 0)


def test_cadence_and_exposure_headers_must_match_the_declaration(tmp_path, monkeypatch):
    install_fake_lightkurve(monkeypatch)
    result = write_result(tmp_path / "result.fits")
    catalog = write_catalog(tmp_path / "system.csv")

    short = write_pdcsap(tmp_path / "short.fits", obsmode="short cadence")
    with pytest.raises(preparation.PreparationInputError, match="do not match declared cadence_mode"):
        prepare_arguments([short], result, catalog)

    incorrect_exposure = write_pdcsap(tmp_path / "incorrect-exposure.fits")
    with pytest.raises(preparation.PreparationInputError, match="does not match PDCSAP TIMEDEL"):
        prepare_arguments([incorrect_exposure], result, catalog, exposure_minutes=20.0)

    first = write_pdcsap(tmp_path / "first.fits")
    second = write_pdcsap(tmp_path / "second.fits", timedel_days=0.0200)
    with pytest.raises(preparation.PreparationInputError, match="inconsistent TIMEDEL"):
        prepare_arguments([first, second], result, catalog)

    contradictory_mode = write_pdcsap(
        tmp_path / "contradictory-mode.fits", obsmode_location="both"
    )
    with fits.open(contradictory_mode, mode="update") as hdul:
        hdul[1].header["OBSMODE"] = "short cadence"
    with pytest.raises(preparation.PreparationInputError, match="contradictory primary and LIGHTCURVE OBSMODE"):
        preparation.read_pdcsap_file(contradictory_mode, 123456, 0)

    contradictory_timedelta = write_pdcsap(
        tmp_path / "contradictory-timedelta.fits", timedel_location="both"
    )
    with fits.open(contradictory_timedelta, mode="update") as hdul:
        hdul[1].header["TIMEDEL"] = 0.0200
    with pytest.raises(preparation.PreparationInputError, match="contradictory primary and LIGHTCURVE TIMEDEL"):
        preparation.read_pdcsap_file(contradictory_timedelta, 123456, 0)


def test_target_and_kic_mismatches_are_rejected(tmp_path, monkeypatch):
    install_fake_lightkurve(monkeypatch)
    pdcsap = write_pdcsap(tmp_path / "q4.fits")
    result = write_result(tmp_path / "result.fits", target="K00001")
    with pytest.raises(preparation.PreparationInputError, match="TARGET=K00001"):
        prepare_arguments([pdcsap], result, write_catalog(tmp_path / "system.csv"))
    correct_result = write_result(tmp_path / "correct-result.fits")
    with pytest.raises(preparation.PreparationInputError, match="has KIC 123456"):
        prepare_arguments([pdcsap], correct_result, write_catalog(tmp_path / "other.csv", kic=999999), kic=999999)


def test_preparation_fails_without_pinned_lightkurve(tmp_path, monkeypatch):
    monkeypatch.setattr(
        preparation,
        "_load_pinned_lightkurve",
        lambda: (_ for _ in ()).throw(preparation.LightkurveDependencyError("missing pinned lightkurve")),
    )
    with pytest.raises(preparation.LightkurveDependencyError, match="missing pinned lightkurve"):
        prepare_arguments(
            [write_pdcsap(tmp_path / "q4.fits")],
            write_result(tmp_path / "result.fits"),
            write_catalog(tmp_path / "system.csv"),
        )


def test_output_immutability_and_preflight_contract_failures(tmp_path, monkeypatch):
    install_fake_lightkurve(monkeypatch)
    arrays, manifest = prepare_arguments(
        [write_pdcsap(
            tmp_path / "q4.fits",
            time=np.array([100.0, 100.05, 100.10, 100.15, 100.20, 102.0]),
            quality=np.array([0, 0, 128, 0, 0, 0], dtype=np.int32),
        )],
        write_result(tmp_path / "result.fits", models=((100.05, 101.05),)),
        write_catalog(tmp_path / "system.csv", npl=1),
    )
    output_dir = tmp_path / "prepared"
    preparation.write_prepared_input(arrays, manifest, output_dir)
    with pytest.raises(preparation.PreparationInputError, match="absent or empty"):
        preparation.write_prepared_input(arrays, manifest, output_dir)

    with pytest.raises(preparation.PreparationInputError, match="strictly increasing"):
        preparation.conservative_transit_mask(np.array([2.0, 1.0]), np.array([1.5]), 0.1)
    with pytest.raises(preparation.PreparationInputError, match="finite and positive"):
        preparation.conservative_transit_mask(np.array([1.0, 2.0]), np.array([1.5]), 0.0)


@pytest.mark.parametrize(
    ("values", "message"),
    [
        ({"INDEX": [-1, 1], "TTIME": [100.0, 101.0], "MODEL": [100.0, 101.0], "OUT_PROB": [0, 0], "OUT_FLAG": [0, 0]}, "INDEX must contain"),
        ({"INDEX": [0.5, 1], "TTIME": [100.0, 101.0], "MODEL": [100.0, 101.0], "OUT_PROB": [0, 0], "OUT_FLAG": [0, 0]}, "INDEX must contain"),
        ({"INDEX": [0, 0], "TTIME": [100.0, 101.0], "MODEL": [100.0, 101.0], "OUT_PROB": [0, 0], "OUT_FLAG": [0, 0]}, "INDEX must be strictly"),
        ({"INDEX": [0, 1], "TTIME": [100.0, 99.0], "MODEL": [100.0, 101.0], "OUT_PROB": [0, 0], "OUT_FLAG": [0, 0]}, "TTIME must be strictly"),
        ({"INDEX": [0, 1], "TTIME": [100.0, 101.0], "MODEL": [100.0, 99.0], "OUT_PROB": [0, 0], "OUT_FLAG": [0, 0]}, "MODEL must be strictly"),
        ({"INDEX": [0, 1], "TTIME": [100.0, 101.0], "MODEL": [100.0, 101.0], "OUT_PROB": [0, 1.1], "OUT_FLAG": [0, 0]}, "OUT_PROB"),
        ({"INDEX": [0, 1], "TTIME": [100.0, 101.0], "MODEL": [100.0, 101.0], "OUT_PROB": [0, 0], "OUT_FLAG": [0, 2]}, "OUT_FLAG"),
    ],
)
def test_timing_validation_matches_readiness_contract(tmp_path, values, message):
    result = tmp_path / "result.fits"
    primary = fits.PrimaryHDU()
    primary.header["TARGET"] = "K00367"
    primary.header["NPL"] = 1
    fits.HDUList([primary, timing_hdu(0, [100.0, 101.0], values=values)]).writeto(result)
    with pytest.raises(preparation.PreparationInputError, match=message):
        preparation.read_result_timings(result, "K00367", "BKJD", "d", "TDB")


def test_catalog_and_time_reference_gates_are_explicit(tmp_path, monkeypatch):
    install_fake_lightkurve(monkeypatch)
    pdcsap = write_pdcsap(tmp_path / "q4.fits")
    result = write_result(tmp_path / "result.fits")
    with pytest.raises(preparation.PreparationInputError, match="has 1 rows"):
        prepare_arguments([pdcsap], result, write_catalog(tmp_path / "incomplete.csv", npl=1))
    with pytest.raises(preparation.PreparationInputError, match="maps K00367 to KIC values"):
        prepare_arguments([pdcsap], result, write_catalog(tmp_path / "wrong-kic.csv", kic=999999))
    ambiguous = write_catalog(tmp_path / "ambiguous.csv")
    with ambiguous.open("a", encoding="utf-8") as stream:
        stream.write("K99999,123456,1,4.0,104.0,6.0\n")
    with pytest.raises(preparation.PreparationInputError, match="maps KIC 123456 to KOI identifiers"):
        prepare_arguments([pdcsap], result, ambiguous)

    with fits.open(pdcsap, mode="update") as hdul:
        del hdul[0].header["TIMEUNIT"]
        del hdul[1].header["TIMEUNIT"]
    with pytest.raises(preparation.PreparationInputError, match="missing standard TIMEUNIT"):
        prepare_arguments([pdcsap], result, write_catalog(tmp_path / "system.csv"))

    pdcsap = write_pdcsap(tmp_path / "contradictory.fits")
    with fits.open(pdcsap, mode="update") as hdul:
        hdul[1].header["BJDREFI"] = 2454834
    with pytest.raises(preparation.PreparationInputError, match="contradictory primary and LIGHTCURVE BJDREF"):
        prepare_arguments([pdcsap], result, write_catalog(tmp_path / "same-system.csv"))
