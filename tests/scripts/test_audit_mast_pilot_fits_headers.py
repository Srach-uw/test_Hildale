from __future__ import annotations

import csv
import hashlib
import io
import sys
from pathlib import Path

import pytest
from astropy.io import fits


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_mast_pilot_fits_headers import (  # noqa: E402
    AuditInputError,
    MAX_RANGE_BYTES,
    base_record,
    fetch_product,
    mast_url,
    parse_content_range,
    read_and_validate_availability,
    read_and_validate_pilot,
    read_and_validate_summary,
)


class FakeResponse:
    def __init__(self, status, headers, body=b""):
        self.status_code = status
        self.headers = headers
        self.body = body
        self.closed = False
        self.iterated = False

    def iter_content(self, chunk_size):
        self.iterated = True
        for index in range(0, len(self.body), chunk_size):
            yield self.body[index:index + chunk_size]

    def close(self):
        self.closed = True


def header_bytes(kepid=123456789, obsmode="long cadence", with_hdu=True):
    primary = fits.PrimaryHDU()
    primary.header["KEPLERID"] = kepid
    primary.header["OBSMODE"] = obsmode
    primary.header["QUARTER"] = 4
    hdus = [primary]
    if with_hdu:
        table = fits.BinTableHDU.from_columns([fits.Column(name="TIME", format="D", array=[1.0])])
        table.header["EXTNAME"] = "LIGHTCURVE"
        table.header["TSTART"] = 1.0
        table.header["TSTOP"] = 2.0
        table.header["TIMEDEL"] = 0.02
        table.header["INT_TIME"] = 1625.35
        table.header["NUM_FRM"] = 270
        hdus.append(table)
    buffer = io.BytesIO()
    fits.HDUList(hdus).writeto(buffer)
    return buffer.getvalue()


def record():
    availability = {"koi_target": "K00001", "kepid": "123456789", "cadence": "long"}
    pilot = {"source_class": "source_original_archive", "candidate_role": "canonical_leverage", "priority_planet": "K00001.01"}
    result = base_record(availability, pilot, "kplr123456789-2009166043257_llc.fits")
    result["url"] = "https://example.test/file"
    return result


def requester(response):
    def call(*args, **kwargs):
        assert kwargs["headers"]["Range"] == "bytes=0-65535"
        assert kwargs["stream"] is True
        return response
    return call


def test_valid_partial_range_parses_headers():
    body = header_bytes()
    response = FakeResponse(206, {"Content-Range": f"bytes 0-{len(body)-1}/{len(body)}", "Content-Length": str(len(body))}, body)
    result = fetch_product(record(), 1, requester(response))
    assert result["status"] == "ok"
    assert result["identity_consistent"] is True
    assert result["cadence_consistent"] is True
    assert result["hdu1_ttype_columns"] == "TIME"
    assert response.closed is True


def test_status_200_is_closed_without_reading_body():
    response = FakeResponse(200, {}, b"must not be read")
    result = fetch_product(record(), 1, requester(response))
    assert result["status"] == "range_not_honored"
    assert response.iterated is False
    assert response.closed is True


@pytest.mark.parametrize("value", ["bytes 0-65536/70000", "bytes 1-10/100", "bytes 0-1/*", ""])
def test_malformed_or_oversized_content_range_is_rejected(value):
    with pytest.raises(ValueError):
        parse_content_range(value)


def test_missing_headers_are_explicitly_incomplete():
    body = header_bytes(with_hdu=False)
    response = FakeResponse(206, {"Content-Range": f"bytes 0-{len(body)-1}/{len(body)}"}, body)
    result = fetch_product(record(), 1, requester(response))
    assert result["status"] == "header_incomplete"
    assert "missing_HDU1_header" in result["error"]


def test_cadence_mismatch_is_explicitly_incomplete():
    body = header_bytes(obsmode="short cadence")
    response = FakeResponse(206, {"Content-Range": f"bytes 0-{len(body)-1}/{len(body)}"}, body)
    result = fetch_product(record(), 1, requester(response))
    assert result["status"] == "header_incomplete"
    assert "primary_OBSMODE_cadence_mismatch" in result["error"]


def test_wrong_kic_or_path_injection_filename_is_rejected():
    with pytest.raises(AuditInputError, match="KIC"):
        mast_url(123456789, "kplr987654321-2009166043257_llc.fits", "long")
    with pytest.raises(AuditInputError, match="Unsafe"):
        mast_url(123456789, "../kplr123456789-2009166043257_llc.fits", "long")


def write_csv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def strict_summary_rows():
    rows = []
    for index in range(2123):
        target = f"K{index % 1583:05d}"
        rows.append({
            "kepoi_name": f"{target}.{index // 1583 + 1:02d}", "koi_target": target,
            "kepid": str(100000000 + index % 1583), "qc_primary_exclude": "False",
            "qc_reasons": "", "impact_mode": "alderaan", "nested_weight_mode": "dynesty",
            "include_transit_prior": "False", "posterior_source": "alderaan_direct_importance",
        })
    for row in rows[-4:]:
        row["qc_primary_exclude"] = "True"
        row["qc_reasons"] = "importance_ess_below_threshold"
    return rows


def test_summary_hash_and_pilot_availability_mapping_reject_duplicates(tmp_path):
    summary = tmp_path / "summary.csv"
    fields = list(strict_summary_rows()[0])
    write_csv(summary, fields, strict_summary_rows())
    expected = hashlib.sha256(summary.read_bytes()).hexdigest()
    mapping = read_and_validate_summary(summary, expected)
    pilot_rows = [
        {"koi_target": f"K{i:05d}", "source_class": "source_original_archive", "candidate_role": "control", "priority_planet": f"K{i:05d}.01"}
        for i in range(24)
    ]
    pilot = tmp_path / "pilot.csv"
    write_csv(pilot, list(pilot_rows[0]), pilot_rows)
    validated_pilot = read_and_validate_pilot(pilot, mapping)
    availability_rows = []
    for target in validated_pilot:
        for cadence in ("long", "short"):
            availability_rows.append({"koi_target": target, "kepid": str(mapping[target]), "cadence": cadence, "query_status": "success", "filenames": ""})
    availability_rows[-1]["cadence"] = "long"
    availability = tmp_path / "availability.csv"
    write_csv(availability, list(availability_rows[0]), availability_rows)
    with pytest.raises(AuditInputError, match="24 pilot targets"):
        read_and_validate_availability(availability, validated_pilot, mapping)
    with pytest.raises(AuditInputError, match="SHA-256"):
        read_and_validate_summary(summary, "0" * 64)


def test_range_cap_is_constant():
    assert MAX_RANGE_BYTES == 65_536
