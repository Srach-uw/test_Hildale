from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from prepare_six_target_timing_readiness import (  # noqa: E402
    ENV_VAR,
    PILOT_TARGETS,
    ReadinessInputError,
    build_readiness,
    main,
    parse_environment_mapping,
    parse_result_specifications,
    write_readiness,
)


def timing_hdu(
    index: int,
    *,
    missing_column: str | None = None,
    ttime: np.ndarray | None = None,
    nonfinite_column: str | None = None,
) -> fits.BinTableHDU:
    values = {
        "INDEX": np.array([0, 1, 3], dtype=np.int64),
        "TTIME": np.array([100.0, 110.0, 130.0]) if ttime is None else ttime,
        "MODEL": np.array([100.1, 110.1, 130.1]),
        "OUT_PROB": np.array([0.1, 0.2, 0.3]),
        "OUT_FLAG": np.array([0, 0, 1], dtype=np.int16),
    }
    if nonfinite_column is not None:
        values[nonfinite_column] = np.asarray(values[nonfinite_column], dtype=float)
        values[nonfinite_column][1] = np.nan
    columns = []
    for name, array in values.items():
        if name == missing_column:
            continue
        fmt = "K" if name in {"INDEX", "OUT_FLAG"} and array.dtype.kind in "iu" else "D"
        columns.append(fits.Column(name=name, format=fmt, array=array))
    return fits.BinTableHDU.from_columns(columns, name=f"TTIMES_{index:02d}")


def write_result(
    path: Path,
    target: str,
    *,
    npl: int = 1,
    extensions: list[fits.BinTableHDU] | None = None,
) -> Path:
    primary = fits.PrimaryHDU()
    primary.header["TARGET"] = target
    primary.header["NPL"] = npl
    if extensions is None:
        extensions = [timing_hdu(index) for index in range(npl)]
    fits.HDUList([primary, *extensions]).writeto(path)
    return path


def valid_six_target_mapping(tmp_path: Path) -> dict[str, Path]:
    return {
        target: write_result(tmp_path / f"{target}-results.fits", target)
        for target in PILOT_TARGETS
    }


def reason_codes(record: dict) -> set[str]:
    return {reason["code"] for reason in record["reasons"]}


def test_valid_six_target_fits_are_ready_and_hashed(tmp_path):
    paths = valid_six_target_mapping(tmp_path)
    payload = build_readiness(paths)

    assert payload["overall_ready"] is True
    assert payload["ready_target_count"] == 6
    assert payload["timing_inference_performed"] is False
    for record in payload["targets"]:
        expected_hash = hashlib.sha256(paths[record["target"]].read_bytes()).hexdigest()
        assert record["input_sha256"] == expected_hash
        assert record["ready"] is True
        assert record["extensions"][0]["index_gap_count"] == 1


def test_missing_extension_is_reported_without_inference(tmp_path):
    paths = valid_six_target_mapping(tmp_path)
    paths["K00367"] = write_result(
        tmp_path / "missing-extension.fits", "K00367", npl=2, extensions=[timing_hdu(0)]
    )

    payload = build_readiness(paths)
    record = payload["targets"][0]

    assert record["ready"] is False
    assert "missing_ttimes_extension" in reason_codes(record)
    assert record["extensions"][1] == {
        "planet_index": 1,
        "extension": "TTIMES_01",
        "present": False,
        "ready": False,
        "rows": 0,
    }
    assert payload["timing_inference_performed"] is False


@pytest.mark.parametrize(
    ("extension", "expected_code"),
    [
        (timing_hdu(0, missing_column="MODEL"), "missing_ttimes_columns"),
        (timing_hdu(0, nonfinite_column="TTIME"), "non_finite_column"),
        (
            timing_hdu(0, ttime=np.array([100.0, 130.0, 120.0])),
            "nonmonotonic_ttime",
        ),
    ],
)
def test_malformed_extension_is_not_ready(tmp_path, extension, expected_code):
    paths = valid_six_target_mapping(tmp_path)
    paths["K00367"] = write_result(
        tmp_path / f"{expected_code}.fits", "K00367", extensions=[extension]
    )

    record = build_readiness(paths)["targets"][0]

    assert record["ready"] is False
    assert expected_code in reason_codes(record)


def test_npl_requires_every_index_and_rejects_unexpected_extension(tmp_path):
    paths = valid_six_target_mapping(tmp_path)
    paths["K00367"] = write_result(
        tmp_path / "coverage.fits",
        "K00367",
        npl=2,
        extensions=[timing_hdu(0), timing_hdu(2)],
    )

    record = build_readiness(paths)["targets"][0]

    assert {"missing_ttimes_extension", "unexpected_ttimes_extension"} <= reason_codes(record)


def test_cli_writes_only_explicit_output_directory(tmp_path):
    paths = valid_six_target_mapping(tmp_path)
    output = tmp_path / "readiness"
    arguments = ["--output-dir", str(output)]
    for target, path in paths.items():
        arguments.extend(["--result", f"{target}={path}"])

    assert main(arguments) == 0
    assert sorted(path.name for path in output.iterdir()) == [
        "timing_readiness.csv",
        "timing_readiness.json",
    ]
    payload = json.loads((output / "timing_readiness.json").read_text(encoding="utf-8"))
    assert payload["overall_ready"] is True


def test_environment_json_and_exact_target_set(monkeypatch, tmp_path):
    paths = valid_six_target_mapping(tmp_path)
    raw = json.dumps({target: str(path) for target, path in paths.items()})
    parsed = parse_environment_mapping(raw)
    assert list(parsed) == list(PILOT_TARGETS)

    monkeypatch.setenv(ENV_VAR, raw)
    output = tmp_path / "environment-output"
    assert main(["--output-dir", str(output)]) == 0

    incomplete = [f"{target}={path}" for target, path in list(paths.items())[:-1]]
    with pytest.raises(ReadinessInputError, match="missing K01852"):
        parse_result_specifications(incomplete)


def test_nonempty_output_directory_is_refused(tmp_path):
    paths = valid_six_target_mapping(tmp_path)
    output = tmp_path / "occupied"
    output.mkdir()
    (output / "keep.txt").write_text("do not replace", encoding="utf-8")

    with pytest.raises(ReadinessInputError, match="absent or empty"):
        write_readiness(build_readiness(paths), output)
    assert (output / "keep.txt").read_text(encoding="utf-8") == "do not replace"
