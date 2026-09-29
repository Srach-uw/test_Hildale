import ast
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from verify_diagnostic_fixture_reader import (
    ReaderAuditError,
    _reader_arrays,
    classify_reader_times,
    compare_reader_provenance,
)


def test_reader_audit_avoids_python39_incompatible_path_write_text_newline():
    source = (SCRIPTS / "verify_diagnostic_fixture_reader.py").read_text(encoding="utf-8")
    calls = [
        node
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "write_text"
        and any(keyword.arg == "newline" for keyword in node.keywords)
    ]
    assert not calls


def test_reader_audit_counts_masked_and_eligible_cadences():
    counts = classify_reader_times(
        np.array([1.0, 2.0, 3.0, 4.0]),
        np.array([True, False, False, True]),
        np.array([2.0, 3.0, 8.0]),
    )
    assert counts == {
        "prepared_total": 4,
        "eligible_total": 2,
        "eligible_retained": 2,
        "eligible_missing": 0,
        "known_transit_total": 2,
        "known_transit_retained": 0,
        "reader_total": 3,
        "reader_extra": 1,
    }


def test_reader_audit_rejects_misaligned_mask():
    with pytest.raises(ReaderAuditError, match="aligned"):
        classify_reader_times(np.array([1.0]), np.array([False, True]), np.array([1.0]))


class _Values:
    def __init__(self, values):
        self.value = np.asarray(values)


class _LightCurve:
    def __init__(self, *, time, flux, flux_err, cadenceno, quarter):
        self.time = _Values(time)
        self.flux = _Values(flux)
        self.flux_err = _Values(flux_err)
        self.cadenceno = _Values(cadenceno)
        self.quarter = quarter


def test_reader_array_collection_preserves_scalar_quarter_and_cadence_identity():
    payload = _reader_arrays([
        _LightCurve(time=[1.0, 2.0], flux=[1.0, 1.1], flux_err=[0.01, 0.01], cadenceno=[10, 11], quarter=3),
        _LightCurve(time=[3.0], flux=[1.0], flux_err=[0.02], cadenceno=[20], quarter=4),
    ])
    assert payload["cadence"].tolist() == [10, 11, 20]
    assert payload["quarter"].tolist() == [3, 3, 4]
    assert payload["flux_err"].tolist() == [0.01, 0.01, 0.02]


def _fixture_provenance():
    return {
        "time": np.array([1.0, 2.0, 3.0, 4.0]),
        "flux": np.array([100.0, 102.0, 200.0, 204.0]),
        "flux_err": np.array([2.0, 2.0, 4.0, 4.0]),
        "cadence": np.array([11, 12, 21, 22]),
        "quarter": np.array([1, 1, 2, 2]),
        "known_transit_mask": np.array([False, True, False, False]),
    }


def test_reader_provenance_accepts_per_quarter_stitch_normalization():
    reader = {
        "time": np.array([1.0, 3.0, 4.0, 8.0]),
        "flux": np.array([1.0, 1.0, 1.02, 9.0]),
        "flux_err": np.array([0.02, 0.02, 0.02, 0.3]),
        "cadence": np.array([11, 21, 22, 99]),
        "quarter": np.array([1, 2, 2, 9]),
    }
    comparison = compare_reader_provenance(_fixture_provenance(), reader)
    assert comparison["pass"] is True
    assert comparison["counts"]["reader_extra"] == 1
    assert comparison["quarter_stitch_scales"] == {"1": 0.01, "2": 0.005}
    assert all(item["pass"] for item in comparison["field_checks"].values() if item["applicable"])


@pytest.mark.parametrize(
    ("field", "replacement"),
    [
        ("cadence", np.array([111, 21, 22])),
        ("quarter", np.array([1, 9, 2])),
        ("flux", np.array([1.0, 1.0, 1.08])),
        ("flux_err", np.array([0.02, 0.03, 0.02])),
    ],
)
def test_reader_provenance_rejects_changed_fixture_fields(field, replacement):
    reader = {
        "time": np.array([1.0, 3.0, 4.0]),
        "flux": np.array([1.0, 1.0, 1.02]),
        "flux_err": np.array([0.02, 0.02, 0.02]),
        "cadence": np.array([11, 21, 22]),
        "quarter": np.array([1, 2, 2]),
    }
    reader[field] = replacement
    comparison = compare_reader_provenance(_fixture_provenance(), reader)
    assert comparison["pass"] is False
    assert comparison["field_checks"][field]["pass"] is False


def test_reader_provenance_requires_reader_metadata_when_fixture_has_it():
    reader = {
        "time": np.array([1.0, 3.0, 4.0]),
        "flux": np.array([1.0, 1.0, 1.02]),
        "cadence": np.array([11, 21, 22]),
        "quarter": np.array([1, 2, 2]),
    }
    comparison = compare_reader_provenance(_fixture_provenance(), reader)
    assert comparison["pass"] is False
    assert comparison["field_checks"]["flux_err"]["reason"] == "reader has no flux_err values"
