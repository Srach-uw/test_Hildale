from pathlib import Path

from alderaan_failure_taxonomy import build_taxonomy, classify, infer_stage


def test_celerite_failure_is_not_mislabeled_as_transit_fit():
    text = "celerite2.backprop.LinAlgError: failed to factorize or solve matrix"
    assert classify(text) == "celerite_linalg"
    assert infer_stage(Path("K00001.stderr.log"), "failed_noise_exit_1") == "noise"


def test_missing_lightcurve_is_distinct():
    assert classify("No long-cadence FITS downloaded") == "missing_lightcurve"


def test_successful_result_overrides_recovered_celerite_warning(tmp_path):
    logs = tmp_path / "logs"
    status = tmp_path / "status"
    results = tmp_path / "Results" / "run" / "K00064"
    logs.mkdir()
    status.mkdir()
    results.mkdir(parents=True)
    (logs / "K00064.stderr.log").write_text(
        "celerite2.backprop.LinAlgError\nfallback succeeded\n", encoding="utf-8"
    )
    (status / "K00064.status").write_text("complete\n", encoding="utf-8")
    (results / "K00064-results.fits").write_bytes(b"FITS")

    frame = build_taxonomy(logs.rglob("*"), status, tmp_path / "Results")
    row = frame.iloc[0]
    assert row["outcome"] == "success"
    assert row["failure_family"] == "none"


def test_failure_uses_terminal_traceback_and_stage_status(tmp_path):
    logs = tmp_path / "logs"
    status = tmp_path / "status"
    logs.mkdir()
    status.mkdir()
    (logs / "K00001.stderr.log").write_text(
        "earlier warning: LD_U1 fallback applied\n"
        "Traceback (most recent call last):\n"
        "celerite2.backprop.LinAlgError: failed to factorize\n",
        encoding="utf-8",
    )
    (status / "K00001.status").write_text("failed_detrend_exit_1\n", encoding="utf-8")

    frame = build_taxonomy(logs.rglob("*"), status, tmp_path / "Results")
    row = frame.iloc[0]
    assert row["outcome"] == "failed"
    assert row["stage"] == "detrend"
    assert row["failure_family"] == "celerite_linalg"
