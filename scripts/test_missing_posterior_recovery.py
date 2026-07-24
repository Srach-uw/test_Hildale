from pathlib import Path

import pandas as pd

from audit_missing_posterior_recovery import main


def test_exact_recovery_inputs_have_expected_counts(tmp_path, monkeypatch, capsys):
    root = Path(__file__).resolve().parents[1]
    output = tmp_path / "audit.csv"
    monkeypatch.setattr(
        "sys.argv",
        [
            "audit_missing_posterior_recovery.py",
            "--root",
            str(root),
            "--manifest",
            str(
                root
                / "metadata"
                / "recovery_preflight_20260724"
                / "eccentricity_posterior_manifest_equal_nested_published_inventory_pre_visual_qc.csv"
            ),
            "--batch",
            str(root / "cloud" / "recovery_142"),
            "--output",
            str(output),
        ],
    )
    main()
    text = capsys.readouterr().out
    assert "missing_population_planets,174" in text
    assert "missing_systems,142" in text
    assert "new_systems,116" in text
    assert "retry_systems,26" in text
    frame = pd.read_csv(output)
    assert len(frame) == 142
    assert frame["missing_population_planets"].sum() == 174
    assert frame["new_system"].sum() == 116
    assert frame["prior_retry"].sum() == 26
