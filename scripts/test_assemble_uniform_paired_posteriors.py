from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from assemble_uniform_paired_posteriors import assemble, validate_summary


def row(name: str, path: Path, e50: float) -> dict[str, object]:
    return {
        "kepoi_name": name,
        "kepid": 1,
        "koi_target": name.split(".")[0],
        "koi_period": 10.0,
        "posterior_file": str(path),
        "posterior_source": "alderaan_direct_importance",
        "impact_mode": "alderaan",
        "formalism": "direct",
        "include_transit_prior": False,
        "density_source": "berger2020_table2",
        "density_sampling_mode": "draw_split_normal",
        "density_error_mode": "symmetric-average",
        "posterior_sampling_mode": "weighted_importance",
        "period_sampling_mode": "paired_alderaan",
        "nested_weight_mode": "dynesty",
        "e_max": 0.95,
        "qc_primary_exclude": False,
        "qc_reasons": "",
        "e50": e50,
    }


def test_archive_wins_overlap_and_labels_come_from_sample(tmp_path: Path) -> None:
    posterior = tmp_path / "posterior.npz"
    posterior.touch()
    archive = pd.DataFrame([row("K00001.01", posterior, 0.1)])
    new = pd.DataFrame([row("K00001.01", posterior, 0.8), row("K00002.01", posterior, 0.2)])
    sample = pd.DataFrame(
        [
            {"kepoi_name": "K00001.01", "kepid": 1, "koi_target": "K00001", "koi_period": 10.0, "disk": "thin", "system": "single"},
            {"kepoi_name": "K00002.01", "kepid": 2, "koi_target": "K00002", "koi_period": 20.0, "disk": "thick", "system": "multi"},
        ]
    )
    merged, overlap = assemble(archive, new, sample)
    assert len(overlap) == 1
    assert merged.set_index("kepoi_name").loc["K00001.01", "e50"] == pytest.approx(0.1)
    assert set(merged["posterior_source"]) == {"uniform_paired_direct_importance"}
    assert merged.set_index("kepoi_name").loc["K00002.01", "system"] == "multi"


def test_geometric_summary_is_rejected(tmp_path: Path) -> None:
    posterior = tmp_path / "posterior.npz"
    posterior.touch()
    summary = pd.DataFrame([row("K00001.01", posterior, 0.1)])
    summary["impact_mode"] = "geometric"
    with pytest.raises(ValueError, match="not uniformly paired-impact"):
        validate_summary(summary, "archive")


def test_cross_source_density_mixture_is_rejected(tmp_path: Path) -> None:
    posterior = tmp_path / "posterior.npz"
    posterior.touch()
    archive = pd.DataFrame([row("K00001.01", posterior, 0.1)])
    new = pd.DataFrame(
        [row("K00002.01", posterior, 0.2) | {"density_source": "berger2018_kg"}]
    )
    sample = pd.DataFrame(
        [
            {"kepoi_name": "K00001.01", "kepid": 1, "koi_target": "K00001", "koi_period": 10.0, "disk": "thin", "system": "single"},
            {"kepoi_name": "K00002.01", "kepid": 2, "koi_target": "K00002", "koi_period": 20.0, "disk": "thin", "system": "single"},
        ]
    )
    with pytest.raises(ValueError, match="uniform density_source"):
        assemble(archive, new, sample)
