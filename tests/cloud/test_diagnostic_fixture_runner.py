import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import textwrap

import pytest


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "cloud" / "diagnostic_injection" / "run_one_fixture.sh"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _bash() -> str:
    candidates = [shutil.which("bash"), r"C:\Program Files\Git\bin\bash.exe"]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return candidate
    pytest.skip("bash is required for the diagnostic runner test")


def _shell_path(path: Path) -> str:
    value = path.resolve().as_posix()
    if os.name == "nt" and len(value) > 2 and value[1] == ":":
        return "/" + value[0].lower() + value[2:]
    return value


def _write_stage(path: Path, name: str, creates_result: bool = False) -> None:
    code = textwrap.dedent(
        f"""\
        import os
        import numpy as np
        from pathlib import Path
        import sys

        arguments = sys.argv[1:]
        project = Path(arguments[arguments.index("--project_dir") + 1])
        timing = project / "Catalogs" / "holczer_2016_kepler_ttvs.txt"
        rows = [line for line in timing.read_text().splitlines() if line and not line.startswith("#")]
        assert all(int(float(line.split()[0])) != 367 for line in rows)
        assert any(int(float(line.split()[0])) == 999 for line in rows)
        assert os.environ["ALDERAAN_SEED"] == "1729"
        with Path(os.environ["STAGE_LOG"]).open("a", encoding="utf-8") as handle:
            handle.write("{name}:" + str(np.random.get_state()[1][0]) + "\\n")
        """
    )
    if creates_result:
        code += textwrap.dedent(
            """\
            target = arguments[arguments.index("--target") + 1]
            run_id = arguments[arguments.index("--run_id") + 1]
            result = project / "Results" / run_id / target / f"{target}-results.fits"
            result.parent.mkdir(parents=True, exist_ok=True)
            result.write_bytes(b"result")
            """
        )
    path.write_text(code, encoding="utf-8")


def _case(tmp_path: Path, *, audit_pass: bool = True, fixture_hash: str | None = None,
          injection_hash: str | None = None, include_audit: bool = True) -> dict[str, Path | dict[str, str]]:
    fixture = tmp_path / "fixture"
    (fixture / "Data").mkdir(parents=True)
    (fixture / "Data" / "kplr000000123-2000000000000_llc.fits").write_bytes(b"fixture")
    (fixture / "diagnostic_catalog.csv").write_text("index,koi_id,kic_id,npl\n0,K00367,123,1\n", encoding="utf-8")
    manifest = fixture / "diagnostic_fixture_manifest.json"
    expected_injection = "a" * 64
    manifest.write_text(json.dumps({"injection": {"sha256": expected_injection}}), encoding="utf-8")

    repo = tmp_path / "alderaan"
    (repo / "bin").mkdir(parents=True)
    (repo / "Catalogs").mkdir()
    holczer = repo / "Catalogs" / "holczer_2016_kepler_ttvs.txt"
    holczer.write_text(
        "# Holczer fixture\n367.01 1 10.0 0.0 1.0\n367.02 2 20.0 0.0 1.0\n999.01 1 30.0 0.0 1.0\n",
        encoding="utf-8",
    )
    _write_stage(repo / "bin" / "detrend_and_estimate_ttvs.py", "detrend")
    _write_stage(repo / "bin" / "analyze_autocorrelated_noise.py", "noise")
    _write_stage(repo / "bin" / "fit_transit_shape_simultaneous_nested.py", "fit", creates_result=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.email", "test@example.invalid"], check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.name", "Test"], check=True)
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-qm", "fixture"], check=True)

    audit = tmp_path / "reader_audit.json"
    if include_audit:
        audit.write_text(
            json.dumps(
                {
                    "pass": audit_pass,
                    "fixture_manifest_sha256": fixture_hash or _sha256(manifest),
                    "injection_sha256": injection_hash or expected_injection,
                }
            ),
            encoding="utf-8",
        )
    project = tmp_path / "project"
    stage_log = tmp_path / "stages.log"
    env = {
        **os.environ,
        "FIXTURE_DIR": _shell_path(fixture),
        "PROJECT_DIR": _shell_path(project),
        "ALDERAAN_REPO": _shell_path(repo),
        "STAGE_LOG": _shell_path(stage_log),
        "HILDALE_DIAGNOSTIC_SEED": "1729",
    }
    if include_audit:
        env["READER_AUDIT"] = _shell_path(audit)
    return {"fixture": fixture, "manifest": manifest, "repo": repo, "holczer": holczer,
            "audit": audit, "project": project, "stage_log": stage_log, "env": env}


def _run(case: dict[str, Path | dict[str, str]]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [_bash(), _shell_path(RUNNER), "K00367", "123"],
        env=case["env"],
        text=True,
        capture_output=True,
        check=False,
    )


def test_fixture_runner_cannot_download_or_reuse_project():
    text = RUNNER.read_text(encoding="utf-8")
    assert "bin/get_kepler_data.py" not in text
    assert "wget" not in text
    assert "curl" not in text
    assert "Refusing to reuse an existing diagnostic project" in text
    assert "diagnostic_fixture_manifest.json" in text
    assert "diagnostic_catalog.csv" in text
    assert "provenance_sha256.txt" in text
    assert "READER_AUDIT" in text


def test_fixture_runner_avoids_python39_incompatible_path_write_text_newline():
    text = RUNNER.read_text(encoding="utf-8")
    assert "write_text(" not in text
    assert 'destination.open("w", encoding="utf-8", newline="\\n")' in text


def test_fixture_runner_uses_public_three_stage_sequence():
    text = RUNNER.read_text(encoding="utf-8")
    for stage in (
        "detrend_and_estimate_ttvs.py",
        "analyze_autocorrelated_noise.py",
        "fit_transit_shape_simultaneous_nested.py",
    ):
        assert stage in text
    assert "--use_sc True" not in text


def test_fixture_runner_refuses_missing_reader_audit_before_project_or_stages(tmp_path):
    case = _case(tmp_path, include_audit=False)
    result = _run(case)
    assert result.returncode == 2
    assert "READER_AUDIT JSON path is required" in result.stderr
    assert not case["project"].exists()
    assert not case["stage_log"].exists()


@pytest.mark.parametrize("seed", [None, "", "-1", "+1", "1.5", "not-a-seed"])
def test_fixture_runner_refuses_missing_or_invalid_seed_before_project_or_stages(tmp_path, seed):
    case = _case(tmp_path)
    if seed is None:
        case["env"].pop("HILDALE_DIAGNOSTIC_SEED")
    else:
        case["env"]["HILDALE_DIAGNOSTIC_SEED"] = seed

    result = _run(case)

    assert result.returncode == 2
    assert "HILDALE_DIAGNOSTIC_SEED" in result.stderr
    assert not case["project"].exists()
    assert not case["stage_log"].exists()


def test_fixture_runner_refuses_empty_reader_audit_before_project_or_stages(tmp_path):
    case = _case(tmp_path)
    case["audit"].write_bytes(b"")
    result = _run(case)
    assert result.returncode == 2
    assert "Reader audit is missing or empty" in result.stderr
    assert not case["project"].exists()
    assert not case["stage_log"].exists()


def test_fixture_runner_refuses_failing_reader_audit_before_project_or_stages(tmp_path):
    case = _case(tmp_path, audit_pass=False)
    result = _run(case)
    assert result.returncode != 0
    assert "pass must be true" in result.stderr
    assert not case["project"].exists()
    assert not case["stage_log"].exists()


@pytest.mark.parametrize(
    ("fixture_hash", "injection_hash", "message"),
    [
        ("0" * 64, None, "fixture manifest SHA-256 mismatch"),
        (None, "0" * 64, "injection SHA-256 mismatch"),
    ],
)
def test_fixture_runner_refuses_mismatched_reader_audit_before_stages(
    tmp_path, fixture_hash, injection_hash, message
):
    case = _case(tmp_path, fixture_hash=fixture_hash, injection_hash=injection_hash)
    result = _run(case)
    assert result.returncode != 0
    assert message in result.stderr
    assert not case["project"].exists()
    assert not case["stage_log"].exists()


def test_fixture_runner_filters_target_timing_rows_and_runs_three_stages(tmp_path):
    case = _case(tmp_path)
    source_before = case["holczer"].read_bytes()
    result = _run(case)
    assert result.returncode == 0, result.stderr
    assert case["stage_log"].read_text(encoding="utf-8").splitlines() == [
        "detrend:1729",
        "noise:1729",
        "fit:1729",
    ]
    assert case["holczer"].read_bytes() == source_before

    local_timing = case["project"] / "Catalogs" / "holczer_2016_kepler_ttvs.txt"
    text = local_timing.read_text(encoding="utf-8")
    assert "target=K00367 removed_rows=2" in text
    data_rows = [line for line in text.splitlines() if line and not line.startswith("#")]
    assert all(int(float(line.split()[0])) != 367 for line in data_rows)
    assert any(int(float(line.split()[0])) == 999 for line in data_rows)
    provenance = (case["project"] / "run_provenance.tsv").read_text(encoding="utf-8")
    assert "diagnostic_seed\t1729" in provenance
    assert "alderaan_seed\t1729" in provenance
    assert "holczer_removed_rows\t2" in provenance
    hashes = (case["project"] / "provenance_sha256.txt").read_text(encoding="utf-8")
    assert "reader_audit.json" in hashes
    assert "holczer_2016_kepler_ttvs.txt" in hashes
