from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATION_HELPERS = REPO_ROOT / "cloud" / "ld_validation"
sys.path.insert(0, str(VALIDATION_HELPERS))

from patch_alderaan_repro import main


def test_repro_patch_seeds_all_stages_and_is_idempotent(
    tmp_path: Path, monkeypatch
) -> None:
    repo = tmp_path / "alderaan"
    target_bin = repo / "bin"
    target_bin.mkdir(parents=True)
    (target_bin / "detrend_and_estimate_ttvs.py").write_text(
        "import os\nimport numpy as np\n"
        "from alderaan.validate import remove_known_transits, inject_synthetic_transits\n",
        encoding="utf-8",
    )
    (target_bin / "analyze_autocorrelated_noise.py").write_text(
        "import os\nimport numpy as np\n",
        encoding="utf-8",
    )
    (target_bin / "fit_transit_shape_simultaneous_nested.py").write_text(
        "import os\nimport numpy as np\n    USE_MULTIPRO = False\n"
        "                sample=\"rwalk\",\n                pool=pool,\n"
        "            sample=\"rwalk\",\n            logl_args=logl_args,\n",
        encoding="utf-8",
    )
    names = ["detrend_and_estimate_ttvs.py", "fit_transit_shape_simultaneous_nested.py"]

    monkeypatch.setattr("sys.argv", ["patch_alderaan_repro.py", str(repo)])
    main()
    first = {name: (target_bin / name).read_text(encoding="utf-8") for name in names}
    main()
    second = {name: (target_bin / name).read_text(encoding="utf-8") for name in names}

    assert first == second
    for text in second.values():
        assert text.count("ALDERAAN_NUMPY_SEED") == 1
        assert 'os.environ.get("ALDERAAN_SEED")' in text
    assert "rstate=RSTATE" in second["fit_transit_shape_simultaneous_nested.py"]
