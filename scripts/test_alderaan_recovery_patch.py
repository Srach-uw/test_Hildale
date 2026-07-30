import importlib.util
from pathlib import Path


def load_patch_module():
    path = Path(__file__).parents[1] / "cloud" / "recovery_142" / "patch_alderaan_repro.py"
    spec = importlib.util.spec_from_file_location("alderaan_recovery_patch", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_control_flow_repairs_are_idempotent(tmp_path: Path) -> None:
    repo = tmp_path / "alderaan"
    script = repo / "bin" / "detrend_and_estimate_ttvs.py"
    module = repo / "alderaan" / "detrend.py"
    script.parent.mkdir(parents=True)
    module.parent.mkdir(parents=True)

    script.write_text(
        "import os\nimport numpy as np\n"
        "from alderaan.validate import remove_known_transits, inject_synthetic_transits\n"
        "    for j, q in enumerate(quarters):\n"
        "        if all_dtype[q] == \"long\":\n"
        "            use = lc.quarter == q\n"
        "            good_cadno_lc.append(lc.cadno[use][~bad[j]])\n\n"
        "        if all_dtype[q] == \"short\":\n"
        "            use = sc.quarter == q\n"
        "            good_cadno_sc.append(sc.cadno[use][~bad[j]])\n",
        encoding="utf-8",
    )
    module.write_text(
        "def flatten_with_gp(*args, **kwargs):\n    return None\n\n"
        "            litecurve = detrend.flatten_with_gp(\n",
        encoding="utf-8",
    )

    patch = load_patch_module()
    patch.patch_detrend_control_flow(script, module)
    patch.patch_detrend_control_flow(script, module)

    assert "ALDERAAN_NO_TRANSIT_QUARTER_GUARD" in script.read_text(encoding="utf-8")
    assert "if bad[j] is None:\n            continue" in script.read_text(encoding="utf-8")
    assert "ALDERAAN_SHOTERM_FALLBACK" in module.read_text(encoding="utf-8")
    assert "detrend.flatten_with_gp" not in module.read_text(encoding="utf-8")
