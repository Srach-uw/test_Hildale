"""Apply the minimal, idempotent ALDERAAN patches required by this project.

The patches do not change the transit model or priors. They make the unused
validation import optional, seed NumPy plus dynesty through the ALDERAAN_SEED
environment variable, and repair two upstream detrending control-flow errors:
quarters with no modelled transits are excluded rather than indexed with an
undefined outlier mask, and the documented SHOTerm fallback calls the local
function rather than an undefined module name.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def replace_once(path: Path, old: str, new: str, marker: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker in text:
        return
    if old not in text:
        raise RuntimeError(f"Expected patch context not found in {path}: {old!r}")
    path.write_text(text.replace(old, new, 1), encoding="utf-8")


def seed_numpy(path: Path) -> None:
    replace_once(
        path,
        "import numpy as np\n",
        "import numpy as np\n\n"
        "# ALDERAAN_NUMPY_SEED: deterministic target-level pipeline seed\n"
        "_alderaan_seed = os.environ.get(\"ALDERAAN_SEED\")\n"
        "if _alderaan_seed is not None:\n"
        "    np.random.seed(int(_alderaan_seed))\n",
        "ALDERAAN_NUMPY_SEED",
    )


def patch_detrend_control_flow(detrend_script: Path, detrend_module: Path) -> None:
    replace_once(
        detrend_script,
        "    for j, q in enumerate(quarters):\n"
        "        if all_dtype[q] == \"long\":\n"
        "            use = lc.quarter == q\n"
        "            good_cadno_lc.append(lc.cadno[use][~bad[j]])\n\n"
        "        if all_dtype[q] == \"short\":\n"
        "            use = sc.quarter == q\n"
        "            good_cadno_sc.append(sc.cadno[use][~bad[j]])\n",
        "    for j, q in enumerate(quarters):\n"
        "        # ALDERAAN_NO_TRANSIT_QUARTER_GUARD: a refined ephemeris can\n"
        "        # leave a covered quarter without a transit model. Exclude that\n"
        "        # quarter rather than applying a nonexistent outlier mask.\n"
        "        if bad[j] is None:\n"
        "            continue\n\n"
        "        if all_dtype[q] == \"long\":\n"
        "            use = lc.quarter == q\n"
        "            good_cadno_lc.append(lc.cadno[use][~bad[j]])\n\n"
        "        if all_dtype[q] == \"short\":\n"
        "            use = sc.quarter == q\n"
        "            good_cadno_sc.append(sc.cadno[use][~bad[j]])\n",
        "ALDERAAN_NO_TRANSIT_QUARTER_GUARD",
    )
    replace_once(
        detrend_module,
        "            litecurve = detrend.flatten_with_gp(\n",
        "            # ALDERAAN_SHOTERM_FALLBACK: this module already owns the helper.\n"
        "            litecurve = flatten_with_gp(\n",
        "ALDERAAN_SHOTERM_FALLBACK",
    )


def _apply_reproducibility_patches(
    detrend: Path,
    analyze: Path,
    fit: Path,
) -> None:
    replace_once(
        detrend,
        "from alderaan.validate import remove_known_transits, inject_synthetic_transits",
        "try:\n    from alderaan.validate import remove_known_transits, inject_synthetic_transits\n"
        "except ImportError:\n    remove_known_transits = inject_synthetic_transits = None",
        "remove_known_transits = inject_synthetic_transits = None",
    )
    for script in (detrend, analyze, fit):
        seed_numpy(script)

    replace_once(
        fit,
        "    USE_MULTIPRO = False\n",
        "    seed_text = os.environ.get(\"ALDERAAN_SEED\")\n"
        "    RSTATE = np.random.default_rng(int(seed_text)) if seed_text is not None else None\n\n"
        "    USE_MULTIPRO = False\n",
        "RSTATE = np.random.default_rng",
    )
    text = fit.read_text(encoding="utf-8")
    constructor = "                sample=\"rwalk\",\n                pool=pool,\n"
    if "pool=pool,\n                rstate=RSTATE," not in text:
        if constructor not in text:
            raise RuntimeError("Multiprocessing dynesty constructor context changed")
        text = text.replace(
            constructor,
            "                sample=\"rwalk\",\n                pool=pool,\n                rstate=RSTATE,\n",
            1,
        )
    constructor = "            sample=\"rwalk\",\n            logl_args=logl_args,\n"
    if "sample=\"rwalk\",\n            rstate=RSTATE," not in text:
        if constructor not in text:
            raise RuntimeError("Serial dynesty constructor context changed")
        text = text.replace(
            constructor,
            "            sample=\"rwalk\",\n            rstate=RSTATE,\n            logl_args=logl_args,\n",
            1,
        )
    fit.write_text(text, encoding="utf-8")


def apply_patches(repo: Path) -> None:
    detrend = repo / "bin" / "detrend_and_estimate_ttvs.py"
    detrend_module = repo / "alderaan" / "detrend.py"
    analyze = repo / "bin" / "analyze_autocorrelated_noise.py"
    fit = repo / "bin" / "fit_transit_shape_simultaneous_nested.py"
    replace_once(
        detrend,
        "from alderaan.validate import remove_known_transits, inject_synthetic_transits",
        "try:\n    from alderaan.validate import remove_known_transits, inject_synthetic_transits\n"
        "except ImportError:\n    remove_known_transits = inject_synthetic_transits = None",
        "remove_known_transits = inject_synthetic_transits = None",
    )
    patch_detrend_control_flow(detrend, detrend_module)
    _apply_reproducibility_patches(detrend, analyze, fit)
    print(f"Patched ALDERAAN reproducibly at {repo}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("repo")
    args = parser.parse_args()
    apply_patches(Path(args.repo).resolve())


if __name__ == "__main__":
    main()
