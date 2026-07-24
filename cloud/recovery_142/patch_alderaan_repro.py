"""Apply the minimal, idempotent ALDERAAN patches required by this project.

The patches do not change the transit model or priors. They make the unused
validation import optional and seed NumPy plus dynesty through the
ALDERAAN_SEED environment variable so the detrending, synthetic-noise, and
transit-fit stages are reproducible.
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("repo")
    args = parser.parse_args()
    repo = Path(args.repo).resolve()

    detrend = repo / "bin" / "detrend_and_estimate_ttvs.py"
    replace_once(
        detrend,
        "from alderaan.validate import remove_known_transits, inject_synthetic_transits",
        "try:\n    from alderaan.validate import remove_known_transits, inject_synthetic_transits\n"
        "except ImportError:\n    remove_known_transits = inject_synthetic_transits = None",
        "remove_known_transits = inject_synthetic_transits = None",
    )
    analyze = repo / "bin" / "analyze_autocorrelated_noise.py"
    fit = repo / "bin" / "fit_transit_shape_simultaneous_nested.py"
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
    print(f"Patched ALDERAAN reproducibly at {repo}")


if __name__ == "__main__":
    main()
