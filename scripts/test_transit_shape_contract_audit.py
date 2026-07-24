from pathlib import Path

import numpy as np
import pandas as pd

from transit_shape_contract_audit import audit


def test_shape_audit_flags_bad_duration_and_density(tmp_path: Path) -> None:
    factorial = pd.DataFrame(
        {
            "arm": ["original_lc"],
            "kepid": [1],
            "koi_period": [5.0],
            "fit_period_days": [5.0],
            "t14_hr_p50": [2.0],
            "rp_over_rs_p50": [0.1],
            "impact_p50": [0.2],
            "e50": [0.4],
        }
    )
    catalog = pd.DataFrame(
        {
            "kepid": [1],
            "koi_period": [5.0],
            "koi_duration": [1.0],
            "koi_impact": [0.2],
            "koi_ror": [0.1],
            "koi_dor": [10.0],
            "koi_snr": [20.0],
            "koi_count": [1],
        }
    )
    berger = pd.DataFrame(
        {"kepid": [1], "rho_log": [0.0], "rho_log_upper": [-1.0], "rho_log_lower": [-1.0]}
    )
    fp = tmp_path / "factorial.csv"
    cp = tmp_path / "catalog.csv"
    bp = tmp_path / "berger.dat.gz"
    factorial.to_csv(fp, index=False)
    catalog.to_csv(cp, index=False)
    import gzip

    with gzip.open(bp, "wt") as handle:
        line = [" "] * 200
        text = "00000001"
        line[0:8] = text
        line[116:122] = " 0.000"
        line[123:129] = "-1.000"
        line[130:136] = "-1.000"
        handle.write("".join(line) + "\n")
    result = audit(fp, cp, bp)
    assert bool(result.loc[0, "duration_flag"])
    assert bool(result.loc[0, "contract_flag"])
    assert np.isfinite(result.loc[0, "rho_circular_from_fit"])
