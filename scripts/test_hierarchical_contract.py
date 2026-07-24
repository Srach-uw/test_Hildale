from __future__ import annotations

import pandas as pd
import pytest

from hierarchical_rayleigh import validate_summary_contract


def valid_summary() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "kepoi_name": ["K00001.01"],
            "disk": ["thin"],
            "system": ["single"],
            "posterior_file": ["posterior.npz"],
            "qc_primary_exclude": [False],
            "qc_reasons": [""],
            "posterior_source": ["alderaan_direct_importance"],
            "impact_mode": ["alderaan"],
            "nested_weight_mode": ["dynesty"],
        }
    )


def test_valid_direct_summary_contract_passes() -> None:
    validate_summary_contract(valid_summary())


def test_mixed_sources_fail_canonical_contract() -> None:
    summary = pd.concat([valid_summary(), valid_summary().assign(posterior_source="legacy_zeta")])
    with pytest.raises(ValueError, match="cannot mix posterior constructions"):
        validate_summary_contract(summary)


def test_nonpaired_impact_fails_canonical_contract() -> None:
    with pytest.raises(ValueError, match="paired ALDERAAN"):
        validate_summary_contract(valid_summary().assign(impact_mode="geometric"))


def test_equal_nested_weights_fail_canonical_contract() -> None:
    with pytest.raises(ValueError, match="LN_WT"):
        validate_summary_contract(valid_summary().assign(nested_weight_mode="equal"))


def test_equal_nested_weights_require_explicit_diagnostic_override() -> None:
    validate_summary_contract(
        valid_summary().assign(nested_weight_mode="equal"),
        allow_non_dynesty_weights=True,
    )
