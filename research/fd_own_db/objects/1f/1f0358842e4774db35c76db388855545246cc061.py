import hashlib

import pandas as pd
import pytest

from fd.functional_source_placebo import (
    validate_predictions, safe_grouped_scores,
)


def rows(n_islands=44, states=None):
    out = []
    for i in range(n_islands):
        for j in range(3):
            name = f"Species_{i:03d}_{j:02d}"
            state = 1 if states is None else states(i, j)
            out.append({
                "scope": "All", "trait": "self_compatibility", "region": "tropical",
                "island_id": f"secret_island_{i}", "taxon_name": name,
                "spatial_block": f"private_block_{i//3}",
                "fold": hashlib.sha256(name.encode("utf-8")).digest()[0] % 5,
                "state": state, "p_actual": .75,
                "p_placebo_1": .50, "p_placebo_2": .45,
                "source_pool_ratio_1": .9, "source_pool_ratio_2": 1.1,
                "source_coverage_diff_1": .01, "source_coverage_diff_2": .015,
                "source_jaccard_1": .35, "source_jaccard_2": .60,
            })
    return pd.DataFrame(out)


def test_source_advantage_positive_and_private():
    s = safe_grouped_scores(rows(), bootstrap_draws=399)
    assert len(s) == 1
    assert s.iloc[0].brier_advantage_equal_block > 0
    assert s.iloc[0].brier_ci025 > 0
    assert s.iloc[0].n_islands == 44
    assert s.iloc[0].n_spatial_blocks == 15
    assert not {"taxon_name", "island_id", "spatial_block", "source_id"} & set(s)


def test_empty_when_geographical_independence_insufficient():
    z = rows()
    z["spatial_block"] = "one_connected_block"
    out = safe_grouped_scores(z, bootstrap_draws=399)
    assert out.empty


def test_no_fake_name_fold_allowed():
    z = rows()
    z.loc[0, "fold"] = (z.loc[0, "fold"] + 1) % 5
    with pytest.raises(ValueError, match="leakage"):
        validate_predictions(z)


def test_missing_traits_and_fallback_zero_are_not_imputed():
    z = rows()
    z.loc[0, "state"] = None
    with pytest.raises(ValueError, match="nonfinite"):
        validate_predictions(z)
    z = rows()
    z.loc[0, "p_actual"] = 1.5
    with pytest.raises(ValueError, match="probability"):
        validate_predictions(z)


def test_control_evidence_caliper_is_mandatory():
    z = rows()
    z.loc[0, "source_coverage_diff_1"] = .040
    with pytest.raises(ValueError, match="donor coverage"):
        validate_predictions(z)
    z = rows()
    z.loc[0, "source_jaccard_2"] = .90
    with pytest.raises(ValueError, match="too similar"):
        validate_predictions(z)


def test_duplicate_island_species_fails():
    z = pd.concat([rows(), rows().iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="Duplicate"):
        validate_predictions(z)
