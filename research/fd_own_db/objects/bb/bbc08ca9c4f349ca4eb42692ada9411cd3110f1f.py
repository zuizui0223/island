import hashlib

import pandas as pd
import pytest

from fd.functional_matched_ledger_audit import (
    validate_paired_ledgers, equal_block_paired_summary,
)


def _private_rows(n_islands=45, per_island=23):
    result = []
    for i in range(n_islands):
        for seed in ("20261110", "20261111"):
            for j in range(per_island):
                name = f"Taxon_{i:03d}_{j:02d}"
                fold = hashlib.sha256(name.encode()).digest()[0] % 5
                for ledger, real in (("snapshot", .70), ("source", .80)):
                    result.append({
                        "scope": "All", "trait": "self_incompatibility",
                        "region": "northern_high_latitude",
                        "island_id": f"secret_island_{i}",
                        "control_seed": seed, "taxon_name": name,
                        "spatial_block": f"secret_block_{i//3}",
                        "ledger": ledger,
                        "source_real": "real_geography",
                        "source_placebo_1": "placebo_a",
                        "source_placebo_2": "placebo_b",
                        "fold": fold, "state": 1, "p_real": real,
                        "p_placebo_1": .52, "p_placebo_2": .48,
                        "max_coverage_difference_1": .01,
                        "max_coverage_difference_2": .015,
                        "pool_ratio_1": .95, "pool_ratio_2": 1.15,
                        "source_jaccard_1": .40, "source_jaccard_2": .50,
                    })
    return pd.DataFrame(result)


def test_same_target_same_placebo_blocked_comparison_is_private_safe():
    out = equal_block_paired_summary(_private_rows(), draws=299)
    assert len(out) == 2
    assert set(out.n_same_target_islands) == {45}
    assert set(out.n_spatial_blocks) == {15}
    assert (out.paired_gain_delta > 0).all()
    assert (out.paired_ci_low > 0).all()
    assert not {"island_id", "taxon_name", "spatial_block",
                "source_real", "source_placebo_1"} & set(out)


def test_same_target_name_support_is_mandatory():
    x = _private_rows()
    row = x.index[(x.ledger == "snapshot")][0]
    with pytest.raises(ValueError, match="sets differ"):
        validate_paired_ledgers(x.drop(index=row))


def test_same_binary_target_state_is_mandatory():
    x = _private_rows()
    row = x.index[(x.ledger == "snapshot")][0]
    x.loc[row, "state"] = 0
    with pytest.raises(ValueError, match="conflicting|Different"):
        validate_paired_ledgers(x)


def test_actual_and_placebo_sources_must_match_between_ledgers():
    x = _private_rows()
    row = x.index[(x.ledger == "snapshot")][0]
    x.loc[row, "source_placebo_1"] = "different_placebo"
    with pytest.raises(ValueError, match="Different|Conflicting"):
        validate_paired_ledgers(x)


def test_name_fold_and_donor_calipers_are_checked():
    x = _private_rows()
    x.loc[0, "fold"] = (int(x.loc[0, "fold"])+1) % 5
    with pytest.raises(ValueError, match="fold"):
        validate_paired_ledgers(x)
    x = _private_rows()
    x.loc[0, "max_coverage_difference_2"] = .05
    with pytest.raises(ValueError, match="coverage"):
        validate_paired_ledgers(x)
    x = _private_rows()
    x.loc[0, "source_jaccard_1"] = .85
    with pytest.raises(ValueError, match="overlap"):
        validate_paired_ledgers(x)


def test_nonindependent_or_low_carrier_groups_are_not_reported():
    x = _private_rows()
    x["spatial_block"] = "same_block"
    assert equal_block_paired_summary(x, draws=299).empty
    x = _private_rows(per_island=19)
    with pytest.raises(ValueError, match="Under-supported"):
        equal_block_paired_summary(x, draws=299)


def test_two_control_seeds_not_independent_replications():
    x = _private_rows()
    summary = equal_block_paired_summary(x, draws=299)
    assert len(summary) == 2
    assert set(summary.control_seed) == {"20261110", "20261111"}
    assert summary.n_same_target_islands.sum() == 90  # repeated 45, not 90 distinct
