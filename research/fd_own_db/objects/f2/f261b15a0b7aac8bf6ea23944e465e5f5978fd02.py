import numpy as np
import pandas as pd
import pytest

from fd.functional_rarefaction import audited_rarefaction


def _counts(n=40, m_species=35, m_genus=25, positive=10):
    rows = []
    for stage, count in [
        ("all", n), ("species_exclusive", m_species),
        ("genus_exclusive", m_genus),
    ]:
        rows.append({
            "scope": "All", "region": "test", "trait_scope": stage,
            "island_a": "a", "island_b": "b",
            "n_self_compatibility": count, "n_self_compatibility_b": count,
            "k_self_compatibility_a": min(positive, count),
            "k_self_compatibility_b": min(positive, count),
        })
    return pd.DataFrame(rows)


def _blocks():
    return pd.DataFrame([
        {"island_id": "a", "spatial_block": "x"},
        {"island_id": "b", "spatial_block": "y"},
    ])


def _small(df, **kwargs):
    return audited_rarefaction(
        df, _blocks(), traits=("self_compatibility",),
        min_carriers=1, min_pairs=5, min_islands=5,
        min_spatial_blocks=5, draws=100, **kwargs,
    )


def test_same_carriers_give_exact_zero_rarefaction_residual():
    rows = []
    for i in range(10):
        for stage in ("all", "species_exclusive", "genus_exclusive"):
            rows.append({
                "scope": "All", "region": "north", "trait_scope": stage,
                "island_a": f"a{i}", "island_b": f"b{i}",
                "n_self_compatibility": 40, "n_self_compatibility_b": 40,
                "k_self_compatibility_a": 12, "k_self_compatibility_b": 18,
            })
    x = pd.DataFrame(rows)
    blocks = pd.DataFrame([
        {"island_id": v, "spatial_block": v}
        for v in list(x.island_a.unique()) + list(x.island_b.unique())
    ])
    z = audited_rarefaction(
        x, blocks, traits=("self_compatibility",), draws=100,
        min_pairs=5, min_islands=5, min_spatial_blocks=5,
    )
    assert len(z) == 2
    assert np.allclose(z.mean_paired_excess_gap, 0, atol=1e-12)
    assert np.allclose(z.median_observed_gap, 0.15)
    assert not {"island_a", "island_b", "spatial_block"} & set(z.columns)


def test_reject_impossible_retained_positives():
    x = _counts()
    x.loc[x.trait_scope == "genus_exclusive", "k_self_compatibility_a"] = 11
    with pytest.raises(ValueError, match="Impossible"):
        _small(x)


def test_spatial_metadata_required():
    with pytest.raises(ValueError, match="Missing blocks"):
        audited_rarefaction(_counts(), _blocks().iloc[:1],
                           traits=("self_compatibility",), draws=100)


def test_reject_repeated_pair_stage():
    x = _counts()
    with pytest.raises(ValueError, match="Duplicate"):
        _small(pd.concat([x, x.iloc[[0]]]))


def test_suppress_small_pair_or_island_groups():
    y = _small(_counts())
    assert y.empty
    assert "n_pair_comparisons" in y.columns


def test_stage_coding_and_subset_consistency():
    x = _counts()
    x.loc[x.trait_scope == "genus_exclusive", "k_self_compatibility_a"] = 0
    assert _small(x).empty
    x.loc[x.trait_scope == "genus_exclusive", "k_self_compatibility_a"] = 11
    with pytest.raises(ValueError, match="Impossible"):
        _small(x)
