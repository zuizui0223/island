from pathlib import Path

import pandas as pd

from island_v2.chapter1_p1c_matched_genus_qualification import (
    make_matched_pseudo_taxonomy,
)


def test_pseudo_taxonomy_preserves_family_group_sizes() -> None:
    taxonomy = pd.DataFrame(
        {
            "accepted_species": [f"sp{i}" for i in range(10)],
            "family": ["F1"] * 6 + ["F2"] * 4,
            "genus": ["A", "A", "A", "B", "B", "C", "D", "D", "E", "E"],
        }
    )
    pseudo, audit = make_matched_pseudo_taxonomy(taxonomy, seed=123)
    assert audit["group_size_multiset_preserved_all_families"] is True
    assert audit["n_pseudo_genera"] == 5
    for family in ("F1", "F2"):
        observed_sizes = sorted(
            taxonomy.loc[taxonomy["family"].eq(family)].groupby("genus").size().tolist()
        )
        pseudo_sizes = sorted(
            pseudo.loc[pseudo["family"].eq(family)].groupby("genus").size().tolist()
        )
        assert observed_sizes == pseudo_sizes


def test_pseudo_taxonomy_is_deterministic() -> None:
    taxonomy = pd.DataFrame(
        {
            "accepted_species": [f"sp{i}" for i in range(8)],
            "family": ["F"] * 8,
            "genus": ["A", "A", "A", "B", "B", "C", "C", "D"],
        }
    )
    one, _ = make_matched_pseudo_taxonomy(taxonomy, seed=2026091501)
    two, _ = make_matched_pseudo_taxonomy(taxonomy, seed=2026091501)
    pd.testing.assert_series_equal(one["genus"], two["genus"])
