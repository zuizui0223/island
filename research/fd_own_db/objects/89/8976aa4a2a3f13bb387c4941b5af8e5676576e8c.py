"""No private flora/GIFT data: synthetic N2 block-fold leakage tests."""
import pandas as pd
import pytest

from fd.n2_source_coverage_audit import (
    balanced_block_folds, retained_donor_names, geographic_fold_coverage,
)


def test_balanced_folds_never_split_source_geographic_block():
    eligible = pd.DataFrame({
        "island_id": [f"i{x}" for x in range(12)],
        "spatial_block": ["A"] * 4 + ["B"] * 3 + ["C"] * 2 +
                         ["D", "E", "F"],
    })
    fold = balanced_block_folds(eligible, folds=3)
    assert len(fold) == 6
    assignment = eligible.spatial_block.map(fold)
    assert all(assignment[eligible.spatial_block == b].nunique() == 1
               for b in eligible.spatial_block.unique())
    sizes = assignment.value_counts()
    assert max(sizes) - min(sizes) <= 1
    assert fold == balanced_block_folds(eligible.sample(frac=1, random_state=19),
                                        folds=3)
    with pytest.raises(ValueError):
        balanced_block_folds(eligible.assign(spatial_block=None), folds=3)


def test_global_holdout_removes_exact_species_from_every_mainland_source():
    species = {"G one", "G two", "H one", "H two"}
    withheld = {"G one", "H two"}
    assert retained_donor_names(species, withheld) == {"G two", "H one"}
    with pytest.raises(ValueError):
        retained_donor_names(species, set())


def test_fold_source_genus_support_uses_no_focal_taxon_names():
    # Two geography-disjoint island folds: A includes G one; B includes G two.
    # Both species would occur in the same GIFT source without fold exclusion.
    membership = pd.DataFrame({
        "island_id": ["i1", "i1", "i2", "i2"],
        "accepted_species": ["G one", "H one", "G two", "J one"],
        "genus": ["G", "H", "G", "J"],
        "spatial_block": ["A", "A", "B", "B"],
    })
    donors = {"G one", "G two", "H one", "J one", "G three"}
    source = {"S1": {"G one", "G two", "H one", "G three"},
              "S2": {"G two", "J one"}}
    assigned = {"i1": {"S1"}, "i2": {"S2"}}
    result = geographic_fold_coverage(
        membership, donors, source, assigned, {"A": 0, "B": 1},
        min_islands=1)
    assert len(result) == 2
    a = next(x for x in result if x["fold"] == 0)
    b = next(x for x in result if x["fold"] == 1)
    assert a["n_global_donor_names_after_holdout"] == 3
    assert a["fraction_focal_carriers_with_local_source_genus_donor"] == .5
    # G two is withheld in fold 1; S2 contains J one, which is ALSO withheld.
    assert b["fraction_focal_carriers_with_local_source_genus_donor"] == 0.
    assert all("accepted_species" not in r and "island_id" not in r and
               "spatial_block" not in r for r in result)


def test_small_cell_suppresses_both_global_and_local_effect_quantities():
    membership = pd.DataFrame({
        "island_id": ["one", "one", "two"],
        "accepted_species": ["G a", "H a", "G b"],
        "genus": ["G", "H", "G"],
        "spatial_block": ["A", "A", "B"],
    })
    result = geographic_fold_coverage(
        membership, {"G a", "G b", "H a", "G c"},
        {"src": {"G a", "G b", "H a", "G c"}},
        {"one": {"src"}, "two": {"src"}}, {"A": 0, "B": 1},
        min_islands=5)
    assert all(r["status"] == "small_cell_suppressed" for r in result)
    assert all("fraction_focal_carriers_with_global_genus_donor" not in r and
               "fraction_focal_carriers_with_local_source_genus_donor" not in r
               for r in result)
