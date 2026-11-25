import pandas as pd
import pytest
from fd.functional_turnover import audit_pairs, safe_summary, validate_trait_table, TRAITS


def _row(island, species, genus, value):
    return {
        "island_id": island, "taxon_name": species, "genus": genus,
        **{t: value for t in TRAITS}
    }


def test_species_removal_and_pairwise_prevalence():
    traits = pd.DataFrame([
        _row("A", "T1", "G1", 1),
        _row("A", "T2", "G2", 0),
        _row("B", "T1", "G1", 1),
        _row("B", "T3", "G3", 1),
    ])
    pairs = pd.DataFrame([{
        "island_a": "A", "island_b": "B", "comparison": "near"
    }])
    scored = audit_pairs(traits, pairs, min_carriers=1)
    sc = scored.loc[scored.trait.eq("self_compatibility")].set_index("trait_scope")
    assert sc.loc["all_taxa", "abs_prevalence_gap"] == 0.5
    assert sc.loc["species_exclusive", "abs_prevalence_gap"] == 1.0
    assert sc.loc["genus_exclusive", "abs_prevalence_gap"] == 1.0
    assert (sc["taxon_jaccard"] == 1 / 3).all()


def test_missing_traits_are_not_zero():
    r = [_row("A", "T1", "G1", 1), _row("B", "T2", "G2", 0)]
    r[0]["self_compatibility"] = None
    traits = pd.DataFrame(r)
    pairs = pd.DataFrame([{
        "island_a": "A", "island_b": "B", "comparison": "far"
    }])
    scored = audit_pairs(traits, pairs, min_carriers=1)
    sc = scored[
        (scored.trait == "self_compatibility") &
        (scored.trait_scope == "all_taxa")
    ].iloc[0]
    assert pd.isna(sc.abs_prevalence_gap)
    assert sc.n_a == 0


def test_conflicting_species_state_fails():
    r = [_row("A", "T1", "G1", 1), _row("B", "T1", "G1", 0)]
    with pytest.raises(ValueError, match="Inconsistent self_compatibility"):
        validate_trait_table(pd.DataFrame(r))


def test_duplicate_island_taxon_fails():
    r = [_row("A", "T1", "G1", 1)] * 2
    with pytest.raises(ValueError, match="Duplicate"):
        validate_trait_table(pd.DataFrame(r))


def test_group_summary_does_not_export_identities_and_suppresses_small_n():
    rows = pd.DataFrame([{
        "comparison": "near", "trait_scope": "all_taxa", "domain": "test",
        "trait": "self_compatibility", "island_a": f"secret{i}",
        "island_b": f"other{i}", "abs_prevalence_gap": .2
    } for i in range(12)])
    assert safe_summary(rows, minimum_pairs=13).empty
    result = safe_summary(rows, minimum_pairs=10)
    assert len(result) == 1
    assert not {"island_a", "island_b"} & set(result.columns)
    assert result.iloc[0].n_pairs == 12
