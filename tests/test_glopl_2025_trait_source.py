from __future__ import annotations

import pandas as pd

from island_v2.glopl_2025_trait_source import build_package


def _source() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Author": "Example",
                "Year": "2020",
                "DOI": "10.example/a",
                "Species_accepted_names": "Alpha_one",
                "Family": "Exampleaceae",
                "AF": "yes",
                "BS1": "sc",
                "BS2": "no",
                "BS3": "hm",
                "BS4": "none",
                "BS5": "yes",
                "flower_symmetry_1z_0a": "1",
                "floral_structure": "s",
                "type_of_reward": "n",
                "accessibility_1hard_0easy": "1",
            },
            {
                "Author": "Example",
                "Year": "2020",
                "DOI": "10.example/b",
                "Species_accepted_names": "Beta_two",
                "Family": "Wrongaceae",
                "AF": "no",
                "BS1": "si",
                "BS2": "",
                "BS3": "none",
                "BS4": "",
                "BS5": "",
                "flower_symmetry_1z_0a": "0",
                "floral_structure": "g",
                "type_of_reward": "p",
                "accessibility_1hard_0easy": "0",
            },
        ]
    )


def _master() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "accepted_species": ["Alpha one", "Beta two"],
            "family": ["Exampleaceae", "Exampleaceae"],
        }
    )


def test_strict_mapping_keeps_reproductive_traits_separate() -> None:
    strict, _, audit = build_package(_source(), _master())
    observed = strict.set_index("trait_name")["normalized_value"].to_dict()
    assert observed == {
        "autonomous_selfing_capacity": "autonomous",
        "floral_symmetry": "zygomorphic",
        "self_incompatibility": "SC",
    }
    assert set(strict["axis"]) == {
        "reproductive_assurance",
        "floral_structural_complexity",
    }
    assert not audit.empty


def test_reward_and_specialisation_never_enter_strict_axis() -> None:
    strict, independent, _ = build_package(_source(), _master())
    assert "reward_type" not in set(strict["trait_name"])
    assert "floral_specialisation_class" not in set(strict["trait_name"])
    assert {"reward_type", "floral_specialisation_class"}.issubset(
        set(independent["trait_name"])
    )
    assert not independent["strict_three_axis_included"].any()


def test_family_conflict_is_fail_closed() -> None:
    strict, independent, audit = build_package(_source(), _master())
    assert "Beta two" not in set(strict["accepted_species"])
    assert "Beta two" not in set(independent["accepted_species"])
    excluded = audit.loc[audit["submitted_species"].eq("Beta two")]
    assert set(excluded["reason"]) == {"family_conflict"}


def test_within_source_conflict_is_excluded_not_order_selected() -> None:
    source = pd.concat([_source().iloc[[0]], _source().iloc[[0]]], ignore_index=True)
    source.loc[1, "AF"] = "no"
    strict, _, audit = build_package(source, _master())
    assert not (
        strict["trait_name"].eq("autonomous_selfing_capacity")
        & strict["accepted_species"].eq("Alpha one")
    ).any()
    affected = audit.loc[
        audit["trait_name"].eq("autonomous_selfing_capacity")
        & audit["submitted_species"].eq("Alpha one")
    ]
    assert set(affected["reason"]) == {"within_source_direct_conflict"}


def test_lineage_is_underlying_compilation_not_provider_url() -> None:
    strict, _, _ = build_package(_source(), _master())
    assert strict["source_lineage"].str.startswith(
        "glopl-burns-trait-compilation:"
    ).all()
    assert not strict["source_lineage"].str.contains("nature.com").any()
