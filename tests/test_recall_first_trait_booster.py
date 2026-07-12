from __future__ import annotations

import pandas as pd

from island_v2.recall_first_trait_booster import boost_table, coverage


def test_direct_values_win_and_family_priors_fill_only_gaps() -> None:
    master = pd.DataFrame(
        {
            "accepted_species": ["Plantus directus", "Plantus brassica"],
            "family": ["Lamiaceae", "Brassicaceae"],
        }
    )
    direct = pd.DataFrame(
        [
            {
                "species": "Plantus directus",
                "flower_color": "red",
                "flower_shape": "unknown",
                "pollination_guild": "flies",
                "pollination_notes": "observed fly visitation",
                "mating_system": "unknown",
                "self_incompatibility": "unknown",
                "evidence_type": "field_study",
                "confidence": "medium",
            },
            {
                "species": "Plantus brassica",
                "flower_color": "unknown",
                "flower_shape": "unknown",
                "pollination_guild": "unknown",
                "pollination_notes": "unknown",
                "mating_system": "unknown",
                "self_incompatibility": "unknown",
                "evidence_type": "inference",
                "confidence": "low",
            },
        ]
    )
    priority = pd.DataFrame(
        {
            "species": ["Plantus directus", "Plantus brassica"],
            "flower_color": ["unknown", "unknown"],
            "flower_shape": ["unknown", "unknown"],
            "pollination_guild": ["bees", "unknown"],
            "mating_system": ["unknown", "unknown"],
            "selfing_status": ["autonomous_selfing", "unknown"],
            "self_incompatibility": ["unknown", "unknown"],
        }
    )

    table, provenance = boost_table(master, direct, priority)
    first = table.loc[table["species"].eq("Plantus directus")].iloc[0]
    second = table.loc[table["species"].eq("Plantus brassica")].iloc[0]

    assert first["flower_color"] == "red"
    assert first["pollination_guild"] == "flies"
    assert first["flower_shape"] == "zygomorphic / bilaterally symmetric"
    assert first["mating_system"] == "mainly_selfing"
    assert "field_study" in first["evidence_type"]
    assert "inference" in first["evidence_type"]
    assert second["pollination_guild"] == "bees"
    assert second["self_incompatibility"] == "likely_SI"
    assert second["confidence"] == "low"
    assert set(provenance["mode"]) == {"direct", "likely"}
    assert {"source", "evidence_scope", "confidence"}.issubset(provenance.columns)


def test_simple_synonym_genus_family_cascade() -> None:
    master = pd.DataFrame(
        {
            "accepted_species": ["Alpha one", "Alpha two", "Beta three"],
            "genus": ["Alpha", "Alpha", "Beta"],
            "family": ["Testaceae", "Testaceae", "Testaceae"],
            "synonyms": ["Oldalpha one", "", ""],
        }
    )
    direct = pd.DataFrame(
        [
            {
                "species": "Oldalpha one",
                "flower_color": "red",
                "flower_shape": "tubular",
                "pollination_guild": "unknown",
                "pollination_notes": "synonym account",
                "mating_system": "unknown",
                "self_incompatibility": "unknown",
                "evidence_type": "flora",
                "confidence": "medium",
            },
            {
                "species": "Beta three",
                "flower_color": "white",
                "flower_shape": "open",
                "pollination_guild": "flies",
                "pollination_notes": "direct account",
                "mating_system": "mixed_mating",
                "self_incompatibility": "SC",
                "evidence_type": "flora",
                "confidence": "medium",
            },
        ]
    )
    priority = pd.DataFrame(
        {
            "species": ["Alpha one", "Alpha two", "Beta three"],
            "flower_color": ["unknown", "unknown", "unknown"],
            "flower_shape": ["unknown", "unknown", "unknown"],
            "pollination_guild": ["unknown", "unknown", "unknown"],
            "mating_system": ["unknown", "unknown", "unknown"],
            "selfing_status": ["unknown", "unknown", "unknown"],
            "self_incompatibility": ["unknown", "unknown", "unknown"],
        }
    )

    table, provenance = boost_table(master, direct, priority)
    one = table.loc[table["species"].eq("Alpha one")].iloc[0]
    two = table.loc[table["species"].eq("Alpha two")].iloc[0]

    assert one["flower_color"] == "red"
    assert two["flower_color"] == "red"
    assert two["pollination_guild"] == "flies"
    scopes = set(provenance["evidence_scope"])
    assert "synonym_direct" in scopes
    assert "genus_fallback" in scopes
    assert "family_fallback" in scopes


def test_coverage_counts_non_unknown_values() -> None:
    table = pd.DataFrame(
        {
            "flower_color": ["red", "unknown"],
            "flower_shape": ["unknown", "tubular"],
            "pollination_guild": ["bees", "wind"],
            "mating_system": ["unknown", "mainly_selfing"],
            "self_incompatibility": ["SC", "likely_SI"],
        }
    )
    assert coverage(table) == {
        "flower_color": 1,
        "flower_shape": 1,
        "pollination_guild": 2,
        "mating_system": 1,
        "self_incompatibility": 2,
    }
