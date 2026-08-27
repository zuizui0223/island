from __future__ import annotations

import pandas as pd

from island_v2.floraweb_trait_source import _colour_states, build_source_package


def _master(*species: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "accepted_species": species,
            "family": ["Exampleaceae"] * len(species),
        }
    )


def test_colour_extraction_keeps_flower_multistate_and_rejects_leaf_colour() -> None:
    states, quote = _colour_states(
        "Blätter blaugrün. Blüten weiß bis rosa, Kronblätter mit gelbem Fleck. "
        "Frucht rot."
    )
    assert states == ["red_pink", "white", "yellow_orange"]
    assert "Blätter blaugrün" not in quote
    assert "Frucht rot" not in quote


def test_colour_extraction_rejects_faded_flower_colour() -> None:
    states, _ = _colour_states("Krone 5 mm lang. Welke Blüten braun, nicht abfallend.")
    assert states == []


def test_controlled_reproductive_traits_are_not_substituted() -> None:
    raw = pd.DataFrame(
        [
            ["1", "Alpha beta", "biologie", "SI-Reaktion", "selbstkompatibel"],
            ["1", "Alpha beta", "biologie", "Befruchtungstyp", "xenogam"],
            [
                "1",
                "Alpha beta",
                "biologie",
                "Bestäubung (Pollenvektoren)",
                "in der Regel Selbstbestäubung",
            ],
            [
                "1",
                "Alpha beta",
                "biologie",
                "Bestäubung (Pollenvektoren)",
                "in der Regel Windbestäubung",
            ],
            [
                "1",
                "Alpha beta",
                "biologie",
                "Belohnung für Bestäuber",
                "vorhanden Nektar",
            ],
        ],
        columns=["name_usage_id", "canonical_name", "page", "label", "value"],
    )
    strict, independent, _, names = build_source_package(raw, _master("Alpha beta"))
    values = dict(zip(strict["trait_name"], strict["normalized_value"], strict=False))
    assert values == {
        "autonomous_selfing_capacity": "autonomous",
        "mating_system": "predominantly_outcrossing",
        "self_incompatibility": "SC",
    }
    independent_values = dict(
        zip(independent["trait_name"], independent["normalized_value"], strict=False)
    )
    assert independent_values == {
        "reward_type": "nectar",
        "pollen_vector_mode": "abiotic_wind",
    }
    assert independent["strict_three_axis_included"].eq("false").all()
    assert names.iloc[0]["accepted"]


def test_unmatched_name_is_audited_and_never_accepted() -> None:
    raw = pd.DataFrame(
        [["2", "Other species", "biologie", "SI-Reaktion", "selbstinkompatibel"]],
        columns=["name_usage_id", "canonical_name", "page", "label", "value"],
    )
    strict, independent, _, names = build_source_package(raw, _master("Alpha beta"))
    assert strict.empty
    assert independent.empty
    assert names.iloc[0]["name_match_method"] == "unmatched_no_synonym_claim"
    assert not names.iloc[0]["accepted"]


def test_structured_flower_type_maps_only_anatomical_classes() -> None:
    raw = pd.DataFrame(
        [
            [
                "3",
                "Alpha beta",
                "biologie",
                "Blumentyp (nach Kugler 1970)",
                "Trichterblumen, großblütig",
            ],
            [
                "3",
                "Alpha beta",
                "biologie",
                "Blumentyp (nach Kugler 1970)",
                "Pollenblumen",
            ],
        ],
        columns=["name_usage_id", "canonical_name", "page", "label", "value"],
    )
    strict, _, _, _ = build_source_package(raw, _master("Alpha beta"))
    assert strict[["trait_name", "normalized_value"]].to_dict("records") == [
        {"trait_name": "floral_form", "normalized_value": "funnel_trumpet"}
    ]
