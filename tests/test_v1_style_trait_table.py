from __future__ import annotations

import pandas as pd

from island_v2.v1_style_trait_table import build_v1_style_table


def test_builds_requested_nine_column_schema_and_direct_si() -> None:
    species = pd.DataFrame([{"accepted_species": "Plantus alba"}])
    candidates = pd.DataFrame(
        [
            {
                "accepted_species": "Plantus alba",
                "trait_name": "flower_primary_color",
                "provisional_candidate_value": "white",
                "candidate_class": "reported",
                "confidence": "medium",
                "evidence_quote": "Plantus alba has white flowers.",
            },
            {
                "accepted_species": "Plantus alba",
                "trait_name": "floral_form",
                "provisional_candidate_value": "tubular",
                "candidate_class": "reported",
                "confidence": "medium",
                "evidence_quote": "Plantus alba has tubular flowers.",
            },
            {
                "accepted_species": "Plantus alba",
                "trait_name": "pollination_functional_guild",
                "provisional_candidate_value": "other_bees",
                "candidate_class": "reported",
                "confidence": "high",
                "evidence_quote": "Pollinated by bees.",
            },
        ]
    )
    search = pd.DataFrame(
        [
            {
                "accepted_species": "Plantus alba",
                "searched_name": "Plantus alba",
                "result_title": "Plantus alba breeding study",
                "result_snippet": (
                    "Plantus alba is self-compatible and shows mixed mating. "
                    "Flowers are pollinated by bees."
                ),
                "result_url": "https://example.org/journal/article",
            }
        ]
    )

    frame, report = build_v1_style_table(species, candidates, search)

    assert list(frame.columns) == [
        "species",
        "flower_color",
        "flower_shape",
        "pollination_guild",
        "pollination_notes",
        "mating_system",
        "self_incompatibility",
        "evidence_type",
        "confidence",
    ]
    row = frame.iloc[0]
    assert row["flower_color"] == "white"
    assert row["flower_shape"] == "tubular"
    assert row["pollination_guild"] == "bees"
    assert row["mating_system"] == "mixed_mating"
    assert row["self_incompatibility"] == "SC"
    assert row["confidence"] == "high"
    assert report["coverage"]["self_incompatibility"] == 1


def test_unknown_when_information_is_insufficient() -> None:
    species = pd.DataFrame([{"accepted_species": "Plantus ignota"}])
    candidates = pd.DataFrame(
        columns=[
            "accepted_species",
            "trait_name",
            "provisional_candidate_value",
            "candidate_class",
            "confidence",
            "evidence_quote",
        ]
    )
    search = pd.DataFrame(
        columns=[
            "accepted_species",
            "searched_name",
            "result_title",
            "result_snippet",
            "result_url",
        ]
    )

    frame, _ = build_v1_style_table(species, candidates, search)
    row = frame.iloc[0]
    assert row["flower_color"] == "unknown"
    assert row["mating_system"] == "unknown"
    assert row["self_incompatibility"] == "unknown"
    assert row["confidence"] == "low"
