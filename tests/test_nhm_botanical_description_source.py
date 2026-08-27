from __future__ import annotations

import pandas as pd

from island_v2.nhm_botanical_description_source import REVIEW_COLUMNS, build_package


def _row(**overrides: str) -> dict[str, str]:
    row = {
        "accepted_species": "Alpha beta",
        "searched_taxon": "Alpha beta",
        "axis": "flower_colour",
        "trait_name": "flower_primary_color",
        "normalized_value": "red_pink|white",
        "provider": "flora_china",
        "resource_id": "resource-one",
        "query_url": "https://data.nhm.ac.uk/api/3/action/datastore_search?q=Alpha+beta",
        "source_name": "efloras.flora_of_china",
        "source_record_ids": "record-one",
        "supporting_excerpt": "flower colour = red, white",
        "decision": "accept",
        "reason": "source_specific_structured_cell_reviewed",
        "name_match_method": "accepted_name_exact",
        "reviewer": "reviewer",
        "reviewed_at_utc": "2026-08-27T15:00:00Z",
    }
    row.update(overrides)
    assert set(row) == set(REVIEW_COLUMNS)
    return row


def test_promotes_reviewed_structured_cell_as_medium() -> None:
    evidence, audit = build_package(pd.DataFrame([_row()]))
    assert len(evidence) == 1
    assert evidence.iloc[0].quality == "medium"
    assert evidence.iloc[0].trait_name == "flower_primary_color"
    assert evidence.iloc[0].source_excerpt == "flower colour = red, white"
    assert audit.iloc[0].package_decision == "accept"


def test_strict_synonym_gate_is_allowed_but_fuzzy_match_is_not() -> None:
    accepted = _row(
        accepted_species="Newname beta",
        searched_taxon="Oldname beta",
        name_match_method="wfo_gbif_two_backbone_strict_synonym",
    )
    fuzzy = _row(
        accepted_species="Other beta",
        name_match_method="fuzzy_name",
    )
    evidence, audit = build_package(pd.DataFrame([accepted, fuzzy]))
    assert set(evidence.accepted_species) == {"Newname beta"}
    rejected = audit.loc[audit.accepted_species.eq("Other beta")].iloc[0]
    assert rejected.package_reason == "strict_species_identity_gate_failed"


def test_wfo_stable_identifier_orthographic_variant_is_allowed() -> None:
    variant = _row(
        accepted_species="Osteomeles schweriniae",
        searched_taxon="Osteomeles schwerinae",
        name_match_method="wfo_stable_identifier_orthographic_variant",
    )
    evidence, _ = build_package(pd.DataFrame([variant]))
    assert set(evidence.accepted_species) == {"Osteomeles schweriniae"}


def test_reward_and_pollination_are_not_strict_structure_traits() -> None:
    rows = [
        _row(trait_name="reward_type", normalized_value="nectar"),
        _row(trait_name="pollen_vector_mode", normalized_value="biotic"),
    ]
    evidence, audit = build_package(pd.DataFrame(rows))
    assert evidence.empty
    assert set(audit.package_reason) == {"trait_not_allowed_in_strict_axes"}
