from __future__ import annotations

import pandas as pd

from island_v2.open_web_network_pilot import (
    AXIS_TRAITS,
    audit_sample,
    build_priority_queue,
)


def test_strict_axis_mapping_excludes_reward_and_pollen_vector() -> None:
    structure = set(AXIS_TRAITS["floral_structural_complexity"])
    assert "reward_type" not in structure
    assert "pollen_vector_mode" not in structure
    assert "floral_form" in structure
    assert "flower_size_class" in structure


def test_priority_queue_prefers_two_agreeing_congeners() -> None:
    coverage = pd.DataFrame(
        [
            {"accepted_species": species, "axis": axis, "after_quality": quality}
            for species, quality in (
                ("Alpha one", "high"),
                ("Alpha two", "medium"),
                ("Alpha target", ""),
                ("Beta target", ""),
            )
            for axis in ("flower_colour", "floral_structural_complexity", "reproductive_assurance")
        ]
    )
    evidence = pd.DataFrame(
        [
            {
                "accepted_species": species,
                "trait_name": "floral_form",
                "normalized_value": "tubular",
                "quality": quality,
                "evidence_scope": "species_direct",
                "source_lineage": f"url:{species}",
            }
            for species, quality in (("Alpha one", "high"), ("Alpha two", "medium"))
        ]
    )
    master = pd.DataFrame(
        [
            {"accepted_species": "Alpha one", "genus": "Alpha", "family": "Testaceae"},
            {"accepted_species": "Alpha two", "genus": "Alpha", "family": "Testaceae"},
            {"accepted_species": "Alpha target", "genus": "Alpha", "family": "Testaceae"},
            {"accepted_species": "Beta target", "genus": "Beta", "family": "Testaceae"},
        ]
    )
    queue = build_priority_queue(
        coverage,
        evidence,
        master,
        species_limit=1,
        traits_per_species=1,
    )
    assert queue.iloc[0]["accepted_species"] == "Alpha target"
    assert queue.iloc[0]["trait_name"] == "floral_form"
    assert queue.iloc[0]["priority_reason"] == "two_agree_one_more_unlocks_rule"


def test_audit_sample_is_deterministic_trait_stratified_and_page_unique() -> None:
    candidates = pd.DataFrame(
        [
            {
                "candidate_id": f"{trait}-{index}",
                "trait_name": trait,
                "source_url": f"https://example{index}.org/{trait}",
            }
            for trait in ("floral_form", "flower_primary_color")
            for index in range(8)
        ]
    )
    first = audit_sample(candidates, limit=10)
    second = audit_sample(candidates, limit=10)
    assert first["candidate_id"].tolist() == second["candidate_id"].tolist()
    assert first["source_url"].nunique() == len(first)
    assert first.groupby("trait_name").size().to_dict() == {
        "floral_form": 5,
        "flower_primary_color": 5,
    }
    assert "domain_quality_decision" in first.columns
