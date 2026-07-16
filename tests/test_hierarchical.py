import json

import pandas as pd

from island_v2.hierarchical import make_candidates


def test_genus_inference_preserves_category_distribution():
    taxa = pd.DataFrame(
        {
            "accepted_species": ["Genus alpha", "Genus beta", "Genus gamma", "Genus delta"],
            "genus": ["Genus", "Genus", "Genus", "Genus"],
            "family": ["Family", "Family", "Family", "Family"],
        }
    )
    evidence = pd.DataFrame(
        {
            "evidence_id": ["e1", "e2", "e3"],
            "accepted_species": ["Genus alpha", "Genus beta", "Genus gamma"],
            "trait_name": ["flower_primary_color"] * 3,
            "standardized_value": ["white", "white", "red_pink"],
            "evidence_scope": ["species_direct"] * 3,
            "source_reliability": ["B_curated_database_or_institution"] * 3,
            "review_status": ["accepted"] * 3,
        }
    )
    ontology = {
        "traits": {
            "flower_primary_color": {
                "allowed_values": ["white", "red_pink", "unresolved"]
            }
        }
    }
    policy = {
        "accepted_direct_scopes": ["species_direct"],
        "accepted_source_reliability": ["B_curated_database_or_institution"],
        "default": {
            "genus": {"minimum_supporting_species": 3, "minimum_modal_probability": 0.5},
            "family": {"minimum_supporting_species": 10, "minimum_modal_probability": 0.5},
        },
        "trait_rules": {},
    }

    output = make_candidates(taxa, evidence, ontology, policy)
    target = output.loc[output["accepted_species"].eq("Genus delta")].iloc[0]

    assert target["evidence_scope"] == "genus_inference"
    assert target["standardized_value"] == "white"
    assert json.loads(target["category_probabilities_json"]) == {"red_pink": 1 / 3, "white": 2 / 3}
    assert target["analysis_track"] == "expanded_genus_distribution"


def test_disabled_taxonomic_level_emits_no_candidate():
    taxa = pd.DataFrame(
        {
            "accepted_species": ["Genus alpha", "Genus beta", "Genus gamma"],
            "genus": ["Genus"] * 3,
            "family": ["Family"] * 3,
        }
    )
    evidence = pd.DataFrame(
        {
            "evidence_id": ["e1", "e2"],
            "accepted_species": ["Genus alpha", "Genus beta"],
            "trait_name": ["flower_primary_color"] * 2,
            "standardized_value": ["white", "white"],
            "evidence_scope": ["species_direct"] * 2,
            "source_reliability": ["B_curated_database_or_institution"] * 2,
            "review_status": ["accepted"] * 2,
        }
    )
    policy = {
        "accepted_direct_scopes": ["species_direct"],
        "accepted_source_reliability": ["B_curated_database_or_institution"],
        "default": {
            "genus": {
                "enabled": True,
                "minimum_supporting_species": 2,
                "minimum_modal_probability": 1.0,
            },
            "family": {
                "enabled": False,
                "minimum_supporting_species": 2,
                "minimum_modal_probability": 1.0,
            },
        },
        "trait_rules": {
            "flower_primary_color": {
                "genus": {
                    "enabled": False,
                    "minimum_supporting_species": 2,
                    "minimum_modal_probability": 1.0,
                }
            }
        },
    }
    ontology = {"traits": {"flower_primary_color": {"allowed_values": ["white"]}}}

    assert make_candidates(taxa, evidence, ontology, policy).empty
