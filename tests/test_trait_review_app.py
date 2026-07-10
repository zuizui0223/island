import pandas as pd
import yaml

from island_v2.trait_review_app import (
    apply_decision,
    final_value_options,
    next_pending_index,
)


def _ontology() -> dict:
    return yaml.safe_load(open("config/trait_ontology.yml", encoding="utf-8"))


def test_next_pending_index_prefers_review_priority_then_species():
    table = pd.DataFrame(
        {
            "review_priority": ["30", "10", "10"],
            "review_group": [
                "P3_proxy_candidate",
                "P1_species_direct_reported",
                "P1_species_direct_reported",
            ],
            "accepted_species": ["C c", "B b", "A a"],
            "trait_name": ["floral_syndrome_proxy", "flower_primary_color", "flower_primary_color"],
            "candidate_value": ["wind_like", "red", "white"],
            "adjudication_decision": ["", "", ""],
        }
    )

    assert next_pending_index(table, ["P1_species_direct_reported", "P3_proxy_candidate"]) == 2


def test_final_value_options_map_raw_colour_to_ontology_value_first():
    row = pd.Series({"trait_name": "flower_primary_color", "candidate_value": "pink"})

    options = final_value_options(row, _ontology())

    assert options[0] == "red_pink"
    assert "pink" not in options


def test_apply_decision_blanks_final_value_for_non_acceptance():
    table = pd.DataFrame(
        {
            "adjudication_decision": [""],
            "final_value": [""],
            "decision_reason": [""],
            "reviewer": [""],
            "review_date": [""],
        }
    )

    result = apply_decision(
        table,
        0,
        decision="rejected",
        final_value="red_pink",
        reason="Fruit colour, not flower colour.",
        reviewer="reviewer",
        review_date="2026-07-10",
    )

    assert result.loc[0, "adjudication_decision"] == "rejected"
    assert result.loc[0, "final_value"] == ""
    assert result.loc[0, "decision_reason"] == "Fruit colour, not flower colour."
