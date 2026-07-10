import pandas as pd
import yaml

from island_v2 import trait_machine_method_eval as evaluator


def _ontology() -> dict:
    return yaml.safe_load(open("config/trait_ontology.yml", encoding="utf-8"))


def test_machine_final_value_maps_raw_colour_to_ontology():
    ontology = _ontology()

    assert evaluator.machine_final_value("flower_primary_color", "pink", ontology) == "red_pink"
    assert evaluator.machine_final_value("flower_primary_color", "yellowish", ontology) == (
        "yellow_orange"
    )
    assert evaluator.machine_final_value("flower_primary_color", "chartreuse", ontology) == ""


def test_machine_selection_collapses_multiple_colour_bins():
    candidates = pd.DataFrame(
        [
            {
                "machine_candidate_id": "a",
                "accepted_species": "Colora example",
                "trait_name": "flower_primary_color",
                "machine_final_value": "white",
                "machine_method": "web_reported_species_direct",
                "source_url": "https://example.test/a",
                "raw_description": "White flowers.",
                "machine_use_tier": "main_text_reported",
            },
            {
                "machine_candidate_id": "b",
                "accepted_species": "Colora example",
                "trait_name": "flower_primary_color",
                "machine_final_value": "red_pink",
                "machine_method": "web_reported_species_direct",
                "source_url": "https://example.test/b",
                "raw_description": "Pink flowers.",
                "machine_use_tier": "main_text_reported",
            },
        ]
    )

    selected = evaluator.select_machine_traits(candidates)

    assert len(selected) == 1
    assert selected.loc[0, "machine_final_value"] == "multicolored_variable"
    assert selected.loc[0, "machine_selection_rule"] == "multiple_explicit_colour_bins"


def test_machine_candidates_from_queue_assigns_method_tiers():
    queue = pd.DataFrame(
        [
            {
                "source_lane": "web_reported_scout",
                "evidence_scope": "species_direct",
                "candidate_class": "reported",
                "accepted_species": "Directa example",
                "trait_layer": "M0",
                "trait_name": "flower_primary_color",
                "candidate_value": "red",
                "source": "wikipedia",
                "source_type": "wikipedia",
                "source_url": "https://example.test/direct",
                "source_citation": "",
                "basis_traits": "",
                "basis_family": "",
                "raw_description": "Red flowers.",
            },
            {
                "source_lane": "trait_proxy_generator",
                "evidence_scope": "proxy_not_reported",
                "candidate_class": "proxy",
                "accepted_species": "Proxya example",
                "trait_layer": "proxy",
                "trait_name": "floral_syndrome_proxy",
                "candidate_value": "open_or_generalist_insect_like",
                "source": "",
                "source_type": "rule_based_proxy",
                "source_url": "",
                "source_citation": "",
                "basis_traits": "flower_primary_color=white",
                "basis_family": "Exampleaceae",
                "raw_description": "Proxy only.",
            },
        ]
    )

    result = evaluator.machine_candidates_from_queue(queue, _ontology())

    assert list(result["machine_method"]) == [
        "web_reported_species_direct",
        "rule_based_trait_proxy",
    ]
    assert list(result["machine_use_tier"]) == [
        "main_text_reported",
        "proxy_sensitivity_only",
    ]
    assert list(result["machine_final_value"]) == [
        "red_pink",
        "open_or_generalist_insect_like",
    ]


def test_method_comparison_includes_globi_and_future_methods():
    species = pd.DataFrame({"accepted_species": ["A a", "B b"]})
    candidates = pd.DataFrame(
        {
            "machine_method": ["web_reported_species_direct"],
            "accepted_species": ["A a"],
            "trait_name": ["flower_primary_color"],
            "machine_final_value": ["white"],
        }
    )
    globi = pd.DataFrame({"query_taxon": ["B b"], "evidence_id": ["g1"]})

    comparison = evaluator.method_comparison(species, candidates, globi)

    assert "web_reported_species_direct" in set(comparison["method"])
    assert "globi_interaction_claims" in set(comparison["method"])
    assert "llm_reported_ecology_excerpt" in set(comparison["method"])


def test_method_comparison_counts_llm_reported_ecology_rows():
    species = pd.DataFrame({"accepted_species": ["A a"]})
    candidates = pd.DataFrame(
        {
            "machine_method": ["llm_reported_ecology_excerpt"],
            "accepted_species": ["A a"],
            "trait_name": ["self_incompatibility"],
            "machine_final_value": ["SC"],
        }
    )

    comparison = evaluator.method_comparison(species, candidates, pd.DataFrame())
    llm = comparison.loc[comparison["method"].eq("llm_reported_ecology_excerpt")].iloc[0]

    assert llm["n_rows"] == 1
    assert llm["n_species"] == 1
    assert llm["n_valid_machine_values"] == 1
