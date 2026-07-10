import json

import pandas as pd
import yaml

from island_v2 import (
    reported_ecology_machine as ecology,
    trait_machine_method_eval as evaluator,
    trait_multimodal_candidates as multimodal,
)


def _ontology() -> dict:
    return yaml.safe_load(open("config/trait_ontology.yml", encoding="utf-8"))


def test_normalize_value_maps_reported_ecology_aliases():
    ontology = _ontology()

    assert ecology.normalize_value("self_incompatibility", "self-compatible", ontology) == "SC"
    assert ecology.normalize_value("autonomous_selfing_capacity", "autogamy", ontology) == (
        "autonomous"
    )
    assert ecology.normalize_value("pollination_functional_guild", "Bombus", ontology) == (
        "bumblebees"
    )
    assert ecology.normalize_value("self_incompatibility", "maybe", ontology) == ""


def test_extract_candidates_from_text_sources_preserves_excerpt_and_scope():
    sources = pd.DataFrame(
        [
            {
                "accepted_species": "Ecologia example",
                "source_text": (
                    "Ecologia example is self-compatible and shows autonomous selfing. "
                    "Its flowers are pollinated by bees."
                ),
                "source_url": "https://example.test/source",
                "source_citation": "Example source",
                "source_type": "paper_abstract",
                "evidence_scope": "species_direct",
            }
        ]
    )

    candidates, holdouts = ecology.extract_candidates_from_text_sources(sources, _ontology())

    assert holdouts.empty
    assert set(candidates["trait_name"]) == {
        "self_incompatibility",
        "autonomous_selfing_capacity",
        "pollination_functional_guild",
    }
    assert set(candidates["candidate_value"]) == {"SC", "autonomous", "other_bees"}
    assert candidates["source_excerpt"].str.contains("self-compatible").any()
    assert set(candidates["evidence_scope"]) == {"species_direct"}


def test_jsonl_import_validates_and_normalizes_candidates(tmp_path):
    path = tmp_path / "llm.jsonl"
    path.write_text(
        "\n".join(
            [
                json.dumps(
                    {
                        "accepted_species": "Jsonia example",
                        "trait_name": "self_incompatibility",
                        "candidate_value": "self-incompatible",
                        "source_excerpt": "The species is self-incompatible.",
                        "source_url": "https://example.test/paper",
                    }
                ),
                json.dumps(
                    {
                        "accepted_species": "Jsonia example",
                        "trait_name": "self_incompatibility",
                        "candidate_value": "unsupported",
                        "source_excerpt": "Unsupported value.",
                    }
                ),
            ]
        ),
        encoding="utf-8",
    )

    candidates, errors = ecology.read_jsonl_candidates(path, _ontology())

    assert len(candidates) == 1
    assert len(errors) == 1
    assert candidates.loc[0, "candidate_value"] == "SI"


def test_reported_ecology_candidates_flow_to_queue_and_machine_eval():
    candidates = pd.DataFrame(
        [
            {
                "accepted_species": "Flowia example",
                "trait_name": "self_incompatibility",
                "candidate_value": "SC",
                "source_url": "https://example.test/paper",
                "source_citation": "Example paper",
                "source_excerpt": "Flowia example is self-compatible.",
                "source_type": "paper_abstract",
                "extraction_model": "llm_or_rule",
                "prompt_id": "reported_ecology_v1",
                "evidence_scope": "species_direct",
            }
        ]
    )

    queue = multimodal.build_multimodal_review_queue(pd.DataFrame(), candidates)
    machine = evaluator.machine_candidates_from_queue(queue, _ontology())
    selected = evaluator.select_machine_traits(machine)

    assert len(queue) == 1
    assert queue.loc[0, "source_lane"] == "llm_reported_ecology"
    assert machine.loc[0, "machine_method"] == "llm_reported_ecology_excerpt"
    assert selected.loc[0, "machine_final_value"] == "SC"


def test_indirect_reported_ecology_stays_sensitivity_only():
    candidates = pd.DataFrame(
        [
            {
                "accepted_species": "Indirecta example",
                "trait_name": "self_incompatibility",
                "candidate_value": "SC",
                "source_url": "https://example.test/paper",
                "source_citation": "Example paper",
                "source_excerpt": "A related taxon is self-compatible.",
                "source_type": "paper_abstract",
                "extraction_model": "llm_or_rule",
                "prompt_id": "reported_ecology_v1",
                "evidence_scope": "species_indirect",
            }
        ]
    )

    queue = multimodal.build_multimodal_review_queue(pd.DataFrame(), candidates)
    machine = evaluator.machine_candidates_from_queue(queue, _ontology())
    selected = evaluator.select_machine_traits(machine)

    assert machine.loc[0, "machine_method"] == "llm_reported_ecology_excerpt"
    assert machine.loc[0, "machine_use_tier"] == "sensitivity_text_reported"
    assert selected.empty
