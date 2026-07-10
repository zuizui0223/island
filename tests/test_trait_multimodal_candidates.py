import pandas as pd
import yaml

from island_v2 import (
    trait_candidate_review_queue as review_queue,
    trait_multimodal_candidates as multimodal,
)


def _ontology() -> dict:
    return yaml.safe_load(open("config/trait_ontology.yml", encoding="utf-8"))


def test_visual_candidates_become_visible_trait_review_rows_only():
    visual = pd.DataFrame(
        [
            {
                "accepted_species": "Visua example",
                "trait_name": "flower_primary_color",
                "candidate_value": "red_pink",
                "image_url": "https://example.test/flower.jpg",
                "model_name": "vision-model",
                "model_confidence": "0.72",
                "evidence_note": "Flower appears red-pink.",
            },
            {
                "accepted_species": "Visua example",
                "trait_name": "self_incompatibility",
                "candidate_value": "SC",
                "image_url": "https://example.test/flower.jpg",
            },
        ]
    )

    result = multimodal.build_multimodal_review_queue(visual, pd.DataFrame())

    assert len(result) == 1
    row = result.iloc[0]
    assert row["source_lane"] == "visual_trait_candidate"
    assert row["review_group"] == "P1_visual_trait_candidate"
    assert row["trait_name"] == "flower_primary_color"
    assert row["source_type"] == "image"
    assert row["source_url"] == "https://example.test/flower.jpg"
    assert "model_name=vision-model" in row["basis_traits"]


def test_reported_ecology_candidates_preserve_llm_excerpt_for_review():
    reported = pd.DataFrame(
        [
            {
                "accepted_species": "Ecologia example",
                "trait_name": "self_incompatibility",
                "candidate_value": "SC",
                "source_url": "https://example.test/paper",
                "source_citation": "Author 2026",
                "source_excerpt": "Controlled crosses showed the species is self-compatible.",
                "extraction_model": "text-model",
                "prompt_id": "reported_ecology_v1",
            }
        ]
    )

    result = multimodal.build_multimodal_review_queue(pd.DataFrame(), reported)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["source_lane"] == "llm_reported_ecology"
    assert row["review_group"] == "P1_reported_ecology_llm_source_check"
    assert row["trait_layer"] == "M1"
    assert row["raw_description"] == "Controlled crosses showed the species is self-compatible."
    assert "extraction_model=text-model" in row["basis_traits"]


def test_multimodal_queue_validates_after_human_acceptance():
    visual = pd.DataFrame(
        [
            {
                "accepted_species": "Accepta visualis",
                "trait_name": "floral_form",
                "candidate_value": "tubular",
                "image_url": "https://example.test/tubular.jpg",
            }
        ]
    )
    result = multimodal.build_multimodal_review_queue(visual, pd.DataFrame())
    result.loc[0, "adjudication_decision"] = "accepted"
    result.loc[0, "final_value"] = "tubular"
    result.loc[0, "decision_reason"] = "Visible corolla form is tubular."
    result.loc[0, "reviewer"] = "reviewer"
    result.loc[0, "review_date"] = "2026-07-10"

    errors = review_queue.validate_review_queue(result, ontology=_ontology())

    assert errors == []


def test_templates_include_expected_input_columns():
    templates = multimodal.template_tables()

    assert "visual_trait_candidates_template.csv" in templates
    assert "reported_ecology_candidates_template.csv" in templates
    assert list(templates["visual_trait_candidates_template.csv"].columns) == (
        multimodal.VISUAL_TEMPLATE_COLUMNS
    )
    assert list(templates["reported_ecology_candidates_template.csv"].columns) == (
        multimodal.REPORTED_ECOLOGY_TEMPLATE_COLUMNS
    )
