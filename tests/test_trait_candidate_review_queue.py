import pandas as pd
import yaml

from island_v2 import trait_candidate_review_queue as queue


def _ontology() -> dict:
    return yaml.safe_load(open("config/trait_ontology.yml", encoding="utf-8"))


def test_review_queue_prioritizes_direct_then_indirect_then_proxy():
    web = pd.DataFrame(
        [
            {
                "accepted_species": "Indirecta example",
                "trait_layer": "M0",
                "trait_name": "flower_primary_color",
                "provisional_candidate_value": "white",
                "matched_term": "white",
                "candidate_class": "reported",
                "source": "wikipedia",
                "source_type": "wikipedia",
                "source_url": "https://example.test/indirect",
                "raw_description": "The flowers are white.",
                "evidence_scope": "species_indirect",
                "review_status": "unreviewed",
            },
            {
                "accepted_species": "Directa example",
                "trait_layer": "M0",
                "trait_name": "flower_primary_color",
                "provisional_candidate_value": "red",
                "matched_term": "red",
                "candidate_class": "reported",
                "source": "wikipedia",
                "source_type": "wikipedia",
                "source_url": "https://example.test/direct",
                "raw_description": "The species has red flowers.",
                "evidence_scope": "species_direct",
                "review_status": "unreviewed",
            },
        ]
    )
    proxies = pd.DataFrame(
        [
            {
                "accepted_species": "Proxya example",
                "trait_name": "floral_syndrome_proxy",
                "proxy_value": "butterfly_or_moth_like",
                "candidate_class": "proxy",
                "basis_traits": "flower_primary_color=red",
                "basis_family": "Exampleaceae",
                "source_type": "rule_based_proxy",
                "raw_description": "Proxy only.",
                "evidence_scope": "proxy_not_reported",
                "review_status": "unreviewed_proxy",
            }
        ]
    )

    result = queue.build_review_queue(web, proxies)

    assert list(result["review_group"]) == [
        "P1_species_direct_reported",
        "P2_species_indirect_source_check",
        "P3_proxy_candidate",
    ]
    assert set(result["adjudication_decision"]) == {""}
    assert set(result["final_value"]) == {""}
    assert result["review_candidate_id"].str.startswith("vpilot_trait:").all()


def test_descriptive_scout_rows_are_kept_after_primary_queue():
    descriptive = pd.DataFrame(
        [
            {
                "accepted_species": "Descriptiva example",
                "trait_layer": "M0",
                "trait_name": "flower_primary_color",
                "provisional_candidate_value": "green",
                "matched_term": "green",
                "candidate_kind": "descriptive_keyword_match",
                "evidence_excerpt": "Leaves are green.",
                "source": "gbif_species_description",
                "description_type": "description",
                "source_url": "https://example.test/gbif",
                "source_citation": "Example source",
                "review_status": "unreviewed",
            }
        ]
    )

    result = queue.build_review_queue(pd.DataFrame(), pd.DataFrame(), descriptive)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["review_group"] == "P4_descriptive_scout_source_check"
    assert row["candidate_kind"] == "descriptive_keyword_match"
    assert row["source_lane"] == "m0_descriptive_scout"


def test_summary_and_accepted_extraction_are_fail_closed():
    result = pd.DataFrame([{column: "" for column in queue.OUTPUT_COLUMNS}])
    result.loc[0, "review_candidate_id"] = "vpilot_trait:test"
    result.loc[0, "accepted_species"] = "Accepta example"
    result.loc[0, "trait_name"] = "flower_primary_color"
    result.loc[0, "candidate_class"] = "reported"
    result.loc[0, "source_lane"] = "web_reported_scout"
    result.loc[0, "review_group"] = "P1_species_direct_reported"

    summary = queue.review_queue_summary(result)
    accepted = queue.accepted_for_curation(result)

    assert summary["n_review_decisions_recorded"] == 0
    assert summary["n_accepted_ready_for_curation"] == 0
    assert accepted.empty

    result.loc[0, "adjudication_decision"] = "accepted"
    result.loc[0, "final_value"] = "red_pink"
    result.loc[0, "decision_reason"] = "Source explicitly states red flowers."

    summary = queue.review_queue_summary(result)
    accepted = queue.accepted_for_curation(result)

    assert summary["n_review_decisions_recorded"] == 1
    assert summary["n_accepted_ready_for_curation"] == 1
    assert list(accepted["final_value"]) == ["red_pink"]
    assert list(accepted["source_type"]) == [""]
    assert list(accepted["raw_description"]) == [""]


def test_review_queue_validation_rejects_invalid_or_incomplete_decisions():
    result = pd.DataFrame([{column: "" for column in queue.OUTPUT_COLUMNS} for _ in range(3)])
    result.loc[0, "review_candidate_id"] = "vpilot_trait:bad_decision"
    result.loc[0, "adjudication_decision"] = "accept"
    result.loc[1, "review_candidate_id"] = "vpilot_trait:missing_fields"
    result.loc[1, "adjudication_decision"] = "accepted"
    result.loc[1, "final_value"] = "white"
    result.loc[2, "review_candidate_id"] = "vpilot_trait:stale_value"
    result.loc[2, "adjudication_decision"] = "rejected"
    result.loc[2, "final_value"] = "red"

    errors = queue.validate_review_queue(result)

    assert any("bad_decision" in error and "accept" in error for error in errors)
    assert any("missing_fields" in error and "decision_reason" in error for error in errors)
    assert any("missing_fields" in error and "reviewer" in error for error in errors)
    assert any("stale_value" in error and "must not carry final_value" in error for error in errors)


def test_accepted_extraction_preserves_provenance_for_curation():
    result = pd.DataFrame([{column: "" for column in queue.OUTPUT_COLUMNS}])
    result.loc[0, "review_candidate_id"] = "vpilot_trait:accepted"
    result.loc[0, "accepted_species"] = "Accepta example"
    result.loc[0, "trait_layer"] = "M0"
    result.loc[0, "trait_name"] = "flower_primary_color"
    result.loc[0, "candidate_value"] = "white"
    result.loc[0, "candidate_class"] = "reported"
    result.loc[0, "candidate_kind"] = "reported_keyword_match"
    result.loc[0, "evidence_scope"] = "species_direct"
    result.loc[0, "source_lane"] = "web_reported_scout"
    result.loc[0, "source"] = "wikipedia"
    result.loc[0, "source_type"] = "wikipedia"
    result.loc[0, "source_url"] = "https://example.test/source"
    result.loc[0, "matched_term"] = "white"
    result.loc[0, "raw_description"] = "Flowers are white."
    result.loc[0, "adjudication_decision"] = "accepted"
    result.loc[0, "final_value"] = "white"
    result.loc[0, "decision_reason"] = "Explicit species-page floral statement."
    result.loc[0, "reviewer"] = "reviewer"
    result.loc[0, "review_date"] = "2026-07-10"

    errors = queue.validate_review_queue(result, ontology=_ontology())
    accepted = queue.accepted_for_curation(result)

    assert errors == []
    assert list(accepted["candidate_value"]) == ["white"]
    assert list(accepted["source_type"]) == ["wikipedia"]
    assert list(accepted["evidence_scope"]) == ["species_direct"]
    assert list(accepted["raw_description"]) == ["Flowers are white."]


def test_review_queue_validation_checks_accepted_final_value_vocabulary():
    result = pd.DataFrame([{column: "" for column in queue.OUTPUT_COLUMNS}])
    result.loc[0, "review_candidate_id"] = "vpilot_trait:bad_value"
    result.loc[0, "trait_name"] = "flower_primary_color"
    result.loc[0, "adjudication_decision"] = "accepted"
    result.loc[0, "final_value"] = "chartreuse"
    result.loc[0, "decision_reason"] = "Reviewer typo."
    result.loc[0, "reviewer"] = "reviewer"
    result.loc[0, "review_date"] = "2026-07-10"

    errors = queue.validate_review_queue(result, ontology=_ontology())

    assert any(
        "bad_value" in error and "outside the allowed vocabulary" in error for error in errors
    )


def test_review_queue_validation_checks_declared_proxy_values():
    result = pd.DataFrame([{column: "" for column in queue.OUTPUT_COLUMNS}])
    result.loc[0, "review_candidate_id"] = "vpilot_trait:proxy"
    result.loc[0, "trait_name"] = "floral_syndrome_proxy"
    result.loc[0, "adjudication_decision"] = "accepted"
    result.loc[0, "final_value"] = "open_or_generalist_insect_like"
    result.loc[0, "decision_reason"] = "Accepted as proxy only."
    result.loc[0, "reviewer"] = "reviewer"
    result.loc[0, "review_date"] = "2026-07-10"

    assert queue.validate_review_queue(result, ontology=_ontology()) == []

    result.loc[0, "final_value"] = "bee_pollinated"
    errors = queue.validate_review_queue(result, ontology=_ontology())

    assert any("outside the allowed vocabulary" in error for error in errors)
