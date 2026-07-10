import pandas as pd

from island_v2 import trait_candidate_review_queue as queue


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
    result = pd.DataFrame(
        [
            {
                column: ""
                for column in queue.OUTPUT_COLUMNS
            }
        ]
    )
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
