from __future__ import annotations

from island_v2 import interaction_evidence as evidence


def test_explicit_flower_visit_is_retained_without_effectiveness_inference():
    record = evidence._record(
        "Plantus alba",
        "target",
        {
            "interaction_type": "visits flowers of",
            "source_taxon_name": "Bee species",
            "target_taxon_name": "Plantus alba",
            "study_source_id": "example-dataset",
            "study_citation": "Example citation",
        },
    )

    assert record is not None
    assert record["interaction_claim_class"] == "flower_visit_claim"
    assert record["partner_taxon_name"] == "Bee species"
    assert record["query_taxon_exact_in_record"] == "true"
    assert record["effectiveness_class"] == "not_assessed"
    assert record["review_status"] == "unreviewed_source_backed"


def test_explicit_pollination_claim_is_retained():
    record = evidence._record(
        "Plantus alba",
        "source",
        {
            "interaction_type": "pollinated by",
            "source_taxon_name": "Plantus alba",
            "target_taxon_name": "Bee species",
        },
    )

    assert record is not None
    assert record["interaction_claim_class"] == "pollination_claim"
    assert record["partner_taxon_name"] == "Bee species"


def test_non_flower_interaction_is_not_exported_as_pollination_evidence():
    record = evidence._record(
        "Plantus alba",
        "source",
        {
            "interaction_type": "eats",
            "source_taxon_name": "Herbivore species",
            "target_taxon_name": "Plantus alba",
        },
    )

    assert record is None
