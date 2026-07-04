import importlib

import pandas as pd

packet = importlib.import_module("island_v2.trait_review_packet")


def test_packet_is_empty_of_biological_claims_but_keeps_staging_gate():
    taxa = pd.DataFrame([{
        "accepted_species": "Begonia annobonensis",
        "genus": "Begonia",
        "family": "Begoniaceae",
        "island_id": "annobon",
        "trait_candidate_status": "staged_not_curated",
        "release_gate": "scope and establishment accepted",
    }])
    ontology = {
        "flower_primary_color": {
            "domain": "floral_signal",
            "allowed_values": ["white", "unresolved"],
        },
        "pollen_vector_mode": {
            "domain": "pollen_vector",
            "allowed_values": ["biotic", "unresolved"],
        },
    }

    result = packet.build_review_packet(taxa, ontology)

    assert len(result) == 2
    assert set(result["evidence_status"]) == {"not_searched"}
    assert set(result["provisional_standardized_value"]) == {""}
    assert set(result["review_status"]) == {"pending"}
    assert set(result["trait_candidate_status"]) == {"staged_not_curated"}
    assert all("Begonia annobonensis" in query for query in result["recommended_search_query"])
