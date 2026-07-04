import importlib

trait = importlib.import_module("island_v2.trait_extraction")
schemas = importlib.import_module("island_v2.schemas")


def _response():
    return schemas.TraitExtractionResponse.model_validate(
        {
            "prompt_version": "v2.1",
            "model_interpretation_note": "Candidate evidence only.",
            "species": [
                {
                    "accepted_species": "Begonia annobonensis",
                    "taxon_id": None,
                    "batch_notes": "Review before use.",
                    "no_evidence_found": False,
                    "trait_candidates": [
                        {
                            "trait_name": "pollen_vector_mode",
                            "raw_value": "Animal-pollinated candidate",
                            "standardized_value": "biotic",
                            "candidate_kind": "source_backed",
                            "evidence_scope": "species_direct",
                            "source_type": "flora_or_monograph",
                            "source_reliability": "A_primary_or_monograph",
                            "source_citation": "Example flora",
                            "source_url": "https://example.org/trait",
                            "evidence_excerpt": "Example paraphrase.",
                            "supporting_taxa": [],
                            "inference_rule": None,
                            "confidence": "medium",
                            "inference_note": "Candidate only.",
                            "needs_human_review": True,
                        }
                    ],
                }
            ],
        }
    )


def test_taxon_context_excludes_identifying_columns_but_preserves_release_gate():
    contexts = trait.taxon_context_by_species(
        [
            {
                "accepted_species": "Begonia annobonensis",
                "genus": "Begonia",
                "family": "Begoniaceae",
                "trait_candidate_status": "staged_not_curated",
                "release_gate": "Accepted flora membership required.",
                "island_id": "annobon",
            }
        ]
    )

    context = contexts["Begonia annobonensis"]
    assert "genus" not in context
    assert context["trait_candidate_status"] == "staged_not_curated"
    assert context["island_id"] == "annobon"


def test_flatten_candidates_carries_staging_context_without_curating_trait():
    rows = trait.flatten_candidates(
        _response(),
        "run_1",
        [{"url": "https://example.org/trait", "title": "Example"}],
        {
            "Begonia annobonensis": {
                "trait_candidate_status": "staged_not_curated",
                "release_gate": "Accepted flora membership required.",
                "island_id": "annobon",
            }
        },
    )

    assert len(rows) == 1
    assert rows[0]["input_trait_candidate_status"] == "staged_not_curated"
    assert rows[0]["input_release_gate"] == "Accepted flora membership required."
    assert "annobon" in rows[0]["input_taxon_context_json"]
    assert rows[0]["review_status"] == "pending"
