import importlib

import pandas as pd

intake = importlib.import_module("island_v2.taxon_intake")


def _taxa() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "accepted_species": "Abrus",
                "genus": "Abrus",
                "family": "Fabaceae",
                "n_islands": 1,
                "n_records": 1,
            },
            {
                "accepted_species": "Asplenium africanum",
                "genus": "Asplenium",
                "family": "Aspleniaceae",
                "n_islands": 1,
                "n_records": 3,
            },
            {
                "accepted_species": "Acacia cyclops",
                "genus": "Acacia",
                "family": "Fabaceae",
                "n_islands": 1,
                "n_records": 8,
            },
        ]
    )


def _matches() -> dict[str, dict[str, object] | None]:
    return {
        "Abrus": {
            "canonicalName": "Abrus",
            "rank": "GENUS",
            "matchType": "EXACT",
            "confidence": 100,
            "class": "Magnoliopsida",
        },
        "Asplenium africanum": {
            "usageKey": 1,
            "scientificName": "Asplenium africanum Desv.",
            "canonicalName": "Asplenium africanum",
            "species": "Asplenium africanum",
            "rank": "SPECIES",
            "status": "ACCEPTED",
            "matchType": "EXACT",
            "confidence": 100,
            "class": "Polypodiopsida",
            "phylum": "Tracheophyta",
        },
        "Acacia cyclops": {
            "usageKey": 2,
            "scientificName": "Acacia cyclops A.Cunn. ex G.Don",
            "canonicalName": "Acacia cyclops",
            "species": "Acacia cyclops",
            "rank": "SPECIES",
            "status": "ACCEPTED",
            "matchType": "EXACT",
            "confidence": 100,
            "class": "Magnoliopsida",
            "phylum": "Tracheophyta",
        },
    }


def test_name_rank_triage_is_conservative():
    assert intake.classify_name_rank("Abrus") == "genus_or_higher_rank"
    assert intake.classify_name_rank("Acacia cyclops") == "binomial_candidate"
    assert intake.classify_name_rank("Acacia cf. cyclops") == "indeterminate_or_placeholder"


def test_taxon_intake_keeps_all_rows_pending_and_separates_fern_from_angiosperm_candidate():
    queue = intake.build_taxon_intake(_taxa(), _matches()).set_index("reported_name")

    assert queue.loc["Abrus", "provisional_match_status"] == "rank_review_required"
    assert queue.loc["Abrus", "taxon_scope_review_status"] == "pending"
    assert queue.loc["Asplenium africanum", "provisional_taxonomic_group"] == "non_seed_vascular"
    assert queue.loc["Asplenium africanum", "review_priority"] == "medium"
    assert queue.loc["Acacia cyclops", "provisional_taxonomic_group"] == "angiosperm"
    assert queue.loc["Acacia cyclops", "provisional_match_status"] == "exact_backbone_candidate"
    assert queue.loc["Acacia cyclops", "taxon_scope_review_status"] == "pending"


def test_establishment_queue_never_infers_native_status_from_gbif_occurrence():
    queue = intake.build_taxon_intake(_taxa(), _matches())
    island_species = pd.DataFrame(
        {
            "island_id": ["island_a", "island_a"],
            "species": ["Acacia cyclops", "Asplenium africanum"],
            "n_records": [8, 3],
            "n_unique_gbif_ids": [8, 3],
            "basis_of_record_set": ["PRESERVED_SPECIMEN", "PRESERVED_SPECIMEN"],
        }
    )
    establishment = intake.build_island_establishment_queue(island_species, queue).set_index("species")

    assert establishment.loc["Acacia cyclops", "island_establishment_status"] == "unresolved"
    assert establishment.loc["Acacia cyclops", "review_status"] == "pending"
    assert establishment.loc["Asplenium africanum", "record_evidence_status"] == "GBIF_occurrence_only"


def test_unmatched_taxon_is_retained_for_review_not_dropped():
    queue = intake.build_taxon_intake(_taxa().iloc[[2]], {"Acacia cyclops": None})

    assert len(queue) == 1
    assert queue.loc[0, "provisional_match_status"] == "unmatched"
    assert queue.loc[0, "provisional_taxonomic_group"] == "unresolved"
    assert queue.loc[0, "taxon_scope_review_status"] == "pending"
