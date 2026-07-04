import importlib

import pandas as pd

pilot = importlib.import_module("island_v2.pilot_curation")


def _registry() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "island_id": "source_ready",
                "canonical_island_name": "Checklist Island",
                "country": "Example",
                "country_iso3": "EXA",
                "source_region_id": "example_region",
                "flora_anchor_status": "bibliographic_anchor_ready",
                "flora_anchor_citation": "Example checklist",
                "flora_anchor_url": "",
                "flora_anchor_scope": "Island checklist",
                "curation_route": "island_checklist_then_taxon_review",
                "identity_review_status": "accepted",
            },
            {
                "island_id": "source_gap",
                "canonical_island_name": "Gap Island",
                "country": "Example",
                "country_iso3": "EXA",
                "source_region_id": "example_region",
                "flora_anchor_status": "source_gap",
                "flora_anchor_citation": "",
                "flora_anchor_url": "",
                "flora_anchor_scope": "",
                "curation_route": "source_gap_then_taxon_review",
                "identity_review_status": "accepted",
            },
        ]
    )


def _manifest() -> pd.DataFrame:
    return pd.DataFrame({"island_id": ["source_ready", "source_gap"], "area_km2": [10, 20]})


def _taxon_intake() -> pd.DataFrame:
    rows = []
    for name, rank, group, match, priority in [
        ("Flower species", "binomial_candidate", "angiosperm", "exact_backbone_candidate", "high"),
        ("Fern species", "binomial_candidate", "non_seed_vascular", "exact_backbone_candidate", "medium"),
        ("GenusOnly", "genus_or_higher_rank", "angiosperm", "rank_review_required", "high"),
    ]:
        rows.append(
            {
                "reported_name": name,
                "reported_genus": name.split()[0],
                "reported_family": "Exampleaceae",
                "n_islands": "1",
                "n_records": "2",
                "name_rank_triage": rank,
                "provisional_accepted_species": name if rank == "binomial_candidate" else "",
                "provisional_taxonomic_group": group,
                "provisional_match_status": match,
                "taxon_scope_review_status": "pending",
                "establishment_review_status": "pending",
                "review_priority": priority,
                "next_review_action": "review",
                "gbif_match_url": "https://example.org/match",
            }
        )
    return pd.DataFrame(rows)


def _establishment() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "island_id": "source_ready",
                "species": "Flower species",
                "provisional_accepted_species": "Flower species",
                "provisional_taxonomic_group": "angiosperm",
                "provisional_match_status": "exact_backbone_candidate",
                "n_records": "8",
                "n_unique_gbif_ids": "8",
                "basis_of_record_set": "PRESERVED_SPECIMEN",
                "island_establishment_status": "unresolved",
                "review_status": "pending",
            },
            {
                "island_id": "source_ready",
                "species": "Fern species",
                "provisional_accepted_species": "Fern species",
                "provisional_taxonomic_group": "non_seed_vascular",
                "provisional_match_status": "exact_backbone_candidate",
                "n_records": "2",
                "n_unique_gbif_ids": "2",
                "basis_of_record_set": "PRESERVED_SPECIMEN",
                "island_establishment_status": "unresolved",
                "review_status": "pending",
            },
            {
                "island_id": "source_gap",
                "species": "Flower species",
                "provisional_accepted_species": "Flower species",
                "provisional_taxonomic_group": "angiosperm",
                "provisional_match_status": "exact_backbone_candidate",
                "n_records": "4",
                "n_unique_gbif_ids": "4",
                "basis_of_record_set": "HUMAN_OBSERVATION",
                "island_establishment_status": "unresolved",
                "review_status": "pending",
            },
            {
                "island_id": "source_ready",
                "species": "GenusOnly",
                "provisional_accepted_species": "",
                "provisional_taxonomic_group": "angiosperm",
                "provisional_match_status": "rank_review_required",
                "n_records": "1",
                "n_unique_gbif_ids": "1",
                "basis_of_record_set": "HUMAN_OBSERVATION",
                "island_establishment_status": "unresolved",
                "review_status": "pending",
            },
        ]
    )


def test_worklist_routes_taxonomy_nonangiosperm_and_source_gap_separately():
    worklist = pilot.build_worklist(_manifest(), _taxon_intake(), _establishment(), _registry())
    stages = {(row.island_id, row.reported_species): row.review_stage for row in worklist.itertuples()}

    assert stages[("source_ready", "Flower species")] == "review_island_establishment"
    assert stages[("source_ready", "Fern species")] == "confirm_non_angiosperm_scope"
    assert stages[("source_ready", "GenusOnly")] == "resolve_taxon_name"
    assert stages[("source_gap", "Flower species")] == "resolve_island_flora_source_gap"


def test_trait_staging_is_blocked_and_does_not_include_source_gap_or_ferns():
    worklist = pilot.build_worklist(_manifest(), _taxon_intake(), _establishment(), _registry())
    staged = pilot.build_trait_staging(worklist, max_per_island=5)

    assert staged["accepted_species"].tolist() == ["Flower species"]
    assert staged.loc[0, "trait_candidate_status"] == "staged_not_curated"
    assert "accepted taxon scope" in staged.loc[0, "release_gate"]


def test_missing_island_context_fails_closed():
    registry = _registry().iloc[[0]]

    try:
        pilot.build_worklist(_manifest(), _taxon_intake(), _establishment(), registry)
    except Exception as exc:
        assert "does not cover island IDs" in str(exc)
    else:
        raise AssertionError("Expected a missing pilot-island context error")
