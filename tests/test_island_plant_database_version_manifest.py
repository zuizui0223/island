from __future__ import annotations

from pathlib import Path

import yaml


VERSION_DIR = Path("config/island_plant_database_versions")
VERSION = VERSION_DIR / "v2.0.0-alpha1.yml"
CURRENT = VERSION_DIR / "current.yml"
TAXONOMY_VERSION = VERSION_DIR / "v2.0.0-alpha2-taxonomy.yml"
CURRENT_TAXONOMY = VERSION_DIR / "current_taxonomy.yml"


def _load(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def test_current_manifest_is_exact_alpha1_version() -> None:
    assert CURRENT.read_bytes() == VERSION.read_bytes()


def test_alpha1_receipt_pins_canonical_main_build() -> None:
    manifest = _load(VERSION)
    database = manifest["database"]
    scope = manifest["scientific_scope"]
    rights = manifest["rights"]
    build = manifest["canonical_build"]

    assert database == {
        "database_id": "global_island_plant_database",
        "version": "2.0.0-alpha1",
        "version_status": "frozen_alpha",
        "primary_object": "island_x_taxon",
    }
    assert scope["island_universe"]["islands"] == 8265
    assert scope["taxa"]["provisional_species_name_identities"] == 115328
    assert scope["candidate_flora"]["island_taxon_rows"] == 1039757
    assert scope["candidate_flora"]["evidence_rows"] == 1039757
    assert scope["candidate_flora"]["islands_with_exact_occurrence_records"] == 4549
    assert scope["candidate_flora"]["islands_with_species_candidate_records"] == 4505
    assert scope["candidate_flora"]["islands_with_occurrence_but_no_species_candidate"] == 44
    assert scope["candidate_flora"]["absence_inferred_from_missing_records"] is False

    assert rights["gbif_downloads"] == 103
    assert rights["candidate_value_release_status"] == "review_required"
    assert rights["redistributable_candidate_values"] is False

    assert build["git_sha"] == "8a6d3e4b984d019ce18f575e58c59ab9bd31cbbf"
    assert build["workflow_run_id"] == 34665438815
    assert build["artifact_id"] == 10288721725
    assert build["artifact_digest"] == (
        "sha256:3a5af33e60f9eb4156b94f4d5e0f91cf36fc417ee55074af2c042a1378b28a40"
    )
    assert build["validation_valid"] is True

    expected_files = {
        "islands.csv": "c2d9bb54db34c005c34ea650c37b2cc01cd6d39beb638a7ea8faeaa6d3d4de1b",
        "taxa.csv.gz": "6680a2f100f8389c672290185c24cd3dd5cc58681d027e512f6d9ae9719e6055",
        "island_taxa.csv.gz": "c5b35d6d684eb1d5aee8d78cc4e5664e3d11e7a655ade1d8cafec7abd4401abf",
        "evidence.csv.gz": "21e73023b25c5092e7dc7e1af61b9097c7279d1eb81a839571cd3346a5e137bd",
        "RIGHTS_LEDGER.csv": "fe676224195fb9d74e68826ce6ded1d25670fe619e1cbc210db56420971bbf67",
        "DATABASE_MANIFEST.json": "5e778b435d813b00efec0e8cf1e45a12fb647ed8163177f9f1ec7e762dd7c63e",
        "VALIDATION_REPORT.json": "9302794cbaad5eb9df083b963fd5deaf8a4328cd74388755e8b8690aa5f7d08e",
    }
    assert {name: entry["sha256"] for name, entry in build["files"].items()} == expected_files


def test_alpha1_is_not_claimed_as_publicly_deposited() -> None:
    publication = _load(VERSION)["publication"]
    assert publication["release_status"] == "not_deposited"
    assert publication["zenodo_doi"] is None
    assert publication["zenodo_record_url"] is None


def test_current_taxonomy_manifest_is_exact_alpha2_taxonomy_version() -> None:
    assert CURRENT_TAXONOMY.read_bytes() == TAXONOMY_VERSION.read_bytes()


def test_alpha2_taxonomy_receipt_pins_reproduced_main_build() -> None:
    manifest = _load(TAXONOMY_VERSION)
    database = manifest["database"]
    scope = manifest["scientific_scope"]
    taxonomy = scope["taxonomy"]
    boundary = scope["scientific_boundary"]
    build = manifest["canonical_build"]
    reproducibility = manifest["reproducibility_check"]

    assert database == {
        "database_id": "global_island_plant_database",
        "version": "2.0.0-alpha2-taxonomy",
        "version_status": "frozen_alpha_layer",
        "layer_type": "taxonomy_crosswalk",
        "base_database_version": "2.0.0-alpha1",
    }
    assert scope["source_taxa"]["provisional_species_name_identities"] == 115328
    assert scope["source_taxa"]["source_taxa_summary_sha256"] == (
        "c0586264d9a877b88c26e866264f2c77b20423cd93380a8eebe269582d53352e"
    )

    assert taxonomy["authority"] == "Catalogue of Life Extended Release"
    assert taxonomy["checklist_key"] == "7ddf754f-d193-4cc9-b351-99906754a03b"
    assert taxonomy["match_type_counts"] == {
        "EXACT": 112406,
        "HIGHERRANK": 466,
        "NONE": 2065,
        "VARIANT": 391,
    }
    assert taxonomy["automatic_resolution_candidates"] == 111561
    assert taxonomy["exact_accepted_candidates"] == 107439
    assert taxonomy["exact_synonym_to_accepted_candidates"] == 4122
    assert taxonomy["review_required"] == 1702
    assert taxonomy["unmatched"] == 2065
    assert taxonomy["review_queue_rows"] == 3767

    assert boundary["alpha1_mutated"] is False
    assert boundary["automatic_results_are_candidates"] is True
    assert boundary["review_required_results_promoted"] is False
    assert boundary["unmatched_results_promoted"] is False
    assert boundary["native_status_inferred"] is False
    assert boundary["introduced_status_inferred"] is False
    assert boundary["establishment_status_inferred"] is False
    assert boundary["endemic_status_inferred"] is False
    assert boundary["absence_inferred"] is False

    assert build["git_sha"] == "20fd40cbdde3a93519c561c19a851cd1a1e1f4a6"
    assert build["workflow_run_id"] == 34668274725
    assert build["artifact_id"] == 10289759727
    assert build["artifact_digest"] == (
        "sha256:b46d9be4701c561aa2e53a35a5aa694c5bb83ac618e75f087f35a9337b9f961b"
    )
    assert build["validation_valid"] is True
    assert build["segment_count"] == 12
    assert build["metadata_canonical_sha256"] == (
        "0d37b729825bb849f0e2a45bbc04748a7906176d270f0b48c2994f5912904b3e"
    )

    expected_files = {
        "taxonomy_crosswalk.csv.gz": "079d9c27b0ca3aa51f0c7310f43ad4901b07a319d6063eab5a7d62b291bdd9e3",
        "taxonomy_review_queue.csv.gz": "66669f60d96e9356e3d7e6ad6b803de9a15b60be608443e7694b0eee600828bf",
        "gbif_colxr_metadata.json": "a80dfd4980c0262e7b7d3ac462dcf2a1978744437202f2b0bbd0fbf15fbf8cd9",
        "TAXONOMY_MANIFEST.json": "d534aab496caf9f9ca230548eb7343c8d5205eb7c36a80386a42c04147c55db1",
    }
    assert {name: entry["sha256"] for name, entry in build["files"].items()} == expected_files

    assert reproducibility["independent_pr_build"]["workflow_run_id"] == 34667780773
    assert reproducibility["independent_pr_build"]["artifact_id"] == 10289900308
    assert reproducibility["scientific_files_byte_identical_to_canonical"] is True
    assert reproducibility["repeated_live_api_summary_identical"] is True


def test_alpha2_taxonomy_is_not_claimed_as_publicly_deposited() -> None:
    publication = _load(TAXONOMY_VERSION)["publication"]
    assert publication["release_status"] == "not_deposited"
    assert publication["zenodo_doi"] is None
    assert publication["zenodo_record_url"] is None
