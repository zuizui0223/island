from __future__ import annotations

from pathlib import Path

import yaml


VERSION_DIR = Path("config/island_plant_database_versions")
VERSION = VERSION_DIR / "v2.0.0-alpha1.yml"
CURRENT = VERSION_DIR / "current.yml"


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
