from __future__ import annotations

import copy
import hashlib
from pathlib import Path

import pandas as pd
import yaml

from island_v2.chapter1_database_manifest import (
    CANONICAL_WORKFLOW,
    Chapter1DatabaseManifest,
    dispatch_command,
    load_manifest,
    validate_species_axis,
)


def test_current_manifest_is_bound_to_frozen_contract() -> None:
    manifest = load_manifest(Path("config/chapter1_database_versions/current.yml"))
    db = manifest.database
    assert db.version == "1.0.0"
    assert db.analysis_contract == "chapter1_progressive_analysis_v1"
    assert db.resolved_cells == 222688
    assert db.axis_resolved_cells["reproductive_assurance"] == 48497
    assert db.source.run_id == 34191508045
    assert db.source.relative_path == "species_axis_coverage.csv.gz"


def test_database_version_swap_changes_data_source_not_analysis_contract() -> None:
    payload = yaml.safe_load(Path("config/chapter1_database_versions/current.yml").read_text())
    v1 = Chapter1DatabaseManifest.model_validate(payload)
    v2_payload = copy.deepcopy(payload)
    v2_payload["database"]["version"] = "2.0.0"
    v2_payload["database"]["source"]["run_id"] = 999999
    v2_payload["database"]["source"]["artifact_name"] = "chapter1-database-v2"
    v2_payload["database"]["source"]["sha256"] = "0" * 64
    v2 = Chapter1DatabaseManifest.model_validate(v2_payload)

    v1_cmd = dispatch_command(v1)
    v2_cmd = dispatch_command(v2)
    assert v1.database.analysis_contract == v2.database.analysis_contract
    assert v1_cmd[3] == CANONICAL_WORKFLOW
    assert v2_cmd[3] == CANONICAL_WORKFLOW
    assert "trait_run_id=34191508045" in v1_cmd
    assert "trait_run_id=999999" in v2_cmd


def test_species_axis_validation_is_versioned_and_hash_locked(tmp_path: Path) -> None:
    axes = ["flower_colour", "floral_structural_complexity", "reproductive_assurance"]
    rows = []
    for species in ("Alpha one", "Beta two"):
        for axis in axes:
            rows.append(
                {
                    "accepted_species": species,
                    "axis": axis,
                    "trait_composition": "state_a",
                    "trait_names": "trait_a",
                    "source_groups": "source_a",
                    "source_lineages": "source_a:1",
                    "quality": "high",
                }
            )
    path = tmp_path / "species_axis.csv.gz"
    pd.DataFrame(rows).to_csv(path, index=False, compression="gzip")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = Chapter1DatabaseManifest.model_validate(
        {
            "schema_version": 1,
            "database": {
                "database_id": "test",
                "version": "2.0.0",
                "analysis_contract": "chapter1_progressive_analysis_v1",
                "denominator_species": 2,
                "axes": axes,
                "resolved_cells": 6,
                "axis_resolved_cells": {axis: 2 for axis in axes},
                "source": {
                    "kind": "github_artifact",
                    "repository": "zuizui0223/island",
                    "run_id": 1,
                    "artifact_name": "mock-v2",
                    "relative_path": "species_axis.csv.gz",
                    "sha256": digest,
                },
            },
        }
    )
    report = validate_species_axis(path, manifest)
    assert report["resolved_cells"] == 6
    assert report["sha256"] == digest
