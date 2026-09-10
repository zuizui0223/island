from __future__ import annotations

import copy
import hashlib
from pathlib import Path

import pandas as pd
import pytest
import yaml

from island_v2.chapter1_database_manifest import (
    CANONICAL_WORKFLOW,
    Chapter1DatabaseManifest,
    dispatch_command,
    load_manifest,
    validate_rights_registry,
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
    assert db.rights_registry is None


def test_database_version_swap_changes_data_source_not_analysis_contract() -> None:
    payload = yaml.safe_load(Path("config/chapter1_database_versions/current.yml").read_text())
    v1 = Chapter1DatabaseManifest.model_validate(payload)
    v2_payload = copy.deepcopy(payload)
    v2_payload["database"]["version"] = "2.0.0"
    v2_payload["database"]["source"]["run_id"] = 999999
    v2_payload["database"]["source"]["artifact_name"] = "chapter1-database-v2"
    v2_payload["database"]["source"]["sha256"] = "0" * 64
    v2_payload["database"]["rights_registry"] = {
        "relative_path": "source_lineage_registry.csv.gz",
        "sha256": "1" * 64,
        "required_for_public_release": True,
    }
    v2 = Chapter1DatabaseManifest.model_validate(v2_payload)

    v1_cmd = dispatch_command(v1)
    v2_cmd = dispatch_command(v2)
    assert v1.database.analysis_contract == v2.database.analysis_contract
    assert v1_cmd[3] == CANONICAL_WORKFLOW
    assert v2_cmd[3] == CANONICAL_WORKFLOW
    assert "trait_run_id=34191508045" in v1_cmd
    assert "trait_run_id=999999" in v2_cmd


def test_database2_manifest_fails_without_rights_registry() -> None:
    payload = yaml.safe_load(Path("config/chapter1_database_versions/current.yml").read_text())
    payload["database"]["version"] = "2.0.0"
    with pytest.raises(ValueError, match="Database 2.x and later require"):
        Chapter1DatabaseManifest.model_validate(payload)


def _mock_species_axis(tmp_path: Path) -> tuple[Path, list[str]]:
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
                    "source_lineages": "source_a:1|source_b:2",
                    "quality": "high",
                }
            )
    path = tmp_path / "species_axis.csv.gz"
    pd.DataFrame(rows).to_csv(path, index=False, compression="gzip")
    return path, axes


def _mock_registry(tmp_path: Path, include_source_b: bool = True) -> Path:
    rows = [
        {
            "source_lineage": "source_a:1",
            "source_family": "dataset:source_a",
            "provider": "Provider A",
            "source_url": "https://example.org/a",
            "dataset_doi": "10.0000/a",
            "license": "CC-BY-4.0",
            "redistribution_status": "redistributable",
            "rights_evidence": "provider license page",
        }
    ]
    if include_source_b:
        rows.append(
            {
                "source_lineage": "source_b:2",
                "source_family": "dataset:source_b",
                "provider": "Provider B",
                "source_url": "https://example.org/b",
                "dataset_doi": "",
                "license": "",
                "redistribution_status": "review_required",
                "rights_evidence": "not yet resolved",
            }
        )
    registry = tmp_path / "source_lineage_registry.csv.gz"
    pd.DataFrame(rows).to_csv(registry, index=False, compression="gzip")
    return registry


def _mock_manifest(
    species_axis: Path,
    axes: list[str],
    registry: Path | None = None,
    version: str = "2.0.0",
) -> Chapter1DatabaseManifest:
    db: dict[str, object] = {
        "database_id": "test",
        "version": version,
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
            "sha256": hashlib.sha256(species_axis.read_bytes()).hexdigest(),
        },
    }
    if registry is not None:
        db["rights_registry"] = {
            "relative_path": "source_lineage_registry.csv.gz",
            "sha256": hashlib.sha256(registry.read_bytes()).hexdigest(),
            "required_for_public_release": True,
        }
    return Chapter1DatabaseManifest.model_validate({"schema_version": 1, "database": db})


def test_species_axis_validation_is_versioned_and_hash_locked(tmp_path: Path) -> None:
    path, axes = _mock_species_axis(tmp_path)
    manifest = _mock_manifest(path, axes, version="1.1.0")
    report = validate_species_axis(path, manifest)
    assert report["resolved_cells"] == 6
    assert report["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()


def test_database2_rights_registry_covers_all_resolved_lineages(tmp_path: Path) -> None:
    species_axis, axes = _mock_species_axis(tmp_path)
    registry = _mock_registry(tmp_path)
    manifest = _mock_manifest(species_axis, axes, registry)
    report = validate_rights_registry(registry, manifest, species_axis)
    assert report["rows"] == 2
    assert report["redistributable_rows"] == 1
    assert report["review_required_rows"] == 1


def test_database2_rights_registry_fails_closed_on_missing_lineage(tmp_path: Path) -> None:
    species_axis, axes = _mock_species_axis(tmp_path)
    registry = _mock_registry(tmp_path, include_source_b=False)
    manifest = _mock_manifest(species_axis, axes, registry)
    with pytest.raises(ValueError, match="does not cover 1 resolved source lineages"):
        validate_rights_registry(registry, manifest, species_axis)
