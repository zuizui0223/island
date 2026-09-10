from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from island_v2.chapter1_database_release import (
    build_inventory,
    family_for_lineage,
    load_lineage_family_map,
    source_family,
)


def test_source_family_normalizes_row_level_lineages() -> None:
    assert source_family("dataset:dryad.ghx3ffbrh:v7:row") == "dataset:dryad"
    assert source_family("doi:10.5061/dryad.cc2fqz6hr") == "dataset:dryad"
    assert source_family("origin:austraits:nhnsw_2014_2") == "dataset:austraits"
    assert source_family("pladias:feature:319:abc") == "database:pladias"
    assert source_family("url:https://en.wikipedia.org/wiki/Carex") == "domain:en.wikipedia.org"
    assert source_family("citation:abc123") == "unresolved:citation_hash"


def test_release_inventory_is_fail_closed_by_source_family() -> None:
    frame = pd.DataFrame(
        {
            "source_lineages": [
                "dataset:dryad.ghx3ffbrh:v7:row|pladias:feature:319:abc",
                "doi:10.5061/dryad.cc2fqz6hr",
            ],
        }
    )
    policy = {
        "default_status": "review_required",
        "default_license": None,
        "rules": [
            {
                "pattern": "^dataset:dryad$",
                "status": "redistributable",
                "license": "CC0-1.0",
                "note": "audited",
            }
        ],
    }
    inventory, lineages = build_inventory(frame, policy)
    inventory = inventory.set_index("source_family")
    assert inventory.loc["dataset:dryad", "redistribution_status"] == "redistributable"
    assert inventory.loc["dataset:dryad", "source_license"] == "CC0-1.0"
    assert inventory.loc["dataset:dryad", "resolved_cell_mentions"] == 2
    assert inventory.loc["database:pladias", "redistribution_status"] == "review_required"
    assert inventory.loc["database:pladias", "source_license"] == ""
    assert len(lineages) == 3
    assert not lineages["family_override_applied"].any()


def test_lineage_family_override_collapses_hashes_without_granting_rights(tmp_path: Path) -> None:
    mapping_path = tmp_path / "map.csv"
    pd.DataFrame(
        {
            "source_lineage": [
                "florml-content:Alpha:hash1",
                "florml-content:Beta:hash2",
            ],
            "source_family": [
                "provider_treatment:flora_malesiana",
                "provider_treatment:flora_malesiana",
            ],
        }
    ).to_csv(mapping_path, index=False)
    overrides = load_lineage_family_map(mapping_path)
    assert family_for_lineage("florml-content:Alpha:hash1", overrides) == (
        "provider_treatment:flora_malesiana"
    )

    frame = pd.DataFrame(
        {
            "source_lineages": [
                "florml-content:Alpha:hash1",
                "florml-content:Beta:hash2",
                "florml-content:Gamma:hash3",
            ]
        }
    )
    policy = {"default_status": "review_required", "default_license": None, "rules": []}
    inventory, lineages = build_inventory(frame, policy, overrides)
    families = set(inventory["source_family"])
    assert "provider_treatment:flora_malesiana" in families
    assert "florml-content:gamma" in families
    assert len(inventory) == 2
    normalized = inventory.set_index("source_family")
    assert normalized.loc[
        "provider_treatment:flora_malesiana", "redistribution_status"
    ] == "review_required"
    assert normalized.loc["provider_treatment:flora_malesiana", "source_license"] == ""
    assert int(lineages["family_override_applied"].sum()) == 2


def test_lineage_family_map_fails_on_conflicting_assignments(tmp_path: Path) -> None:
    mapping_path = tmp_path / "bad.csv"
    pd.DataFrame(
        {
            "source_lineage": ["x:1", "x:1"],
            "source_family": ["provider:a", "provider:b"],
        }
    ).to_csv(mapping_path, index=False)
    with pytest.raises(ValueError, match="conflicting lineages"):
        load_lineage_family_map(mapping_path)
