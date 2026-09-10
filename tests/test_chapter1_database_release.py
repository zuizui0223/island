from __future__ import annotations

import pandas as pd

from island_v2.chapter1_database_release import build_inventory, source_family


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
