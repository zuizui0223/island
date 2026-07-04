import importlib

import pandas as pd
import pytest
import typer

guarded = importlib.import_module("island_v2.trait_source_discovery_guarded")


def _staged(island_id: str = "island_core") -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "accepted_species": "Example species",
                "island_id": island_id,
                "genus": "Example",
                "family": "Exampleaceae",
                "trait_candidate_status": "staged_not_curated",
                "release_gate": "accepted taxon scope + establishment review",
            }
        ]
    )


def test_guard_requires_all_staged_rows_to_match_requested_island():
    result = guarded.validate_staged_taxa_context(_staged(), "island_core")

    assert result.iloc[0]["island_id"] == "island_core"

    with pytest.raises(typer.BadParameter, match="exactly island_id"):
        guarded.validate_staged_taxa_context(_staged("wrong_island"), "island_core")


def test_guard_requires_curation_gate_columns_and_rejects_finalized_values():
    missing_gate = _staged().drop(columns=["release_gate"])
    with pytest.raises(typer.BadParameter, match="release_gate"):
        guarded.validate_staged_taxa_context(missing_gate, "island_core")

    finalized = _staged().assign(trait_value="white")
    with pytest.raises(typer.BadParameter, match="prohibited"):
        guarded.validate_staged_taxa_context(finalized, "island_core")
