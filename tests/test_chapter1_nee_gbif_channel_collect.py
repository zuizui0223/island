from __future__ import annotations

import pandas as pd
import pytest

from island_v2.chapter1_nee_gbif_channel_collect import (
    exact_island_records,
    summarize_background,
)


def test_exact_island_records_drops_buffer_only_rows_and_sorts() -> None:
    assigned = pd.DataFrame(
        [
            {"island_id": "i2", "gbif_id": "3", "block_id": "b2", "species": "Species c"},
            {"island_id": pd.NA, "gbif_id": "2", "block_id": "b1", "species": "Mainland buffer"},
            {"island_id": "i1", "gbif_id": "1", "block_id": "b1", "species": "Species a"},
        ]
    )
    result = exact_island_records(assigned)
    assert list(result["island_id"]) == ["i1", "i2"]
    assert set(result["gbif_id"]) == {"1", "3"}
    assert "Mainland buffer" not in set(result["species"])


def test_exact_island_records_requires_assignment_column() -> None:
    with pytest.raises(ValueError, match="requires island_id"):
        exact_island_records(pd.DataFrame({"gbif_id": ["1"]}))


def test_summarize_background_preserves_auditable_counts() -> None:
    records = pd.DataFrame(
        [
            {"island_id": "i1", "species": "Species a", "dataset_key": "d1", "year": "2020"},
            {"island_id": "i1", "species": "Species a", "dataset_key": "d2", "year": "2022"},
            {"island_id": "i1", "species": "Species b", "dataset_key": "d2", "year": "2021"},
            {"island_id": "i2", "species": "", "dataset_key": "d3", "year": ""},
        ]
    )
    result = summarize_background(records, "Lepidoptera")
    i1 = result.loc[result["island_id"].eq("i1")].iloc[0]
    i2 = result.loc[result["island_id"].eq("i2")].iloc[0]
    assert i1["campaign_taxon_name"] == "Lepidoptera"
    assert i1["n_records"] == 3
    assert i1["n_species"] == 2
    assert i1["n_datasets"] == 2
    assert i1["year_min"] == 2020
    assert i1["year_max"] == 2022
    assert i2["n_records"] == 1
    assert i2["n_species"] == 0
    assert pd.isna(i2["year_min"])
    assert pd.isna(i2["year_max"])


def test_summarize_background_empty_is_schema_stable() -> None:
    result = summarize_background(pd.DataFrame(), "Diptera")
    assert result.empty
    assert list(result.columns) == [
        "island_id",
        "campaign_taxon_name",
        "n_records",
        "n_species",
        "n_datasets",
        "year_min",
        "year_max",
    ]
