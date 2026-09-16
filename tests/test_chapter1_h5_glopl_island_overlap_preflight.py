from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import pandas as pd
import yaml
from shapely.geometry import Polygon

from island_v2.chapter1_h5_glopl_island_overlap_preflight import (
    build_study_key,
    evaluate_admission_gate,
    match_glopl_to_islands,
    read_glopl_metadata,
    summarize_overlap,
)


def _config() -> dict:
    with open(
        "config/chapter1_h5_glopl_island_overlap_preflight_v1.yml",
        encoding="utf-8",
    ) as handle:
        return yaml.safe_load(handle)


def test_reader_never_materializes_forbidden_effect_size_columns(tmp_path: Path) -> None:
    path = tmp_path / "GloPL.csv"
    pd.DataFrame(
        {
            "Latitude": [1.0],
            "Longitude": [2.0],
            "DOI": ["10.1/example"],
            "Author": ["A"],
            "Year": ["2000"],
            "Species_accepted_names": ["Plant a"],
            "PL_Effect_Size": [99.0],
            "Natural_X_FS": [123.0],
        }
    ).to_csv(path, index=False)

    got = read_glopl_metadata(path, _config())
    assert list(got.columns) == _config()["allowed_columns"]
    assert "PL_Effect_Size" not in got.columns
    assert "Natural_X_FS" not in got.columns


def test_study_key_prefers_doi_and_falls_back_to_author_year() -> None:
    frame = pd.DataFrame(
        {
            "DOI": [" 10.1/ABC ", "", ""],
            "Author": ["X", "Smith", ""],
            "Year": ["1999", "2005", ""],
        }
    )
    key = build_study_key(frame)
    assert key.tolist() == ["doi:10.1/abc", "author_year:Smith|2005", ""]


def test_matching_uses_unique_within_and_keeps_boundary_as_sensitivity_only() -> None:
    islands = gpd.GeoDataFrame(
        [
            {"island_id": "i1", "geometry": Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])},
            {"island_id": "i2", "geometry": Polygon([(2, 0), (3, 0), (3, 1), (2, 1)])},
        ],
        crs="EPSG:4326",
    )
    glopl = pd.DataFrame(
        {
            "Latitude": [0.5, 0.5, 0.5],
            "Longitude": [0.5, 1.0, 2.5],
            "DOI": ["a", "b", "c"],
            "Author": ["A", "B", "C"],
            "Year": ["2000", "2001", "2002"],
            "Species_accepted_names": ["P a", "P b", "P c"],
        }
    )
    matched, audit = match_glopl_to_islands(glopl, islands, _config())
    assert matched.loc[matched["DOI"].eq("a"), "island_id"].iloc[0] == "i1"
    assert matched.loc[matched["DOI"].eq("c"), "island_id"].iloc[0] == "i2"
    boundary = matched.loc[matched["DOI"].eq("b")].iloc[0]
    assert pd.isna(boundary["island_id"])
    assert bool(boundary["boundary_only_match"])
    assert audit["boundary_only_match_count"] == 1


def test_summary_deduplicates_sites_and_gate_requires_both_contexts() -> None:
    matched = pd.DataFrame(
        [
            {"row_id": 0, "island_id": "n1", "analysis_regime": "northern_midlatitude", "site_key": "1|1", "study_key": "doi:a", "Species_accepted_names": "P a"},
            {"row_id": 1, "island_id": "n1", "analysis_regime": "northern_midlatitude", "site_key": "1|1", "study_key": "doi:a", "Species_accepted_names": "P b"},
            {"row_id": 2, "island_id": "t1", "analysis_regime": "tropical", "site_key": "2|2", "study_key": "doi:b", "Species_accepted_names": "P c"},
        ]
    )
    summary = summarize_overlap(matched, _config())
    north = summary["by_analysis_regime"]["northern_midlatitude"]
    assert north["n_rows"] == 2
    assert north["n_sites"] == 1
    assert north["n_islands"] == 1
    assert north["n_studies"] == 1

    decision = evaluate_admission_gate(summary, _config())
    assert not decision["admitted"]
    assert not decision["context_pass"]["northern_midlatitude"]
    assert not decision["context_pass"]["tropical"]
