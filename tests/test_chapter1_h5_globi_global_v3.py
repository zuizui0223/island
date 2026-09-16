from __future__ import annotations

import pandas as pd

from island_v2.chapter1_h5_globi_global_v3 import (
    _all_observed_count_matrix,
    fit_within_context,
)


def test_all_observed_count_matrix_ignores_origin_status() -> None:
    flora = pd.DataFrame(
        {
            "island_id": ["a", "a", "b"],
            "accepted_species": ["GenusA one", "GenusA two", "GenusB one"],
            "origin_status": ["native", "unresolved", "introduced"],
        }
    )
    matrix = _all_observed_count_matrix(
        flora,
        {"a": 0, "b": 1},
        {"GenusA": 0, "GenusB": 1},
    ).toarray()
    assert matrix.tolist() == [[2.0, 0.0], [0.0, 1.0]]


def test_within_context_returns_distance_and_area_interaction() -> None:
    rows = []
    for index in range(60):
        distance = float(index) / 10.0
        area = float((index % 7) + 1)
        rows.append(
            {
                "island_id": f"i{index}",
                "stratum": "all_observed",
                "source_mode": "geo_k5",
                "source_matching": "prevalence_richness_effort",
                "min_independent_references": 3,
                "analysis_regime": "tropical",
                "spatial_block": f"b{index % 12}",
                "entry_enrichment": 0.2 * distance - 0.05 * distance * area,
                "log_distance_to_continent_km": distance,
                "log_island_area_km2": area,
                "climate_pc1": float(index % 5),
                "climate_pc2": float(index % 6),
                "climate_pc3": float(index % 4),
                "climate_pc4": float(index % 3),
            }
        )
    config = {
        "model": {"pilot_min_islands": 30, "confirmatory_min_islands": 50},
    }
    result = fit_within_context(pd.DataFrame(rows), config)
    assert len(result) == 1
    row = result.iloc[0]
    assert row["status"] == "fit"
    assert row["support_class"] == "confirmatory"
    assert row["distance_by_area"] < 0
