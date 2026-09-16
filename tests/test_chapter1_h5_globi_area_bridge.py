from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_h5_globi_area_bridge import run_bridge


def _config() -> dict:
    return {
        "contract": "chapter1_h5_pollination_triangulation_v2",
        "independent_globi_area_bridge": {
            "response_metric": "effective_channel_number",
            "primary_min_independent_references": 3,
            "sensitivity_min_independent_references": [2, 5],
            "primary_source_matching": "prevalence_richness_effort",
            "sensitivity_source_matching": "prevalence_richness",
            "source_modes": ["geo_k5", "geo_k10", "geo_k20", "geo50_climate10"],
            "strata": ["all_native", "native_nonendemic"],
            "contexts": ["northern_midlatitude", "tropical"],
            "support": {"minimum_islands": 20, "alpha": 0.05},
        },
    }


def test_bridge_recovers_negative_northern_distance_by_area() -> None:
    rng = np.random.default_rng(3)
    cov_rows = []
    enrichment_rows = []
    island_counter = 0
    for context in ["northern_midlatitude", "tropical"]:
        for i in range(80):
            island_counter += 1
            island = f"i{island_counter}"
            distance = rng.normal()
            area = rng.normal()
            controls = rng.normal(size=4)
            cov_rows.append(
                {
                    "island_id": island,
                    "analysis_regime": context,
                    "spatial_block": f"b{i // 4}_{context}",
                    "log_distance_to_continent_km": distance,
                    "log_island_area_km2": area,
                    "climate_pc1": controls[0],
                    "climate_pc2": controls[1],
                    "climate_pc3": controls[2],
                    "climate_pc4": controls[3],
                }
            )
            moderation = -0.20 if context == "northern_midlatitude" else 0.03
            value = 0.10 * distance + 0.02 * area + moderation * distance * area + rng.normal(0, 0.04)
            for stratum in ["all_native", "native_nonendemic"]:
                for mode in ["geo_k5", "geo_k10", "geo_k20", "geo50_climate10"]:
                    for effort in [2, 3, 5]:
                        enrichment_rows.append(
                            {
                                "island_id": island,
                                "metric": "effective_channel_number",
                                "min_independent_references": effort,
                                "stratum": stratum,
                                "source_mode": mode,
                                "source_matching": "prevalence_richness_effort",
                                "entry_enrichment": value,
                            }
                        )
    within, between, classification = run_bridge(
        pd.DataFrame(enrichment_rows), pd.DataFrame(cov_rows), _config()
    )
    north = within.loc[
        within["analysis_regime"].eq("northern_midlatitude")
        & within["min_independent_references"].eq(3)
    ]
    assert not north.empty
    assert north["distance_by_area"].lt(0).all()
    assert classification["promoted"].all()
    assert between.loc[between["min_independent_references"].eq(3), "status"].eq("fit").all()
