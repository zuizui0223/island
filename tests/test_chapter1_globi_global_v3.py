from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_globi_global_v3 import (
    _prepare_model_data,
    _scope_mask,
    classify_primary,
    fit_global_context_tests,
)


def _config() -> dict:
    return {
        "contract": "chapter1_globi_global_v3",
        "predictor": {
            "metric": "effective_channel_number",
            "primary_min_independent_references": 3,
            "sensitivity_min_independent_references": [2, 5],
        },
        "source_matching": {
            "source_modes": ["m1", "m2", "m3", "m4"],
            "primary": "prevalence_richness_effort",
            "sensitivity": "prevalence_richness",
            "minimum_represented_genera": 5,
        },
        "flora_scopes": {
            "primary": "all_observed",
            "sensitivities": ["all_native", "native_nonendemic"],
        },
        "contexts": [
            "northern_midlatitude",
            "northern_high_latitude",
            "tropical",
            "southern_extratropical",
        ],
        "reference_context": "northern_midlatitude",
        "primary_pair": ["northern_midlatitude", "tropical"],
        "model": {
            "response": "entry_enrichment",
            "distance": "log_distance_to_continent_km",
            "area": "log_island_area_km2",
            "controls": ["climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"],
            "cluster": "spatial_block",
        },
        "support": {"confirmatory_min_islands_per_context": 20, "alpha": 0.05},
        "claim_ceiling": {},
    }


def test_scope_mask_includes_unresolved_in_all_observed() -> None:
    frame = pd.DataFrame(
        {
            "origin_status": ["native", "unresolved", "introduced"],
            "floristic_status": ["native_nonendemic", "unresolved", "introduced"],
        }
    )
    assert int(_scope_mask(frame, "all_observed").sum()) == 3
    assert int(_scope_mask(frame, "all_native").sum()) == 1
    assert int(_scope_mask(frame, "native_nonendemic").sum()) == 1


def test_global_context_model_detects_heterogeneous_distance_response() -> None:
    rng = np.random.default_rng(44)
    config = _config()
    enrichment_rows = []
    covariate_rows = []
    slopes = {
        "northern_midlatitude": 0.18,
        "northern_high_latitude": 0.02,
        "tropical": -0.15,
        "southern_extratropical": 0.10,
    }
    island_counter = 0
    for context in config["contexts"]:
        for i in range(70):
            island_counter += 1
            island = f"i{island_counter}"
            distance = rng.normal()
            area = rng.normal()
            controls = rng.normal(size=4)
            covariate_rows.append(
                {
                    "island_id": island,
                    "analysis_regime": context,
                    "spatial_block": f"{context}_{i // 5}",
                    "log_distance_to_continent_km": distance,
                    "log_island_area_km2": area,
                    "climate_pc1": controls[0],
                    "climate_pc2": controls[1],
                    "climate_pc3": controls[2],
                    "climate_pc4": controls[3],
                }
            )
            for mode in config["source_matching"]["source_modes"]:
                enrichment_rows.append(
                    {
                        "island_id": island,
                        "flora_scope": "all_observed",
                        "source_mode": mode,
                        "source_matching": "prevalence_richness_effort",
                        "min_independent_references": 3,
                        "entry_enrichment": (
                            slopes[context] * distance + 0.03 * area + rng.normal(0, 0.04)
                        ),
                    }
                )
    data = _prepare_model_data(
        pd.DataFrame(enrichment_rows), pd.DataFrame(covariate_rows), config
    )
    tests = fit_global_context_tests(data, config)
    primary = tests.loc[
        tests["flora_scope"].eq("all_observed")
        & tests["source_matching"].eq("prevalence_richness_effort")
        & tests["min_independent_references"].eq(3)
    ]
    assert len(primary) == 4
    assert primary["distance_heterogeneity_q"].le(0.05).all()
    classes = classify_primary(tests, config)
    hit = classes.loc[classes["family"].eq("distance_heterogeneity")].iloc[0]
    assert hit["n_source_modes_supported"] == 4
    assert hit["classification"] == "robust_global_context_heterogeneity"
