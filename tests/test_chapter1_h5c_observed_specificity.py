from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_h5c_observed_specificity import fit_observed


def test_fit_classifies_positive_biotic_specificity():
    rng = np.random.default_rng(14)
    rows = []
    for i in range(100):
        distance = -1 + 2 * i / 99
        for mode in ("abiotic_wind", "biotic"):
            biotic = mode == "biotic"
            response = 0.05 * distance + (0.35 * distance if biotic else 0.0) + rng.normal(0, 0.03)
            rows.append(
                {
                    "island_id": f"i{i}",
                    "pollen_vector_mode": mode,
                    "response": response,
                    "n_species": 10,
                    "log_distance_to_continent_km": distance,
                    "spatial_block": f"b{i // 4}",
                    "log_island_area_km2": np.sin(i / 7),
                    "climate_pc1": np.cos(i / 9),
                    "climate_pc2": np.sin(i / 11 + 0.2),
                    "climate_pc3": np.cos(i / 13 + 0.4),
                    "climate_pc4": np.sin(i / 17 + 0.6),
                }
            )
    result = fit_observed(pd.DataFrame(rows))
    assert result["classification"] == "expected_biotic_specificity_supported"
    assert result["interaction_estimate"] > 0
    assert result["interaction_p_value"] <= 0.05
    assert result["biotic_distance_slope"] > result["wind_distance_slope"]
