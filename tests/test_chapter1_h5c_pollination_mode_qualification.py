from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_h5c_pollination_mode_qualification import (
    collapse_pollination_mode,
    simulate_cell,
)


def test_collapse_pollination_mode_excludes_mixed_and_conflicting_species():
    gift = pd.DataFrame(
        [
            {"scientific_name": "Alpha one", "trait_name": "pollen_vector_mode", "trait_value": "biotic"},
            {"scientific_name": "Beta two", "trait_name": "pollen_vector_mode", "trait_value": "abiotic_wind"},
            {"scientific_name": "Gamma three", "trait_name": "pollen_vector_mode", "trait_value": "mixed"},
            {"scientific_name": "Delta four", "trait_name": "pollen_vector_mode", "trait_value": "biotic"},
            {"scientific_name": "Delta four", "trait_name": "pollen_vector_mode", "trait_value": "abiotic_wind"},
        ]
    )
    out = collapse_pollination_mode(gift).set_index("accepted_species")
    assert set(out.index) == {"Alpha one", "Beta two"}
    assert out.loc["Alpha one", "pollen_vector_mode"] == "biotic"
    assert out.loc["Beta two", "pollen_vector_mode"] == "abiotic_wind"


def _design() -> pd.DataFrame:
    rows = []
    for i in range(80):
        island = f"i{i}"
        distance = -1 + 2 * i / 79
        for mode in ("biotic", "abiotic_wind"):
            rows.append(
                {
                    "island_id": island,
                    "pollen_vector_mode": mode,
                    "n_species": 8 + i % 4,
                    "log_distance_to_continent_km": distance,
                    "analysis_regime": "tropical",
                    "spatial_block": f"b{i // 4}",
                    "log_island_area_km2": np.sin(i / 8),
                    "climate_pc1": np.cos(i / 9),
                    "climate_pc2": np.sin(i / 11 + 0.3),
                    "climate_pc3": np.cos(i / 13 + 0.2),
                    "climate_pc4": np.sin(i / 17 + 0.7),
                }
            )
    return pd.DataFrame(rows)


def test_simulation_detects_prespecified_interaction_without_inflated_null():
    out = simulate_cell(
        _design(),
        effects=[0.0, 0.75],
        replicates=120,
        seed=5,
        cluster_sd=0.2,
        residual_sd=1.0,
        alpha=0.05,
    )
    null = out.loc[np.isclose(out["interaction_effect_sd"], 0.0)].iloc[0]
    strong = out.loc[np.isclose(out["interaction_effect_sd"], 0.75)].iloc[0]
    assert null["detection_rate"] <= 0.15
    assert strong["detection_rate"] >= 0.80
