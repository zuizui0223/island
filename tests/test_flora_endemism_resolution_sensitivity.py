import importlib

import pandas as pd

mod = importlib.import_module("island_v2.flora_endemism_resolution_sensitivity")


def test_stricter_resolution_threshold_reduces_support_without_changing_input():
    rows = []
    for i in range(20):
        n_endemic = 5 + (i % 3)
        n_nonendemic = 25
        n_unresolved = 0 if i < 10 else 10
        rows.append(
            {
                "island_id": f"i{i}",
                "analysis_regime": "tropical",
                "spatial_block": f"b{i}",
                "n_endemic_species": n_endemic,
                "n_native_nonendemic_species": n_nonendemic,
                "n_native_species": n_endemic + n_nonendemic + n_unresolved,
                "distance_to_continent_km": 10 + i,
                "log_distance_to_continent_km": 1 + i / 10,
                "log_island_area_km2": 2 + i / 100,
                "climate_pc1": i / 20,
                "climate_pc2": (i % 5) / 5,
                "climate_pc3": (i % 4) / 4,
                "climate_pc4": (i % 3) / 3,
            }
        )
    raw = pd.DataFrame(rows)
    config = {
        "regime_column": "analysis_regime",
        "spatial_cluster": "spatial_block",
        "baseline_covariates": [
            "log_distance_to_continent_km",
            "log_island_area_km2",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
        ],
        "nominal_z_threshold": 1.96,
        "regimes": ["tropical"],
        "min_resolved_native_species": 1,
    }
    _, support = mod.run_resolution_sensitivity(raw, config, thresholds=(0.0, 0.95))
    indexed = support.set_index("minimum_endemism_resolution_fraction")
    assert indexed.loc[0.95, "n_islands"] < indexed.loc[0.0, "n_islands"]
    assert len(raw) == 20
