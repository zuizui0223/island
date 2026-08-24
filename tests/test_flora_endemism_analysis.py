import importlib

import numpy as np
import pandas as pd

mod = importlib.import_module("island_v2.flora_endemism_analysis")

CONFIG = {
    "regime_column": "analysis_regime",
    "regimes": ["northern_midlatitude", "tropical", "southern_extratropical"],
    "spatial_cluster": "spatial_block",
    "min_resolved_native_species": 30,
    "nominal_z_threshold": 1.96,
    "baseline_covariates": [
        "log_distance_to_continent_km",
        "log_island_area_km2",
        "climate_pc1",
    ],
}


def test_unresolved_native_endemism_is_not_counted_as_nonendemic():
    frame = pd.DataFrame(
        {
            "island_id": ["i1"],
            "n_native_species": [100],
            "n_native_nonendemic_species": [30],
            "n_endemic_species": [20],
            "analysis_regime": ["northern_midlatitude"],
            "spatial_block": ["b1"],
            "log_distance_to_continent_km": [2.0],
            "log_island_area_km2": [3.0],
            "climate_pc1": [0.0],
        }
    )
    out = mod.prepare_endemism_response(frame, CONFIG).iloc[0]
    assert out["endemism_trials"] == 50
    assert out["endemic_share_resolved_native"] == 0.4
    assert out["endemism_resolution_fraction"] == 0.5
    assert bool(out["endemism_support_eligible"])


def test_support_threshold_uses_resolved_native_count():
    frame = pd.DataFrame(
        {
            "island_id": ["i1", "i2"],
            "n_native_species": [100, 100],
            "n_native_nonendemic_species": [10, 25],
            "n_endemic_species": [5, 5],
            "analysis_regime": ["northern_midlatitude", "northern_midlatitude"],
            "spatial_block": ["b1", "b2"],
            "log_distance_to_continent_km": [1.0, 2.0],
            "log_island_area_km2": [3.0, 4.0],
            "climate_pc1": [-1.0, 1.0],
        }
    )
    out = mod.prepare_endemism_response(frame, CONFIG).set_index("island_id")
    assert not bool(out.loc["i1", "endemism_support_eligible"])
    assert bool(out.loc["i2", "endemism_support_eligible"])


def test_model_support_reports_resolution_fraction():
    rows = []
    for i in range(60):
        regime = "northern_midlatitude" if i < 20 else ("tropical" if i < 40 else "southern_extratropical")
        rows.append(
            {
                "island_id": f"i{i}",
                "n_native_species": 50,
                "n_native_nonendemic_species": 30 + (i % 3),
                "n_endemic_species": 10 + (i % 5),
                "analysis_regime": regime,
                "spatial_block": f"b{i % 8}",
                "log_distance_to_continent_km": np.log1p(i + 1),
                "log_island_area_km2": 2.0 + (i % 7) / 10,
                "climate_pc1": float((i % 9) - 4),
            }
        )
    prepared = mod.prepare_endemism_response(pd.DataFrame(rows), CONFIG)
    _, _, support = mod.fit_endemism_models(prepared, CONFIG)
    assert set(support["regime"]) == set(CONFIG["regimes"])
    assert (support["n_islands"] == 20).all()
    assert support["median_endemism_resolution_fraction"].between(0, 1).all()
