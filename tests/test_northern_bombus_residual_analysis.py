import importlib

import numpy as np
import pandas as pd
import pytest

mod = importlib.import_module("island_v2.northern_bombus_residual_analysis")


def _data(n=90):
    null_rows = []
    cov_rows = []
    for i in range(n):
        distance = 1.0 + (i % 30) / 10
        compatibility = 0.2 + (i % 17) / 25
        climate1 = np.sin(i / 7)
        # Compatibility carries signal not reducible to the baseline synthetic climate.
        response = 0.45 * compatibility + 0.02 * distance + (i % 3) * 0.002
        null_rows.append(
            {
                "island_id": f"i{i}",
                "outcome": "generalized_form",
                "stratum": "native_nonendemic",
                "observed_n_species": 20 + (i % 7),
                "deviation_observed_minus_null": response,
            }
        )
        cov_rows.append(
            {
                "island_id": f"i{i}",
                "analysis_regime": "northern_midlatitude",
                "spatial_block": f"b{i // 3}",
                "log_distance_to_continent_km": distance,
                "log_island_area_km2": 2.0 + (i % 11) / 10,
                "climate_pc1": climate1,
                "climate_pc2": np.cos(i / 5),
                "climate_pc3": (i % 4) / 4,
                "climate_pc4": (i % 5) / 5,
                "environmental_compatibility_max": compatibility,
            }
        )
    return pd.DataFrame(null_rows), pd.DataFrame(cov_rows)


def test_bombus_compatibility_adds_incremental_cv_signal_in_synthetic_data():
    null, cov = _data()
    coef, cv, support = mod.fit_northern_bombus_models(
        null, cov, strata=("native_nonendemic",)
    )
    assert not coef.empty
    assert support.iloc[0]["n_islands"] == 90
    scored = cv.set_index("model")
    assert scored.loc["M1_add_bombus_compatibility", "status"] == "scored"
    assert scored.loc["M1_add_bombus_compatibility", "rmse_improvement_vs_m0"] > 0
    focal = coef.loc[
        (coef["model"] == "M1_add_bombus_compatibility")
        & (coef["predictor"] == "environmental_compatibility_max")
    ].iloc[0]
    assert focal["estimate"] > 0


def test_non_northern_rows_are_excluded():
    null, cov = _data(40)
    cov["analysis_regime"] = "tropical"
    coef, cv, support = mod.fit_northern_bombus_models(
        null, cov, strata=("native_nonendemic",)
    )
    assert coef.empty
    assert cv.empty
    assert support.empty


def test_missing_compatibility_fails_closed():
    null, cov = _data(40)
    with pytest.raises(Exception):
        mod.fit_northern_bombus_models(
            null,
            cov.drop(columns=["environmental_compatibility_max"]),
            strata=("native_nonendemic",),
        )
