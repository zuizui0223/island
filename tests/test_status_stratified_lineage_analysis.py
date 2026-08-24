import importlib

import numpy as np
import pandas as pd
import pytest

mod = importlib.import_module("island_v2.status_stratified_lineage_analysis")


def _synthetic_null(n=60):
    rows = []
    for i in range(n):
        distance = 1.0 + i / 10
        rows.append(
            {
                "island_id": f"i{i}",
                "outcome": "plain_colour",
                "stratum": "native_nonendemic",
                "observed_n_species": 10 + (i % 5),
                "deviation_observed_minus_null": 0.4 * distance + (i % 3) * 0.01,
            }
        )
    return pd.DataFrame(rows)


def _covariates(n=60):
    return pd.DataFrame(
        {
            "island_id": [f"i{i}" for i in range(n)],
            "analysis_regime": ["tropical"] * n,
            "spatial_block": [f"b{i // 3}" for i in range(n)],
            "log_distance_to_continent_km": [1.0 + i / 10 for i in range(n)],
            "log_island_area_km2": [2.0 + (i % 7) / 10 for i in range(n)],
            "climate_pc1": [np.sin(i) for i in range(n)],
            "climate_pc2": [np.cos(i) for i in range(n)],
            "climate_pc3": [(i % 4) / 4 for i in range(n)],
            "climate_pc4": [(i % 5) / 5 for i in range(n)],
        }
    )


def test_weighted_clustered_linear_recovers_positive_distance_signal():
    data = _synthetic_null().merge(_covariates(), on="island_id")
    coef, fit = mod.fit_weighted_linear_clustered(
        data,
        response_column="deviation_observed_minus_null",
        weight_column="observed_n_species",
        predictors=mod.DEFAULT_BASELINE,
        cluster_column="spatial_block",
    )
    assert fit["status"] == "fit"
    distance = coef.loc[coef["predictor"].eq("log_distance_to_continent_km")].iloc[0]
    assert distance["estimate"] > 0
    assert distance["n_islands"] == 60
    assert distance["n_clusters"] == 20


def test_status_models_keep_strata_and_support_class_separate():
    null = _synthetic_null()
    cov = _covariates()
    coef, support = mod.fit_status_stratified_lineage_models(
        null,
        cov,
        strata=("native_nonendemic",),
        regimes=("tropical",),
        pilot_min_islands=30,
        confirmatory_min_islands=50,
    )
    assert not coef.empty
    assert support.iloc[0]["support_class"] == "confirmatory_count_met"
    assert set(coef["stratum"]) == {"native_nonendemic"}
    assert set(coef["regime"]) == {"tropical"}


def test_30_to_49_islands_is_targeted_zone_even_if_model_can_fit():
    null = _synthetic_null(40)
    cov = _covariates(40)
    _, support = mod.fit_status_stratified_lineage_models(
        null,
        cov,
        strata=("native_nonendemic",),
        regimes=("tropical",),
        pilot_min_islands=30,
        confirmatory_min_islands=50,
    )
    assert support.iloc[0]["n_islands"] == 40
    assert support.iloc[0]["support_class"] == "targeted_acquisition_zone"


def test_missing_covariate_fails_closed():
    with pytest.raises(Exception):
        mod.fit_status_stratified_lineage_models(
            _synthetic_null(),
            _covariates().drop(columns=["climate_pc4"]),
            strata=("native_nonendemic",),
            regimes=("tropical",),
        )
