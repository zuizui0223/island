import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from island_v2.chapter1_context_analysis import _load_config, fit_chapter1_context_models
from island_v2.chapter1_context_input import build_island_trait_composition
from island_v2.flora_status_support import attach_floristic_status
from island_v2.genus_fixed_trait_null import run_genus_fixed_null
from island_v2.status_stratified_lineage_analysis import fit_status_stratified_lineage_models


def test_status_conflict_fails_closed_to_unresolved():
    flora = pd.DataFrame({"island_id": ["i1"], "species": ["A a"]})
    status = pd.DataFrame(
        {
            "island_id": ["i1", "i1"],
            "accepted_species": ["A a", "A a"],
            "origin_status": ["native", "introduced"],
            "endemic_status": ["nonendemic", "nonendemic"],
        }
    )
    out = attach_floristic_status(flora, status)
    assert out.iloc[0]["origin_status"] == "unresolved"
    assert bool(out.iloc[0]["status_conflict"])


def test_context_input_excludes_conflicting_direct_trait_states():
    flora = pd.DataFrame(
        {"island_id": ["i1", "i1", "i2"], "species": ["A a", "A b", "A a"]}
    )
    status = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i2"],
            "accepted_species": ["A a", "A b", "A a"],
            "origin_status": ["native"] * 3,
            "endemic_status": ["nonendemic"] * 3,
        }
    )
    evidence = pd.DataFrame(
        {
            "accepted_species": ["A a", "A a", "A b"],
            "trait_name": ["flower_primary_color"] * 3,
            "normalized_value": ["red_pink", "white", "red_pink"],
            "evidence_scope": ["species_direct"] * 3,
            "resolution_status": ["resolved"] * 3,
        }
    )
    composition, _, audit = build_island_trait_composition(flora, status, evidence)
    assert not audit.loc[audit["accepted_species"].eq("A a"), "resolved_for_primary"].iloc[0]
    assert set(composition["trait_state"]) == {"red_pink"}
    assert set(composition["island_id"]) == {"i1"}


def _synthetic_context_data():
    rng = np.random.default_rng(5)
    contexts = ["northern_midlatitude", "tropical", "southern_extratropical"]
    rows = []
    cov = []
    idx = 0
    for context in contexts:
        for j in range(30):
            island_id = f"i{idx}"
            isolation = -1.5 + 3.0 * j / 29
            if context == "northern_midlatitude":
                slope = 1.2
            elif context == "tropical":
                slope = -1.1
            else:
                slope = 0.1
            # Give contexts different intercepts as well, so the test ensures
            # M1->M2 is not merely detecting regional mean composition.
            intercept = {
                "northern_midlatitude": -0.8,
                "tropical": 0.5,
                "southern_extratropical": 1.0,
            }[context]
            eta = intercept + slope * isolation + 0.08 * math.sin(j)
            p = 1.0 / (1.0 + math.exp(-eta))
            trials = 120
            successes = int(rng.binomial(trials, p))
            rows.append(
                {
                    "island_id": island_id,
                    "stratum": "all_native",
                    "trait_name": "floral_form",
                    "trait_state": "open_radial",
                    "successes": successes,
                    "trials": trials,
                }
            )
            cov.append(
                {
                    "island_id": island_id,
                    "analysis_regime": context,
                    "spatial_block": f"b{idx // 3}",
                    "log_distance_to_continent_km": isolation,
                    "log_island_area_km2": 2.0 + rng.normal(0, 0.3),
                    "climate_pc1": rng.normal(),
                    "climate_pc2": rng.normal(),
                    "climate_pc3": rng.normal(),
                    "climate_pc4": rng.normal(),
                }
            )
            idx += 1
    return pd.DataFrame(rows), pd.DataFrame(cov)


def test_m2_recovers_context_dependent_isolation_direction_with_clean_nesting():
    composition, covariates = _synthetic_context_data()
    config = {
        "baseline_covariates": [
            "log_island_area_km2",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
        ],
        "isolation_column": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "reference_context": "northern_midlatitude",
        "cluster_column": "spatial_block",
        "min_islands_per_fit": 30,
        "min_islands_per_context": 30,
    }
    coefficients, fits, slopes, support = fit_chapter1_context_models(
        composition, covariates, config
    )
    assert not coefficients.empty
    assert set(fits["model"]) == {"M0", "M1", "M2"}
    assert support.iloc[0]["status"] == "fit"

    m0 = set(coefficients.loc[coefficients["model"].eq("M0"), "predictor"])
    m1 = set(coefficients.loc[coefficients["model"].eq("M1"), "predictor"])
    m2 = set(coefficients.loc[coefficients["model"].eq("M2"), "predictor"])
    assert any(x.startswith("context[") for x in m0)
    assert "z_log_distance_to_continent_km" not in m0
    assert "z_log_distance_to_continent_km" in m1
    assert not any(x.startswith("z_isolation:context[") for x in m1)
    assert any(x.startswith("z_isolation:context[") for x in m2)

    north = slopes.loc[slopes["context"].eq("northern_midlatitude")].iloc[0]
    tropical = slopes.loc[slopes["context"].eq("tropical")].iloc[0]
    assert north["isolation_slope_log_odds_per_sd"] > 0
    assert tropical["isolation_slope_log_odds_per_sd"] < 0
    assert set(slopes["support_class"]) == {"pilot_count_met"}
    assert not coefficients["predictor"].str.contains("bombus", case=False).any()
    fit_one = fits.set_index("model")
    assert fit_one.loc["M2", "interaction_improvement_m1_to_m2"] > 0


def test_bombus_predictor_is_rejected_from_canonical_config(tmp_path: Path):
    path = tmp_path / "bad.yml"
    path.write_text(
        "\n".join(
            [
                "baseline_covariates: [bombus_deficit]",
                "isolation_column: log_distance_to_continent_km",
                "context_column: analysis_regime",
                "cluster_column: spatial_block",
                "reference_context: northern_midlatitude",
            ]
        ),
        encoding="utf-8",
    )
    with pytest.raises(Exception):
        _load_config(path)


def test_genus_fixed_null_retains_status_strata():
    flora = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i2", "i2", "i3", "i3"],
            "species": ["A a", "A b", "A a", "A c", "A b", "A c"],
        }
    )
    status = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i2", "i2", "i3", "i3"],
            "accepted_species": ["A a", "A b", "A a", "A c", "A b", "A c"],
            "origin_status": ["native"] * 6,
            "endemic_status": ["nonendemic"] * 6,
        }
    )
    evidence = pd.DataFrame(
        {
            "accepted_species": ["A a", "A b", "A c"],
            "trait_name": ["self_incompatibility"] * 3,
            "normalized_value": ["SC", "SI", "SI"],
        }
    )
    taxa = pd.DataFrame(
        {"accepted_species": ["A a", "A b", "A c"], "genus": ["A", "A", "A"]}
    )
    outcomes = {
        "self_compatibility": {
            "trait_name": "self_incompatibility",
            "positive_states": ["SC"],
            "negative_states": ["SI"],
        }
    }
    result, audit = run_genus_fixed_null(
        flora, status, evidence, taxa, outcomes, n_draws=50, seed=3, min_species_per_island=1
    )
    assert len(result.loc[result["stratum"].eq("all_native")]) == 3
    assert audit.iloc[0]["n_permutable_genera"] == 1


def test_status_lineage_layer_recovers_positive_residual_slope():
    n = 60
    null = pd.DataFrame(
        {
            "island_id": [f"i{i}" for i in range(n)],
            "outcome": ["plain_colour"] * n,
            "stratum": ["native_nonendemic"] * n,
            "observed_n_species": [10 + i % 4 for i in range(n)],
            "deviation_observed_minus_null": [
                0.35 * (1 + i / 10) + (i % 3) * 0.01 for i in range(n)
            ],
        }
    )
    cov = pd.DataFrame(
        {
            "island_id": [f"i{i}" for i in range(n)],
            "analysis_regime": ["tropical"] * n,
            "spatial_block": [f"b{i // 3}" for i in range(n)],
            "log_distance_to_continent_km": [1 + i / 10 for i in range(n)],
            "log_island_area_km2": [2 + (i % 7) / 10 for i in range(n)],
            "climate_pc1": [np.sin(i) for i in range(n)],
            "climate_pc2": [np.cos(i) for i in range(n)],
            "climate_pc3": [(i % 4) / 4 for i in range(n)],
            "climate_pc4": [(i % 5) / 5 for i in range(n)],
        }
    )
    coef, support = fit_status_stratified_lineage_models(
        null,
        cov,
        strata=("native_nonendemic",),
        pilot_min_islands=30,
        confirmatory_min_islands=50,
    )
    distance = coef.loc[coef["predictor"].eq("log_distance_to_continent_km")].iloc[0]
    assert distance["estimate"] > 0
    assert support.iloc[0]["support_class"] == "confirmatory_count_met"
