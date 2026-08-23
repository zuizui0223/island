import importlib

import numpy as np
import pandas as pd

mod = importlib.import_module("island_v2.status_category_decomposition")


def test_direct_category_mapping_excludes_cross_category_values():
    evidence = pd.DataFrame(
        {
            "accepted_species": ["A a", "B b", "C c"],
            "trait_name": ["flower_primary_color"] * 3,
            "normalized_value": ["white", "red_pink", "white|red_pink"],
            "evidence_scope": ["species_direct"] * 3,
            "resolution_status": ["resolved"] * 3,
        }
    )
    master = pd.DataFrame(
        {"accepted_species": ["A a", "B b", "C c"], "genus": ["G", "G", "H"]}
    )
    out = mod.direct_categorical_species(evidence, master, mod.CATEGORY_SPECS["colour"])
    assert set(out["accepted_species"]) == {"A a", "B b"}
    assert out.set_index("accepted_species").loc["A a", "category"] == "plain"


def test_genus_expected_category_share_is_exact():
    status = pd.DataFrame(
        {
            "island_id": ["i1", "i2"],
            "accepted_species": ["A a", "B b"],
            "origin_status": ["native", "native"],
            "endemic_status": ["nonendemic", "nonendemic"],
        }
    )
    evidence = pd.DataFrame(
        {
            "accepted_species": ["A a", "B b"],
            "trait_name": ["flower_primary_color", "flower_primary_color"],
            "normalized_value": ["white", "red_pink"],
            "evidence_scope": ["species_direct", "species_direct"],
            "resolution_status": ["resolved", "resolved"],
        }
    )
    master = pd.DataFrame(
        {"accepted_species": ["A a", "B b"], "genus": ["G", "G"]}
    )
    residuals, _ = mod.build_category_residuals(
        status,
        evidence,
        master,
        specs={"colour": mod.CATEGORY_SPECS["colour"]},
        strata=("native_nonendemic",),
    )
    plain = residuals.loc[residuals["category"].eq("plain")].set_index("island_id")
    assert plain.loc["i1", "expected_genus_share"] == 0.5
    assert plain.loc["i1", "deviation_observed_minus_genus"] == 0.5
    assert plain.loc["i2", "deviation_observed_minus_genus"] == -0.5


def test_category_isolation_model_keeps_categories_separate():
    residual_rows = []
    cov_rows = []
    for i in range(70):
        distance = 1 + i / 10
        for category, sign in (("plain", 1.0), ("red_pink", -1.0)):
            residual_rows.append(
                {
                    "island_id": f"i{i}",
                    "stratum": "native_nonendemic",
                    "domain": "colour",
                    "category": category,
                    "trials": 20,
                    "deviation_observed_minus_genus": sign * 0.02 * distance,
                }
            )
        cov_rows.append(
            {
                "island_id": f"i{i}",
                "analysis_regime": "northern_midlatitude",
                "spatial_block": f"b{i // 3}",
                "log_distance_to_continent_km": distance,
                "log_island_area_km2": 2 + (i % 7) / 10,
                "climate_pc1": np.sin(i),
                "climate_pc2": np.cos(i),
                "climate_pc3": (i % 4) / 4,
                "climate_pc4": (i % 5) / 5,
            }
        )
    coef, support = mod.fit_category_isolation_models(
        pd.DataFrame(residual_rows), pd.DataFrame(cov_rows), regimes=("northern_midlatitude",)
    )
    isolation = coef.loc[coef["predictor"].eq("log_distance_to_continent_km")]
    estimates = isolation.set_index("category")["estimate"]
    assert estimates["plain"] > 0
    assert estimates["red_pink"] < 0
    assert set(support["category"]) == {"plain", "red_pink"}
