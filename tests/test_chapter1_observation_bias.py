import math

import numpy as np
import pandas as pd

from island_v2.chapter1_observation_bias import run_observation_bias_audit


def test_observation_bias_keeps_zero_record_islands_and_separates_stages():
    rng = np.random.default_rng(41)
    islands = []
    species_rows = []
    comp_rows = []
    for i in range(160):
        context = "northern_midlatitude" if i < 80 else "tropical"
        distance = float(i % 80) * 10.0
        area = 1.5 + (i % 20) / 10
        islands.append(
            {
                "island_id": f"i{i}",
                "analysis_regime": context,
                "spatial_block": f"b{i // 4}",
                "log_distance_to_continent_km": math.log1p(distance),
                "log_island_area_km2": area,
                "climate_pc1": rng.normal(),
                "climate_pc2": rng.normal(),
                "climate_pc3": rng.normal(),
                "climate_pc4": rng.normal(),
            }
        )
        # Deliberately make small/far islands less observable while retaining them.
        flora = (i % 5 != 0) and (area > 1.7 or distance < 300)
        if flora:
            n_species = 4 + (i % 7)
            for s in range(n_species):
                species_rows.append({"island_id": f"i{i}", "species": f"sp{s}"})
            direct = i % 3 != 0
            if direct:
                for trait in ["flower_primary_color", "floral_form", "self_incompatibility"]:
                    comp_rows.append(
                        {
                            "island_id": f"i{i}",
                            "trait_name": trait,
                            "trait_state": "state",
                            "trials": 5,
                            "successes": 2,
                            "share": 0.4,
                            "stratum": "all_native",
                        }
                    )

    config = {
        "baseline_covariates": [
            "log_island_area_km2",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
        ],
        "distance_gradient": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "reference_context": "northern_midlatitude",
        "contexts": ["northern_midlatitude", "tropical"],
        "trait_domains": {
            "flower_primary_color": "flower_primary_color",
            "floral_form": "floral_form",
            "self_incompatibility": "self_incompatibility",
        },
        "selection_models": [
            "flora_recorded",
            "direct_trait_any",
            "direct_trait_any_given_flora",
            "direct_trait_all_three_given_flora",
        ],
    }
    support, context, coef, fit = run_observation_bias_audit(
        pd.DataFrame(islands),
        pd.DataFrame(species_rows),
        pd.DataFrame(comp_rows),
        config,
    )
    assert support["island_id"].nunique() == 160
    assert (~support["flora_recorded"]).any()
    assert support.loc[~support["flora_recorded"], "observation_stage"].eq("no_flora_record").all()
    assert support.loc[~support["flora_recorded"], "direct_trait_any"].eq(False).all()
    assert context["n_all_islands"].sum() == 160
    assert set(fit["endpoint"]) == set(config["selection_models"])
    assert not coef.empty
    assert not coef["predictor"].str.contains("bombus", case=False).any()
