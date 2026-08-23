import importlib

import numpy as np
import pandas as pd

mod = importlib.import_module("island_v2.status_category_coverage_sensitivity")


def _inputs(n=70):
    residuals = []
    status = []
    cov = []
    for i in range(n):
        island = f"i{i}"
        distance = 1 + i / 10
        total = 20
        trials = 8 + (i % 5)
        residuals.append(
            {
                "island_id": island,
                "stratum": "native_nonendemic",
                "domain": "colour",
                "category": "red_pink",
                "trials": trials,
                "deviation_observed_minus_genus": -0.02 * distance + 0.001 * (i % 3),
            }
        )
        for j in range(total):
            status.append(
                {
                    "island_id": island,
                    "accepted_species": f"S{i}_{j}",
                    "origin_status": "native",
                    "endemic_status": "nonendemic",
                }
            )
        cov.append(
            {
                "island_id": island,
                "analysis_regime": "tropical",
                "spatial_block": f"b{i // 3}",
                "log_distance_to_continent_km": distance,
                "log_island_area_km2": 2 + (i % 7) / 10,
                "climate_pc1": np.sin(i),
                "climate_pc2": np.cos(i),
                "climate_pc3": (i % 4) / 4,
                "climate_pc4": (i % 5) / 5,
            }
        )
    fdr = pd.DataFrame(
        {
            "stratum": ["native_nonendemic"],
            "regime": ["tropical"],
            "domain": ["colour"],
            "category": ["red_pink"],
            "fdr_supported": [True],
            "n_islands": [n],
        }
    )
    return pd.DataFrame(residuals), fdr, pd.DataFrame(status), pd.DataFrame(cov)


def test_supported_category_gets_continuous_and_threshold_sensitivities():
    residuals, fdr, status, cov = _inputs()
    out = mod.run_coverage_sensitivity(residuals, fdr, status, cov)
    assert {"baseline", "continuous_coverage_adjusted"}.issubset(set(out["model"]))
    baseline = out.loc[out["model"].eq("baseline")].iloc[0]
    adjusted = out.loc[out["model"].eq("continuous_coverage_adjusted")].iloc[0]
    assert baseline["distance_beta"] < 0
    assert adjusted["distance_beta"] < 0
    assert 0 < baseline["median_direct_fraction"] <= 1


def test_non_fdr_category_is_not_reanalysed():
    residuals, fdr, status, cov = _inputs()
    fdr["fdr_supported"] = False
    out = mod.run_coverage_sensitivity(residuals, fdr, status, cov)
    assert out.empty
