import numpy as np
import pandas as pd

from island_v2.chapter1_area_support_artifact import (
    add_capped_information_weights,
    common_support_area_overlap,
    heteroskedastic_null_simulation,
)


def test_capped_information_weights_are_positive_capped_and_mean_one() -> None:
    data = pd.DataFrame(
        {
            "response": ["a", "a", "a", "b", "b"],
            "n_response_species": [1, 25, 100, 10, 50],
        }
    )
    weighted = add_capped_information_weights(data, cap_species=50)
    assert weighted["analysis_weight"].gt(0).all()
    assert np.allclose(
        weighted.groupby("response")["analysis_weight"].mean().to_numpy(),
        1.0,
    )
    assert np.isclose(
        weighted.loc[2, "analysis_weight"],
        50 / np.mean([1, 25, 50]),
    )


def test_common_support_overlap_uses_support_and_area_not_response() -> None:
    rows = []
    for index in range(80):
        support = 10 if index % 2 == 0 else 50
        area = index / 10
        for response, score in [("a", index * 100.0), ("b", -index * 100.0)]:
            rows.append(
                {
                    "island_id": f"i{index}",
                    "response": response,
                    "response_score": score,
                    "n_response_species": support + (15 if response == "b" else 0),
                    "area": area,
                }
            )
    data = pd.DataFrame(rows)
    retained, diagnostic = common_support_area_overlap(
        data,
        area_column="area",
        minimum_per_half=30,
    )
    changed = data.copy()
    changed["response_score"] *= -999
    changed = changed.iloc[::-1].reset_index(drop=True)
    retained_changed, diagnostic_changed = common_support_area_overlap(
        changed,
        area_column="area",
        minimum_per_half=30,
    )
    assert diagnostic["status"] == "retained"
    assert diagnostic == diagnostic_changed
    assert set(retained["island_id"]) == set(retained_changed["island_id"])


def test_heteroskedastic_null_simulation_retains_true_null_as_artifact_risk() -> None:
    rng = np.random.default_rng(99)
    n = 160
    distance = np.linspace(-2, 2, n)
    area = np.sin(np.arange(n) / 11)
    support = np.where(np.arange(n) < n / 2, 8.0, 40.0)
    y = 0.4 * distance + 0.1 * area + rng.normal(0, 0.18, n) * np.sqrt(20 / support)
    families = pd.DataFrame(
        {
            "island_id": [f"i{x}" for x in range(n)],
            "stratum": "all_native",
            "source_mode": "geo_k5",
            "evidence_scope": "all_analysis_eligible",
            "family": "source_lineage_broad",
            "response": "entry_enrichment",
            "response_score": y,
            "n_response_species": support,
        }
    )
    covariates = pd.DataFrame(
        {
            "island_id": families["island_id"],
            "distance": distance,
            "area": area,
            "climate": np.cos(np.arange(n) / 7),
            "cluster": [f"b{x // 5}" for x in range(n)],
        }
    )
    realm = pd.DataFrame({"island_id": families["island_id"], "realm": "Palearctic"})
    area_config = {
        "distance_column": "distance",
        "area_column": "area",
        "cluster_column": "cluster",
        "control_columns": ["climate"],
        "context_layers": {"biogeographic_realm": {"column": "realm"}},
    }
    v3 = {
        "heteroskedastic_null_simulation": {"draws": 300, "seed": 20260828},
        "inference": {"simulation_false_positive_threshold": 0.05},
    }
    result = heteroskedastic_null_simulation(families, covariates, realm, area_config, v3)
    assert len(result) == 1
    assert result.loc[0, "status"] == "simulated"
    assert result.loc[0, "null_exceedance_probability"] > 0.05
    assert not bool(result.loc[0, "passes_false_positive_threshold"])
