import numpy as np
import pandas as pd

from island_v2.chapter1_pr136_source_exposure import run_source_exposure_test


def _config(confirmatory=50):
    return {
        "geography_column": "log_distance_to_continent_km",
        "cluster_column": "spatial_block",
        "baseline_covariates": [
            "log_island_area_km2",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
        ],
        "strata": ["all_native"],
        "support_tiers": {"pilot": 30, "confirmatory": confirmatory},
        "minimum_outcomes_per_vector": 2,
        "alpha": 0.05,
    }


def _synthetic(differential=True, n_each=70):
    rng = np.random.default_rng(136)
    genus_rows = []
    app_rows = []
    cov_rows = []
    outcomes = {
        "generalized_form": 0.90,
        "self_compatibility": 0.70,
        "plain_colour": 0.15,
    }
    island_index = 0
    for exposed in (0, 1):
        for j in range(n_each):
            island_id = f"i{island_index}"
            distance = -1.6 + 3.2 * j / (n_each - 1)
            app_rows.append(
                {
                    "island_id": island_id,
                    "applicability": "applicable"
                    if exposed
                    else "structurally_not_applicable",
                }
            )
            cov_rows.append(
                {
                    "island_id": island_id,
                    "log_distance_to_continent_km": distance,
                    "log_island_area_km2": rng.normal(),
                    "climate_pc1": rng.normal(),
                    "climate_pc2": rng.normal(),
                    "climate_pc3": rng.normal(),
                    "climate_pc4": rng.normal(),
                    "spatial_block": f"b{island_index // 4}",
                }
            )
            for outcome, exposed_extra in outcomes.items():
                shared_slope = 0.15
                slope = shared_slope + (exposed * exposed_extra if differential else 0.0)
                value = slope * distance + rng.normal(0, 0.20)
                genus_rows.append(
                    {
                        "island_id": island_id,
                        "outcome": outcome,
                        "stratum": "all_native",
                        "observed_n_species": 40,
                        "deviation_observed_minus_null": value,
                    }
                )
            island_index += 1
    return pd.DataFrame(genus_rows), pd.DataFrame(app_rows), pd.DataFrame(cov_rows)


def test_source_exposure_joint_vector_recovers_pr136_contingency():
    genus_null, applicability, covariates = _synthetic(differential=True)
    coefficients, support, joint, slopes = run_source_exposure_test(
        genus_null, applicability, covariates, _config()
    )

    assert not coefficients.empty
    assert support.loc[
        support["support_tier"].eq("confirmatory"), "retained_in_joint_vector"
    ].all()

    row = joint.loc[joint["support_tier"].eq("confirmatory")].iloc[0]
    assert row["status"] == "fit"
    assert row["n_retained_outcomes"] == 3
    assert row["joint_exposure_interaction_p"] < 0.05
    assert bool(row["source_exposure_contingency_supported"])

    wide = slopes.loc[slopes["support_tier"].eq("confirmatory")].pivot(
        index="outcome", columns="channel_regime", values="geography_slope"
    )
    assert wide.loc["generalized_form", "source_exposed"] > wide.loc[
        "generalized_form", "structurally_absent"
    ]
    assert wide.loc["self_compatibility", "source_exposed"] > wide.loc[
        "self_compatibility", "structurally_absent"
    ]


def test_same_slopes_do_not_manufacture_source_exposure_contingency():
    genus_null, applicability, covariates = _synthetic(differential=False)
    _, _, joint, _ = run_source_exposure_test(
        genus_null, applicability, covariates, _config()
    )
    row = joint.loc[joint["support_tier"].eq("confirmatory")].iloc[0]
    assert row["status"] == "fit"
    assert row["joint_exposure_interaction_p"] > 0.01


def test_source_exposure_test_fails_closed_when_structural_group_is_under_supported():
    genus_null, applicability, covariates = _synthetic(differential=True, n_each=35)
    _, support, joint, _ = run_source_exposure_test(
        genus_null, applicability, covariates, _config(confirmatory=50)
    )
    row = joint.loc[joint["support_tier"].eq("confirmatory")].iloc[0]
    assert row["status"] == "not_testable"
    assert row["n_retained_outcomes"] == 0
    confirmatory = support.loc[support["support_tier"].eq("confirmatory")]
    assert not confirmatory["retained_in_joint_vector"].any()


def test_unresolved_applicability_is_never_recoded_as_structural_absence():
    genus_null, applicability, covariates = _synthetic(differential=True)
    unresolved_ids = applicability.loc[applicability["applicability"].eq("applicable"), "island_id"].head(10)
    applicability.loc[applicability["island_id"].isin(unresolved_ids), "applicability"] = "unresolved"
    _, support, _, _ = run_source_exposure_test(
        genus_null, applicability, covariates, _config()
    )
    row = support.loc[
        support["support_tier"].eq("confirmatory")
        & support["outcome"].eq("generalized_form")
    ].iloc[0]
    assert row["n_source_exposed_islands"] == 60
    assert row["n_structurally_absent_islands"] == 70
