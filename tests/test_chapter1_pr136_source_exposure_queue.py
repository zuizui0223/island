import pandas as pd

from island_v2.chapter1_pr136_source_exposure_queue import (
    build_source_exposure_queue,
)


def _inputs():
    genus_rows = []
    app_rows = []
    cov_rows = []
    for i in range(48):
        island_id = f"i{i:02d}"
        context = "northern_midlatitude" if i < 24 else "tropical"
        app_rows.append(
            {
                "island_id": island_id,
                "applicability": "applicable" if i == 0 else "unresolved",
            }
        )
        cov_rows.append(
            {
                "island_id": island_id,
                "analysis_regime": context,
                "log_distance_to_continent_km": float(i % 24),
            }
        )
        for stratum in ["all_native", "native_nonendemic"]:
            for outcome_index, outcome in enumerate(
                ["plain_colour", "generalized_form", "self_compatibility"]
            ):
                genus_rows.append(
                    {
                        "island_id": island_id,
                        "stratum": stratum,
                        "outcome": outcome,
                        "observed_n_species": 10 + outcome_index + (i % 5),
                        "deviation_observed_minus_null": (-1) ** i * (i + outcome_index),
                    }
                )
    return pd.DataFrame(genus_rows), pd.DataFrame(app_rows), pd.DataFrame(cov_rows)


def _config():
    return {
        "priority_strata": ["all_native", "native_nonendemic"],
        "minimum_outcomes_per_stratum": 2,
        "context_column": "analysis_regime",
        "geography_column": "log_distance_to_continent_km",
        "distance_bins": 4,
        "wave_size_per_context_distance_bin": 2,
    }


def test_queue_is_outcome_blind_to_residual_values_and_signs():
    genus_null, applicability, covariates = _inputs()
    queue_a, manifest_a = build_source_exposure_queue(
        genus_null, applicability, covariates, _config()
    )

    altered = genus_null.copy()
    altered["deviation_observed_minus_null"] = (
        -1000 * altered["deviation_observed_minus_null"] + 77
    )
    queue_b, manifest_b = build_source_exposure_queue(
        altered, applicability, covariates, _config()
    )

    pd.testing.assert_frame_equal(queue_a, queue_b)
    assert manifest_a == manifest_b
    assert manifest_a["outcome_blind"] is True


def test_resolved_islands_are_not_put_back_into_curation_queue():
    genus_null, applicability, covariates = _inputs()
    queue, _ = build_source_exposure_queue(
        genus_null, applicability, covariates, _config()
    )
    assert "i00" not in set(queue["island_id"])
    assert queue["applicability"].eq("unresolved").all()


def test_wave_one_balances_context_and_distance_cells():
    genus_null, applicability, covariates = _inputs()
    queue, manifest = build_source_exposure_queue(
        genus_null, applicability, covariates, _config()
    )
    wave1 = queue.loc[queue["wave1_priority"]]
    cell_counts = wave1.groupby(
        ["analysis_regime", "distance_support_bin"]
    ).size()
    assert (cell_counts <= 2).all()
    assert set(manifest["eligible_by_context"]) == {
        "northern_midlatitude",
        "tropical",
    }
    assert manifest["n_wave1_islands"] == len(wave1)


def test_queue_requires_support_in_both_priority_strata():
    genus_null, applicability, covariates = _inputs()
    genus_null = genus_null.loc[
        ~(
            genus_null["island_id"].eq("i01")
            & genus_null["stratum"].eq("native_nonendemic")
            & genus_null["outcome"].isin(["generalized_form", "self_compatibility"])
        )
    ]
    queue, _ = build_source_exposure_queue(
        genus_null, applicability, covariates, _config()
    )
    assert "i01" not in set(queue["island_id"])
