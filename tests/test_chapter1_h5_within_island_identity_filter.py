from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import typer
import yaml

from island_v2.chapter1_h5_within_island_identity_filter import prepare_rows, run_analysis


def _config() -> dict:
    return yaml.safe_load(
        open("config/chapter1_h5_within_island_identity_filter_v1.yml", encoding="utf-8")
    )


def _synthetic(n_per_context: int = 24):
    channels = ["bombus", "lepidoptera", "flower_visiting_birds"]
    axes = {
        "bombus": "large_bee_like",
        "lepidoptera": "butterfly_like",
        "flower_visiting_birds": "bird_like",
    }
    scores = []
    observations = []
    covariates = []
    rng = np.random.default_rng(17)
    counter = 0
    for context in ["northern_midlatitude", "tropical"]:
        for j in range(n_per_context):
            island = f"{context}_{j}"
            disrupted_channel = channels[j % len(channels)]
            covariates.append(
                {
                    "island_id": island,
                    "analysis_regime": context,
                    "spatial_block": f"block_{counter // 2}",
                }
            )
            counter += 1
            for channel in channels:
                disrupted = channel == disrupted_channel
                observations.append(
                    {
                        "island_id": island,
                        "channel_id": channel,
                        "observation_state": "adequate_non_detection" if disrupted else "detected",
                        "background_record_count": 80,
                        "background_spatial_units": 5,
                        "background_temporal_units": 4,
                        "distinct_dataset_count": 3,
                        "latest_background_year": 2024,
                    }
                )
            island_shift = rng.normal(0.0, 0.15)
            channel_baseline = {
                "bombus": 0.15,
                "lepidoptera": -0.05,
                "flower_visiting_birds": 0.05,
            }
            for channel, axis in axes.items():
                value = island_shift + channel_baseline[channel]
                if channel == disrupted_channel:
                    value -= 0.45
                value += rng.normal(0.0, 0.02)
                scores.append(
                    {
                        "island_id": island,
                        "stratum": "all_observed",
                        "syndrome": axis,
                        "syndrome_score": value,
                    }
                )
    return pd.DataFrame(scores), pd.DataFrame(observations), pd.DataFrame(covariates)


def test_within_island_matched_disruption_recovers_negative_effects() -> None:
    scores, observations, covariates = _synthetic()
    models, omnibus, audit, decision = run_analysis(
        scores,
        observations,
        covariates,
        _config(),
        evidence_scope="all_analysis_eligible",
    )
    matched = models.loc[models["mapping"].eq("matched")].set_index("context")
    assert matched.loc["northern_midlatitude", "estimate"] < -0.2
    assert matched.loc["tropical", "estimate"] < -0.2
    primary = omnibus.loc[omnibus["mapping"].eq("matched")].iloc[0]
    assert bool(primary["support_passed"])
    assert bool(primary["estimable"])
    assert primary["joint_df"] == 2
    assert primary["joint_p_value"] < 0.05
    assert decision["n_eligible_mixed_islands"] == 48
    assert decision["estimability_gate_passed"]
    assert audit["eligible_mixed_island"].sum() == 48


def test_support_gate_and_estimability_gate_are_separate() -> None:
    scores, observations, covariates = _synthetic(n_per_context=6)
    models, omnibus, _, decision = run_analysis(
        scores,
        observations,
        covariates,
        _config(),
        evidence_scope="all_analysis_eligible",
    )
    matched = omnibus.loc[omnibus["mapping"].eq("matched")].iloc[0]
    assert not bool(matched["support_passed"])
    assert bool(matched["estimable"])
    assert matched["target_rank_increment"] == 2
    assert pd.isna(matched["joint_p_value"])
    assert not decision["support_gate_passed"]
    assert decision["estimability_gate_passed"]
    assert models.loc[models["mapping"].eq("matched"), "estimate"].isna().all()


def test_estimability_gate_fails_when_disruption_is_confounded_with_channel_identity() -> None:
    scores, observations, covariates = _synthetic(n_per_context=12)
    observations = observations.copy()
    observations["observation_state"] = np.where(
        observations["channel_id"].eq("bombus"),
        "adequate_non_detection",
        "detected",
    )

    models, omnibus, _, decision = run_analysis(
        scores,
        observations,
        covariates,
        _config(),
        evidence_scope="all_analysis_eligible",
    )
    matched = omnibus.loc[omnibus["mapping"].eq("matched")].iloc[0]

    assert bool(matched["support_passed"])
    assert not bool(matched["estimable"])
    assert matched["target_rank_increment"] == 0
    assert pd.isna(matched["joint_p_value"])
    assert not decision["estimability_gate_passed"]
    assert models.loc[models["mapping"].eq("matched"), "estimate"].isna().all()


def test_prepare_rows_rejects_conflicting_island_context_or_block_assignments() -> None:
    scores, observations, covariates = _synthetic(n_per_context=10)
    duplicated = pd.concat(
        [
            covariates,
            pd.DataFrame(
                [
                    {
                        "island_id": "northern_midlatitude_0",
                        "analysis_regime": "tropical",
                        "spatial_block": "conflicting_block",
                    }
                ]
            ),
        ],
        ignore_index=True,
    )

    with pytest.raises(typer.BadParameter, match="conflicting covariate rows"):
        prepare_rows(scores, observations, duplicated, _config())
