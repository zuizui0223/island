from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from island_v2.chapter1_species_detection_tipping import (
    adjust_species_detection,
    build_recorded_stratum_counts,
    build_tipping_surface,
    completeness_from_distance,
    distance_standardization,
)


def _status_flora() -> pd.DataFrame:
    rows = []
    for island in ("near", "far"):
        for i in range(20):
            rows.append(
                {
                    "island_id": island,
                    "accepted_species": f"{island}_sp_{i}",
                    "origin_status": "native",
                    "endemic_status": "nonendemic",
                    "floristic_status": "native_nonendemic",
                }
            )
    return pd.DataFrame(rows)


def _scores() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "island_id": island,
                "stratum": stratum,
                "syndrome": syndrome,
                "syndrome_score": score,
                "n_species": n_species,
            }
            for island, score in (("near", -0.2), ("far", 0.2))
            for stratum in ("all_native", "native_nonendemic")
            for syndrome, n_species in (
                ("generalized_accessible", 10),
                ("selfing_core", 8),
            )
        ]
    )


def _covariates() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "island_id": ["near", "far"],
            "log_distance_to_continent_km": [1.0, 5.0],
        }
    )


def test_completeness_declines_with_distance_when_or_c_below_one() -> None:
    cov = _covariates()
    mean, sd = distance_standardization(cov, "log_distance_to_continent_km")
    c = completeness_from_distance(
        cov["log_distance_to_continent_km"].to_numpy(float),
        distance_mean=mean,
        distance_sd=sd,
        median_completeness=0.8,
        distance_completeness_odds_ratio=0.5,
    )
    assert c[0] > 0.8 > c[1]


def test_or_d_one_exactly_preserves_score_and_original_information() -> None:
    scores = _scores()
    counts = build_recorded_stratum_counts(
        _status_flora(), ["all_native", "native_nonendemic"]
    )
    cov = _covariates()
    mean, sd = distance_standardization(cov, "log_distance_to_continent_km")
    adjusted, diagnostic = adjust_species_detection(
        scores,
        counts,
        cov,
        geography="log_distance_to_continent_km",
        distance_mean=mean,
        distance_sd=sd,
        median_completeness=0.5,
        distance_completeness_odds_ratio=0.25,
        state_recording_odds_ratio=1.0,
        affected_syndrome="generalized_accessible",
    )
    merged = scores.merge(
        adjusted,
        on=["island_id", "stratum", "syndrome"],
        suffixes=("_before", "_after"),
    )
    assert np.max(
        np.abs(
            merged["syndrome_score_before"].to_numpy(float)
            - merged["syndrome_score_after"].to_numpy(float)
        )
    ) <= 1e-12
    assert (merged["n_species_before"] == merged["n_species_after"]).all()
    assert diagnostic["assumed_unrecorded_species"].gt(0).all()


def test_underrecorded_positive_state_correction_raises_completed_prevalence() -> None:
    scores = _scores()
    counts = build_recorded_stratum_counts(
        _status_flora(), ["all_native", "native_nonendemic"]
    )
    cov = _covariates()
    mean, sd = distance_standardization(cov, "log_distance_to_continent_km")
    _, diagnostic = adjust_species_detection(
        scores,
        counts,
        cov,
        geography="log_distance_to_continent_km",
        distance_mean=mean,
        distance_sd=sd,
        median_completeness=0.8,
        distance_completeness_odds_ratio=0.5,
        state_recording_odds_ratio=0.5,
        affected_syndrome="generalized_accessible",
    )
    assert (
        diagnostic["completed_soft_membership"]
        >= diagnostic["observed_soft_membership"] - 1e-12
    ).all()
    far = diagnostic.loc[diagnostic["island_id"].eq("far")].iloc[0]
    near = diagnostic.loc[diagnostic["island_id"].eq("near")].iloc[0]
    assert far["assumed_list_completeness"] < near["assumed_list_completeness"]


def test_tipping_surface_distinguishes_baseline_not_supported_and_break() -> None:
    rows = []
    for target, baseline_supported in (
        ("Palearctic_accessibility", True),
        ("north_tropical_vector_difference", False),
    ):
        for or_d, supported, sign_ok in (
            (1.0, baseline_supported, True),
            (0.8, baseline_supported, True),
            (0.5, False, True),
        ):
            rows.append(
                {
                    "evidence_scope": "all_analysis_eligible",
                    "target": target,
                    "stratum": "all_native",
                    "median_distance_completeness": 0.8,
                    "distance_completeness_odds_ratio": 0.5,
                    "state_recording_odds_ratio": or_d,
                    "supported": supported,
                    "sign_ok": sign_ok,
                }
            )
    headline = pd.DataFrame(rows)
    config = {
        "primary_bias_direction": {
            "positive_state_recording_odds_ratio_grid": [1.0, 0.8, 0.5]
        }
    }
    tipping = build_tipping_surface(headline, config)
    pal = tipping.loc[tipping["target"].eq("Palearctic_accessibility")].iloc[0]
    assert pal["verdict"] == "break_detected"
    assert pal["tipping_state_recording_odds_ratio"] == pytest.approx(0.5)
    assert pal["recording_advantage_nonfocal_over_focal"] == "2"
    vector = tipping.loc[
        tipping["target"].eq("north_tropical_vector_difference")
    ].iloc[0]
    assert vector["verdict"] == "not_applicable"
    assert vector["tipping_event"] == "baseline_not_supported"
