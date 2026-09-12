from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import yaml

from island_v2.chapter1_nee_n1_run import (
    prepare_canonical_covariates,
    run_n1,
    write_outputs,
)


def _config() -> dict:
    return yaml.safe_load(open("config/chapter1_nee_n1_model.yml", encoding="utf-8"))


def _raw_covariates(n: int) -> pd.DataFrame:
    rng = np.random.default_rng(101)
    log_distance = np.linspace(0.0, 6.0, n)
    log_area = np.linspace(1.0, 5.0, n)
    return pd.DataFrame(
        {
            "island_id": [f"isl-{i:03d}" for i in range(n)],
            "distance_to_continent_km": np.expm1(log_distance),
            "log_distance_to_continent_km": log_distance,
            "area_km2": np.expm1(log_area),
            "log_island_area_km2": log_area,
            "climate_pc1": rng.normal(size=n),
            "climate_pc2": rng.normal(size=n),
            "climate_pc3": rng.normal(size=n),
            "climate_pc4": rng.normal(size=n),
            "spatial_block": [f"block-{i % 12:02d}" for i in range(n)],
        }
    )


def _qualification(channels: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "channel_id": channels,
            "support_tier": ["confirmatory"] * len(channels),
            "N1_gate_eligible": [True] * len(channels),
        }
    )


def _projected(channels: list[str], n: int) -> pd.DataFrame:
    rng = np.random.default_rng(202)
    x = np.linspace(-2.5, 2.5, n)
    slopes = {
        "bombus": -1.8,
        "non_bombus_bees": -0.5,
        "lepidoptera": 0.7,
        "flower_visiting_birds": 1.4,
        "diptera": 0.1,
    }
    rows: list[dict[str, object]] = []
    for i in range(n):
        for channel in channels:
            p = 1.0 / (1.0 + np.exp(-slopes[channel] * x[i]))
            retained = bool(rng.random() < p)
            rows.append(
                {
                    "island_id": f"isl-{i:03d}",
                    "source_region_id": f"source-{i % 3}",
                    "channel_id": channel,
                    "source_state": "available",
                    "channel_state": "retained" if retained else "disrupted",
                    "background_record_count": 100 + (i % 30),
                    "background_spatial_units": 5 + (i % 5),
                    "background_temporal_units": 4 + (i % 4),
                    "distinct_dataset_count": 3 + (i % 3),
                    "latest_background_year": 2025 - (i % 4),
                }
            )
    return pd.DataFrame(rows)


def test_canonical_covariate_adapter_verifies_and_renames_frozen_logs() -> None:
    raw = _raw_covariates(12)
    result = prepare_canonical_covariates(raw, expected_islands=12)
    assert list(result.columns) == [
        "island_id",
        "log1p_distance_to_continent_km",
        "log_area",
        "climate_pc1",
        "climate_pc2",
        "climate_pc3",
        "climate_pc4",
        "spatial_block",
    ]
    assert np.allclose(
        result["log1p_distance_to_continent_km"], raw["log_distance_to_continent_km"]
    )
    assert np.allclose(result["log_area"], raw["log_island_area_km2"])


def test_canonical_covariate_adapter_hard_stops_changed_log_definition() -> None:
    raw = _raw_covariates(8)
    raw.loc[3, "log_distance_to_continent_km"] += 0.01
    with pytest.raises(ValueError, match="not log1p"):
        prepare_canonical_covariates(raw, expected_islands=8)


def test_fewer_than_three_confirmatory_channels_stops_without_fit() -> None:
    channels = ["bombus", "lepidoptera"]
    result = run_n1(
        _projected(channels, 40),
        _qualification(channels),
        _raw_covariates(40),
        _config(),
        expected_islands=40,
    )
    assert result["frame"].empty
    assert result["gate"]["N1_pass"] is False
    assert result["gate"]["status"] == "N1_not_evaluable_for_NEE_gate"
    assert result["gate"]["failure_action"] == "stop_before_N2_and_keep_frozen_Chapter1"
    assert result["coefficients"].empty


def test_qualified_runner_executes_both_required_deletion_robustness_checks(tmp_path) -> None:
    channels = ["bombus", "non_bombus_bees", "lepidoptera", "flower_visiting_birds"]
    result = run_n1(
        _projected(channels, 90),
        _qualification(channels),
        _raw_covariates(90),
        _config(),
        expected_islands=90,
    )
    assert result["metadata"]["status"] == "qualified"
    assert not result["coefficients"].empty
    assert set(result["source_deletions"]["deleted"]) == {"source-0", "source-1", "source-2"}
    assert set(result["block_deletions"]["deleted"]) == {
        f"block-{i:02d}" for i in range(12)
    }
    assert result["gate"]["failure_action"] in {
        "N2_may_open",
        "stop_before_N2_and_keep_frozen_Chapter1",
    }

    write_outputs(result, tmp_path)
    expected = {
        "N1_model_support.csv",
        "N1_coefficients.csv",
        "N1_channel_slopes.csv",
        "N1_global_heterogeneity_test.json",
        "N1_model_comparison.json",
        "N1_leave_one_spatial_block_out.csv",
        "N1_leave_one_source_region_out.csv",
        "N1_gate_receipt.json",
        "N1_run_metadata.json",
    }
    assert expected.issubset({path.name for path in tmp_path.iterdir()})
