from __future__ import annotations

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_nee_n1_model import (
    build_primary_frame,
    deletion_robustness,
    fit_n1_primary,
    gate_receipt,
    global_heterogeneity_test,
    qualified_channels,
)


def _config() -> dict:
    return yaml.safe_load(open("config/chapter1_nee_n1_model.yml", encoding="utf-8"))


def _qualification(channels: list[str], tier: str = "confirmatory") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "channel_id": channels,
            "support_tier": [tier] * len(channels),
            "N1_gate_eligible": [tier == "confirmatory"] * len(channels),
        }
    )


def _covariates(n: int) -> pd.DataFrame:
    rng = np.random.default_rng(123)
    return pd.DataFrame(
        {
            "island_id": [f"isl-{i:03d}" for i in range(n)],
            "log1p_distance_to_continent_km": np.linspace(0.1, 6.0, n),
            "log_area": rng.normal(4, 1, n),
            "climate_pc1": rng.normal(size=n),
            "climate_pc2": rng.normal(size=n),
            "climate_pc3": rng.normal(size=n),
            "climate_pc4": rng.normal(size=n),
            "spatial_block": [f"block-{i % 18:02d}" for i in range(n)],
        }
    )


def _projected(channels: list[str], n: int, heterogeneous: bool = False) -> pd.DataFrame:
    rng = np.random.default_rng(456)
    distance = np.linspace(-2.5, 2.5, n)
    slopes = {
        "bombus": -2.2,
        "non_bombus_bees": -0.6,
        "lepidoptera": 0.9,
        "flower_visiting_birds": 1.8,
        "diptera": 0.2,
    }
    rows = []
    for i in range(n):
        for channel in channels:
            slope = slopes[channel] if heterogeneous else -0.7
            eta = -0.1 + slope * distance[i]
            p = 1.0 / (1.0 + np.exp(-eta))
            retained = bool(rng.random() < p)
            rows.append(
                {
                    "island_id": f"isl-{i:03d}",
                    "source_region_id": f"source-{i % 3}",
                    "channel_id": channel,
                    "source_state": "available",
                    "channel_state": "retained" if retained else "disrupted",
                    "background_record_count": 120 + (i % 20),
                    "background_spatial_units": 6 + (i % 4),
                    "background_temporal_units": 5 + (i % 3),
                    "distinct_dataset_count": 4 + (i % 2),
                    "latest_background_year": 2024 - (i % 3),
                }
            )
    return pd.DataFrame(rows)


def test_qualified_channels_require_confirmatory_gate() -> None:
    config = _config()
    q = pd.DataFrame(
        {
            "channel_id": ["bombus", "lepidoptera", "diptera"],
            "support_tier": ["confirmatory", "pilot", "confirmatory"],
            "N1_gate_eligible": [True, False, True],
        }
    )
    assert qualified_channels(q, config) == ["bombus", "diptera"]


def test_primary_frame_excludes_structural_and_unresolved_rows() -> None:
    config = _config()
    channels = ["bombus", "non_bombus_bees", "lepidoptera"]
    projected = _projected(channels, 60)
    projected.loc[0, ["source_state", "channel_state"]] = ["structurally_absent", "structurally_absent"]
    projected.loc[1, "channel_state"] = "unresolved"
    frame, metadata = build_primary_frame(
        projected, _qualification(channels), _covariates(60), config
    )
    assert metadata["status"] == "qualified"
    assert len(frame) == len(projected) - 2
    assert set(frame["channel_state"]) == {"retained", "disrupted"}
    assert set(frame["source_state"]) == {"available"}


def test_fewer_than_three_channels_hard_stops_primary_N1() -> None:
    config = _config()
    channels = ["bombus", "lepidoptera"]
    frame, metadata = build_primary_frame(
        _projected(channels, 80), _qualification(channels), _covariates(80), config
    )
    assert frame.empty
    assert metadata["status"] == "N1_not_evaluable_for_NEE_gate"
    assert metadata["required_confirmatory_channels"] == 3


def test_global_wald_math_uses_joint_interaction_covariance() -> None:
    fit = {
        "names": ["intercept", "z_isolation", "a", "b"],
        "beta": np.array([0.0, -0.5, 1.2, -1.1]),
        "covariance": np.diag([1.0, 1.0, 0.04, 0.04]),
    }
    result = global_heterogeneity_test(fit, ["a", "b"])
    assert result["df"] == 2
    assert result["Wald_statistic"] > 60
    assert result["p_value"] < 1e-10


def test_synthetic_channel_specific_retention_curves_pass_global_test() -> None:
    config = _config()
    channels = ["bombus", "non_bombus_bees", "lepidoptera", "flower_visiting_birds"]
    frame, coefficients, global_test, result = fit_n1_primary(
        _projected(channels, 180, heterogeneous=True),
        _qualification(channels),
        _covariates(180),
        config,
    )
    assert len(frame) == 180 * len(channels)
    assert not coefficients.empty
    assert global_test["df"] == len(channels) - 1
    assert global_test["p_value"] < 1e-5
    slopes = result["slopes"].set_index("channel_id")["isolation_slope_log_odds_per_sd"]
    assert slopes["bombus"] < slopes["non_bombus_bees"] < slopes["lepidoptera"]
    assert slopes["flower_visiting_birds"] > slopes["lepidoptera"]
    assert result["model_comparison"]["heterogeneous_minus_common_log_likelihood"] > 0


def test_gate_requires_deletion_stability() -> None:
    config = _config()
    good = pd.DataFrame(
        {"deleted": ["a", "b", "c", "d", "e"], "status": ["fitted"] * 5, "p_value": [0.001, 0.01, 0.02, 0.03, 0.04]}
    )
    bad = good.copy()
    bad["p_value"] = [0.001, 0.3, 0.4, 0.5, 0.6]
    global_test = {"p_value": 0.001}
    comparison = {"heterogeneous_minus_common_log_likelihood": 10.0}
    passed = gate_receipt(global_test, comparison, good, good, 4, config)
    failed = gate_receipt(global_test, comparison, bad, good, 4, config)
    assert passed["N1_pass"] is True
    assert passed["failure_action"] == "N2_may_open"
    assert failed["N1_pass"] is False
    assert failed["failure_action"] == "stop_before_N2_and_keep_frozen_Chapter1"


def test_deletion_robustness_returns_one_row_per_level() -> None:
    config = _config()
    channels = ["bombus", "non_bombus_bees", "lepidoptera", "flower_visiting_birds"]
    frame, metadata = build_primary_frame(
        _projected(channels, 90, heterogeneous=True),
        _qualification(channels),
        _covariates(90),
        config,
    )
    result = deletion_robustness(
        frame,
        metadata["confirmatory_channels"],
        metadata["reference_channel"],
        "source_region_id",
    )
    assert set(result["deleted"]) == {"source-0", "source-1", "source-2"}
    assert set(result["status"]).issubset({"fitted", "not_evaluable", "fit_failed"})
