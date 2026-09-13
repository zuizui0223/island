import numpy as np
import pandas as pd

from island_v2.chapter1_response_geometry_observed import (
    classify_cross_scope,
    fit_scope_cell,
)


def _v2_config() -> dict:
    return {
        "primary_exposure": {"column": "log_distance_to_continent_km"},
        "baseline_covariates": [
            "log_island_area_km2",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
        ],
        "cluster_column": "spatial_block",
        "candidate_geometries": {
            "G0_flat": {},
            "G1_cline": {},
            "G2_step": {
                "breakpoint_search_quantiles": [0.25, 0.35, 0.50, 0.65, 0.75],
                "breakpoint_parameter_penalty": 1,
            },
            "G3_hinge": {
                "breakpoint_search_quantiles": [0.25, 0.35, 0.50, 0.65, 0.75],
                "breakpoint_parameter_penalty": 1,
            },
            "G4_reversal": {
                "turning_point_required_between_quantiles": [0.20, 0.80]
            },
        },
    }


def _observed_config() -> dict:
    return {
        "primary_exposure": "log_distance_to_continent_km",
        "baseline_covariates": [
            "log_island_area_km2",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
        ],
        "cluster_column": "spatial_block",
        "uncertainty_reporting": {
            "breakpoint_consensus": {
                "maximum_scope_difference_quantile": 0.15
            }
        },
    }


def _step_frame(n: int = 120) -> pd.DataFrame:
    x = np.linspace(0.0, 7.0, n)
    successes = np.where(x <= np.median(x), 8, 42)
    return pd.DataFrame(
        {
            "island_id": [f"i{i}" for i in range(n)],
            "successes": successes,
            "trials": np.full(n, 50),
            "log_distance_to_continent_km": x,
            "log_island_area_km2": np.log1p(np.arange(n) + 1),
            "climate_pc1": np.sin(0.73 * x),
            "climate_pc2": np.cos(1.11 * x),
            "climate_pc3": np.sin(1.67 * x + 0.3),
            "climate_pc4": np.cos(2.31 * x - 0.2),
            "spatial_block": [f"b{i // 4}" for i in range(n)],
        }
    )


def _scope_row(scope: str, critical_D: float = 5.0) -> dict:
    return fit_scope_cell(
        _step_frame(),
        evidence_scope=scope,
        measurement_axis="flower_colour",
        outcome="plain_colour",
        context="tropical",
        stratum="native_nonendemic",
        critical_D=critical_D,
        observed_config=_observed_config(),
        v2_config=_v2_config(),
    )


def test_strong_step_passes_frozen_scope_gate() -> None:
    row = _scope_row("all_analysis_eligible")
    assert row["nonlinear_gate_pass"] is True
    assert row["best_nonlinear_shape"] == "G2_step"
    assert row["robust_fit_status"] == "fit"
    assert row["effect_ci95_low"] > 0


def test_cross_scope_step_requires_v2_qualification() -> None:
    scope = pd.DataFrame(
        [_scope_row("all_analysis_eligible"), _scope_row("direct_only")]
    )
    qualified = pd.DataFrame(
        [
            {
                "measurement_axis": "flower_colour",
                "outcome": "plain_colour",
                "context": "tropical",
                "stratum": "native_nonendemic",
                "shape": "G2_step",
                "headline_shape_qualified": True,
            }
        ]
    )
    result = classify_cross_scope(scope, qualified, _observed_config())
    assert result.iloc[0]["classification"] == "identified_step"

    qualified.loc[0, "headline_shape_qualified"] = False
    result = classify_cross_scope(scope, qualified, _observed_config())
    assert result.iloc[0]["classification"] == "nonlinear_shape_unresolved"


def test_both_scopes_must_cross_nonlinear_gate() -> None:
    scope = pd.DataFrame(
        [_scope_row("all_analysis_eligible"), _scope_row("direct_only", 1e9)]
    )
    qualified = pd.DataFrame(
        [
            {
                "measurement_axis": "flower_colour",
                "outcome": "plain_colour",
                "context": "tropical",
                "stratum": "native_nonendemic",
                "shape": "G2_step",
                "headline_shape_qualified": True,
            }
        ]
    )
    result = classify_cross_scope(scope, qualified, _observed_config())
    assert result.iloc[0]["classification"] == "monotonic_or_unresolved"
