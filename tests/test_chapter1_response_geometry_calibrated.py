from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_response_geometry_calibrated import (
    calibrate_cell,
    run_audit,
    validate_cell,
)


def _config(replicates: int = 48) -> dict:
    return {
        "contract": "chapter1_response_geometry_identifiability_v2",
        "parent_failure": {"workflow_run_id": 34764480158},
        "measurement_domain_representatives": {
            "flower_colour": {"response": "plain_colour"},
            "floral_structural_complexity": {"response": "generalized_form"},
            "reproductive_assurance": {"response": "self_compatibility"},
        },
        "contexts": ["northern_midlatitude", "tropical"],
        "floristic_strata": ["all_native", "native_nonendemic"],
        "minimum_islands_per_design_cell": 20,
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
            "G0_flat": {"family": "monotonic_null_family"},
            "G1_cline": {"family": "monotonic_null_family"},
            "G2_step": {
                "family": "nonlinear_family",
                "breakpoint_search_quantiles": [0.25, 0.35, 0.50, 0.65, 0.75],
                "breakpoint_parameter_penalty": 1,
            },
            "G3_hinge": {
                "family": "nonlinear_family",
                "breakpoint_search_quantiles": [0.25, 0.35, 0.50, 0.65, 0.75],
                "breakpoint_parameter_penalty": 1,
            },
            "G4_reversal": {
                "family": "nonlinear_family",
                "turning_point_required_between_quantiles": [0.20, 0.80],
                "endpoint_derivatives_must_have_opposite_signs": True,
            },
        },
        "calibration": {
            "seed": 20260914,
            "replicates_per_scenario": replicates,
            "monotonic_truths": {
                "G0_flat": {"effect_logit_sd": [0.0]},
                "G1_cline": {"effect_logit_sd": [0.25, 0.50, 0.80]},
            },
            "cluster_random_intercept_sd": 0.25,
            "nonlinear_gate_quantile": 0.95,
        },
        "validation": {
            "seed": 20260915,
            "replicates_per_scenario": replicates,
            "cluster_random_intercept_sd": 0.25,
            "monotonic_truths": {
                "G0_flat": {"effect_logit_sd": [0.0]},
                "G1_cline": {"effect_logit_sd": [0.25, 0.50, 0.80]},
            },
            "nonlinear_target_effect_logit_sd": 0.80,
            "nonlinear_truth_breakpoint_quantiles": [0.35, 0.50, 0.65],
        },
        "qualification": {
            "maximum_validation_false_nonlinear_rate": 0.20,
            "minimum_shape_recovery_rate": 0.60,
            "maximum_breakpoint_quantile_mae": 0.20,
        },
    }


def _frame(n: int = 120) -> pd.DataFrame:
    x = np.linspace(0.0, 7.0, n)
    return pd.DataFrame(
        {
            "island_id": [f"i{i}" for i in range(n)],
            "successes": np.full(n, 20),
            "trials": np.full(n, 40),
            "log_distance_to_continent_km": x,
            "log_island_area_km2": np.log1p(np.arange(n) + 1),
            "climate_pc1": np.sin(0.73 * x),
            "climate_pc2": np.cos(1.11 * x),
            "climate_pc3": np.sin(1.67 * x + 0.3),
            "climate_pc4": np.cos(2.31 * x - 0.2),
            "spatial_block": [f"b{i // 4}" for i in range(n)],
        }
    )


def test_calibration_and_validation_are_held_out() -> None:
    config = _config(replicates=64)
    frame = _frame()
    critical, calibration, z, weights, candidates = calibrate_cell(
        frame,
        evidence_scope="direct_only",
        measurement_axis="flower_colour",
        outcome="plain_colour",
        context="tropical",
        stratum="native_nonendemic",
        config=config,
    )
    detail, qualification = validate_cell(
        frame,
        evidence_scope="direct_only",
        measurement_axis="flower_colour",
        outcome="plain_colour",
        context="tropical",
        stratum="native_nonendemic",
        critical_D=critical,
        z=z,
        weights=weights,
        candidates=candidates,
        config=config,
    )
    assert np.isfinite(critical)
    assert len(calibration) == 4
    assert set(qualification["shape"]) == {"G2_step", "G3_hinge", "G4_reversal"}
    assert detail.loc[detail["truth_family"].eq("monotonic"), "false_nonlinear_rate"].notna().all()


def test_strong_step_survives_calibrated_gate() -> None:
    config = _config(replicates=96)
    frame = _frame()
    critical, _, z, weights, candidates = calibrate_cell(
        frame,
        evidence_scope="all_analysis_eligible",
        measurement_axis="flower_colour",
        outcome="plain_colour",
        context="tropical",
        stratum="native_nonendemic",
        config=config,
    )
    _, qualification = validate_cell(
        frame,
        evidence_scope="all_analysis_eligible",
        measurement_axis="flower_colour",
        outcome="plain_colour",
        context="tropical",
        stratum="native_nonendemic",
        critical_D=critical,
        z=z,
        weights=weights,
        candidates=candidates,
        config=config,
    )
    step = qualification.loc[qualification["shape"].eq("G2_step")].iloc[0]
    assert float(step["minimum_shape_recovery_rate"]) >= 0.80
    assert float(step["maximum_breakpoint_quantile_mae"]) <= 0.15


def test_full_v2_audit_keeps_observed_geometry_closed(tmp_path: Path) -> None:
    config = _config(replicates=8)
    config_path = tmp_path / "config.yml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    cov_rows = []
    count_rows = []
    for context_index, context in enumerate(config["contexts"]):
        for i in range(28):
            island = f"{context}_{i}"
            x = i / 4 + context_index / 10
            cov_rows.append(
                {
                    "island_id": island,
                    "analysis_regime": context,
                    "spatial_block": f"{context}_b{i // 4}",
                    "log_distance_to_continent_km": x,
                    "log_island_area_km2": np.log1p(i + 2),
                    "climate_pc1": np.sin(0.73 * x),
                    "climate_pc2": np.cos(1.11 * x),
                    "climate_pc3": np.sin(1.67 * x + 0.3),
                    "climate_pc4": np.cos(2.31 * x - 0.2),
                }
            )
            for stratum in config["floristic_strata"]:
                for outcome in ["plain_colour", "generalized_form", "self_compatibility"]:
                    count_rows.append(
                        {
                            "island_id": island,
                            "successes": 20,
                            "trials": 40,
                            "share": 0.5,
                            "outcome": outcome,
                            "stratum": stratum,
                        }
                    )
    covariates = pd.DataFrame(cov_rows)
    counts = pd.DataFrame(count_rows)
    cov_path = tmp_path / "cov.csv"
    all_path = tmp_path / "all.csv"
    direct_path = tmp_path / "direct.csv"
    covariates.to_csv(cov_path, index=False)
    counts.to_csv(all_path, index=False)
    counts.to_csv(direct_path, index=False)
    out = tmp_path / "out"

    manifest = run_audit(
        all_counts=all_path,
        direct_counts=direct_path,
        covariates_csv=cov_path,
        config_path=config_path,
        output_dir=out,
    )
    assert manifest["observed_geometry_fitted"] is False
    assert manifest["observed_breakpoint_reported"] is False
    assert manifest["calibration_and_validation_seeds_differ"] is True
    assert (out / "geometry_v2_cross_scope_shape_qualification.csv").exists()
