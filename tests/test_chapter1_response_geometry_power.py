from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_response_geometry_power import (
    GEOMETRIES,
    build_candidates,
    run_audit,
    simulate_cell,
)


def _config(replicates: int = 32) -> dict:
    return {
        "contract": "chapter1_response_geometry_identifiability_v1",
        "measurement_domain_representatives": {
            "flower_colour": {"response": "plain_colour"},
            "floral_structural_complexity": {"response": "generalized_form"},
            "reproductive_assurance": {"response": "self_compatibility"},
        },
        "contexts": ["northern_midlatitude", "tropical"],
        "floristic_strata": ["all_native", "native_nonendemic"],
        "minimum_islands_per_design_cell": 20,
        "primary_exposure": {"column": "log_distance_to_continent_km"},
        "future_physical_exposure": {
            "source_specific_water_gap": {
                "status": "not_materialized_do_not_substitute_post_hoc"
            }
        },
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
        "selection_diagnostic": {"strong_selection_delta_AICc": 2.0},
        "simulation": {
            "seed": 20260913,
            "replicates": replicates,
            "effect_logit_sd_grid": [0.50, 0.80],
            "cluster_random_intercept_sd_grid": [0.0, 0.25],
            "truth_breakpoint_quantiles": [0.35, 0.50, 0.65],
        },
        "qualification": {
            "target_effect_logit_sd": 0.50,
            "clustered_stress_sd": 0.25,
            "minimum_true_shape_recovery": 0.80,
            "maximum_false_nonmonotonic_selection_under_true_cline": 0.10,
            "maximum_false_step_or_hinge_selection_under_true_cline": 0.10,
            "maximum_breakpoint_quantile_mae_when_shape_recovered": 0.15,
        },
    }


def _frame(n: int = 80) -> pd.DataFrame:
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


def test_candidate_family_contains_all_five_geometries() -> None:
    config = _config()
    _, _, candidates = build_candidates(_frame(), config)
    assert set(GEOMETRIES) == {candidate.shape for candidate in candidates}


def test_strong_step_is_recoverable_without_cluster_stress() -> None:
    config = _config(replicates=64)
    config["simulation"]["effect_logit_sd_grid"] = [0.80]
    config["simulation"]["cluster_random_intercept_sd_grid"] = [0.0]
    result = simulate_cell(
        _frame(),
        evidence_scope="direct_only",
        measurement_axis="flower_colour",
        outcome="plain_colour",
        context="tropical",
        stratum="native_nonendemic",
        config=config,
    )
    step = result.loc[result["truth_shape"].eq("G2_step")]
    assert step["shape_recovered"].mean() >= 0.90
    recovered = step.loc[step["shape_recovered"]]
    assert recovered["breakpoint_abs_error"].mean() <= 0.15


def test_full_audit_never_opens_observed_geometry(tmp_path: Path) -> None:
    config = _config(replicates=4)
    config["simulation"]["effect_logit_sd_grid"] = [0.50]
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
                for outcome in [
                    "plain_colour",
                    "generalized_form",
                    "self_compatibility",
                ]:
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
    assert (out / "geometry_recovery_summary.csv").exists()
    assert (out / "geometry_cross_scope_qualification.csv").exists()
