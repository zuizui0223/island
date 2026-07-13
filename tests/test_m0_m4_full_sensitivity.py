from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.m0_m4_full_sensitivity import (
    collinearity_diagnostics,
    fit_sensitivity_snapshot,
    residualize_bombus,
    spatial_block_bootstrap,
)


def _island_data(n: int = 96) -> tuple[pd.DataFrame, list[str]]:
    rng = np.random.default_rng(123)
    rows = []
    baseline = [
        "log_island_area_km2",
        "abs_latitude",
        "climate_pc1",
        "climate_pc2",
        "climate_pc3",
        "climate_pc4",
        "log_trait_species_richness",
    ]
    for i in range(n):
        climate1 = np.sin(i / 9)
        climate2 = np.cos(i / 11)
        latitude = 5 + i % 65
        area = np.log1p(5 + i * 1.5)
        richness = np.log1p(70 + i % 50)
        bombus = (
            0.45
            - 0.12 * climate1
            + 0.003 * latitude
            + 0.02 * area
            + rng.normal(0, 0.08)
        )
        mixed_share = np.clip(0.35 - 0.10 * bombus + rng.normal(0, 0.04), 0.03, 0.95)
        sc_share = np.clip(0.50 - 0.08 * bombus + rng.normal(0, 0.04), 0.03, 0.95)
        alt_share = np.clip(0.30 - 0.05 * bombus + rng.normal(0, 0.04), 0.03, 0.95)
        floral_share = np.clip(
            0.12 + 0.12 * bombus - 0.03 * mixed_share - 0.02 * sc_share + rng.normal(0, 0.025),
            0.01,
            0.95,
        )
        trials = 90 + i % 30
        rows.append(
            {
                "island_id": f"i{i:03d}",
                "spatial_block": f"b{i // 8:02d}",
                "log_island_area_km2": area,
                "abs_latitude": latitude,
                "climate_pc1": climate1,
                "climate_pc2": climate2,
                "climate_pc3": np.sin(i / 6),
                "climate_pc4": np.cos(i / 7),
                "log_trait_species_richness": richness,
                "bombus_channel_state": bombus,
                "mixed_mating_share": mixed_share,
                "self_compatibility_share": sc_share,
                "alternative_pollen_vector": alt_share,
                "bombus_x_alternative": bombus * alt_share,
                "floral_phenotype_trials": trials,
                "floral_phenotype_successes": round(floral_share * trials),
                "mixed_mating_trials": trials,
                "mixed_mating_successes": round(mixed_share * trials),
                "self_compatibility_trials": trials,
                "self_compatibility_successes": round(sc_share * trials),
            }
        )
    return pd.DataFrame(rows), baseline


def test_residualization_removes_linear_baseline_component() -> None:
    data, baseline = _island_data()
    residualized, metadata = residualize_bombus(data, baseline)
    assert metadata["baseline_r2_for_bombus"] > 0.1
    complete = residualized.dropna(subset=["bombus_residualized", *baseline])
    for column in baseline:
        correlation = complete["bombus_residualized"].corr(complete[column])
        assert abs(float(correlation)) < 1e-8


def test_collinearity_and_sensitivity_snapshot_are_finite() -> None:
    data, baseline = _island_data()
    residualized, _ = residualize_bombus(data, baseline)
    correlations, vif, diagnostics = collinearity_diagnostics(
        residualized, [*baseline, "bombus_channel_state"]
    )
    assert not correlations.empty
    assert np.isfinite(vif["vif"]).all()
    assert diagnostics["condition_number_standardized_design"] > 1

    coefficients, fits = fit_sensitivity_snapshot(residualized, baseline, 1.96)
    assert set(fits["status"]) <= {"fitted", "max_iter_reached"}
    assert {
        "M1_raw_bombus",
        "M1_residualized_bombus",
        "M3_raw_direct_bombus",
        "M3_residualized_direct_bombus",
        "M4_raw_interaction",
        "M4_residualized_interaction",
    } <= set(coefficients["analysis"])


def test_spatial_block_bootstrap_returns_percentile_intervals() -> None:
    data, baseline = _island_data()
    residualized, _ = residualize_bombus(data, baseline)
    draws, summary = spatial_block_bootstrap(
        residualized,
        baseline,
        replicates=20,
        seed=20260713,
        z_threshold=1.96,
    )
    assert draws["replicate"].nunique() == 20
    assert {
        "M1_raw_bombus",
        "M1_residualized_bombus",
        "M3_raw_direct_bombus",
        "M3_residualized_direct_bombus",
        "M3_raw_parallel_total_indirect",
        "M3_residualized_parallel_total_indirect",
        "M4_raw_interaction",
        "M4_residualized_interaction",
    } <= set(summary["metric"])
    assert (summary["n_valid_replicates"] > 0).all()
    assert (summary["ci_2_5"] <= summary["ci_97_5"]).all()
