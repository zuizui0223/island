from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.m0_m4_spatial_cv import (
    assign_balanced_block_folds,
    run_spatial_cv,
)


def _synthetic_cv_data(n: int = 120) -> tuple[pd.DataFrame, list[str]]:
    rng = np.random.default_rng(99)
    baseline = [
        "log_island_area_km2",
        "abs_latitude",
        "climate_pc1",
        "climate_pc2",
        "climate_pc3",
        "climate_pc4",
        "log_trait_species_richness",
    ]
    rows = []
    for i in range(n):
        # Generate Bombus independently of the baseline so the contract test asks
        # whether spatial CV can recover genuinely incremental held-out signal.
        bombus = float(rng.uniform(0.0, 1.0))
        mixed = np.clip(0.35 - 0.05 * bombus + rng.normal(0, 0.03), 0.05, 0.9)
        sc = np.clip(0.45 - 0.04 * bombus + rng.normal(0, 0.03), 0.05, 0.9)
        alt = np.clip(0.25 + rng.normal(0, 0.04), 0.03, 0.8)
        floral = np.clip(
            0.08 + 0.35 * bombus - 0.03 * mixed + rng.normal(0, 0.01),
            0.01,
            0.95,
        )
        trials = 100 + i % 20
        rows.append(
            {
                "island_id": f"i{i:03d}",
                "spatial_block": f"block_{i // 6:02d}",
                "log_island_area_km2": np.log1p(5 + i),
                "abs_latitude": 5 + i % 70,
                "climate_pc1": np.sin(i / 9),
                "climate_pc2": np.cos(i / 10),
                "climate_pc3": np.sin(i / 6),
                "climate_pc4": np.cos(i / 7),
                "log_trait_species_richness": np.log1p(70 + i % 50),
                "bombus_channel_state": bombus,
                "mixed_mating_share": mixed,
                "self_compatibility_share": sc,
                "alternative_pollen_vector": alt,
                "bombus_x_alternative": bombus * alt,
                "floral_phenotype_trials": trials,
                "floral_phenotype_successes": round(floral * trials),
            }
        )
    return pd.DataFrame(rows), baseline


def test_spatial_blocks_are_never_split_across_folds() -> None:
    data, _ = _synthetic_cv_data()
    mapping = assign_balanced_block_folds(data, n_folds=5, seed=20260713)
    assigned = data["spatial_block"].map(mapping)
    assert assigned.notna().all()
    for _, group in data.assign(fold=assigned).groupby("spatial_block"):
        assert group["fold"].nunique() == 1
    assert assigned.nunique() == 5


def test_repeated_spatial_cv_detects_incremental_bombus_signal() -> None:
    data, baseline = _synthetic_cv_data()
    folds, repeats, summary = run_spatial_cv(
        data,
        baseline,
        n_folds=5,
        repeats=2,
        seed=20260713,
    )
    assert not folds.empty
    assert repeats["repeat"].nunique() == 2
    assert set(summary["comparison"]) == {
        "M0_vs_M1_full_baseline",
        "M0_vs_M1_no_abs_latitude",
        "M2_vs_M3_full_baseline",
        "M2_vs_M3_no_abs_latitude",
        "M4_additive_vs_interaction",
    }
    m1 = summary.loc[summary["comparison"].eq("M0_vs_M1_full_baseline")].iloc[0]
    assert m1["proportion_repeats_full_model_better_loglik"] >= 0.5
