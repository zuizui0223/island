from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_h5d_distributed_threshold_identifiability import (
    _distributed_probability,
    _smooth_probability,
    simulate_design,
    summarize_simulation,
)


def _config() -> dict:
    return {
        "simulation": {
            "replicates": 12,
            "seed": 7,
            "train_fraction": 0.7,
            "latent_genera": 200,
            "distributed_threshold_sigma": [0.25, 0.5],
            "smooth_cline_midpoint_sigma": 0.5,
            "smooth_cline_log_slope_mean": np.log(2.0),
            "smooth_cline_log_slope_sd": 0.35,
            "probability_clip": 0.0001,
        },
        "qualification": {
            "minimum_classification_accuracy": 0.8,
            "maximum_false_distributed_selection_under_smooth_clines": 0.1,
            "maximum_relative_sigma_error": 0.25,
            "minimum_sigma_recovery_fraction": 0.8,
        },
    }


def _design() -> pd.DataFrame:
    x = np.linspace(-2.0, 2.0, 80)
    return pd.DataFrame(
        {
            "island_id": [f"i{i}" for i in range(len(x))],
            "n_genera": np.repeat(40, len(x)),
            "z_distance": x,
            "spatial_block": [f"b{i // 4}" for i in range(len(x))],
        }
    )


def test_generators_are_monotone_probability_surfaces():
    x = np.linspace(-2, 2, 25)
    rng = np.random.default_rng(11)
    distributed = _distributed_probability(x, rng, sigma=0.5, latent_genera=500, eps=1e-4)
    smooth = _smooth_probability(
        x,
        rng,
        midpoint_sigma=0.5,
        log_slope_mean=np.log(2),
        log_slope_sd=0.35,
        latent_genera=500,
        eps=1e-4,
    )
    assert np.all(np.diff(distributed) >= 0)
    assert np.all(np.diff(smooth) >= 0)
    assert np.all((distributed > 0) & (distributed < 1))
    assert np.all((smooth > 0) & (smooth < 1))


def test_simulation_materializes_both_generators_without_observed_outcomes():
    config = _config()
    result = simulate_design(_design(), config, seed_offset=0)
    assert set(result["generator"]) == {"distributed_threshold", "smooth_cline"}
    assert len(result) == 36
    summary = summarize_simulation(result, config)
    assert 0 <= summary["classification_accuracy"] <= 1
    assert 0 <= summary["false_distributed_selection_under_smooth"] <= 1
    assert 0 <= summary["sigma_recovery_fraction"] <= 1
