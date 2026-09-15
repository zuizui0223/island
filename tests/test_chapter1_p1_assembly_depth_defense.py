from __future__ import annotations

import math

import numpy as np
import pandas as pd

from island_v2.chapter1_p1_assembly_depth_defense import (
    _attenuation_from_norms,
    _compare_reconstructed_to_frozen,
    _weighted_standardized_slope,
)


def test_attenuation_from_norms() -> None:
    result = _attenuation_from_norms(1.0, 0.8, 0.2, 1.0e-12)
    assert math.isclose(result["family_attenuation"], 0.2)
    assert math.isclose(result["genus_attenuation"], 0.8)
    assert math.isclose(result["family_to_genus_extra_attenuation"], 0.6)
    assert math.isclose(result["conditional_genus_attenuation"], 0.75)


def test_attenuation_fails_closed_near_zero() -> None:
    result = _attenuation_from_norms(0.0, 0.8, 0.2, 1.0e-12)
    assert all(math.isnan(value) for value in result.values())


def test_weighted_standardized_slope_matches_known_signal() -> None:
    n = 30
    x = np.linspace(-2.0, 2.0, n)
    z1 = np.linspace(-1.0, 1.0, n) ** 2
    z2 = np.sin(np.linspace(0.0, 2.0 * np.pi, n))
    z3 = np.cos(np.linspace(0.0, 2.0 * np.pi, n))
    z4 = np.linspace(2.0, -1.0, n) ** 3
    frame = pd.DataFrame(
        {
            "spatial_block": [f"b{i % 6}" for i in range(n)],
            "area": z1,
            "pc1": z1 + 0.1 * z2,
            "pc2": z2,
            "pc3": z3,
            "pc4": z4,
            "distance": x,
        }
    )
    # Construct the outcome from the standardized distance plus nuisance effects.
    distance_z = (x - x.mean()) / x.std(ddof=0)
    frame["response"] = 0.7 * distance_z + 0.1 * z2 - 0.05 * z3
    slope = _weighted_standardized_slope(
        frame,
        "response",
        {f"b{i}": 1 for i in range(6)},
        ["area", "pc1", "pc2", "pc3", "pc4"],
        "distance",
    )
    assert math.isclose(slope, 0.7, rel_tol=0.0, abs_tol=1.0e-10)


def test_compare_reconstructed_to_frozen() -> None:
    frame = pd.DataFrame(
        {
            "island_id": ["i1"],
            "syndrome": ["selfing_core"],
            "source_mode": ["geo_k5"],
            "stratum": ["all_native"],
            "observed_score": [0.2],
            "family_expected": [0.1],
            "genus_expected": [0.15],
            "n_species": [10],
            "n_families": [3],
            "n_genera": [5],
            "after_family_residual": [0.1],
            "family_to_genus_increment": [0.05],
            "after_genus_residual": [0.05],
        }
    )
    max_diff, passed = _compare_reconstructed_to_frozen(frame, frame.copy())
    assert passed is True
    assert max_diff == 0.0

    altered = frame.copy()
    altered.loc[0, "genus_expected"] = 0.151
    _, passed = _compare_reconstructed_to_frozen(frame, altered)
    assert passed is False
