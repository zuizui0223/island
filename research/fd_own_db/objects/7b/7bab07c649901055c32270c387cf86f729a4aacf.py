"""Synthetic tests for origin-paired profile/name-support rarefaction."""
from __future__ import annotations

import numpy as np
import pytest

from fd.core import ProfileSpace
from fd.coverage_rarefaction import (
    _draw_one_carrier_support, _model_one_scope, _sign_swap_scores,
    SCREENS, ENDPOINTS,
)


def _profiles(n: int = 44):
    return {
        i: (frozenset({f"color_{i}"}), frozenset({"bell"} if i % 2 else {"tubular"}))
        for i in range(n)
    }


def _original_carriers(n: int = 44):
    return {i: 1 + (i % 4) for i in range(n)}


def test_rarefy_profile_counts_and_name_support_without_replacement():
    observed = _original_carriers()
    rng = np.random.default_rng(1986)
    assert len(_draw_one_carrier_support(
        observed, screen="original_full", rng=rng)) == 44
    rare = _draw_one_carrier_support(
        observed, screen="rarefy_30_profiles", rng=rng)
    assert len(rare) == 30 and len(set(rare)) == 30
    assert set(rare).issubset(observed)
    subset = _draw_one_carrier_support(
        observed, screen="rarefy_to_pair_min_profiles",
        rng=rng, target_profiles=33)
    assert len(subset) == 33
    name_sample = _draw_one_carrier_support(
        observed, screen="rarefy_40_name_units", rng=rng)
    assert 1 <= len(name_sample) <= 40
    assert set(name_sample).issubset(observed)
    with pytest.raises(ValueError, match="30 profiles"):
        _draw_one_carrier_support({1: 10}, screen="rarefy_30_profiles", rng=rng)
    with pytest.raises(ValueError, match="40 names"):
        _draw_one_carrier_support({1: 6}, screen="rarefy_40_name_units", rng=rng)
    with pytest.raises(ValueError, match="pair-common"):
        _draw_one_carrier_support(observed, screen="rarefy_to_pair_min_profiles",
                                  rng=rng, target_profiles=200)


def test_signswap_metrics_and_stochastic_output_stable():
    space = ProfileSpace(_profiles())
    kernels = {
        "exact_profile_L2": np.eye(len(space.ids)),
        "categorical_similarity_kernel": 1 - space.matrix,
    }
    rows = [{
        "continental_id": f"c{i}", "oceanic_id": f"o{i}",
        "region": "northern_midlatitude" if i < 6 else "tropical",
    } for i in range(8)]
    carriers = {
        **{f"c{i}": {p: 1+(p%3) for p in range(35)} for i in range(8)},
        **{f"o{i}": {p: 1+(p%2) for p in range(9, 44)} for i in range(8)},
    }
    a = _model_one_scope(
        "rarefy_30_profiles", "all", rows, carriers, space, kernels,
        replications=20, permutations=999, seed=33)
    b = _model_one_scope(
        "rarefy_30_profiles", "all", rows, carriers, space, kernels,
        replications=20, permutations=999, seed=33)
    assert a == b
    assert len(a) == len(ENDPOINTS)
    assert all(row["n_resamples"] == 20 for row in a)
    assert all(0 < row["p_median"] <= 1 for row in a)
    assert all(0 <= row["fraction_resamples_p_below_005"] <= 1 for row in a)
    single = _model_one_scope(
        "single_state_20_profiles", "all", rows, carriers, space, kernels,
        replications=20, permutations=999, seed=33)
    assert all(x["n_resamples"] == 1 and x["n_permutations_each"] == 9999
               for x in single)
    small = _model_one_scope(
        "single_state_20_profiles", "all", rows[:5], carriers, space, kernels,
        replications=20, permutations=999, seed=33)
    assert all(x["status"] == "insufficient_pairs" for x in small)


def test_rarefaction_stat_has_expected_zero_for_identical_repertoires():
    space = ProfileSpace(_profiles(34))
    same = {i: 2 for i in range(33)}
    pairs = [{
        "continental_id": f"c{i}", "oceanic_id": f"o{i}",
        "region": "northern_midlatitude",
    } for i in range(8)]
    c = {**{f"c{i}": same for i in range(8)},
         **{f"o{i}": same for i in range(8)}}
    kernels = {"exact_profile_L2": np.eye(len(space.ids)),
               "categorical_similarity_kernel": 1 - space.matrix}
    test = _model_one_scope("original_full", "all", pairs, c, space, kernels,
                            replications=20, permutations=999, seed=33)
    assert all(r["observed_statistic_median"] == pytest.approx(0) for r in test)
    assert all(r["p_median"] == 1 for r in test)
