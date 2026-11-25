"""Synthetic validations for distance-native floral packing and SES inference."""
import numpy as np
import pytest

from fd.distance_ses import (
    ENDPOINTS, _paired_results, _score, _seed, _summary, distance_statistics,
    exact_sign_flip, regional_leave_one_out_pool, regional_profile_counts,
    sample_null,
)


def test_distance_statistics_exclude_diagonal_and_use_original_d():
    d = np.array([
        [0, .2, .4, .8],
        [.2, 0, .6, .7],
        [.4, .6, 0, .9],
        [.8, .7, .9, 0],
    ])
    got = distance_statistics(d, np.array([0, 1, 2]))
    assert got[0] == pytest.approx((.2 + .4 + .6) / 3)
    assert got[1] == pytest.approx((.2 + .2 + .4) / 3)
    assert d[0, 0] == 0  # must not modify original matrix
    with pytest.raises(ValueError):
        distance_statistics(d, np.array([1, 1]))
    with pytest.raises(ValueError):
        distance_statistics(d, np.array([4, 2]))


def test_n0_null_draws_are_fixed_n_no_replacement_and_deterministic():
    d = np.abs(np.arange(6)[:, None] - np.arange(6)[None, :]) / 5
    a = sample_null(d, np.arange(6), 3, draws=100, seed=14)
    b = sample_null(d, np.arange(6), 3, draws=100, seed=14)
    assert a.shape == (100, 2)
    assert np.array_equal(a, b)
    assert np.isfinite(a).all()
    assert np.all(a[:, 0] > 0)
    assert sample_null(d, np.arange(2), 3, draws=30, seed=1) is None
    assert _seed(11, "all", "N0") != _seed(11, "direct", "N0")


def test_weighted_regional_reference_reduces_uniform_rare_profile_bias():
    # Profile 0 occurs on many peer islands and lies near all other profiles.
    # Weighted without-replacement draws MUST favor this common profile.
    d = np.ones((6, 6))
    np.fill_diagonal(d, 0)
    d[0, 1:] = .01
    d[1:, 0] = .01
    uniform = sample_null(d, np.arange(6), 2, draws=999, seed=10)
    weighted = sample_null(d, np.arange(6), 2, draws=999, seed=10,
                           weights=np.array([100, 1, 1, 1, 1, 1]))
    assert weighted[:, 0].mean() < uniform[:, 0].mean() - .25
    assert weighted[:, 1].mean() < uniform[:, 1].mean() - .25
    assert np.array_equal(
        weighted, sample_null(d, np.arange(6), 2, draws=999, seed=10,
                              weights=np.array([100, 1, 1, 1, 1, 1])))
    with pytest.raises(ValueError):
        sample_null(d, np.arange(6), 2, draws=40, seed=1,
                    weights=np.array([1, 1, 1, 0, 1, 1]))


def test_regional_observed_pool_excludes_focal_only_profiles():
    rows = [
        {"island_id": "a", "analysis_regime": "north"},
        {"island_id": "b", "analysis_regime": "north"},
        {"island_id": "c", "analysis_regime": "south"},
    ]
    carriers = {"a": {1: 2, 2: 3}, "b": {2: 9, 3: 1}, "c": {4: 2}}
    counts = regional_profile_counts(rows, carriers)
    pool = regional_leave_one_out_pool(counts, "north", {1, 2},
                                        {1: 0, 2: 1, 3: 2, 4: 3})
    assert pool.tolist() == [1, 2]  # other island has profiles 2 and 3
    assert regional_leave_one_out_pool(counts, "south", {4},
                                       {1: 0, 2: 1, 3: 2, 4: 3}).size == 0


def test_ses_matches_sample_mean_reference_and_no_fake_zero_variance():
    ref = np.array([[0., 1.], [2., 3.], [4., 5.], [6., 7.]])
    scores = _score(np.array([4., 4.]), ref)
    assert scores[0]["ses"] == pytest.approx(
        (4. - 3.) / np.std([0, 2, 4, 6], ddof=1))
    assert scores[0]["delta"] == pytest.approx(1.)
    assert scores[1]["ses"] == pytest.approx(
        (4. - 4.) / np.std([1, 3, 5, 7], ddof=1))
    assert _score(np.array([.5, .4]), np.ones((30, 2)))[0] is None


def test_exact_component_flips_preserve_shared_geographic_sign():
    diffs = np.array([1., 1., -1., -1., 1., 1., -1., -1.])
    # Two blocks of four pairs: only four configurations, NOT 2^8.
    r = exact_sign_flip(diffs, [list(range(4)), list(range(4, 8))])
    assert r["n_sign_configurations"] == 4
    assert r["exact_two_sided_sign_p"] == pytest.approx(1.)
    with pytest.raises(ValueError):
        exact_sign_flip(diffs, [[0], [1]])


def test_paired_family_is_fixed_32_slots_with_spatial_groups():
    pairs = [
        {"continental_id": f"c{i}", "oceanic_id": f"o{i}"}
        for i in range(10)
    ]
    geography = {key: {"spatial_block": f"b{i}{side}"}
                 for i in range(10)
                 for side, key in (("c", f"c{i}"), ("o", f"o{i}"))}
    scores = []
    for scope in ("all", "direct"):
        for null in ("N0_global496", "N1_regional_observed_LOO"):
            for mode in ("full", "equal30"):
                for endpoint in ENDPOINTS:
                    for pair in pairs:
                        for side, value in (("continental_id", 0.), ("oceanic_id", 1.)):
                            scores.append({
                                "scope": scope, "null": null, "mode": mode,
                                "endpoint": endpoint, "island_id": pair[side],
                                "ses": value,
                            })
    result = _paired_results(scores, pairs, pairs, geography)
    assert len(result) == 32
    assert set(r["n_available_pairs"] for r in result) == {10}
    assert set(r["n_independent_sign_groups"] for r in result) == {10}
    assert all(r["exact_two_sided_sign_p"] == pytest.approx(2 / 1024) for r in result)
    assert all(r["q_bh_posthoc32"] == pytest.approx(2 / 1024) for r in result)


def test_small_group_aggregates_are_suppressed_without_island_identity():
    base = dict(scope="all", null="N0_global496", mode="full",
                endpoint="mean_pairwise_distance", region="north",
                geological_origin="oceanic", n_profiles=35,
                observed=.5, null_mean=.6, delta=-.1, ses=-.5,
                null_sd=.2, thinning_q025=-.5, thinning_q975=-.5)
    rows = [{**base, "island_id": f"is{i}"} for i in range(4)]
    report = _summary(rows)
    assert len(report) == 2  # all and north
    assert all(r["status"] == "small_cell_suppressed" for r in report)
    assert not any("median_SES" in r for r in report)
    assert all("island_id" not in r for r in report)
