"""Exact finite pair nulls, taxon exclusion, and coherent regional draws."""

import importlib
import itertools
from collections import Counter

import numpy as np
import pandas as pd
import pytest

from island_v2.floral_space_null import (
    build_membership,
    draw_profile_assignment,
    prepare_region,
)


def pair_module():
    return importlib.import_module("island_v2.floral_space_pairs")


def example_region(*, coherent=False, exclusive_carriers=False):
    """G={0,0,0,1}, H={2}; eligibility 2 keeps pool-only d in K."""
    profiles = pd.DataFrame({
        "accepted_species": ["a", "b", "c", "d", "s"],
        "genus": ["G", "G", "G", "G", "H"],
        "profile_id": [0, 0, 0, 1, 2],
    })
    members = {"I": ["s", "a", "b"], "J": ["s", "c"]}
    members["L" if coherent else "K"] = ["s", "d"] if coherent else ["d"]
    if exclusive_carriers:
        profiles.loc[len(profiles)] = ["x", "H", 2]
        profiles.loc[len(profiles)] = ["y", "H", 2]
        members["I"].append("x")
        members["J"].append("y")
    flora = pd.DataFrame([(island, taxon) for island, taxa in members.items() for taxon in taxa],
                         columns=["island_id", "accepted_species"])
    geography = pd.DataFrame({"island_id": list(members), "analysis_regime": "tropical"})
    metrics = pd.DataFrame({
        "island_id": list(members), "analysis_regime": "tropical",
        "n_species": [len(taxa) for taxa in members.values()],
        "n_profiles": [2, 2, 2 if coherent else 1],
        "expected_profiles_genus_pool": [2.5, 2, 2 if coherent else 1],
    })
    return prepare_region(build_membership(profiles, flora, geography), metrics,
                          "tropical", min_species=2)


def prepare_example(*, coherent=False, exclusive_carriers=False):
    region = example_region(coherent=coherent, exclusive_carriers=exclusive_carriers)
    pairs = pd.DataFrame({"island_a": ["I", "I"] if coherent else ["I"],
                          "island_b": ["J", "L"] if coherent else ["J"]})
    return pair_module().prepare_pair_design(region, pairs, min_unique=1)


def test_exact_labelled_permutations_keep_duplicate_profile_assignments_and_pool_only_taxa():
    design = prepare_example()
    assert design.region.taxon_names == ("a", "b", "c", "d", "s")
    assert design.region.island_ids == ("I", "J")
    values = []
    for labels in itertools.permutations([0, 0, 0, 1]):
        metrics = pair_module().pair_overlap_metrics(design, np.array([*labels, 2]))
        values.append(float(metrics["profile_jaccard"][0]))
    assert Counter(values) == {0.0: 6, 0.5: 12, 1.0: 6}


def test_shared_taxa_are_removed_before_profiles_even_with_exclusive_same_profile_carriers():
    design = prepare_example(exclusive_carriers=True)
    result = pair_module().pair_overlap_metrics(design, design.region.profile_ids)
    for metric in ("n_profiles_a", "n_profiles_b", "profile_intersection", "profile_union"):
        np.testing.assert_array_equal(result[metric], [2])
    np.testing.assert_array_equal(result["profile_jaccard"], [1.0])
    assert "s" not in {design.region.taxon_names[i] for i in design.member_taxon_indices}
    assert design.pairs.loc[0, "n_unique_a"] == 3
    assert design.pairs.loc[0, "n_unique_b"] == 2


def test_overlapping_pairs_use_one_assignment_and_have_exactly_constant_joint_median():
    design = prepare_example(coherent=True)
    individual_values = set()
    for labels in itertools.permutations([0, 0, 0, 1]):
        result = pair_module().pair_overlap_metrics(design, np.array([*labels, 2]))
        assert np.median(result["profile_jaccard"]) == 0.5
        individual_values.update(result["profile_jaccard"])
    assert individual_values == {0.0, 0.5, 1.0}
    summary, draws, _ = pair_module().run_pair_comparison(design, permutations=37, seed=13)
    np.testing.assert_array_equal(draws.median_profile_jaccard, np.full(37, 0.5))
    assert summary["null_sd"] == 0
    assert summary["null_min"] == summary["null_max"] == 0.5
    assert summary["null_n_distinct_sampled_values"] == 1
    assert summary["p_lower"] == summary["p_upper"] == summary["p_two_sided"] == 1


def test_preparation_canonicalizes_endpoints_preserves_row_order_and_counts():
    region = example_region(coherent=True)
    pairs = pd.DataFrame({"pair_id": ["second", "first"], "island_a": ["L", "J"],
                          "island_b": ["I", "I"]}, index=[8, 2])
    original = pairs.copy(deep=True)
    design = pair_module().prepare_pair_design(region, pairs, min_unique=1)
    pd.testing.assert_frame_equal(pairs, original)
    assert design.n_pairs == 2
    assert design.pairs.pair_id.tolist() == ["second", "first"]
    assert list(design.pairs[["island_a", "island_b"]].itertuples(index=False, name=None)) == [
        ("I", "L"), ("I", "J")]
    np.testing.assert_array_equal(design.pairs.n_taxa_a, [3, 3])
    np.testing.assert_array_equal(design.pairs.n_taxa_b, [2, 2])
    np.testing.assert_array_equal(design.pairs.n_shared_taxa, [1, 1])
    np.testing.assert_array_equal(design.pairs.n_unique_a, [2, 2])
    np.testing.assert_array_equal(design.pairs.n_unique_b, [1, 1])
    assert set(design.member_side_indices) == {0, 1, 2, 3}


@pytest.mark.parametrize("pairs, message", [
    (pd.DataFrame(columns=["island_a", "island_b"]), "empty"),
    (pd.DataFrame({"island_a": ["I"]}), "island_a.*island_b|columns"),
    (pd.DataFrame({"island_a": ["I"], "island_b": ["I"]}), "self"),
    (pd.DataFrame({"island_a": ["I"], "island_b": ["K"]}), "eligible|unknown"),
    (pd.DataFrame({"island_a": ["I"], "island_b": [None]}), "missing"),
    (pd.DataFrame({"island_a": ["I", "J"], "island_b": ["J", "I"]}), "repeat|duplicate"),
])
def test_invalid_pairs_cannot_silently_change_the_selected_comparison(pairs, message):
    with pytest.raises(ValueError, match=message):
        pair_module().prepare_pair_design(example_region(), pairs, min_unique=1)


def test_minimum_support_uses_both_nonshared_taxon_counts_and_default_twenty():
    pairs = pd.DataFrame({"island_a": ["I"], "island_b": ["J"]})
    for arguments in ({}, {"min_unique": 2}):
        with pytest.raises(ValueError, match="support|nonshared|unique"):
            pair_module().prepare_pair_design(example_region(), pairs, **arguments)


@pytest.mark.parametrize("minimum", [0, -1, 1.5, True])
def test_invalid_minimum_support_rejected(minimum):
    with pytest.raises(ValueError, match="positive integer"):
        pair_module().prepare_pair_design(example_region(),
                                         pd.DataFrame({"island_a": ["I"], "island_b": ["J"]}),
                                         min_unique=minimum)


@pytest.mark.parametrize("assignment", [
    np.array([0, 0]), np.array([0., 0., 0., 1., 2.]),
    np.array([0, 0, 0, -1, 2]), np.array([0, 0, 0, 1, 3]),
])
def test_assignment_requires_one_valid_original_code_per_regional_taxon(assignment):
    with pytest.raises(ValueError, match="Assignment|profile"):
        pair_module().pair_overlap_metrics(prepare_example(), assignment)


def test_run_stream_is_repeatable_and_matches_one_regional_assignment_per_draw():
    design = prepare_example()
    first, draws, pair_rows = pair_module().run_pair_comparison(design, permutations=31, seed=123)
    second, repeated, repeated_pairs = pair_module().run_pair_comparison(
        design, permutations=31, seed=123)
    assert first == second
    pd.testing.assert_frame_equal(draws, repeated)
    pd.testing.assert_frame_equal(pair_rows, repeated_pairs)
    assert list(draws.columns) == ["permutation", "median_profile_jaccard",
                                  "median_profile_intersection", "median_profile_union",
                                  "median_n_profiles_a", "median_n_profiles_b"]
    np.testing.assert_array_equal(draws.permutation, np.arange(1, 32))
    rng = np.random.default_rng(123)
    values = []
    for _ in range(31):
        assignment = draw_profile_assignment(design.region, rng)
        for genus in design.region.genus_indices:
            np.testing.assert_array_equal(np.sort(assignment[genus]),
                                          np.sort(design.region.profile_ids[genus]))
        values.append(pair_module().pair_overlap_metrics(design, assignment)["profile_jaccard"][0])
    np.testing.assert_array_equal(draws.median_profile_jaccard, values)
    assert 1.0 in values  # Identity-equivalent assignments remain in the reference.
    assert pair_rows.loc[0, "null_mean_jaccard"] == pytest.approx(np.mean(values))
    assert pair_rows.loc[0, "null_sd_jaccard"] == pytest.approx(np.std(values))
    assert pair_rows.loc[0, "profile_jaccard"] == 1
    assert first["observed_median_jaccard"] == 1
    assert first["observed_minus_null_median"] == 1 - np.median(values)
    assert first["observed_to_null_median"] == 1 / np.median(values)
    assert first["n_pairs"] == 1
    assert first["n_islands"] == 2
    assert first["n_taxa_pool"] == 5
    assert first["n_genera"] == 2
    assert first["n_permutable_genera"] == 1
    assert first["n_permutations"] == 31
    assert first["status"] == "ok"
    assert first["null_min"] == min(values)
    assert first["null_max"] == max(values)
    assert first["null_n_distinct_sampled_values"] == len(set(values))
    assert first["null_q025"] == np.quantile(values, 0.025)
    assert first["null_q975"] == np.quantile(values, 0.975)
    assert first["p_lower"] == 1
    assert first["p_upper"] == (1 + values.count(1)) / 32
    assert first["p_two_sided"] == min(1, 2 * first["p_upper"])
    for metric, observed in [("profile_intersection", 1), ("profile_union", 1),
                             ("n_profiles_a", 1), ("n_profiles_b", 1)]:
        assert first[f"observed_median_{metric}"] == observed
        assert first[f"null_median_{metric}"] == draws[f"median_{metric}"].median()


def test_zero_null_median_preserves_absolute_difference_and_flags_undefined_ratio():
    design = prepare_example()
    # A single sampled null draw with 1 at c makes Jaccard zero (seed 12).
    summary, draws, _ = pair_module().run_pair_comparison(design, permutations=1, seed=12)
    assert draws.median_profile_jaccard.tolist() == [0.0]
    assert summary["null_median"] == 0
    assert summary["observed_minus_null_median"] == 1
    assert np.isnan(summary["observed_to_null_median"])
    assert summary["null_median_is_zero"]
    assert summary["status"] == "ok"


@pytest.mark.parametrize("permutations", [0, -1, 1.5, True])
def test_run_requires_positive_integer_draw_count(permutations):
    with pytest.raises(ValueError, match="positive integer"):
        pair_module().run_pair_comparison(prepare_example(), permutations=permutations)
