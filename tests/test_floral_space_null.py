"""Conditional randomization checks using a fully enumerable assemblage."""

import importlib
import itertools
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest


def null_module():
    return importlib.import_module("island_v2.floral_space_null")


def example_frames(region="northern_midlatitude"):
    """Two shared-taxon target islands, plus two one-taxon pool islands.

    Alpha's regional profile multiset is {0, 0, 1}; Beta's is {0, 2}.
    Drawing one of each gives E[profiles] = 2 - (2/3)(1/2) = 5/3.
    """
    profiles = pd.DataFrame({
        "accepted_species": ["Alpha a", "Alpha b", "Alpha c", "Beta d", "Beta e"],
        "genus": ["Alpha", "Alpha", "Alpha", "Beta", "Beta"],
        "profile_id": [0, 0, 1, 0, 2],
    })
    prefix = region + ":"
    flora = pd.DataFrame({
        "island_id": [prefix + s for s in ["one", "one", "two", "two", "pool_b", "pool_c"]],
        "accepted_species": ["Alpha a", "Beta d", "Alpha a", "Beta e", "Alpha b", "Alpha c"],
    })
    geography = pd.DataFrame({
        "island_id": [prefix + s for s in ["one", "two", "pool_b", "pool_c"]],
        "analysis_regime": [region] * 4,
    })
    metrics = pd.DataFrame({
        "island_id": [prefix + s for s in ["one", "two", "pool_b", "pool_c"]],
        "analysis_regime": [region] * 4,
        "n_species": [2, 2, 1, 1],
        "n_profiles": [1, 2, 1, 1],
        "expected_profiles_genus_pool": [5 / 3, 5 / 3, 1.0, 1.0],
    })
    return profiles, flora, geography, metrics


def example_design():
    profiles, flora, geography, metrics = example_frames()
    module = null_module()
    membership = module.build_membership(profiles, flora, geography)
    return module.prepare_region(membership, metrics, "northern_midlatitude", min_species=2)


def test_pool_keeps_taxa_from_ineligible_islands_and_excludes_other_regions():
    profiles, flora, geography, metrics = example_frames()
    profiles.loc[len(profiles)] = ["Other x", "Other", 9]
    flora.loc[len(flora)] = ["elsewhere", "Other x"]
    geography.loc[len(geography)] = ["elsewhere", "tropical"]
    module = null_module()
    design = module.prepare_region(
        module.build_membership(profiles, flora, geography), metrics,
        "northern_midlatitude", min_species=2,
    )
    assert set(design.taxon_names) == {"Alpha a", "Alpha b", "Alpha c", "Beta d", "Beta e"}
    assert len(design.island_ids) == 2
    np.testing.assert_allclose(design.expected_profiles, [5 / 3, 5 / 3])


def test_all_enumerated_permutations_match_hand_calculated_null_distribution():
    module, design = null_module(), example_design()
    statistics, counts = [], []
    # Duplicate permutations of Alpha's two identical zeros correctly remain
    # in the label-permutation distribution: 3! * 2! = 12 equally likely draws.
    for alpha, beta in itertools.product(itertools.permutations([0, 0, 1]),
                                         itertools.permutations([0, 2])):
        assignment = np.array([*alpha, *beta])
        island_counts = module.island_profile_counts(design, assignment)
        counts.append(island_counts)
        statistics.append(module.median_profile_ratio(island_counts, design.expected_profiles))
    np.testing.assert_allclose(np.mean(counts, axis=0), [5 / 3, 5 / 3])
    np.testing.assert_allclose(sorted(statistics), [0.9] * 8 + [1.2] * 4)
    # The median-ratio statistic's null median is 0.9 here, not an assumed 1.
    assert np.median(statistics) == pytest.approx(0.9)


def test_one_taxon_has_one_shuffled_state_across_both_islands():
    module, design = null_module(), example_design()
    # The shared Alpha a is 1 in both islands: each then has two profiles.
    np.testing.assert_array_equal(
        module.island_profile_counts(design, np.array([1, 0, 0, 0, 2])), [2, 2]
    )
    # With shared Alpha a = 0, only the island containing Beta's zero can have
    # one profile. Independent island-level relabelling could violate this.
    for beta in ([0, 2], [2, 0]):
        counts = module.island_profile_counts(design, np.array([0, 0, 1, *beta]))
        assert sorted(counts.tolist()) == [1, 2]


def test_random_assignments_preserve_each_genus_profile_multiset():
    module, design = null_module(), example_design()
    rng = np.random.default_rng(20261106)
    for _ in range(12):
        assigned = module.draw_profile_assignment(design, rng)
        np.testing.assert_array_equal(np.sort(assigned[:3]), [0, 0, 1])
        np.testing.assert_array_equal(np.sort(assigned[3:]), [0, 2])
    np.testing.assert_array_equal(design.profile_ids, [0, 0, 1, 0, 2])


def test_statistic_is_median_of_ratios_not_mean_or_pooled_ratio():
    assert null_module().median_profile_ratio([1, 1, 10], [2, 2, 10]) == 0.5


def test_stale_analytic_expectation_stops_the_comparison():
    profiles, flora, geography, metrics = example_frames()
    metrics.loc[0, "expected_profiles_genus_pool"] = 2.0
    module = null_module()
    with pytest.raises(ValueError, match="expectation"):
        module.prepare_region(module.build_membership(profiles, flora, geography), metrics,
                              "northern_midlatitude", min_species=2)


def test_duplicate_membership_does_not_silently_increase_species_count():
    profiles, flora, geography, _ = example_frames()
    flora = pd.concat([flora, flora.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="unique"):
        null_module().build_membership(profiles, flora, geography)


def test_inconsistent_species_profiles_cannot_generate_island_specific_states():
    profiles, flora, geography, _ = example_frames()
    profiles = pd.concat([profiles, profiles.iloc[[0]].assign(profile_id=2)], ignore_index=True)
    with pytest.raises(ValueError, match="unique"):
        null_module().build_membership(profiles, flora, geography)


def test_unknown_regime_label_cannot_silently_drop_a_planned_comparison():
    profiles, flora, geography, _ = example_frames()
    geography.loc[0, "analysis_regime"] = "northern_highlatitude"
    with pytest.raises(ValueError, match="regime"):
        null_module().build_membership(profiles, flora, geography)


def test_two_sided_p_uses_finite_inclusive_tails():
    result = null_module().finite_two_sided_p(1.2, np.array([0.8, 0.9, 1.0, 1.1]))
    assert result["p_lower"] == 1.0
    assert result["p_upper"] == pytest.approx(0.2)
    assert result["p_two_sided"] == pytest.approx(0.4)


def test_degenerate_null_with_ties_has_p_one_and_no_fictitious_variation():
    result = null_module().finite_two_sided_p(1.0, np.ones(9))
    assert result["p_lower"] == 1.0
    assert result["p_upper"] == 1.0
    assert result["p_two_sided"] == 1.0
    assert result["null_q025"] == result["null_median"] == result["null_q975"] == 1.0


def test_roundoff_ties_are_included_in_both_tails():
    result = null_module().finite_two_sided_p(1.0, np.array([1 - 1e-15, 1 + 1e-15, 2.0]))
    assert result["p_lower"] == 0.75
    assert result["p_upper"] == 1.0
    assert result["p_two_sided"] == 1.0


def test_bh_keeps_missing_comparisons_in_the_fixed_family():
    got = null_module().bh_adjust_fixed([0.01, 0.04, 0.03, np.nan])
    np.testing.assert_allclose(got[:3], [0.04, 0.16 / 3, 0.16 / 3])
    assert np.isnan(got[3])


def test_comparison_is_seed_reproducible_and_does_not_filter_identity_draws():
    module, design = null_module(), example_design()
    first, draws1 = module.run_comparison(design, permutations=17, seed=123)
    second, draws2 = module.run_comparison(design, permutations=17, seed=123)
    np.testing.assert_array_equal(draws1, draws2)
    assert len(draws1) == 17
    assert first == second
    assert first["observed_ratio"] == pytest.approx(0.9)
    assert first["n_taxa"] == 5
    assert first["n_genera"] == 2
    assert first["n_islands"] == 2
    assert first["n_pool_only_taxa"] == 2


def test_cli_preserves_fixed_sixteen_comparisons_and_aggregate_draws(tmp_path):
    regions = ["northern_midlatitude", "northern_high_latitude", "tropical",
               "southern_extratropical"]
    scopes = ["all", "direct", "all_single_state", "direct_single_state"]
    frames = [example_frames(r) for r in regions]
    pd.concat([f[1] for f in frames]).to_csv(tmp_path / "flora.csv", index=False)
    pd.concat([f[2] for f in frames]).to_csv(tmp_path / "geography.csv", index=False)
    for scope in scopes:
        directory = tmp_path / scope
        directory.mkdir()
        frames[0][0].to_csv(directory / "species_profiles.csv.gz", index=False)
        pd.concat([f[3] for f in frames]).to_csv(directory / "island_metrics.csv", index=False)
    process = subprocess.run(
        [sys.executable, "-m", "island_v2.floral_space_null", "--results", str(tmp_path),
         "--flora", str(tmp_path / "flora.csv"), "--geography", str(tmp_path / "geography.csv"),
         "--permutations", "7", "--min-species", "2", "--seed", "123"],
        capture_output=True, text=True, env={**os.environ, "PYTHONPATH": "src"},
    )
    assert process.returncode == 0, process.stderr
    summary = pd.read_csv(tmp_path / "permutation_summary.csv")
    draws = pd.read_csv(tmp_path / "permutation_null_draws.csv.gz")
    manifest = json.loads((tmp_path / "permutation_manifest.json").read_text())
    assert len(summary) == 16
    assert not summary.duplicated(["evidence_scope", "analysis_regime"]).any()
    assert set(summary.evidence_scope) == set(scopes)
    assert set(summary.analysis_regime) == set(regions)
    assert set(summary.n_permutations) == {7}
    assert set(summary.status) == {"ok"}
    np.testing.assert_allclose(summary.q_bh_16, null_module().bh_adjust_fixed(summary.p_two_sided))
    assert len(draws) == 16 * 7
    assert set(draws.columns) == {"evidence_scope", "analysis_regime", "permutation", "median_profile_ratio"}
    assert manifest["fixed_comparison_family_size"] == 16
    assert manifest["status"] == "complete"
