"""Synthetic-only N2 lineage/source floral-distance safeguards.

No restricted GIFT, private flora, or taxon-to-trait mappings enter public CI.
"""
import numpy as np
import pandas as pd
import pytest

from fd.distance_ses import distance_statistics, sample_null
from fd.n2_source_distance import (
    _calibration_summaries, _coverage_summaries, _validate_private_against_public,
    heldout_name_calibration,
    donor_training, expected_occupancy_weights,
    paired_primary, source_coverage, simulate_genus_traits,
)


class ToyDonorModel:
    n_profiles = 3

    def probabilities(self, island, genus):
        if island != "focal":
            raise ValueError("Unknown source island")
        if genus == "GenA":
            global_genus = [.8, .2, 0.]
            local_genus = [.2, .8, 0.]
        elif genus == "GenB":
            global_genus = [0., .2, .8]
            local_genus = [0., .4, .6]
        else:
            global_genus = local_genus = [1/3]*3
        return np.array([
            [1/3]*3, [1/3]*3, global_genus, local_genus,
        ], dtype=float)

    def genus_support(self, island, genus):
        return {"GenA": (True, True), "GenB": (True, False)}[genus]


def test_global_name_holdout_excludes_all_mainland_entities():
    source = pd.DataFrame({
        "entity_ID": ["E1", "E2", "E2", "E3", "E3"],
        "work_species": ["GenA target", "GenA target", "GenA donor",
                         "GenB donor", "GenB target"],
    })
    profile = pd.DataFrame({
        "accepted_species": ["GenA target", "GenA donor",
                              "GenB donor", "GenB target"],
        "genus": ["GenA", "GenA", "GenB", "GenB"],
        "profile_id": [0, 1, 2, 0],
    })
    train = donor_training(source, profile,
                           {"GenA target", "GenB target"})
    assert set(train.accepted_species) == {"GenA donor", "GenB donor"}
    assert len(train) == 2
    assert train.duplicated(["entity_ID", "accepted_species"]).sum() == 0
    with pytest.raises(ValueError):
        donor_training(source, profile, {"GenA target",
                                        "GenB target", "GenA donor", "GenB donor"})
    with pytest.raises(ValueError):
        donor_training(source, profile, set())


def test_genus_conditioned_expected_profile_occupancy_not_genus_shuffle():
    model = ToyDonorModel()
    # Two independent GenA names: P(p present) = 1-(1-p)^2.
    weights, n_unique = expected_occupancy_weights(
        {"GenA": 2}, model, "focal", "N2_global_genus")
    expected = np.array([.96, .36, 0.])
    assert n_unique == pytest.approx(expected.sum())
    assert np.allclose(weights, expected / expected.sum())
    local, local_unique = expected_occupancy_weights(
        {"GenA": 2}, model, "focal", "N2_local_genus")
    assert local_unique == pytest.approx(expected.sum())
    assert local[1] > weights[1]  # local donor shifts form/color affinity
    with pytest.raises(ValueError):
        expected_occupancy_weights({}, model, "focal", "N2_global_genus")
    with pytest.raises(ValueError):
        expected_occupancy_weights({"GenA": 0}, model, "focal", "N2_local_genus")
    with pytest.raises(ValueError):
        expected_occupancy_weights({"GenA": 2}, model, "focal", "N1b")


def test_n2_sampled_support_keeps_fixed_richness_and_respects_weights():
    model = ToyDonorModel()
    weights, _ = expected_occupancy_weights(
        {"GenA": 4, "GenB": 4}, model, "focal", "N2_local_genus")
    d = np.array([[0., .2, .7], [.2, 0., .9], [.7, .9, 0.]])
    a = sample_null(d, np.arange(3), 2, draws=100, seed=8, weights=weights)
    b = sample_null(d, np.arange(3), 2, draws=100, seed=8, weights=weights)
    assert np.array_equal(a, b)
    assert np.isfinite(a).all()
    assert a.shape == (100, 2)
    assert np.all(a[:, 0] > 0)
    assert np.isclose(weights.sum(), 1.)


def test_name_level_holdout_calibration_rewards_correct_source_donor():
    from types import SimpleNamespace
    model = ToyDonorModel()
    fake_space = SimpleNamespace(
        lookup={0: 0, 1: 1, 2: 2}, ids=[0, 1, 2],
        matrix=np.array([[0, .5, 1.], [.5, 0, 1.], [1., 1., 0]]))
    targets = pd.DataFrame({
        "accepted_species": ["GenA test1", "GenA test2"],
        "genus": ["GenA", "GenA"], "profile_id": [1, 1],
    })
    got = heldout_name_calibration(targets, model, "focal", fake_space)
    assert got["log_loss_local_genus"] < got["log_loss_global_genus"]
    assert got["trait_distance_local_genus"] < got["trait_distance_global_genus"]
    assert got["log_loss_global_prior"] == pytest.approx(-np.log(1/3))
    with pytest.raises(ValueError):
        heldout_name_calibration(
            pd.concat([targets, targets.iloc[[0]]], ignore_index=True),
            model, "focal", fake_space)


def test_calibration_summary_suppresses_identity_and_tiny_cells():
    base = dict(scope="all", geological_origin="oceanic",
                log_loss_global_prior=1.0, log_loss_global_genus=.8,
                log_loss_local_genus=.7,
                trait_distance_global_prior=.5,
                trait_distance_global_genus=.4,
                trait_distance_local_genus=.3)
    samples = [{**base, "island_id": f"secret_{i}"} for i in range(4)]
    small = _calibration_summaries(samples)
    assert len(small) == 1
    assert small[0]["status"] == "small_cell_suppressed"
    assert "mean_local_genus_name_log_loss" not in small[0]
    reports = _calibration_summaries(samples + [{**base, "island_id": "secret_5"}])
    assert reports[0]["mean_local_vs_global_genus_logloss_improvement"] == pytest.approx(.1)
    assert reports[0]["fraction_islands_local_genus_better_than_global_genus"] == 1.
    assert reports[0]["mean_local_vs_global_genus_trait_distance_improvement"] == pytest.approx(.1)
    assert all("island_id" not in row for row in reports)


def test_generative_genus_reference_preserves_carrier_counts_and_richness():
    model = ToyDonorModel()
    d = np.array([[0., .5, .8], [.5, 0., .6], [.8, .6, 0.]])
    abundance = {"GenA": 2, "GenB": 2}
    sim, thin = simulate_genus_traits(
        abundance, model, "focal", d, variant="N2_local_genus",
        draws=2000, equal_n=2, seed=5)
    analytic = expected_occupancy_weights(
        abundance, model, "focal", "N2_local_genus")[1]
    assert abs(float(sim.mean()) - analytic) < 0.04
    assert len(sim) == 2000 and thin.shape[1] == 2
    assert len(thin) > 50
    other, alt = simulate_genus_traits(
        abundance, model, "focal", d, variant="N2_local_genus",
        draws=2000, equal_n=2, seed=5)
    assert np.array_equal(sim, other)
    assert np.array_equal(thin, alt)
    with pytest.raises(ValueError):
        simulate_genus_traits(abundance, model, "focal", d,
                              variant="N1b", draws=30, equal_n=2, seed=1)


def test_coverage_uses_taxon_carrier_not_genus_weight():
    m = ToyDonorModel()
    a, b = source_coverage(m, "focal", {"GenA": 1, "GenB": 9})
    assert a == 1.
    assert b == pytest.approx(.1)


def test_private_public_composition_audit_rejects_frequency_mismatch():
    actual = pd.DataFrame({
        "island_id": ["i1", "i1", "i1", "i2", "i2"],
        "accepted_species": ["A one", "A two", "B one", "A three", "B two"],
        "genus": ["A", "A", "B", "A", "B"],
        "profile_id": [0, 0, 1, 0, 1],
    })
    g0 = pd.DataFrame({
        "island_id": ["i1", "i2"], "n_names": [3, 2], "n_profiles": [2, 2],
    })
    _validate_private_against_public(
        actual, g0, {"i1": {0: 2, 1: 1}, "i2": {0: 1, 1: 1}}, "all")
    with pytest.raises(ValueError):
        _validate_private_against_public(
            actual, g0, {"i1": {0: 1, 1: 2}, "i2": {0: 1, 1: 1}}, "all")


def test_four_primary_geographic_tests_use_only_local_genus_equal30_nn():
    pair = [{"continental_id": f"c{i}", "oceanic_id": f"o{i}"}
            for i in range(10)]
    geo = {f"{side}{i}": {"spatial_block": f"{side}-{i}"}
           for i in range(10) for side in ("c", "o")}
    rows = []
    for scope in ("all", "direct"):
        for p in pair:
            for side, ses in (("continental_id", 0.), ("oceanic_id", -1.)):
                rows.append({
                    "scope": scope, "island_id": p[side],
                    "null": "N2_local_genus",
                    "mode": "equal30",
                    "endpoint": "mean_nearest_profile_distance",
                    "ses": ses,
                })
    results = paired_primary(
        rows, {"original_block_components": pair,
               "geography_both_verified": pair}, geo)
    assert len(results) == 4
    assert all(x["exact_two_sided_sign_p"] == pytest.approx(2/1024)
               for x in results)
    assert all(x["q_bh_posthoc4"] == pytest.approx(2/1024)
               for x in results)
    # Unknown/missing test islands must not be silently declared independent.
    incomplete = paired_primary(
        rows[:7], {"original_block_components": pair,
                   "geography_both_verified": pair}, geo)
    assert len(incomplete) == 4
    assert all(x["status"] == "insufficient_complete_matched_pairs"
               for x in incomplete)
    assert all("exact_two_sided_sign_p" not in x for x in incomplete)


def test_donor_coverage_aggregation_emits_no_focal_identity():
    base = dict(scope="all", geological_origin="oceanic", n_names=50,
                n_profiles=30, global_fraction=.8, local_fraction=.3,
                global_expected=24., local_expected=20., expected_ratio=2/3,
                generated_richness_median=22.,
                generated_equal30_valid_fraction=.4,
                generated_obs_outside95=False)
    data = [{**base, "island_id": f"private{i}"} for i in range(4)]
    out = _coverage_summaries(data)
    assert len(out) == 1 and out[0]["status"] == "small_cell_suppressed"
    assert "mean_local_genus_donor_fraction" not in out[0]
    data.append({**base, "island_id": "private5"})
    out = _coverage_summaries(data)
    assert out[0]["mean_local_genus_donor_fraction"] == pytest.approx(.3)
    assert all("island_id" not in row for row in out)
