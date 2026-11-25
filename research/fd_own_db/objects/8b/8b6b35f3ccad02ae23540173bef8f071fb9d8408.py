"""Synthetic-only tests for held-out mainland floral/genus predictive diagnostics."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from fd.core import ProfileSpace
from fd.heldout_mainland_prediction import (
    DonorModel, MODELS, _paired_uncertainty, predict_one_island,
)

PROFILES = {
    0: (frozenset({"white"}), frozenset({"bell"})),
    1: (frozenset({"red"}), frozenset({"bell"})),
    2: (frozenset({"white"}), frozenset({"tubular"})),
    3: (frozenset({"red"}), frozenset({"tubular"})),
}


def synthetic_data():
    donors = pd.DataFrame([
        {"entity_ID": "E1", "accepted_species": "A species1", "genus": "A", "profile_id": 0},
        {"entity_ID": "E1", "accepted_species": "A species2", "genus": "A", "profile_id": 0},
        {"entity_ID": "E2", "accepted_species": "A species1", "genus": "A", "profile_id": 0},
        {"entity_ID": "E2", "accepted_species": "B species1", "genus": "B", "profile_id": 3},
        {"entity_ID": "E2", "accepted_species": "B species2", "genus": "B", "profile_id": 3},
    ])
    assignments = pd.DataFrame([
        {"island_id": "I1", "entity_ID": "E1"},
        {"island_id": "I2", "entity_ID": "E2"},
    ])
    return donors, assignments


def test_donor_model_preserves_normalized_probabilities_and_fallback():
    donors, assignments = synthetic_data()
    model = DonorModel(donors, assignments, 4, shrink=8)
    assert model.n_donor_names == 4   # duplicate A species1 does not count twice globally
    assert model.n_local_source_islands == 2
    a = model.probabilities("I1", "A")
    assert a.shape == (len(MODELS), 4)
    assert np.all(a > 0)
    assert np.allclose(a.sum(axis=1), 1)
    assert a[2, 0] > a[0, 0]  # held-out genus A benefits from genus-specific reference
    assert np.allclose(model.probabilities("I1", "B")[3],
                       model.probabilities("I1", "B")[2])  # local genus unavailable
    assert model.genus_support("I1", "A") == (True, True)
    assert model.genus_support("I1", "B") == (True, False)
    with pytest.raises(ValueError, match="Matched island"):
        model.probabilities("missing_island", "A")


def test_predict_island_holds_out_target_and_returns_occupancy():
    donors, assignments = synthetic_data()
    model = DonorModel(donors, assignments, 4)
    observed = pd.DataFrame([
        {"accepted_species": "A heldout1", "genus": "A", "profile_id": 0},
        {"accepted_species": "A heldout2", "genus": "A", "profile_id": 0},
        {"accepted_species": "B heldout", "genus": "B", "profile_id": 3},
    ])
    assert set(observed.accepted_species).isdisjoint(set(donors.accepted_species))
    stats, actual, predicted = predict_one_island(
        observed, model, "I1", ProfileSpace(PROFILES))
    assert stats["audit"]["n_names"] == 3
    assert stats["audit"]["n_profiles"] == 2
    assert stats["audit"]["global_genus_support"] == 1
    assert stats["audit"]["local_genus_support"] == pytest.approx(2 / 3)
    assert actual.sum() == pytest.approx(1)
    assert np.allclose(predicted.sum(axis=1), 1)
    assert stats["metrics"]["global_genus"]["log_loss"] < (
        stats["metrics"]["global_mainland"]["log_loss"])
    assert stats["metrics"]["local_genus"]["repertoire_L2_error"] >= 0


def test_invalid_training_or_empty_evaluation_rejected():
    donors, assignments = synthetic_data()
    with pytest.raises(ValueError, match="pseudo-count"):
        DonorModel(donors, assignments, 4, shrink=-1)
    with pytest.raises(ValueError, match="distinct complete"):
        predict_one_island(pd.DataFrame(), DonorModel(donors, assignments, 4),
                           "I1", ProfileSpace(PROFILES))
    repeated = pd.DataFrame([
        {"accepted_species": "A heldout", "genus": "A", "profile_id": 0},
        {"accepted_species": "A heldout", "genus": "A", "profile_id": 1},
    ])
    with pytest.raises(ValueError, match="distinct complete"):
        predict_one_island(repeated, DonorModel(donors, assignments, 4),
                           "I1", ProfileSpace(PROFILES))


def test_deterministic_cluster_bootstrap_keeps_matched_unit():
    values = np.array([0.5, 0.75, 1.0, -0.25, 0.2])
    a = _paired_uncertainty(values, seed=174, draws=500)
    b = _paired_uncertainty(values, seed=174, draws=500)
    assert a == b
    assert a[0] < values.mean() < a[1]
    with pytest.raises(ValueError):
        _paired_uncertainty(values, seed=174, draws=5)
