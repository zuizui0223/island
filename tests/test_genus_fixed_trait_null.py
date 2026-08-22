import importlib

import numpy as np
import pandas as pd

mod = importlib.import_module("island_v2.genus_fixed_trait_null")


def test_binary_classifier_excludes_cross_class_multistate():
    positive = {"white", "green"}
    negative = {"red", "blue"}
    assert mod.classify_binary_value("white|green", positive, negative) == 1.0
    assert mod.classify_binary_value("red|blue", positive, negative) == 0.0
    assert np.isnan(mod.classify_binary_value("white|red", positive, negative))


def test_permutation_stays_within_genus():
    states = pd.DataFrame(
        {
            "accepted_species": ["A a", "A b", "B a", "B b"],
            "genus": ["A", "A", "B", "B"],
            "state": [0.0, 1.0, 0.0, 0.0],
        }
    )
    out = mod._permuted_state_map(states, np.random.default_rng(1))
    assert sorted(out.loc[["A a", "A b"]].tolist()) == [0.0, 1.0]
    assert sorted(out.loc[["B a", "B b"]].tolist()) == [0.0, 0.0]


def test_null_runs_by_status_stratum_and_keeps_species_state_shared():
    flora = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i2", "i2", "i3", "i3"],
            "species": ["A a", "A b", "A a", "A c", "A b", "A c"],
        }
    )
    status = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i2", "i2", "i3", "i3"],
            "accepted_species": ["A a", "A b", "A a", "A c", "A b", "A c"],
            "origin_status": ["native"] * 6,
            "endemic_status": ["nonendemic"] * 6,
        }
    )
    evidence = pd.DataFrame(
        {
            "accepted_species": ["A a", "A b", "A c"],
            "trait_name": ["self_incompatibility"] * 3,
            "normalized_value": ["SC", "SI", "SI"],
        }
    )
    master = pd.DataFrame(
        {"accepted_species": ["A a", "A b", "A c"], "genus": ["A", "A", "A"]}
    )
    outcomes = {
        "self_compatibility": {
            "trait_name": "self_incompatibility",
            "positive_states": ["SC"],
            "negative_states": ["SI"],
        }
    }
    result, audit = mod.run_genus_fixed_null(
        flora,
        status,
        evidence,
        master,
        outcomes,
        n_draws=100,
        seed=5,
        min_species_per_island=1,
    )
    native = result.loc[result["stratum"].eq("all_native")]
    assert len(native) == 3
    assert set(native["outcome"]) == {"self_compatibility"}
    assert audit.iloc[0]["n_permutable_genera"] == 1
