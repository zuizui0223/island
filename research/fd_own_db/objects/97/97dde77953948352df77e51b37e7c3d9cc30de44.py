import numpy as np
import pandas as pd
import pytest

from fd.functional_genus_reference import (
    expected_pair_gap_normal,
    name_heldout_genus_probabilities,
    score_private_predictions,
)


def fixture(n=100):
    return pd.DataFrame({
        "taxon_name": [f"Species {i:04d}" for i in range(n)],
        "genus": [f"Genus{i%5}" for i in range(n)],
        "state": [i % 2 for i in range(n)],
    })


def test_fold_names_do_not_train_their_own_prediction():
    x = fixture()
    p = name_heldout_genus_probabilities(x, tau=1.0)
    assert p.taxon_name.is_unique
    assert set(p.fold) == set(range(5))
    for row in p.itertuples(index=False):
        training = x.loc[
            (~x.taxon_name.isin(p.loc[p.fold == row.fold, "taxon_name"]))
            & (x.genus == row.genus), "state"]
        assert row.n_donors == len(training)
        assert row.p_genus == pytest.approx(
            (training.sum() + row.p_global) / (len(training) + 1.0)
        )


def test_singleton_training_genus_falls_back_to_global():
    x = fixture()
    x.loc[0, "genus"] = "SingletonGenus"
    p = name_heldout_genus_probabilities(x)
    r = p[p.taxon_name == x.loc[0, "taxon_name"]].iloc[0]
    assert r.n_donors == 0
    assert r.p_genus == pytest.approx(r.p_global)


def test_binary_and_duplicate_guards():
    x = fixture()
    x.loc[0, "state"] = np.nan
    with pytest.raises(ValueError, match="binary"):
        name_heldout_genus_probabilities(x)
    x = fixture()
    x.loc[1, "taxon_name"] = x.loc[0, "taxon_name"]
    with pytest.raises(ValueError, match="one row"):
        name_heldout_genus_probabilities(x)


def test_genus_heldout_improvement_is_private_safe():
    x = fixture()
    p = name_heldout_genus_probabilities(x, tau=1.0)
    z = score_private_predictions(p)
    assert z["n_species"] == 100
    assert "taxon_name" not in z
    assert "genus" not in z
    assert "fold" not in z
    with pytest.raises(ValueError, match="Privacy gate"):
        score_private_predictions(p, minimum_species=2)


def test_abs_normal_gap_approximation_and_deterministic_edge():
    assert expected_pair_gap_normal([1,1,1], [0,0]) == pytest.approx(1.0)
    a=np.full(80, .3)
    b=np.full(80, .7)
    approx=expected_pair_gap_normal(a,b)
    rng=np.random.default_rng(20261110)
    obs=np.abs(
        rng.binomial(1,a,size=(5000,len(a))).mean(axis=1)
        - rng.binomial(1,b,size=(5000,len(b))).mean(axis=1))
    assert approx == pytest.approx(float(obs.mean()), abs=.012)
    with pytest.raises(ValueError):
        expected_pair_gap_normal([], [0.5])
