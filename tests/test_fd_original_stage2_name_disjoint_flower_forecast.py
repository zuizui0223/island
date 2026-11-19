"""Only synthetic guardrails. No botanical original-data claims here."""
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from fd_original_stage2_name_disjoint_flower_forecast import (
    filtered_genus_matrix, weights_for_flower, _distance_predictions,
    run_analysis,
)
from fd.two_stage_floral_forecast import (
    UNKNOWN_GENUS, predict_flower_distribution, proper_losses,
)


def test_source_donor_name_recurrence_is_removed_even_on_far_islands():
    donors = ["far_a", "far_b"]
    plants = {
        "far_a": {"Alpha focal", "Beta donor", "Gamma donor"},
        "far_b": {"Alpha focal", "Delta donor", "Epsilon donor"},
    }
    names = ["Alpha", "Beta", "Delta", "Epsilon", "Gamma"]
    ids, q, kept, lost = filtered_genus_matrix(
        donors, plants, {"Alpha focal"},
        {name: i for i, name in enumerate(names)},
        min_species=2, min_genera=2)
    assert ids == donors and lost == 0
    assert "Alpha focal" not in kept
    np.testing.assert_allclose(q.toarray()[:, 0], [0, 0])
    np.testing.assert_allclose(q.toarray().sum(axis=1), [1, 1])


def test_whole_donor_island_removed_after_taxon_name_holdout():
    donors = {"island": {"A s1", "B s2", "C s3"}}
    names = {x: j for j, x in enumerate(["A", "B", "C"])}
    ids, q, used, lost = filtered_genus_matrix(
        ["island"], donors, {"A s1", "B s2"}, names,
        min_species=2, min_genera=2)
    assert ids == [] and q.shape == (0, 3) and lost == 1
    assert used == set()


def test_genus_mixture_preserves_exact_unknown_bin_and_only_donor_genus():
    p = weights_for_flower(np.array([.25, .75, .0]), ["A", "B", "C"])
    assert p[UNKNOWN_GENUS] == pytest.approx(.02)
    assert "C" not in p
    assert sum(p.values()) == pytest.approx(1)
    predicted = predict_flower_distribution(
        p, {"A": [.8, .2], "B": [.4, .6]}, [.5, .5])
    np.testing.assert_allclose(predicted, .98*(.25*np.array([.8, .2])
                                            +.75*np.array([.4, .6]))
                                      +.02*np.array([.5, .5]))
    with pytest.raises(ValueError):
        weights_for_flower(np.array([.5,.4,.0]), ["A", "B", "C"])


def test_geo_models_same_exact_source_rows_and_direct_categories():
    Q = csr_matrix([[.9,.1], [.5,.5], [.8,.2], [.2,.8]])
    geo = {
        "features": np.array([[.0]*6, [.1]*6, [2.]*6, [2.1]*6, [.05]*6]),
        "xyz": np.array([[1,0,0],[.99,.1,0],[0,1,0],[.1,.99,0],[.999,.02,0]]),
        "topfive": [{"a"}, {"a"}, {"b"}, {"b"}, {"a"}]
    }
    m = _distance_predictions(4, np.arange(4), Q, geo)
    assert set(m) == {"regional_genus_prior","geographic_nearest",
                      "environment_nearest","environment_source_nearest"}
    for name, p in m.items():
        assert np.isclose(p.sum(),1) and np.isfinite(p).all()
    # The same numeric original-source prior is shared by the four models.
    np.testing.assert_allclose(m["regional_genus_prior"], [.6,.4])
    assert m["environment_source_nearest"][0] > .6
    scores = proper_losses([.6,.4], [0]*5)
    assert scores[0] >= 0 and scores[1] >= 0


def test_guard_prevents_claim_on_changed_4505_full_parent():
    small = pd.DataFrame([
        dict(island_id="a", accepted_species="Alpha sp", origin_status="native")])
    with pytest.raises(ValueError, match="4505"):
        run_analysis(small, pd.DataFrame(), pd.DataFrame(), {})


def test_original_categorical_raoq_models_same_direct_carriers_and_proper_scores():
    from fd.two_stage_floral_forecast import (
        categorical_simpson_diversity_diagnostic,
        predict_flower_distribution, score_focal_flora, REQUIRED_MODELS,
        proper_losses,
    )
    weights={name: {"Alpha":.6,"Beta":.3,UNKNOWN_GENUS:.1}
             for name in REQUIRED_MODELS}
    profiles={"Alpha":[.8,.2],"Beta":[.3,.7]}
    blank=[.5,.5]; local=[.7,.3]
    codes=[0,1,0,0,1,0]
    pred={"regional_blind":local}
    for name in REQUIRED_MODELS:
        pred[name]=predict_flower_distribution(weights[name],profiles,blank)
    q=categorical_simpson_diversity_diagnostic(codes,pred)
    skill=score_focal_flora(observed_codes=codes,weights_by_model=weights,
                            external_genus_profiles=profiles,
                            external_blind=blank,regional_blind=local)
    assert q["n_direct_carriers"]==skill["n_original_direct_plant_carriers"]==6
    assert q["observed_categorical_simpson_unbiased"] > q["observed_categorical_simpson_plugin"]
    for name in REQUIRED_MODELS:
        b,l=proper_losses(pred[name],codes)
        assert skill[f"{name}_Brier"]==pytest.approx(b,abs=1e-12)
        assert skill[f"{name}_logloss"]==pytest.approx(l,abs=1e-12)
        assert 0<=q[f"{name}_predicted_categorical_simpson"]<=1
    assert q["regional_blind_predicted_categorical_simpson"]==pytest.approx(.42)


def test_original_secondary_nominal_q_requires_both_independent_geography_groups():
    from fd_original_stage2_name_disjoint_flower_forecast import (
        is_diversity_error_improvement_block_and_source_supported,
        DIVERSITY_COLUMNS)
    model="environment_source_nearest"
    method="unbiased"
    col=f"{model}_gain_vs_regional_blind_categorical_simpson_{method}_abs_error"
    assert col in DIVERSITY_COLUMNS
    strong={"supported_for_inference":True,col:{"bootstrap_95":[.05,.10]}}
    crossing={"supported_for_inference":True,col:{"bootstrap_95":[-.01,.10]}}
    weak={"supported_for_inference":False,col:{"bootstrap_95":[.01,.10]}}
    assert is_diversity_error_improvement_block_and_source_supported(
        strong,strong,model,method)
    assert not is_diversity_error_improvement_block_and_source_supported(
        strong,crossing,model,method)
    assert not is_diversity_error_improvement_block_and_source_supported(
        strong,weak,model,method)
