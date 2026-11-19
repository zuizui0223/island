"""Synthetic-only guards for original genuine morphology×reproductive assurance pairs."""
import numpy as np
import pandas as pd
import pytest

from fd_original_stage2_name_disjoint_flower_forecast import (
    ARCH_REPRO_PAIRS, PILOT_ENDPOINTS,
    direct_architecture_assurance_pairs,
    full_original_joint_category_ontology,
    product_of_own_marginals, joint_dependency_both_losses_groups,
    run_analysis,
)
from fd.two_stage_floral_forecast import proper_losses, UNKNOWN_GENUS


def source():
    return {
        "floral_form": {"Alpha a": "open", "Beta b": "bell", "Gamma c": "open"},
        "floral_symmetry": {"Alpha a": "radial", "Beta b": "bilateral",
                            "Gamma c": "radial"},
        "tube_depth_class": {"Alpha a": "shallow", "Beta b": "deep",
                             "Gamma c": "deep"},
        "self_incompatibility": {"Alpha a": "SC", "Beta b": "SI",
                                 "Delta d": "SC"},
        "autonomous_selfing_capacity": {"Alpha a": "yes", "Beta b": "no",
                                        "Delta d": "yes"},
    }


def test_exact_six_pair_registration_and_real_same_accepted_species_only():
    d = source()
    out = direct_architecture_assurance_pairs(d)
    assert len(ARCH_REPRO_PAIRS) == len(PILOT_ENDPOINTS) == len(out) == 6
    assert set(out) == set(PILOT_ENDPOINTS)
    for name in PILOT_ENDPOINTS:
        assert set(out[name]) == {"Alpha a", "Beta b"}
        assert "Gamma c" not in out[name] and "Delta d" not in out[name]
    assert out["floral_form__autonomous_selfing_capacity"] == {
        "Alpha a": ("open", "yes"), "Beta b": ("bell", "no"),
    }


def test_source_joint_cartesian_ontology_has_unobserved_combinations():
    d = source()
    ontology = full_original_joint_category_ontology(
        d, "floral_form", "self_incompatibility")
    assert len(ontology) == 4
    assert ("open", "SI") in ontology  # Not a false observed same-species state.
    assert ("bell", "SC") in ontology
    assert len(set(ontology.values())) == 4


def test_independence_projection_retains_marginals_but_not_coupling():
    idx = {("open", "SC"): 0, ("open", "SI"): 1,
           ("bell", "SC"): 2, ("bell", "SI"): 3}
    p = np.array([.4,.1,.1,.4])
    product = product_of_own_marginals(p, idx)
    np.testing.assert_allclose(product, [.25]*4)
    assert np.isclose(product.sum(),1.0)
    b_full, l_full = proper_losses(p,[0,3]*10)
    b_prod, l_prod = proper_losses(product,[0,3]*10)
    assert b_full < b_prod and l_full < l_prod
    with pytest.raises(ValueError,match="cartesian"):
        product_of_own_marginals([.5,.5], {("open","SC"):0,("bell","SI"):1})


def test_signed_direct_joint_dependence_two_scores_two_independent_source_groups():
    def source_g(ci_brier,ci_log):
        return {
            "supported_for_inference": True,
            "predicted_joint_gain_vs_own_independence_Brier":
                {"bootstrap_95":ci_brier},
            "predicted_joint_gain_vs_own_independence_logloss":
                {"bootstrap_95":ci_log},
        }
    a=source_g([.01,.08],[.02,.12])
    b=source_g([.005,.06],[.01,.09])
    assert joint_dependency_both_losses_groups(
        a,b,"predicted_joint_gain_vs_own_independence")
    assert not joint_dependency_both_losses_groups(
        a,source_g([-.03,.06],[.01,.09]),
        "predicted_joint_gain_vs_own_independence")
    assert joint_dependency_both_losses_groups(
        source_g([-.2,-.01],[-.4,-.03]),
        source_g([-.3,-.005],[-.5,-.02]),
        "predicted_joint_gain_vs_own_independence",direction="negative")
    with pytest.raises(ValueError):
        joint_dependency_both_losses_groups(a,b,"x",direction="neither")


def test_original_full_4505_parent_not_skippable_by_joint_endpoint():
    flora=pd.DataFrame([{"island_id":"i", "accepted_species":"Alpha a",
                        "origin_status":"native"}])
    with pytest.raises(ValueError,match="4505"):
        run_analysis(flora,pd.DataFrame(),pd.DataFrame(),source(),nboot=51)



def test_genus_conditional_independence_removes_only_within_genus_pair_dependence():
    from fd_original_stage2_name_disjoint_flower_forecast import (
        genus_conditional_independence_prediction,
        product_of_own_marginals,
    )
    idx={("A0","B0"):0,("A0","B1"):1,
         ("A1","B0"):2,("A1","B1"):3}
    # Two genera have deterministic, *different* joint states. There is no
    # within-genus association to detect, but pooling creates association!
    weights={"Genus0": .49, "Genus1": .49, UNKNOWN_GENUS:.02}
    profiles={"Genus0":np.array([1.,0.,0.,0.]),
              "Genus1":np.array([0.,0.,0.,1.])}
    blind=np.array([.25]*4)
    expected_joint=weights["Genus0"]*profiles["Genus0"]+weights["Genus1"]*profiles["Genus1"]+.02*blind
    stronger=genus_conditional_independence_prediction(weights,profiles,blind,idx)
    np.testing.assert_allclose(expected_joint,stronger)
    pooled_product=product_of_own_marginals(expected_joint,idx)
    assert not np.allclose(pooled_product,stronger)
    # A source "joint" win vs global P(A)P(B) here cannot be claimed as
    # extra within-genus architecture×selfing coupling.
    assert proper_losses(stronger,[0,3]*10) == pytest.approx(
        proper_losses(expected_joint,[0,3]*10))


def test_true_within_genus_coupling_can_add_heldout_joint_skill():
    from fd_original_stage2_name_disjoint_flower_forecast import (
        genus_conditional_independence_prediction,
    )
    idx={("A0","B0"):0,("A0","B1"):1,
         ("A1","B0"):2,("A1","B1"):3}
    prof=np.array([.45,.05,.05,.45])
    blind=np.array([.25]*4)
    weights={"DonorGenus":.8, UNKNOWN_GENUS:.2}
    true_pred=.8*prof+.2*blind
    null=genus_conditional_independence_prediction(
        weights,{"DonorGenus":prof},blind,idx)
    np.testing.assert_allclose(null,blind)
    obs=[0,3]*10
    btrue,ltrue=proper_losses(true_pred,obs)
    bnull,lnull=proper_losses(null,obs)
    assert btrue<bnull and ltrue<lnull


def test_genus_independence_requires_full_original_joint_category_ontology():
    from fd_original_stage2_name_disjoint_flower_forecast import (
        genus_conditional_independence_prediction,
    )
    with pytest.raises(ValueError,match="cartesian"):
        genus_conditional_independence_prediction(
            {"DonorGenus":1.0},
            {"DonorGenus":np.array([.5,.5])},
            np.array([.5,.5]),
            {("A0","B0"):0,("A1","B1"):1})
