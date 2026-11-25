import numpy as np
import pandas as pd
import pytest

from fd.two_stage_floral_forecast import (
    UNKNOWN_GENUS, REQUIRED_MODELS, audit_no_name_or_group_leakage,
    mixture_from_stage1, predict_flower_distribution, proper_losses,
    score_focal_flora, paired_cluster_summary, dual_loss_dual_group_positive,
)


def weights():
    return {x: {"A": .6, "B": .3, UNKNOWN_GENUS: .1} for x in REQUIRED_MODELS}


def test_unknown_and_unsupported_genera_fall_back_to_blind():
    result = predict_flower_distribution(
        {"A": .6, "unobserved": .3, UNKNOWN_GENUS: .1},
        {"A": [.9, .1]}, [.4, .6])
    np.testing.assert_allclose(result, [.7, .3])


def test_stage1_unknown_mass_does_not_use_focal_genus():
    got = mixture_from_stage1(["A", "B"], [.7, .3])
    assert got[UNKNOWN_GENUS] == .02
    assert sum(got.values()) == pytest.approx(1)
    with pytest.raises(ValueError):
        mixture_from_stage1(["A", "A"], [.7, .3])


def test_both_full_focal_and_group_accepted_names_blacklisted():
    x = dict(focal_names={"focal Taxon"},
             forbidden_full_flora_names={"neighbor Species"},
             genus_donor_names={"remote A"}, flower_donor_names={"remote B"},
             focal_block="10N_100E", focal_first_source="entity_1",
             genus_donor_blocks={"other_block"},
             genus_donor_first_sources={"entity_2"},
             flower_donor_blocks={"other_flower_block"},
             flower_donor_first_sources={"entity_3"})
    audit_no_name_or_group_leakage(**x)
    for key, offender in (("genus_donor_names", "focal Taxon"),
                          ("flower_donor_names", "neighbor Species")):
        with pytest.raises(ValueError, match="leakage"):
            audit_no_name_or_group_leakage(**(x | {key: {offender}}))
    with pytest.raises(ValueError, match="group leakage"):
        audit_no_name_or_group_leakage(**(x | {"genus_donor_blocks": {"10N_100E"}}))
    with pytest.raises(ValueError, match="group leakage"):
        audit_no_name_or_group_leakage(**(x | {"flower_donor_first_sources": {"entity_1"}}))


def test_scores_are_proper_and_negative_gains_retained():
    p = [.9, .1]
    perfect = proper_losses(p, [0]*5)
    wrong = proper_losses(p, [1]*5)
    assert perfect[0] < wrong[0] and perfect[1] < wrong[1]
    models = weights()
    out = score_focal_flora(
        observed_codes=[0]*5, weights_by_model=models,
        external_genus_profiles={"A": [.9, .1], "B": [.7, .3]},
        external_blind=[.4, .6], regional_blind=[.99, .01],
        observed_focal_genera_for_ORACLE_only=["A"]*5)
    assert out["n_original_direct_plant_carriers"] == 5
    assert out["environment_source_nearest_gain_vs_external_blind_Brier"] > 0
    assert out["environment_source_nearest_gain_vs_regional_blind_Brier"] < 0
    assert out["oracle_gain_vs_environment_source_nearest_Brier"] > 0


def test_focal_true_genus_identity_never_enters_nonoracle():
    kw = dict(observed_codes=[0]*5, weights_by_model=weights(),
              external_genus_profiles={"A": [.9, .1], "B": [.7, .3]},
              external_blind=[.4, .6], regional_blind=[.5, .5])
    a = score_focal_flora(**kw, observed_focal_genera_for_ORACLE_only=["A"]*5)
    b = score_focal_flora(**kw, observed_focal_genera_for_ORACLE_only=["B"]*5)
    nonoracle = [k for k in a if not k.startswith("oracle")]
    assert all(a[k] == b[k] for k in nonoracle)
    assert a["oracle_known_genus_Brier"] != b["oracle_known_genus_Brier"]


def test_unqualified_grouping_and_dual_support_gate():
    col1 = "environment_source_nearest_gain_vs_regional_blind_Brier"
    col2 = "environment_source_nearest_gain_vs_regional_blind_logloss"
    frame = pd.DataFrame({"physical_block": [f"b{i//3}" for i in range(60)],
                          "first_source": [f"s{i//2}" for i in range(60)],
                          col1: [.05]*60, col2: [.08]*60})
    a = paired_cluster_summary(frame, group="physical_block", columns=[col1,col2], nboot=75)
    b = paired_cluster_summary(frame, group="first_source", columns=[col1,col2], nboot=75)
    assert dual_loss_dual_group_positive(a, b, "environment_source_nearest", "regional_blind")
    neg = frame.copy(); neg[col2] = -.02
    c = paired_cluster_summary(neg, group="first_source", columns=[col1,col2], nboot=75)
    assert not dual_loss_dual_group_positive(a, c, "environment_source_nearest", "regional_blind")
    one_group = frame.assign(physical_block="same")
    d = paired_cluster_summary(one_group, group="physical_block", columns=[col1,col2], nboot=75)
    assert d["supported_for_inference"] is False
    assert not dual_loss_dual_group_positive(d, b, "environment_source_nearest", "regional_blind")


def test_evaluation_only_residual_exact_split_with_direct_coverage_bias():
    from fd.two_stage_floral_forecast import diagnostic_flower_residual_split
    pieces = diagnostic_flower_residual_split(
        predicted_genus_weights={"A": .4, "B": .6},
        observed_whole_flora_genera=["A", "B"],
        observed_direct_carrier_genera=["A"]*5,
        observed_direct_codes=[0]*5,
        external_genus_profiles={"A": [.8, .2], "B": [.2, .8]},
        external_blind=[.5, .5])
    np.testing.assert_allclose(pieces["taxonomic_genus_forecast_residual"], [.06,-.06])
    np.testing.assert_allclose(pieces["direct_trait_documentation_residual"], [.3,-.3])
    np.testing.assert_allclose(pieces["within_genus_documented_trait_residual"], [.2,-.2])
    np.testing.assert_allclose(pieces["observed_minus_predicted_direct_trait_frequency"], [.56,-.56])



def test_nominal_simpson_paired_diversity_diagnostic_and_unbiased_n():
    from fd.two_stage_floral_forecast import categorical_simpson_diversity_diagnostic
    q = categorical_simpson_diversity_diagnostic(
        [0, 0, 1, 1],
        {"regional_blind": [1., 0.], "geographic_nearest": [.5, .5]})
    assert q["n_direct_carriers"] == 4
    assert q["observed_categorical_simpson_plugin"] == pytest.approx(.5)
    assert q["observed_categorical_simpson_unbiased"] == pytest.approx(2/3)
    assert q["geographic_nearest_predicted_categorical_simpson"] == pytest.approx(.5)
    assert q["geographic_nearest_gain_vs_regional_blind_categorical_simpson_abs_error"] == pytest.approx(.5)


def test_simpson_diagnostic_does_not_conflate_ontology_or_single_carrier():
    from fd.two_stage_floral_forecast import categorical_simpson_diversity_diagnostic
    with pytest.raises(ValueError, match=">=2"):
        categorical_simpson_diversity_diagnostic(
            [0], {"regional_blind": [.5, .5]})
    with pytest.raises(ValueError, match="ontology"):
        categorical_simpson_diversity_diagnostic(
            [0, 1], {"regional_blind": [.5, .5],
                     "other": [.3, .3, .4]})
