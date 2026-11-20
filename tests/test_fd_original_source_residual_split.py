"""Original source-private genus vs Direct-documentation vs trait-transfer residuals."""
import numpy as np
import pytest

from fd.two_stage_floral_forecast import (
    UNKNOWN_GENUS, predict_flower_distribution, proper_losses
)
from fd_original_stage2_name_disjoint_flower_forecast import (
    FLORAL_FORECAST_RESIDUAL_COLUMNS,
    FLORAL_FORECAST_RESIDUAL_NAMES,
    source_original_noncategorical_residual_diagnostic,
)


def fixture():
    # Each actual whole-flora and Direct name refers to an individual
    # botanical accepted species in this SYNTHETIC fixture only.
    weights={"Alpha":.65,"Beta":.25,UNKNOWN_GENUS:.10}
    full={"Alpha a", "Alpha b", "Alpha c", "Beta d", "Gamma e"}
    direct=["Alpha a","Beta d","Gamma e"]
    codes=np.array([0,1,0], dtype=int)
    external_prior=np.array([.55,.45])
    external_profiles={
        "Alpha":np.array([.80,.20]),
        "Beta":np.array([.10,.90]),
    }
    pred=predict_flower_distribution(weights,external_profiles,external_prior)
    brier,_=proper_losses(pred,codes)
    kwargs=dict(
        predicted_genus_weights=weights,
        observed_whole_flora_names=full,
        observed_direct_species_names=direct,
        observed_direct_codes=codes,
        congener_profiles=external_profiles,
        external_category_prior=external_prior,
        exact_original_model_Brier=brier,
    )
    return kwargs


def test_three_source_components_exactly_sum_to_observed_flora_category_error():
    result=source_original_noncategorical_residual_diagnostic(**fixture())
    assert set(result)==set(FLORAL_FORECAST_RESIDUAL_COLUMNS)
    squared=result["observed_minus_predicted_direct_trait_squared_norm"]
    proj=sum(result[f"{x}_signed_projection_on_total_error"]
             for x in FLORAL_FORECAST_RESIDUAL_NAMES)
    assert squared==pytest.approx(proj,abs=1e-12)
    assert result["observed_direct_trait_plugin_Q"] == pytest.approx(4/9)
    assert all(result[f"{x}_squared_norm"]>=0
               for x in FLORAL_FORECAST_RESIDUAL_NAMES)


def test_residual_is_a_signed_projection_not_an_additive_variance_partition():
    result=source_original_noncategorical_residual_diagnostic(**fixture())
    squared=result["observed_minus_predicted_direct_trait_squared_norm"]
    components=sum(result[f"{x}_squared_norm"]
                   for x in FLORAL_FORECAST_RESIDUAL_NAMES)
    # Residual norms cannot be naively summed: three vector components may
    # align or cancel. Signed projections, NOT norm sums, are additive.
    assert abs(components-squared)>1e-5


def test_direct_carrier_species_must_really_occur_in_same_whole_flora():
    d=fixture()
    d["observed_direct_species_names"]=["Alpha a","Beta d","Absent species"]
    with pytest.raises(ValueError,match="absent"):
        source_original_noncategorical_residual_diagnostic(**d)


def test_original_source_proper_score_unchanged_under_audit():
    d=fixture()
    d["exact_original_model_Brier"]+=.1
    with pytest.raises(AssertionError,match="Original proper Brier"):
        source_original_noncategorical_residual_diagnostic(**d)


def test_only_after_prediction_genus_truth_may_inform_residuals():
    d=fixture()
    base=source_original_noncategorical_residual_diagnostic(**d)
    new=dict(d)
    new["observed_whole_flora_names"]={"Alpha a","Beta d","Gamma e"}
    alternate=source_original_noncategorical_residual_diagnostic(**new)
    # The audit should respond to a different source whole-flora genus
    # composition even while the precomputed model weights and Brier stay
    # identical. This is *not* an oracle-prediction model change.
    assert alternate["taxonomic_genus_forecast_residual_squared_norm"] != pytest.approx(
        base["taxonomic_genus_forecast_residual_squared_norm"])
    assert alternate["observed_minus_predicted_direct_trait_squared_norm"] == pytest.approx(
        base["observed_minus_predicted_direct_trait_squared_norm"])
