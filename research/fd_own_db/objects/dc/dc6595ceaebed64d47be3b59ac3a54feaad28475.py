"""Exact nominal-Q source interpretation tests; all synthetic, no ecology claims."""
import numpy as np
import pytest

from fd.nominal_diversity_identity import (
    nominal_simpson, category_forecast_identity,
    compare_category_and_nominal_diversity,
)


def test_brier_and_simpson_projection_identity_on_heldout_category_frequencies():
    q=[.8,.2]
    p=[.6,.4]
    z=category_forecast_identity(p,q,n_direct=10)
    assert z["observed_plugin_nominal_Q"]==pytest.approx(.32)
    assert z["predicted_nominal_Q"]==pytest.approx(.48)
    assert z["mean_multiclass_Brier"]==pytest.approx(.4)
    assert z["squared_category_frequency_error"]==pytest.approx(.08)
    assert z["nominal_Q_error_linear_projection"]==pytest.approx(.24)
    assert z["nominal_Q_error_negative_quadratic"]==pytest.approx(-.08)
    assert z["plugin_Q_signed_forecast_error"]==pytest.approx(.16)
    assert z["observed_iid_unbiased_nominal_Q"]==pytest.approx(.32*10/9)
    assert nominal_simpson(q)==pytest.approx(.32)


def test_category_brier_improves_even_as_nominal_rao_q_worsens_same_plants():
    # One exact 10-species island Direct trait roster: eight category-0, two category-1.
    # Moving .6 -> .95 improves empirical Brier but overshoots concentration.
    z=compare_category_and_nominal_diversity(
        empirical_frequency=[.8,.2], blind_reference=[.6,.4],
        improved_model=[.95,.05], n_direct=10)
    assert z["Brier_gain_model_over_blind"]==pytest.approx(.035)
    assert z["squared_category_frequency_error_gain"]==pytest.approx(.035)
    assert z["plugin_Q_absolute_error_gain_model_over_blind"]==pytest.approx(-.065)
    assert z["iid_unbiased_Q_absolute_error_gain_model_over_blind"] < 0
    assert z["same_matched_species_ontology"]
    assert z["categorical_Q_is_not_multivariate_morphological_volume"]


def test_perfect_category_forecast_zero_concentration_and_diversity_error():
    q=np.array([.2,.3,.5])
    z=category_forecast_identity(q,q)
    assert z["squared_category_frequency_error"]==pytest.approx(0,abs=1e-15)
    assert z["plugin_Q_signed_forecast_error"]==pytest.approx(0,abs=1e-15)
    assert z["nominal_Q_error_linear_projection"]==pytest.approx(0,abs=1e-15)
    assert z["mean_multiclass_Brier"]==pytest.approx(nominal_simpson(q))
    cmp=compare_category_and_nominal_diversity(q,q,q)
    assert cmp["Brier_gain_model_over_blind"]==pytest.approx(0)
    assert cmp["plugin_Q_absolute_error_gain_model_over_blind"]==pytest.approx(0)


def test_observed_empirical_sample_count_finite_sample_guard():
    q=[.8,.2]
    with pytest.raises(ValueError,match="integer >=2"):
        category_forecast_identity(q,q,n_direct=1)
    with pytest.raises(ValueError,match="integer >=2"):
        category_forecast_identity(q,q,n_direct=True)
    with pytest.raises(ValueError,match="inconsistent"):
        category_forecast_identity(q,q,n_direct=6)
    assert "observed_iid_unbiased_nominal_Q" not in category_forecast_identity(q,q)
    assert category_forecast_identity(q,q,n_direct=10)["n_direct_species_carriers"]==10


@pytest.mark.parametrize("p,q",[
    ([-.1,1.1],[.5,.5]),
    ([.2,.2],[.5,.5]),
    ([.5,.5,.0],[.5,.5]),
    ([np.nan,.5],[.5,.5]),
    ([.5,.5],[.1,.2,.7]),
    ([1],[1]),
])
def test_invalid_category_ontologies_and_simplexes_fail_closed(p,q):
    with pytest.raises(ValueError):
        category_forecast_identity(p,q)


def test_multicategory_exact_identity_random_nonzero_probabilities():
    rng=np.random.default_rng(20261120)
    for k in [2,3,13,40]:
        for _ in range(16):
            q=rng.dirichlet(np.ones(k))
            base=rng.dirichlet(np.ones(k))
            pred=rng.dirichlet(np.ones(k))
            z=compare_category_and_nominal_diversity(q,base,pred)
            ref=category_forecast_identity(base,q)
            model=category_forecast_identity(pred,q)
            assert z["Brier_gain_model_over_blind"]==pytest.approx(
                np.square(base-q).sum()-np.square(pred-q).sum(), abs=1e-11)
            assert model["plugin_Q_signed_forecast_error"]==pytest.approx(
                -(pred-q)@(pred+q),abs=1e-11)
            assert ref["observed_plugin_nominal_Q"]==pytest.approx(
                model["observed_plugin_nominal_Q"])



def test_general_rao_fixed_distance_matrix_recovers_nominal_simpson_identity():
    from fd.nominal_diversity_identity import rao_distance_decomposition
    q=[.8,.2]
    p=[.95,.05]
    dist=np.ones((2,2))-np.eye(2)
    r=rao_distance_decomposition(p,q,dist)
    z=category_forecast_identity(p,q)
    assert r["observed_Rao_Q_fixed_D"]==pytest.approx(z["observed_plugin_nominal_Q"])
    assert r["predicted_Rao_Q_fixed_D"]==pytest.approx(z["predicted_nominal_Q"])
    assert r["Rao_Q_signed_forecast_error"]==pytest.approx(-.225)
    assert r["Rao_Q_linear_error_from_trait_distances"]==pytest.approx(-.18)
    assert r["Rao_Q_quadratic_error_from_trait_distances"]==pytest.approx(-.045)


def test_generalized_rao_fixed_mixed_trait_dissimilarity_exact_expansion():
    from fd.nominal_diversity_identity import rao_distance_decomposition
    D=np.array([[0,.2,1],[.2,0,.8],[1,.8,0]])
    p=np.array([.2,.5,.3])
    q=np.array([.3,.2,.5])
    r=rao_distance_decomposition(p,q,D)
    assert r["observed_Rao_Q_fixed_D"]==pytest.approx(q@D@q)
    assert r["predicted_Rao_Q_fixed_D"]==pytest.approx(p@D@p)
    assert (r["Rao_Q_linear_error_from_trait_distances"]
            +r["Rao_Q_quadratic_error_from_trait_distances"]
            ==pytest.approx(p@D@p-q@D@q,abs=1e-12))


@pytest.mark.parametrize("D",[
    [[0,1],[.5,0]],
    [[0,-.1],[-.1,0]],
    [[.1,.3],[.3,0]],
    [[0,float("inf")],[float("inf"),0]],
    [[0,1,2],[1,0,3]],
])
def test_rao_quadratic_distance_must_be_finite_symmetric_hollow(D):
    from fd.nominal_diversity_identity import rao_distance_decomposition
    with pytest.raises(ValueError,match="distance"):
        rao_distance_decomposition([.7,.3],[.4,.6],D)
