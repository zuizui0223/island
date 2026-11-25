"""Exact categorical Brier-to-nominal-Rao-Q diagnostic; no ecological mechanism.

If q is the empirical species-Direct categorical frequency on a held-out
island and p is a model's predicted category simplex, the multiclass mean
Brier score with one-hot species outcomes is exactly

    Brier(p) = Q(q) + ||p-q||^2,

where Q(q)=1-||q||^2 is the observed *plugin* nominal Rao/Gini–Simpson
diversity. Yet predicted diversity has signed error

    Q(p)-Q(q) = -2*q.dot(p-q) - ||p-q||^2.

Thus improved proper category prediction does NOT imply improved nominal
functional-diversity prediction. The projection against q determines the
direction of the diversity error. No trait-distance matrix, visitor-guild,
ancestral mainland, or real within-species island adaptation is inferred.

Sample-size correction n/(n-1)Q_plugin is supplied as a separate diagnostic,
and is unbiased only under an appropriate iid random-draw categorical
sampling model. Nonrandom botanical trait documentation may violate it.
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np


def _simplex(x: Sequence[float], label: str) -> np.ndarray:
    p = np.asarray(x, dtype=float)
    if (p.ndim != 1 or len(p) < 2 or not np.isfinite(p).all()
            or (p < 0).any() or not np.isclose(float(p.sum()), 1, atol=1e-10)):
        raise ValueError(f"{label} needs >=2 finite nonnegative categories summing to one")
    return p


def nominal_simpson(p: Sequence[float]) -> float:
    """Nominal Rao Q / Gini–Simpson diversity with unit category distances."""
    p = _simplex(p, "Categorical probability")
    return float(1.0 - p @ p)


def category_forecast_identity(
    probability: Sequence[float],
    empirical_frequency: Sequence[float],
    *,
    n_direct: int | None = None,
) -> dict[str, float]:
    """Exact Brier and diversity algebra for ONE matched empirical Direct flora.

    n_direct is optional and used ONLY for the iid categorical sample-size
    correction. It never changes p, q, sample eligibility, or proper scores.
    """
    p = _simplex(probability, "Predicted category probability")
    q = _simplex(empirical_frequency, "Observed empirical category frequency")
    if len(p) != len(q):
        raise ValueError("Predicted and observed category ontology not aligned")
    if n_direct is not None:
        if isinstance(n_direct, bool) or not isinstance(n_direct, (int, np.integer)) or n_direct < 2:
            raise ValueError("Direct observation n must be integer >=2")
        if not np.allclose(n_direct * q, np.rint(n_direct * q), atol=1e-8):
            raise ValueError("Empirical category frequencies inconsistent with n Direct carriers")

    e = p - q
    concentration = float(e @ e)
    orientation = float(-2.0 * q @ e)
    brier = float(1.0 + p @ p - 2.0 * p @ q)
    q_obs = float(1.0 - q @ q)
    q_pred = float(1.0 - p @ p)
    signed_diversity_error = q_pred - q_obs
    if not np.isclose(brier, q_obs + concentration, atol=1e-12):
        raise AssertionError("Proper Brier/excess squared category error identity failed")
    if not np.isclose(signed_diversity_error, orientation-concentration, atol=1e-12):
        raise AssertionError("Nominal Rao signed error projection identity failed")

    result = {
        "observed_plugin_nominal_Q": q_obs,
        "predicted_nominal_Q": q_pred,
        "plugin_Q_signed_forecast_error": signed_diversity_error,
        "plugin_Q_absolute_forecast_error": abs(signed_diversity_error),
        "mean_multiclass_Brier": brier,
        "squared_category_frequency_error": concentration,
        "nominal_Q_error_linear_projection": orientation,
        "nominal_Q_error_negative_quadratic": -concentration,
    }
    if n_direct is not None:
        q_corrected = q_obs * n_direct / (n_direct - 1)
        result.update(
            n_direct_species_carriers=int(n_direct),
            observed_iid_unbiased_nominal_Q=float(q_corrected),
            iid_unbiased_Q_signed_forecast_error=float(q_pred-q_corrected),
            iid_unbiased_Q_absolute_forecast_error=float(abs(q_pred-q_corrected)),
        )
    return result


def compare_category_and_nominal_diversity(
    empirical_frequency: Sequence[float],
    blind_reference: Sequence[float],
    improved_model: Sequence[float],
    *,
    n_direct: int | None = None,
) -> dict[str, float | bool]:
    """Signed exact matched-model gains, without choosing favourable baselines.

    A positive Brier gain means the improved_model better forecasts category
    frequency. Positive absolute-Q gain means its predicted Q is nearer the
    observed Q. Their directions CAN DISAGREE even for the exact same plants.
    """
    q = _simplex(empirical_frequency, "Observed empirical category frequency")
    base = category_forecast_identity(blind_reference, q, n_direct=n_direct)
    model = category_forecast_identity(improved_model, q, n_direct=n_direct)
    brier_gain = base["mean_multiclass_Brier"]-model["mean_multiclass_Brier"]
    l2_gain = base["squared_category_frequency_error"] - model["squared_category_frequency_error"]
    if not np.isclose(brier_gain, l2_gain, atol=1e-12):
        raise AssertionError("Common observed Brier benchmark does not cancel")
    result: dict[str, float | bool] = {
        "Brier_gain_model_over_blind": float(brier_gain),
        "squared_category_frequency_error_gain": float(l2_gain),
        "plugin_Q_absolute_error_gain_model_over_blind": float(
            base["plugin_Q_absolute_forecast_error"]
            - model["plugin_Q_absolute_forecast_error"]),
        "same_matched_species_ontology": True,
        "categorical_Q_is_not_multivariate_morphological_volume": True,
    }
    if n_direct is not None:
        result["iid_unbiased_Q_absolute_error_gain_model_over_blind"] = float(
            base["iid_unbiased_Q_absolute_forecast_error"]
            - model["iid_unbiased_Q_absolute_forecast_error"])
    return result



def rao_distance_decomposition(
    probability: Sequence[float],
    empirical_frequency: Sequence[float],
    dissimilarity: Sequence[Sequence[float]],
) -> dict[str, float]:
    """Exact general Rao Q(p)=p' D p, not restricted to nominal unit distance.

    D must be a fixed, finite, symmetric, nonnegative categorical distance
    table with zero diagonal, shared by model and source observation. Unlike
    ||p-q||², the quadratic term e' D e may have EITHER sign because ordinary
    dissimilarity matrices need not be positive-semidefinite. This diagnostic
    is a mathematical identity only; it does not construct or estimate a
    flower morphology/trait distance matrix from source data.
    """
    p = _simplex(probability, "Predicted category probability")
    q = _simplex(empirical_frequency, "Observed empirical category frequency")
    if len(p) != len(q):
        raise ValueError("Predicted and observed Rao Q category ontology not aligned")
    d = np.asarray(dissimilarity, dtype=float)
    if (d.shape != (len(p), len(p)) or not np.isfinite(d).all()
            or np.any(d < 0) or not np.allclose(d, d.T, atol=1e-10)
            or not np.allclose(np.diag(d), 0, atol=1e-10)):
        raise ValueError("Shared category distance D must be symmetric, finite and hollow")
    e = p - q
    q_obs = float(q @ d @ q)
    q_pred = float(p @ d @ p)
    linear = float(2 * e @ d @ q)
    quadratic = float(e @ d @ e)
    if not np.isclose(q_pred-q_obs, linear+quadratic, atol=1e-12):
        raise AssertionError("General pairwise Rao Q expansion failed")
    return {
        "observed_Rao_Q_fixed_D": q_obs,
        "predicted_Rao_Q_fixed_D": q_pred,
        "Rao_Q_signed_forecast_error": float(q_pred-q_obs),
        "Rao_Q_linear_error_from_trait_distances": linear,
        "Rao_Q_quadratic_error_from_trait_distances": quadratic,
        "category_squared_error_is_not_Rao_Q_error": True,
    }
