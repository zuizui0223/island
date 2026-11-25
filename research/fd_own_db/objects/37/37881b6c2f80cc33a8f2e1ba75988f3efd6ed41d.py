"""Leakage-audited bridge from geography-predicted genera to floral categories.

This module is a numeric adapter, NOT a source-data analysis. It does not
estimate donor island distances or flower profiles. Those must be fitted on
out-of-block, out-of-first-source and accepted-name-disjoint source material
before calling these functions. No island trait adaptation is inferred.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

UNKNOWN_GENUS = "__UNKNOWN_GENUS__"
REQUIRED_MODELS = ("regional_genus_prior", "geographic_nearest",
                   "environment_nearest", "environment_source_nearest")
COMPARATORS = ("external_blind", "regional_blind")


def _probability(p: Sequence[float], label: str) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    if (p.ndim != 1 or p.size < 2 or not np.isfinite(p).all()
            or (p < 0).any() or not np.isclose(p.sum(), 1.0, atol=1e-8)):
        raise ValueError(f"{label} must be a finite normalized categorical distribution")
    return p


def audit_no_name_or_group_leakage(*, focal_names: set[str],
                                   forbidden_full_flora_names: set[str],
                                   genus_donor_names: set[str],
                                   flower_donor_names: set[str],
                                   focal_block: str,
                                   focal_first_source: str,
                                   genus_donor_blocks: set[str],
                                   genus_donor_first_sources: set[str],
                                   flower_donor_blocks: set[str],
                                   flower_donor_first_sources: set[str]) -> None:
    """Require full accepted-name and geographic deletion *before* model fitting.

    The original environmental genus stage excluded focal-group ISLANDS,
    but other islands can contain the SAME botanical species as the focal.
    A new species-name-disjoint run is required to pass this diagnostic.
    """
    forbidden = focal_names | forbidden_full_flora_names
    if not focal_names or not focal_block or not focal_first_source:
        raise ValueError("Missing focal full-flora identity or spatial/source grouping")
    if genus_donor_names & forbidden or flower_donor_names & forbidden:
        raise ValueError("Accepted-species leakage into genus or flower donor training")
    if focal_block in genus_donor_blocks or focal_first_source in genus_donor_first_sources:
        raise ValueError("Geographical/source donor group leakage into genus training")
    if focal_block in flower_donor_blocks or focal_first_source in flower_donor_first_sources:
        raise ValueError("Geographical/source donor group leakage into flower training")


def mixture_from_stage1(genus_names: Sequence[str], probabilities: Sequence[float],
                        *, unknown_mass: float = 0.02) -> dict[str, float]:
    """Convert old stage-1 conditional known-genus p to an explicit unknown bin.

    Do not use focal accepted species or observed genus composition here.
    """
    p = _probability(probabilities, "stage1 genus weights")
    if (len(genus_names) != len(p) or len(set(genus_names)) != len(p)
            or UNKNOWN_GENUS in genus_names):
        raise ValueError("Stage-1 genus vocabulary must be unique and fully aligned")
    if not (0 <= unknown_mass < 1):
        raise ValueError("Invalid unknown-genus probability mass")
    result = {name: float((1-unknown_mass)*mass)
              for name, mass in zip(genus_names, p)}
    result[UNKNOWN_GENUS] = float(unknown_mass)
    return result


def predict_flower_distribution(genus_weights: Mapping[str, float],
                                external_genus_profiles: Mapping[str, Sequence[float]],
                                external_blind: Sequence[float]) -> np.ndarray:
    """P(flower category) = sum_g P(g | geography) P(category | g, external).

    Unsupported and unknown genera use only the external blind source prior.
    The focal's true genus roster is not an argument to this function.
    """
    blind = _probability(external_blind, "external-blind flower prior")
    weights = np.asarray(list(genus_weights.values()), dtype=float)
    if (weights.size == 0 or not np.isfinite(weights).all() or
            (weights < 0).any() or not np.isclose(weights.sum(), 1, atol=1e-8)):
        raise ValueError("Geography-predicted genus weights must sum to one")
    out = np.zeros_like(blind)
    for genus, weight in genus_weights.items():
        p = blind if genus == UNKNOWN_GENUS or genus not in external_genus_profiles \
                  else _probability(external_genus_profiles[genus], "external congener")
        if len(p) != len(blind):
            raise ValueError("Inconsistent flower-category ontology")
        out += weight * p
    return _probability(out, "composed flower-category prior")


def proper_losses(probability: Sequence[float], observed_codes: Sequence[int]) -> tuple[float, float]:
    """Multiclass Brier and log loss; category frequency is counted by plant."""
    p = _probability(probability, "model flower-category probabilities")
    codes = np.asarray(observed_codes)
    if (codes.ndim != 1 or codes.size == 0 or not np.issubdtype(codes.dtype, np.integer)
            or (codes < 0).any() or (codes >= len(p)).any()):
        raise ValueError("Observed Direct states must index the same fixed ontology")
    return (float(1. + p @ p - 2*p[codes].mean()),
            float(-np.log(np.maximum(p[codes], 1e-15)).mean()))


def score_focal_flora(*, observed_codes: Sequence[int],
                      weights_by_model: Mapping[str, Mapping[str, float]],
                      external_genus_profiles: Mapping[str, Sequence[float]],
                      external_blind: Sequence[float],
                      regional_blind: Sequence[float],
                      observed_focal_genera_for_ORACLE_only: Sequence[str] | None = None,
                      min_direct: int = 5) -> dict[str, float]:
    """Score ALL same-species Direct target carriers on one held-out island.

    Every non-oracle model is evaluated on the same target plant carrier set.
    Pass factual genus identities ONLY for the explicitly named oracle ceiling.
    """
    observed = np.asarray(observed_codes, dtype=int)
    if len(observed) < min_direct:
        raise ValueError("Below independently predeclared focal Direct n")
    missing = set(REQUIRED_MODELS) - set(weights_by_model)
    if missing:
        raise ValueError(f"Missing preregistered genus predictors: {sorted(missing)}")
    ref_external = _probability(external_blind, "external blind")
    ref_local = _probability(regional_blind, "regional blind")
    if len(ref_external) != len(ref_local):
        raise ValueError("Incompatible regional/external flower ontology")
    probs = {"external_blind": ref_external, "regional_blind": ref_local}
    probs.update({label: predict_flower_distribution(
        weights_by_model[label], external_genus_profiles, ref_external)
        for label in REQUIRED_MODELS})
    raw = {name: proper_losses(p, observed) for name, p in probs.items()}
    out = {f"{name}_{metric}": value
           for name, (brier, logloss) in raw.items()
           for metric, value in (("Brier", brier), ("logloss", logloss))}
    for name in REQUIRED_MODELS:
        for base in (*COMPARATORS, "regional_genus_prior"):
            if name == base:
                continue
            for j, metric in enumerate(("Brier", "logloss")):
                out[f"{name}_gain_vs_{base}_{metric}"] = raw[base][j]-raw[name][j]
    if observed_focal_genera_for_ORACLE_only is not None:
        if len(observed_focal_genera_for_ORACLE_only) != len(observed):
            raise ValueError("Explicit oracle genus and Direct carrier counts differ")
        b = l = 0.
        for g, state in zip(observed_focal_genera_for_ORACLE_only, observed):
            pred = external_genus_profiles.get(g, ref_external)
            bi, li = proper_losses(pred, [int(state)])
            b += bi
            l += li
        oracle = (b/len(observed), l/len(observed))
        for j, metric in enumerate(("Brier", "logloss")):
            out[f"oracle_known_genus_{metric}"] = oracle[j]
            for name in REQUIRED_MODELS:
                out[f"oracle_gain_vs_{name}_{metric}"] = raw[name][j]-oracle[j]
    out["n_original_direct_plant_carriers"] = len(observed)
    return out


def paired_cluster_summary(frame: pd.DataFrame, *, group: str,
                           columns: Sequence[str], nboot: int = 499,
                           min_islands: int = 30,
                           min_groups: int = 10, seed: int = 20261119) -> dict:
    """Paired island gains, equal weight per distinct block OR source group."""
    needed = {group, *columns}
    if not needed <= set(frame):
        raise ValueError("Missing required paired-scoring/group columns")
    if frame[list(needed)].isna().any().any():
        raise ValueError("Missing paired score/group is not a zero-valued effect")
    if nboot < 50:
        raise ValueError("Not enough exploratory bootstrap resamples")
    n = len(frame)
    ng = frame[group].nunique()
    if n < min_islands or ng < min_groups:
        return {"supported_for_inference": False, "n_islands": n, "n_groups": ng}
    matrix = frame.groupby(group, sort=True)[list(columns)].mean().to_numpy(dtype=float)
    if not np.isfinite(matrix).all():
        raise ValueError("Nonfinite paired score")
    rng = np.random.default_rng(seed)
    draws = matrix[rng.integers(ng, size=(nboot, ng))].mean(axis=1)
    result = {"supported_for_inference": True, "n_islands": n, "n_groups": ng,
              "resampling_unit": group, "nboot": nboot,
              "exploratory_not_multiplicity_adjusted": True}
    for j, col in enumerate(columns):
        result[col] = {"group_equal_mean": float(matrix[:, j].mean()),
                       "bootstrap_95": [float(x) for x in np.quantile(draws[:, j], [.025, .975])]}
    return result


def dual_loss_dual_group_positive(block_summary: Mapping,
                                  source_summary: Mapping,
                                  model: str, baseline: str) -> bool:
    """Never call a positive stage-2 finding if either loss/group crosses zero."""
    if not block_summary.get("supported_for_inference") or not source_summary.get("supported_for_inference"):
        return False
    for summary in (block_summary, source_summary):
        for loss in ("Brier", "logloss"):
            name = f"{model}_gain_vs_{baseline}_{loss}"
            if name not in summary or summary[name]["bootstrap_95"][0] <= 0:
                return False
    return True


def diagnostic_flower_residual_split(*, predicted_genus_weights: Mapping[str, float],
                                     observed_whole_flora_genera: Sequence[str],
                                     observed_direct_carrier_genera: Sequence[str],
                                     observed_direct_codes: Sequence[int],
                                     external_genus_profiles: Mapping[str, Sequence[float]],
                                     external_blind: Sequence[float]) -> dict[str, np.ndarray]:
    """Exact post-prediction decomposition; focal true genera used ONLY to audit.

    total observed-minus-predicted trait frequency =
      (whole flora genera - predicted genera) projected through donor flowers
      + (Direct scored carriers - whole flora genera) ascertainment component
      + (actual Direct state frequencies - donor-predicted Direct composition).
    These are signed categorical residual vectors, NOT causal selection effects.
    Their squared norms are NOT an additive explained-variance partition.
    """
    from collections import Counter
    whole = list(observed_whole_flora_genera)
    direct = list(observed_direct_carrier_genera)
    if not whole or not direct or len(direct) != len(observed_direct_codes):
        raise ValueError("Whole-flora and Direct-carrier genus audit inputs incomplete")
    blind = _probability(external_blind, "external genus-blind trait prior")
    codes = np.asarray(observed_direct_codes)
    if (codes.ndim != 1 or not np.issubdtype(codes.dtype, np.integer)
            or (codes < 0).any() or (codes >= len(blind)).any()):
        raise ValueError("Direct-carrier phenotype ontology mismatch")
    def freq(genera):
        count = Counter(genera)
        return {g: n/len(genera) for g, n in count.items()}
    p_geo = predict_flower_distribution(predicted_genus_weights, external_genus_profiles, blind)
    p_whole = predict_flower_distribution(freq(whole), external_genus_profiles, blind)
    p_direct = predict_flower_distribution(freq(direct), external_genus_profiles, blind)
    q_direct = np.bincount(codes, minlength=len(blind)).astype(float)/len(codes)
    pieces = {
        "taxonomic_genus_forecast_residual": p_whole - p_geo,
        "direct_trait_documentation_residual": p_direct - p_whole,
        "within_genus_documented_trait_residual": q_direct - p_direct,
        "observed_minus_predicted_direct_trait_frequency": q_direct - p_geo,
    }
    if not np.allclose(sum((pieces[k] for k in list(pieces)[:3]), np.zeros_like(blind)),
                       pieces["observed_minus_predicted_direct_trait_frequency"], atol=1e-10):
        raise AssertionError("Flower residual identity failed")
    return pieces


def categorical_simpson_diversity_diagnostic(
    observed_codes: Sequence[int],
    model_probabilities: Mapping[str, Sequence[float]],
    *, local_blind_model: str = "regional_blind",
) -> dict[str, float]:
    """Descriptive Direct-flower categorical diversity on matched carriers.

    Here all categories use unit between-category distance; Q=1-sum(p_k^2).
    This is categorical Gini-Simpson / nominal Rao Q, NOT a multivariate
    morphology-distance hypervolume. For observed finite n, also report the
    unbiased U-statistic n/(n-1) * Q_plugin, defined only for n>=2.
    Compare all model predictions on the identical original Direct carriers
    and label errors descriptive, never direct tests of ecological selection.
    """
    if local_blind_model not in model_probabilities:
        raise ValueError("Strong local-blind category baseline must be present")
    if not model_probabilities:
        raise ValueError("No flower category probabilities")
    p = {name: _probability(prob, name)
         for name, prob in model_probabilities.items()}
    ncat = len(p[local_blind_model])
    if any(len(v) != ncat for v in p.values()):
        raise ValueError("Flower-category ontology differs between models")
    codes = np.asarray(observed_codes)
    if (codes.ndim != 1 or len(codes) < 2
            or not np.issubdtype(codes.dtype, np.integer)
            or (codes < 0).any() or (codes >= ncat).any()):
        raise ValueError("Observed Direct categories must contain >=2 valid carriers")
    counts = np.bincount(codes, minlength=ncat).astype(float)
    n = len(codes)
    observed = float(1.0 - np.sum((counts/n)**2))
    unbiased = float((n/(n-1))*observed)
    result = {"n_direct_carriers": int(n),
              "observed_categorical_simpson_plugin": observed,
              "observed_categorical_simpson_unbiased": unbiased,
              "observed_simpson_unbiased_is_a_noisy_finite_sample_estimate": True}
    abs_err = {}
    for name, prob in p.items():
        pred = float(1.0 - prob @ prob)
        err = pred - unbiased
        result[f"{name}_predicted_categorical_simpson"] = pred
        result[f"{name}_categorical_simpson_signed_error"] = err
        abs_err[name] = abs(err)
        result[f"{name}_categorical_simpson_absolute_error"] = abs_err[name]
    for name in p:
        if name != local_blind_model:
            result[f"{name}_gain_vs_{local_blind_model}_categorical_simpson_abs_error"] = (
                abs_err[local_blind_model]-abs_err[name])
    return result
