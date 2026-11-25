"""Out-of-species binary-trait genus predictions for island floral-function audits.

Inputs/returned species rows are PRIVATE; publish only explicitly grouped diagnostics.
This predicts recorded trait states, not pollinator performance or evolutionary history.
"""
from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd
from scipy.special import erf

FIELDS = ("taxon_name", "genus", "state")


def name_heldout_genus_probabilities(
    data: pd.DataFrame, *, tau: float = 0.05, folds: int = 5
) -> pd.DataFrame:
    """Make one global, reproducible species-name fold split.

    Exclude the ENTIRE test fold from each candidate genus and global donor pool.
    A species appearing on many islands still has one fold and one prediction.
    Result has private taxon/genus identity columns and MUST NOT be published.
    """
    missing = set(FIELDS) - set(data.columns)
    if missing:
        raise ValueError(f"Missing fields: {sorted(missing)}")
    if not (np.isfinite(tau) and tau > 0 and isinstance(folds, int) and 2 <= folds <= 32):
        raise ValueError("tau must be positive and folds must be in 2..32")
    d = data[list(FIELDS)].copy()
    for col in ("taxon_name", "genus"):
        if d[col].isna().any() or d[col].astype(str).str.strip().eq("").any():
            raise ValueError(f"Missing or blank {col}")
        d[col] = d[col].astype(str).str.strip()
    if d.taxon_name.duplicated().any():
        raise ValueError("A taxon_name must have one row across ALL islands")
    d["state"] = pd.to_numeric(d.state, errors="coerce")
    if d.state.isna().any() or not d.state.isin((0, 1)).all():
        raise ValueError("state must be binary 0/1, never imputed missing")
    if len(d) < folds * 2:
        raise ValueError("Too few species for independent folds")
    d["fold"] = d.taxon_name.map(
        lambda name: hashlib.sha256(name.encode("utf-8")).digest()[0] % folds
    )
    by_genus = d.groupby("genus").state.agg(["count", "sum"])
    by_genus.columns = ["all_n", "all_k"]
    by_fold = d.groupby(["genus", "fold"]).state.agg(["count", "sum"])
    by_fold.columns = ["fold_n", "fold_k"]
    d = d.join(by_genus, on="genus").join(by_fold, on=["genus", "fold"])
    global_fold = d.groupby("fold").state.agg(["count", "sum"])
    global_n = len(d) - d["fold"].map(global_fold["count"]).to_numpy()
    global_k = d["state"].sum() - d["fold"].map(global_fold["sum"]).to_numpy()
    if (global_n == 0).any():
        raise ValueError("A fold lacks training donors")
    d["p_global"] = global_k / global_n
    d["n_donors"] = d["all_n"] - d["fold_n"]
    d["k_donors"] = d["all_k"] - d["fold_k"]
    d["p_genus"] = (
        d["k_donors"] + tau * d["p_global"]
    ) / (d["n_donors"] + tau)
    return d[["taxon_name", "genus", "state", "fold",
              "n_donors", "p_global", "p_genus"]].reset_index(drop=True)


def score_private_predictions(
    predictions: pd.DataFrame, *, minimum_species: int = 50
) -> dict:
    """Return only name-aggregated proper predictive metrics, no identities."""
    if minimum_species < 50:
        raise ValueError("Privacy gate requires >=50 distinct taxa")
    if len(predictions) < minimum_species:
        raise ValueError("Not enough independent taxon names for report")
    if predictions.taxon_name.duplicated().any():
        raise ValueError("Prediction is not a species-level holdout")
    required = {"state", "n_donors", "p_global", "p_genus"}
    if not required.issubset(predictions):
        raise ValueError("Missing prediction columns")
    y = predictions.state.to_numpy(float)
    a = predictions.p_global.to_numpy(float)
    b = predictions.p_genus.to_numpy(float)
    if not (np.isfinite(y).all() and np.isfinite(a).all() and
            np.isfinite(b).all() and np.isin(y, [0, 1]).all() and
            ((a >= 0) & (a <= 1)).all() and ((b >= 0) & (b <= 1)).all()):
        raise ValueError("Invalid scored outcome/probability")
    def brier(p):
        return float(np.mean((p - y) ** 2))
    def log_loss(p):
        p = np.clip(p, 1e-8, 1 - 1e-8)
        return float(np.mean(-y * np.log(p) - (1-y) * np.log1p(-p)))
    return {
        "n_species": int(len(y)),
        "donor_coverage": float(np.mean(predictions.n_donors.to_numpy() > 0)),
        "global_brier": brier(a),
        "genus_brier": brier(b),
        "brier_improvement": brier(a)-brier(b),
        "global_logloss": log_loss(a),
        "genus_logloss": log_loss(b),
        "logloss_improvement": log_loss(a)-log_loss(b),
    }


def expected_pair_gap_normal(p_a, p_b) -> float:
    """Normal approximation to E|prevalence_A - prevalence_B|.

    The two evaluated islands must contain DISJOINT species identities;
    the intended use is the genus-exclusive cohort. q's derive from other
    species (never the evaluated name itself). This approximation is NOT
    a source-/phylogeny-conditioned permutation.
    """
    a, b = np.asarray(p_a, dtype=float), np.asarray(p_b, dtype=float)
    if (a.ndim != 1 or b.ndim != 1 or len(a) == 0 or len(b) == 0 or
        not np.isfinite(a).all() or not np.isfinite(b).all() or
        ((a < 0) | (a > 1)).any() or ((b < 0) | (b > 1)).any()):
        raise ValueError("Expected nonempty 1-D valid probabilities")
    mean = float(a.mean() - b.mean())
    var = float(np.sum(a * (1 - a)) / (len(a)**2) +
                np.sum(b * (1 - b)) / (len(b)**2))
    if var <= 1e-24:
        return abs(mean)
    sigma = np.sqrt(var)
    z = abs(mean) / sigma
    return float(sigma * np.sqrt(2 / np.pi) * np.exp(-z*z/2)
                 + abs(mean) * erf(z / np.sqrt(2)))
