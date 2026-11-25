"""Falsification controls for global GloBI genus breadth and floral morphology.

The GloBI response is DOCUMENTED global interaction-channel breadth, not
local pollinator occurrence or effective service. Both helpers are post-hoc
diagnostics that keep genus identities, botanical traits, and source refs
in memory; outputs are nonidentifying summaries only.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.sparse import issparse

BASE = ("log_refs", "log_matched_plant_species", "log_n_morph")
FOLD_COUNT = 5


def effort_stratified_genus_permutation(
    island_to_genus_weights,
    genus_breadth, genus_ref_counts, genus_plant_counts,
    trait_association, nuisance_design, *,
    repeats: int = 1999, seed: int = 20261111
) -> dict:
    """Shuffle GloBI genus breadth inside reference-effort-matched strata.

    Island×genus weight rows sum to one over source-backed known genera.
    Nuisance must already contain a constant, area/climate/region effects,
    floral Rao Q, trait coverage and GloBI evidence effort terms.
    This tests LABEL ASSOCIATION, not phylogenetic independence.
    """
    W = island_to_genus_weights
    if repeats < 199 or not issparse(W):
        raise ValueError("Insufficient randomization or non-sparse island-genus matrix")
    b, refs, plants = [np.asarray(x, dtype=float) for x in
                       (genus_breadth, genus_ref_counts, genus_plant_counts)]
    y = np.asarray(trait_association, dtype=float)
    Q = np.asarray(nuisance_design, dtype=float)
    n, g = W.shape
    if (n < 40 or g < 50 or len(b) != g or len(refs) != g or len(plants) != g
        or len(y) != n or Q.ndim != 2 or len(Q) != n):
        raise ValueError("Insufficient or misaligned independent groups")
    if (not all(np.isfinite(v).all() for v in (b, refs, plants, y, Q))
        or (refs < 1).any() or (plants < 1).any()
        or not np.allclose(np.asarray(W.sum(axis=1)).ravel(), 1)):
        raise ValueError("Nonfinite, invalid reference evidence or unmatched genus weights")
    if np.linalg.matrix_rank(Q) != Q.shape[1]:
        raise ValueError("Rank-deficient nuisance design")
    proj = np.linalg.pinv(Q)
    resid_y = y - Q @ (proj @ y)

    def slope(v):
        a = np.asarray(W @ v).ravel()
        r = a-Q @ (proj @ a)
        den = float(r @ r)
        if den <= 1e-14:
            raise ValueError("GloBI breadth not separately identifiable")
        return float((r @ resid_y)/den)

    # Six quantile bins of reference sample depth, by three source-plant
    # support classes. Do NOT choose bins based on floral or model outcome.
    logref = np.log1p(refs)
    bins = pd.qcut(logref, q=6, duplicates="drop", labels=False).astype(int)
    plant_bins = np.where(plants == 1, 0, np.where(plants <= 3, 1, 2))
    strata = bins*3 + plant_bins
    groups = [np.flatnonzero(strata == q) for q in np.unique(strata)
              if (strata == q).sum() >= 2]
    rng = np.random.default_rng(seed)
    observed = slope(b)
    null = np.empty(repeats)
    for k in range(repeats):
        shuffled = b.copy()
        for ix in groups:
            shuffled[ix] = rng.permutation(shuffled[ix])
        null[k] = slope(shuffled)
    return {
        "n_islands": int(n), "n_globi_genera": int(g),
        "n_effort_strata_with_shuffle": int(len(groups)),
        "n_shuffle_draws": int(repeats), "observed_partial_slope": observed,
        "null_median": float(np.median(null)),
        "null_025": float(np.quantile(null, .025)),
        "null_975": float(np.quantile(null, .975)),
        "two_sided_posthoc_tail_fraction":
            float((1 + (np.abs(null) >= abs(observed)).sum())/(repeats+1)),
        "biological_phylogenetic_independence_tested": False,
    }


def genus_heldout_guild_breadth(
    genus_frame: pd.DataFrame, *,
    base_features: list[str],
    marginal_features: list[str],
    pairing_features: list[str],
    folds: int = 5, ridge_alpha: float = 10., seed: int = 20261111
) -> pd.DataFrame:
    """Forecast documented GloBI breadth for entirely unseen plant genera.

    Each genus occurs ONCE; the fold may not share the same genus. Evidence-
    only, marginal traits, and full trait-pairing models share exact test rows.
    Uses train-only standardization and fixed alpha; no post-outcome tuning.
    """
    if folds != 5 or ridge_alpha <= 0 or len(genus_frame) < 50:
        raise ValueError("Require 5 folds, positive ridge shrinkage and >=50 genera")
    required = {"genus", "effective_channel_number"} | set(
        base_features+marginal_features+pairing_features)
    if not required.issubset(genus_frame):
        raise ValueError("Missing predictor/genus/interaction response")
    if genus_frame.genus.isna().any() or genus_frame.genus.duplicated().any():
        raise ValueError("Each recorded genus must have exactly one row")
    cols = {
        "evidence_only": base_features,
        "morph_marginal": base_features+marginal_features,
        "morph_marginal_plus_joint": base_features+marginal_features+pairing_features,
    }
    y = genus_frame.effective_channel_number.to_numpy(float)
    if not np.isfinite(y).all() or (y <= 0).any():
        raise ValueError("Unknown/invalid source-based guild breadth")
    n = len(y)
    # Same shuffled 5 test folds for every model, equivalent to sklearn KFold
    # shuffle=True with legacy seeded RandomState, but no sklearn dependency.
    perm = np.arange(n)
    np.random.RandomState(seed).shuffle(perm)
    sizes = np.full(folds, n//folds, dtype=int); sizes[:n % folds] += 1
    indices=[];start=0
    for size in sizes:
        indices.append(perm[start:start+size]);start+=size
    records=[]
    errors={}
    for model, fields in cols.items():
        if len(fields) != len(set(fields)) or len(fields) < 1:
            raise ValueError("Duplicate or empty feature family")
        X = genus_frame[fields].to_numpy(float)
        if not np.isfinite(X).all():
            raise ValueError("Unknown source effort / floral category, never impute 0")
        predicted=np.empty(n)
        for test in indices:
            train=np.setdiff1d(np.arange(n), test, assume_unique=True)
            A=X[train]; B=X[test]
            means=A.mean(axis=0);std=A.std(axis=0);std[std < 1e-12] = 1
            A=(A-means)/std;B=(B-means)/std
            ym=y[train].mean()
            coeff=np.linalg.solve(
                A.T@A + ridge_alpha*np.eye(A.shape[1]),
                A.T@(y[train]-ym),
            )
            predicted[test]=ym + B@coeff
        err=(y-predicted)**2
        errors[model]=err
        records.append({
            "model":model, "n_heldout_genera":n,
            "n_features":len(fields), "folds":folds, "fixed_ridge_alpha":ridge_alpha,
            "RMSE":float(np.sqrt(err.mean())),
            "MAE":float(np.abs(y-predicted).mean()),
            "out_of_genus_R2":float(
                1 - err.sum()/np.sum((y-y.mean())**2)),
        })
    for record in records:
        model = record["model"]
        record["MSE_gain_over_effort_only"] = float(
            (errors["evidence_only"]-errors[model]).mean())
    return pd.DataFrame(records)
