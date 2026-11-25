"""Family-leave-one-out crossvalidation of GloBI DOCUMENTED channel shares.

Source-rights-preserving: input genus/family identities never appear in
returned group-level scores. A documented visit is not pollen effectiveness.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import rel_entr


def _simplex(q: np.ndarray) -> np.ndarray:
    u = np.sort(q, axis=1)[:, ::-1]
    c = np.cumsum(u, axis=1) - 1
    j = np.arange(1, q.shape[1] + 1)
    rho = (u - c / j > 0).sum(axis=1)
    th = c[np.arange(len(q)), rho - 1] / rho
    return np.maximum(q - th[:, None], 0)


def leave_one_family_out_scores(
    features: pd.DataFrame, families: pd.Series,
    documented_five_channels: np.ndarray, models: dict[str, list[str]], *,
    alpha: float = 100., bootstrap_draws: int = 4999,
    seed: int = 20261111,
) -> pd.DataFrame:
    """Fit on all other families; score each genus exactly once.

    Family-equal bootstrap is conditional on the sampled family/record sources.
    It is not proof of phylogenetically independent evolutionary replication.
    """
    if not features.index.is_unique or features.index.hasnans or len(features) < 50:
        raise ValueError("Requires >=50 distinct botanical genera")
    if not families.index.equals(features.index) or families.isna().any():
        raise ValueError("Incomplete or misaligned family taxonomy")
    labels = families.astype(str)
    if (labels.str.lower().isin(["unknown", "unresolved", "nan"]).any()
            or labels.str.startswith("AMBIG:").any()
            or labels.str.strip().eq("").any()):
        raise ValueError("Unresolved or ambiguous family identities")
    fam = labels.to_numpy()
    unique = np.unique(fam)
    if len(unique) < 10 or bootstrap_draws < 199 or not np.isfinite(alpha) or alpha <= 0:
        raise ValueError("Insufficient families or uncertainty support")
    y = np.asarray(documented_five_channels, dtype=float)
    if (y.shape != (len(features), 5) or not np.isfinite(y).all()
            or (y < 0).any() or not np.allclose(y.sum(axis=1), 1)):
        raise ValueError("Invalid GloBI documented five-channel shares")
    if "evidence_only" not in models:
        raise ValueError("Must include frozen reference-effort baseline")
    for cols in models.values():
        if not cols or len(cols) != len(set(cols)) or not set(cols).issubset(features):
            raise ValueError("Invalid fixed model feature family")

    # Use the same ridge normal equations and TRAIN-fold-only z-scoring
    # as fd's existing family-disjoint implementation. No sklearn runtime
    # dependency is needed by the lightweight public CI.
    errors = {}
    for name, cols in models.items():
        X = features[cols].to_numpy(float)
        if not np.isfinite(X).all():
            raise ValueError("Missing Direct trait or reference-effort predictor")
        pred = np.empty_like(y)
        for family in unique:
            test = np.flatnonzero(fam == family)
            train = np.flatnonzero(fam != family)
            mu = X[train].mean(axis=0)
            sd = X[train].std(axis=0)
            sd[sd < 1e-12] = 1
            a_train = (X[train] - mu) / sd
            a_test = (X[test] - mu) / sd
            intercept = y[train].mean(axis=0)
            beta = np.linalg.solve(
                a_train.T @ a_train + alpha * np.eye(a_train.shape[1]),
                a_train.T @ (y[train] - intercept),
            )
            pred[test] = intercept + a_test @ beta
        pred = _simplex(pred)
        mse = np.square(pred - y).mean(axis=1)
        mae = np.abs(pred - y).mean(axis=1)
        mid = np.maximum((pred + y) / 2, 1e-12)
        jsd = (.5 * np.sum(rel_entr(y, mid), axis=1)
               + .5 * np.sum(rel_entr(pred, mid), axis=1))
        errors[name] = (mse, mae, jsd)
    rng = np.random.default_rng(seed)
    results = []
    for name in models:
        mse, mae, jsd = errors[name]
        gain = errors["evidence_only"][0] - mse
        family_gain = pd.DataFrame({"family": fam, "gain": gain}).groupby("family").gain.mean().to_numpy()
        sample = rng.integers(len(family_gain), size=(bootstrap_draws, len(family_gain)))
        ci = np.quantile(family_gain[sample].mean(axis=1), [.025, .975])
        results.append({
            "model": name, "n_genera": len(features), "n_families": len(unique),
            "ridge_alpha": float(alpha), "n_features": len(models[name]),
            "genus_equal_MSE": float(mse.mean()), "genus_equal_MAE": float(mae.mean()),
            "genus_equal_JSD": float(jsd.mean()),
            "genus_equal_MSE_gain": float(gain.mean()),
            "family_equal_MSE_gain": float(family_gain.mean()),
            "n_families_positive": int((family_gain > 0).sum()),
            "family_equal_bootstrap_lo": float(ci[0]),
            "family_equal_bootstrap_hi": float(ci[1]),
            "actual_island_pollination_tested": False,
            "original_provider_heldout_tested": False,
        })
    return pd.DataFrame.from_records(results)
