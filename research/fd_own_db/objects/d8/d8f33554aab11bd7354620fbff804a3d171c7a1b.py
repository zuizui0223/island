"""Family-disjoint GloBI five-channel documentation composition validation.

The input genus-to-family mapping is an independently sourced taxonomy table.
No original taxon, genus, family, GloBI provider, or island IDs are exported.
GloBI documented channel fractions are NOT observed local pollinator services.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .globi_channel_composition_cv import simplex_project

MODEL_ORDER = ("evidence_only", "simple_morphology", "all_marginals",
               "marginals_plus_joints")


def _validate(x: pd.DataFrame, models: dict[str, list[str]], y: np.ndarray,
              families: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    if tuple(models) != MODEL_ORDER:
        raise ValueError("Fixed reference-only through joint morphology models required")
    if not x.index.is_unique or x.index.hasnans or len(x) < 50:
        raise ValueError("One unique genus per test row, >=50 genera")
    if not families.index.equals(x.index) or families.isna().any():
        raise ValueError("Missing or misordered botanical family mapping")
    labels = families.astype(str).str.strip()
    if (labels.eq("").any() or labels.str.lower().isin(("unknown", "unresolved", "nan")).any()
            or labels.str.startswith("AMBIG:").any()):
        raise ValueError("Unresolved or conflicting family identification cannot pass a family holdout")
    f = labels.to_numpy()
    if len(np.unique(f)) < 10:
        raise ValueError("Insufficient botanical family blocks")
    arr = np.asarray(y, float)
    if arr.shape != (len(x), 5) or not np.isfinite(arr).all() or (arr < 0).any():
        raise ValueError("Five observed nonnegative source-reference channel fractions required")
    if not np.allclose(arr.sum(axis=1), 1):
        raise ValueError("GloBI documentation shares must sum to one")
    for name, cols in models.items():
        if len(cols) != len(set(cols)) or not cols or not set(cols).issubset(x):
            raise ValueError("Invalid fixed model feature set")
        if not np.isfinite(x[cols].to_numpy(float)).all():
            raise ValueError("Unknown direct morphology or source effort cannot be zero-imputed")
    return f, arr


def family_folds(families: np.ndarray, *, seed: int, folds: int = 5):
    """Outcome-blind balanced assignment of WHOLE botanical families to test folds."""
    if folds != 5:
        raise ValueError("Frozen five test-fold design")
    f = np.asarray(families, dtype=str)
    counts = pd.Series(f).value_counts()
    if len(counts) < folds * 2:
        raise ValueError("Too few distinct families for blocked CV")
    rng = np.random.default_rng(seed)
    ordered = sorted(counts.index, key=lambda name: (-int(counts[name]), rng.random()))
    size = np.zeros(folds, int)
    test_groups = [[] for _ in range(folds)]
    for family in ordered:
        candidates = np.flatnonzero(size == size.min())
        k = int(rng.choice(candidates))
        ids = np.flatnonzero(f == family)
        test_groups[k].extend(ids)
        size[k] += len(ids)
    out = []
    all_rows = np.arange(len(f))
    for test_rows in test_groups:
        test = np.asarray(test_rows, int)
        train = np.setdiff1d(all_rows, test)
        if not len(test) or set(f[train]) & set(f[test]):
            raise AssertionError("Overlapping botanical families between train and test")
        out.append((train, test))
    if not np.array_equal(np.sort(np.concatenate([x[1] for x in out])), all_rows):
        raise AssertionError("Folds do not partition evaluated genera")
    return out


def _predict(x: np.ndarray, y: np.ndarray, groups, alpha: float):
    result = np.empty_like(y, dtype=float)
    for train, test in groups:
        a, b = x[train], x[test]
        mu, sd = a.mean(axis=0), a.std(axis=0)
        sd[sd < 1e-12] = 1
        a, b = (a - mu) / sd, (b - mu) / sd
        ym = y[train].mean(axis=0)
        w = np.linalg.solve(a.T @ a + alpha * np.eye(a.shape[1]),
                            a.T @ (y[train] - ym))
        result[test] = ym + b @ w
    return simplex_project(result)


def family_disjoint_scores(x: pd.DataFrame, models: dict[str, list[str]],
                           y: np.ndarray, families: pd.Series, *,
                           seeds=(20261111, 20261112, 20261113, 20261114, 20261115),
                           alphas=(10.0, 100.0), bootstrap_draws: int = 999,
                           seed_boot: int = 20261111) -> pd.DataFrame:
    """Genus-weighted predictive scores and conditional family-block intervals.

    Shared evolutionary ancestors and GloBI references beyond family labels
    are not independent; these are NOT confirmatory evolutionary p-values.
    """
    if bootstrap_draws < 199 or len(set(seeds)) != len(seeds):
        raise ValueError("Insufficient family bootstrap / repeated fold diagnostics")
    fam, target = _validate(x, models, y, families)
    uniques = np.unique(fam)
    family_rows = [np.flatnonzero(fam == group) for group in uniques]
    rng = np.random.default_rng(seed_boot)
    output = []
    for seed in seeds:
        split = family_folds(fam, seed=int(seed))
        for alpha in alphas:
            if not np.isfinite(alpha) or alpha <= 0:
                raise ValueError("Ridge regularization must be positive")
            pred = {label: _predict(x[cols].to_numpy(float), target, split, alpha)
                    for label, cols in models.items()}
            reference = np.square(target - pred["evidence_only"]).mean(axis=1)
            for name in MODEL_ORDER:
                error = np.square(target - pred[name]).mean(axis=1)
                gain = reference - error
                fam_sum = np.array([gain[idx].sum() for idx in family_rows])
                fam_n = np.array([len(idx) for idx in family_rows])
                fam_mean = fam_sum / fam_n
                boots = rng.integers(len(uniques), size=(bootstrap_draws, len(uniques)))
                # Two different estimands: prior interval weights sampled
                # families by their GENUS counts; the new interval gives each
                # family one vote. Do not interchange them.
                resamples = fam_sum[boots].sum(axis=1) / fam_n[boots].sum(axis=1)
                family_equal_resamples = fam_mean[boots].mean(axis=1)
                output.append({
                    "model": name, "fold_seed": int(seed), "ridge_alpha": float(alpha),
                    "n_genera": len(x), "n_family_blocks": len(uniques),
                    "n_features": len(models[name]),
                    "MSE_five_documented_channels": float(error.mean()),
                    "MSE_gain_over_reference_only": float(gain.mean()),
                    "family_equal_gain": float(fam_mean.mean()),
                    "family_equal_bootstrap_lower": float(np.quantile(family_equal_resamples, .025)),
                    "family_equal_bootstrap_upper": float(np.quantile(family_equal_resamples, .975)),
                    "n_positive_families": int((fam_mean > 0).sum()),
                    "conditional_family_bootstrap_lower": float(np.quantile(resamples, .025)),
                    "conditional_family_bootstrap_upper": float(np.quantile(resamples, .975)),
                    "no_actual_island_pollination_assessed": True,
                })
    return pd.DataFrame(output)


def within_family_morphology_null(
        x: pd.DataFrame, y: np.ndarray, families: pd.Series,
        *, base_features: list[str], morph_features: list[str],
        seed_fold: int = 20261111, alpha: float = 100.,
        permutations: int = 999, seed: int = 20261111) -> dict:
    """Within-family shuffle of the intact three floral-morphology features."""
    if permutations < 99 or len(morph_features) != 3:
        raise ValueError("Requires >=99 within-family triplet shuffles")
    cols = base_features + morph_features
    model_stub = {
        "evidence_only": base_features,
        "simple_morphology": cols,
        "all_marginals": cols,
        "marginals_plus_joints": cols,
    }
    fam, target = _validate(x, model_stub, y, families)
    splits = family_folds(fam, seed=seed_fold)
    base = x[base_features].to_numpy(float)
    morph = x[morph_features].to_numpy(float)
    ref = np.square(target-_predict(base, target, splits, alpha)).mean(axis=1)

    def gain_of(m):
        p = _predict(np.column_stack([base, m]), target, splits, alpha)
        return float(np.mean(ref - np.square(target - p).mean(axis=1)))

    observed = gain_of(morph)
    groups = [np.flatnonzero(fam == group) for group in np.unique(fam)]
    groups = [idx for idx in groups if len(idx) > 1]
    rng = np.random.default_rng(seed)
    null = np.empty(permutations)
    for k in range(permutations):
        shuffle = morph.copy()
        for idx in groups:
            shuffle[idx] = morph[rng.permutation(idx)]
        null[k] = gain_of(shuffle)
    return {
        "n_genera": len(x), "n_families": len(np.unique(fam)),
        "n_genera_in_shufflable_families": int(sum(map(len, groups))),
        "observed_MSE_gain": observed, "n_within_family_shuffles": permutations,
        "null_median": float(np.median(null)),
        "null_q025": float(np.quantile(null, .025)),
        "null_q975": float(np.quantile(null, .975)),
        "posthoc_upper_tail_fraction": float((1+(null >= observed).sum())/(permutations+1)),
        "actual_island_visitors_or_pollination_effectiveness_tested": False,
    }


def coarse_morphology_ablation(
    x: pd.DataFrame, y: np.ndarray, families: pd.Series, *,
    base_features: list[str], symmetry: str, form: str, colour: str,
    seeds=(20261111, 20261112, 20261113, 20261114, 20261115),
    alpha: float = 100.,
) -> pd.DataFrame:
    """Outcome-blind whole-family ablation of three DIRECT structural/signal axes.

    All feature choices and test genus rows stay fixed across comparisons.
    The output reports documentation composition scores, not local guilds.
    """
    if len(set(base_features+[symmetry,form,colour])) != len(base_features)+3:
        raise ValueError("Duplicate or overlapping evidence/morphology predictors")
    designs = {
        "evidence_only": [],
        "symmetry_only": [symmetry],
        "form_only": [form],
        "colour_only": [colour],
        "symmetry_form": [symmetry,form],
        "symmetry_colour": [symmetry,colour],
        "form_colour": [form,colour],
        "all_three": [symmetry,form,colour],
    }
    complete = base_features + designs["all_three"]
    stub = {name:complete for name in MODEL_ORDER}
    fam, target = _validate(x, stub, y, families)
    if not np.isfinite(alpha) or alpha <= 0:
        raise ValueError("Regularization must be finite and positive")
    result = []
    for seed in seeds:
        split = family_folds(fam, seed=int(seed))
        prediction = {
            name: _predict(x[base_features+features].to_numpy(float),
                           target, split, alpha)
            for name, features in designs.items()
        }
        loss = {
            name: np.square(target-p).mean(axis=1)
            for name, p in prediction.items()
        }
        for name, features in designs.items():
            gain = loss["evidence_only"]-loss[name]
            result.append({
                "fold_seed":int(seed), "ridge_alpha":float(alpha),
                "model":name, "n_genera":len(x),
                "n_family_blocks":len(np.unique(fam)),
                "n_morph_predictors":len(features),
                "five_channel_MSE":float(loss[name].mean()),
                "MSE_gain_over_evidence_only":float(gain.mean()),
                "no_island_visitor_or_pollen_service_observed":True,
            })
    return pd.DataFrame(result)
