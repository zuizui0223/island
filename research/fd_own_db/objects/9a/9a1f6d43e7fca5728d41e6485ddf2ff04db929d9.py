"""Out-of-genus prediction of GloBI DOCUMENTED five-channel reference composition.

Predicts source-reporting fractions, NOT island pollinator abundance, effective
pollen transfer or a validated bee/butterfly/bird syndrome. Original species/
genus identities stay in private input; only aggregate scores are published.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

TRAITS = ("flower_primary_color", "floral_form", "floral_symmetry")
CHANNELS = ("bombus", "non_bombus_bees", "lepidoptera",
            "flower_visiting_birds", "diptera")
REFS = tuple("refs__" + c for c in CHANNELS)
COLORS = ("blue_purple", "green_brown_inconspicuous", "multicolored_variable",
          "other_described", "red_pink", "white", "yellow_orange")
FORMS = ("bell_campanulate", "bilabiate", "brush_puff", "composite_head",
         "funnel_trumpet", "open_radial", "other_described", "papilionaceous",
         "salverform", "spurred", "tubular", "urn_urceolate")
ALLOWED = dict(zip(TRAITS, (COLORS, FORMS,
                           ("actinomorphic", "asymmetric", "zygomorphic"))))
MODELS = ("evidence_only", "simple_morphology",
          "all_marginals", "marginals_plus_joints")


def direct_complete_traits(original: pd.DataFrame) -> pd.DataFrame:
    """Retain only original high/medium species-direct single-state traits."""
    required = {"accepted_species", "trait_name", "normalized_value",
                "evidence_scope", "evidence_quality"}
    if not required.issubset(original):
        raise ValueError("Missing original trait provenance columns")
    x = original[list(required)].copy()
    x = x[x.trait_name.isin(TRAITS) &
          x.evidence_scope.eq("species_direct") &
          x.evidence_quality.isin(("high", "medium"))]
    if x.duplicated(["accepted_species", "trait_name"]).any():
        raise ValueError("Duplicate species-direct trait evidence")
    x = x.dropna(subset=["accepted_species", "normalized_value"])
    x = x[~x.normalized_value.astype(str).str.contains("|", regex=False)]
    x = x[x.apply(lambda row: row.normalized_value in
                  ALLOWED[row.trait_name], axis=1)]
    a = x.pivot(index="accepted_species", columns="trait_name",
                values="normalized_value")
    for trait in TRAITS:
        if trait not in a:
            a[trait] = pd.NA
    a = a.dropna(subset=list(TRAITS)).reset_index()
    a["genus"] = a.accepted_species.astype(str).str.split().str[0]
    if a.accepted_species.duplicated().any():
        raise ValueError("Repeated original accepted species")
    return a[["accepted_species", "genus", *TRAITS]]


def prepare_genus_frame(species: pd.DataFrame, source: pd.DataFrame,
                        *, min_independent_references: int = 3):
    """Build PRIVATE genus-level marginal/joint morphology and channel fractions."""
    required = {"genus", "n_independent_references",
                "n_matched_plant_species", "n_reference_channel_units", *REFS}
    if not required.issubset(source):
        raise ValueError("Missing frozen GloBI genus source evidence columns")
    if min_independent_references < 3:
        raise ValueError("Independent-reference threshold too permissive")
    if source.genus.isna().any() or source.genus.duplicated().any():
        raise ValueError("Missing or repeated GloBI genus identity")
    g = source.copy()
    for col in required - {"genus"}:
        g[col] = pd.to_numeric(g[col], errors="coerce")
        if not np.isfinite(g[col].to_numpy()).all() or (g[col] < 0).any():
            raise ValueError("Invalid GloBI reference or support count")
    if not np.allclose(g[list(REFS)].sum(axis=1),
                       g.n_reference_channel_units):
        raise ValueError("Five GloBI channel units do not sum to total")
    if species.empty or species.accepted_species.duplicated().any():
        raise ValueError("Empty or duplicated Direct source species")
    if not set(TRAITS).issubset(species) or species[list(TRAITS)].isna().any().any():
        raise ValueError("Missing directly measured complete morphology")
    g = g.loc[g.n_independent_references >= min_independent_references]
    s = species.loc[species.genus.isin(g.genus)]
    g = g.loc[g.genus.isin(s.genus)].set_index("genus")
    if len(g) < 50:
        raise ValueError("At least 50 distinct GloBI genera required")
    counts = s.genus.value_counts().reindex(g.index)
    g["log_refs"] = np.log1p(g.n_independent_references)
    g["log_plant"] = np.log1p(g.n_matched_plant_species)
    g["log_n_morph"] = np.log1p(counts)
    margins = []
    for trait in TRAITS:
        f = pd.crosstab(s.genus, s[trait])
        f = f.div(s.groupby("genus").size(), axis=0)
        f.columns = [trait + "::" + str(c) for c in f.columns]
        margins.append(f)
    joint = pd.crosstab(s.genus,
                        s[list(TRAITS)].astype(str).agg("|||".join, axis=1))
    joint = joint.div(s.groupby("genus").size(), axis=0)
    joint.columns = ["joint::" + str(c) for c in joint.columns]
    M = pd.concat(margins + [joint], axis=1).reindex(g.index).fillna(0)
    simple = pd.DataFrame({
        "prop_actinomorphic": M.get("floral_symmetry::actinomorphic",
                                   pd.Series(0., index=g.index)),
        "prop_open_radial": M.get("floral_form::open_radial",
                                pd.Series(0., index=g.index)),
        "prop_inconspicuous": M.get(
            "flower_primary_color::green_brown_inconspicuous",
            pd.Series(0., index=g.index)),
    }, index=g.index)
    base = ["log_refs", "log_plant", "log_n_morph"]
    X = pd.concat([g[base], simple, M], axis=1)
    X = X.loc[:, ~X.columns.duplicated()]
    maps = {
        "evidence_only": base,
        "simple_morphology": base + list(simple),
        "all_marginals": base + [v for v in M if not v.startswith("joint::")],
        "marginals_plus_joints": base + list(M),
    }
    y = g[list(REFS)].to_numpy(float)
    if not (y.sum(axis=1) > 0).all():
        raise ValueError("Source genus with zero documented channel evidence")
    y = y / y.sum(axis=1)[:, None]
    return X, maps, y, {"n_source_direct_species": len(s), "n_genera": len(g)}


def simplex_project(arr: np.ndarray) -> np.ndarray:
    """Euclidean projection onto a valid nonnegative, unit-sum composition."""
    u = np.sort(arr, axis=1)[:, ::-1]
    c = np.cumsum(u, axis=1) - 1.
    j = np.arange(1, arr.shape[1] + 1)
    rho = (u - c/j > 0).sum(axis=1)
    theta = c[np.arange(len(arr)), rho-1] / rho
    return np.maximum(arr - theta[:, None], 0.)


def _folds(n: int, seed: int):
    if n < 50:
        raise ValueError("Too few independent test genera")
    idx = np.arange(n)
    np.random.RandomState(seed).shuffle(idx)
    size = np.full(5, n//5, dtype=int)
    size[:n % 5] += 1
    start = 0
    for k in size:
        test = idx[start:start+k]
        start += k
        train = np.setdiff1d(np.arange(n), test, assume_unique=True)
        yield train, test


def _predict(x: np.ndarray, y: np.ndarray, split, alpha: float) -> np.ndarray:
    result = np.empty_like(y)
    for train, test in split:
        a, b = x[train], x[test]
        mu, sd = a.mean(axis=0), a.std(axis=0)
        sd[sd < 1e-12] = 1
        a, b = (a-mu)/sd, (b-mu)/sd
        y0 = y[train].mean(axis=0)
        beta = np.linalg.solve(a.T@a + alpha*np.eye(a.shape[1]),
                               a.T@(y[train]-y0))
        result[test] = y0 + b@beta
    return simplex_project(result)


def heldout_channel_scores(X: pd.DataFrame, models: dict, Y: np.ndarray,
                           *, seeds=(20261111, 20261112, 20261113),
                           ridge_penalties=(10., 100.)) -> pd.DataFrame:
    """Report aggregate fold-separated CV, one test unit per plant genus."""
    if tuple(models) != MODELS or Y.shape != (len(X), len(CHANNELS)):
        raise ValueError("Invalid five-channel test model family")
    if not np.isfinite(Y).all() or (Y < 0).any() or not np.allclose(Y.sum(axis=1), 1):
        raise ValueError("Invalid reference-share target")
    out = []
    for seed in seeds:
        split = list(_folds(len(Y), int(seed)))
        for alpha in ridge_penalties:
            if not np.isfinite(alpha) or alpha <= 0:
                raise ValueError("Invalid ridge penalty")
            fitted = {name: _predict(X[cols].to_numpy(float), Y, split, alpha)
                      for name, cols in models.items()}
            base = (Y - fitted["evidence_only"])**2
            for name, pred in fitted.items():
                err = (Y-pred)**2
                for j, ch in enumerate((*CHANNELS, "all_five_channels")):
                    b = base[:, j] if j < 5 else base.mean(axis=1)
                    e = err[:, j] if j < 5 else err.mean(axis=1)
                    out.append({
                        "model": name, "channel": ch, "seed": seed,
                        "ridge_alpha": float(alpha), "n_genera": len(Y),
                        "outofgenus_MSE": float(e.mean()),
                        "outofgenus_RMSE": float(np.sqrt(e.mean())),
                        "MSE_gain_vs_reference_only": float((b-e).mean()),
                        "local_guild_identity_observed": False,
                    })
    return pd.DataFrame.from_records(out)
