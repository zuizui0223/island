"""Same-cohort comparison of documented GloBI visitor shares before and after
deleting the largest original archive-level provider group.

No true island visitor census or effective pollen transfer is estimated.
All taxon, genus, source and provider IDs remain local to the runner.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

TRAITS = ("flower_primary_color", "floral_form", "floral_symmetry")
CHANNELS = ("bombus", "non_bombus_bees", "lepidoptera", "flower_visiting_birds", "diptera")


def _project(x):
    u = np.sort(x, axis=1)[:, ::-1]
    css = np.cumsum(u, axis=1) - 1
    denom = np.arange(1, x.shape[1] + 1)
    rho = np.sum(u - css / denom > 0, axis=1)
    theta = css[np.arange(len(x)), rho-1] / rho
    return np.maximum(x - theta[:, None], 0)


def direct_morphology(ledger: pd.DataFrame) -> pd.DataFrame:
    needed = {"accepted_species", "trait_name", "normalized_value",
              "evidence_scope", "evidence_quality"}
    if not needed.issubset(ledger):
        raise ValueError("Missing original direct trait evidence fields")
    d = ledger.loc[ledger.trait_name.isin(TRAITS) &
                   ledger.evidence_scope.eq("species_direct") &
                   ledger.evidence_quality.isin(("high", "medium"))].copy()
    if d.duplicated(["accepted_species", "trait_name"]).any():
        raise ValueError("Duplicated original species x trait evidence")
    d = d.dropna(subset=["accepted_species", "normalized_value"])
    d = d[~d.normalized_value.astype(str).str.contains("|", regex=False)]
    p = d.pivot(index="accepted_species", columns="trait_name",
                values="normalized_value")
    if not set(TRAITS).issubset(p):
        raise ValueError("Missing directly measured morphology")
    p = p.dropna(subset=list(TRAITS)).reset_index()
    p["genus"] = p.accepted_species.astype(str).str.split().str[0]
    return p


def family_crosswalk(taxa: pd.DataFrame) -> tuple[pd.Series, int]:
    if not {"genus", "family"}.issubset(taxa):
        raise ValueError("Original island taxonomy lacks genus and family")
    x = taxa[["genus", "family"]].dropna().copy()
    x = x.loc[~x.family.astype(str).str.lower().isin(("unknown", "nan", "unresolved"))]
    g = x.groupby("genus").family.agg(lambda col: sorted(set(col)))
    conflicting = int(g.map(len).gt(1).sum())
    return g.loc[g.map(len).eq(1)].map(lambda v: v[0]), conflicting


def family_transfer(
    before: pd.DataFrame, after: pd.DataFrame,
    direct: pd.DataFrame, frozen_breadth: pd.DataFrame, taxa: pd.DataFrame, *,
    alpha: float = 100., min_references: int = 3,
    min_genera: int = 50, min_families: int = 15,
    bootstrap: int = 4999, seed: int = 20261111
) -> dict:
    """Compare each matched genus before/after provider omission.

    The outcome is recomputed from independent citation-backed provider units.
    Baseline reference-effort covariates are refit on EACH target dataset;
    don't treat changed response plus stale baseline as a matched test.
    """
    if (alpha <= 0 or min_references < 3 or min_genera < 50 or
            min_families < 15 or bootstrap < 199):
        raise ValueError("Unsafe provider/source support rule")
    if not set(CHANNELS).issubset(before) or not set(CHANNELS).issubset(after):
        raise ValueError("Missing five visitor reference channels")
    if not before.index.is_unique or not after.index.is_unique:
        raise ValueError("Duplicate provider-target genus")
    if not {"genus", "n_matched_plant_species"}.issubset(frozen_breadth):
        raise ValueError("Missing frozen source genus metadata")
    if frozen_breadth.genus.isna().any() or frozen_breadth.genus.duplicated().any():
        raise ValueError("Duplicate frozen plant genus")
    morph = direct_morphology(direct)
    counts = morph.groupby("genus").size()
    open_form = morph.floral_form.eq("open_radial").groupby(morph.genus).mean()
    families, conflicts = family_crosswalk(taxa)
    support = (before[list(CHANNELS)].sum(axis=1).ge(min_references) &
               after.reindex(before.index)[list(CHANNELS)].sum(axis=1).ge(min_references))
    eligible = before.index[support]
    eligible = eligible.intersection(counts.index)
    eligible = eligible.intersection(families.index)
    eligible = eligible.intersection(frozen_breadth.genus)
    fam = families.reindex(eligible)
    # Never use the target to decide which family or which genus is eligible.
    receipt = {
        "n_frozen_genera": int(len(frozen_breadth)),
        "n_source_genera_before_matched": int(before[list(CHANNELS)].sum(axis=1).ge(min_references).sum()),
        "n_source_genera_after_matched": int(after[list(CHANNELS)].sum(axis=1).ge(min_references).sum()),
        "n_eligible_direct_morphology_species": int(len(morph)),
        "n_taxonomy_conflicting_genera": int(conflicts),
        "n_same_genera_before_after_with_direct_traits_and_unambiguous_family": int(len(eligible)),
        "n_same_families": int(fam.nunique()),
        "heldout_original_provider_prediction_completed": False,
    }
    if len(eligible) < min_genera or fam.nunique() < min_families:
        return {"status": "NOT_IDENTIFIABLE_AFTER_PROVIDER_OMISSION", **receipt}
    frame = frozen_breadth.set_index("genus").reindex(eligible)
    genus_n = counts.reindex(eligible).to_numpy(float)
    same = eligible
    X_common = np.column_stack((
        np.log1p(frame.n_matched_plant_species.to_numpy(float)),
        np.log1p(genus_n),
    ))
    form = open_form.reindex(eligible).to_numpy(float)
    if not np.isfinite(X_common).all() or not np.isfinite(form).all():
        raise ValueError("Missing source support / direct morphology information")
    family_labels = fam.astype(str).to_numpy()
    groups = [np.flatnonzero(family_labels == name)
              for name in np.unique(family_labels)]
    rng = np.random.default_rng(seed)
    results = {}
    for label, ref in (("pre_provider_omission", before),
                       ("post_provider_omission", after)):
        count = ref.reindex(same)[list(CHANNELS)].to_numpy(float)
        total = count.sum(axis=1)
        if (total < min_references).any():
            raise ValueError("Matched source reference threshold broken")
        Y = count / total[:, None]
        base = np.column_stack((np.log1p(total), X_common))
        pred = {}
        for model, X in (
            ("source_support_only", base),
            ("source_support_plus_open_form", np.column_stack((base, form))),
        ):
            result = np.empty_like(Y)
            for test in groups:
                train = np.setdiff1d(np.arange(len(Y)), test, assume_unique=True)
                x_train = X[train]
                means = x_train.mean(axis=0)
                sd = x_train.std(axis=0)
                sd[sd < 1e-12] = 1
                a = (x_train - means) / sd
                b = (X[test] - means) / sd
                offset = Y[train].mean(axis=0)
                beta = np.linalg.solve(
                    a.T @ a + alpha * np.eye(a.shape[1]),
                    a.T @ (Y[train] - offset)
                )
                result[test] = offset + b @ beta
            pred[model] = _project(result)
        err0 = (pred["source_support_only"] - Y)**2
        err1 = (pred["source_support_plus_open_form"] - Y)**2
        gain = (err0 - err1).mean(axis=1)
        familyscore = np.array([gain[test].mean() for test in groups])
        draws = familyscore[rng.integers(
            len(groups), size=(bootstrap, len(groups))
        )].mean(axis=1)
        results[label] = {
            "genus_equal_baseline_MSE": float(err0.mean()),
            "genus_equal_open_form_MSE": float(err1.mean()),
            "genus_equal_open_form_gain": float(gain.mean()),
            "family_equal_open_form_gain": float(familyscore.mean()),
            "family_equal_ci025": float(np.quantile(draws, .025)),
            "family_equal_ci975": float(np.quantile(draws, .975)),
            "n_positive_families": int((familyscore > 0).sum()),
        }
    return {
        "status": "COMPLETED_EXPLORATORY_PROVIDER_OMITTED_FAMILY_TRANSFER",
        **receipt,
        "heldout_original_provider_prediction_completed": True,
        "dominant_provider_grouping_is_archive_URI_and_must_be_audited": True,
        "pre": results["pre_provider_omission"],
        "after": results["post_provider_omission"],
        "local_pollinator_presence_or_effectiveness_measured": False,
        "original_taxon_or_provider_IDs_released": False,
    }
