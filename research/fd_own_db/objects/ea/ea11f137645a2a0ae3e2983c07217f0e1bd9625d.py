"""Two-trait donor-opportunity matched predictive audit (licensed inputs stay local).

Caller first excludes ALL island-flora species and heldout mainland source
regions from the donor list. Both states must be actually recorded for
EVERY donor and evaluated taxon. Only coarse aggregate rows are returned.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

TRAITS = ("self_incompatibility", "floral_symmetry")
FIELDS = ("accepted_species", "genus", *TRAITS)


def _checked(frame: pd.DataFrame, title: str) -> pd.DataFrame:
    missing = set(FIELDS) - set(frame.columns)
    if missing:
        raise ValueError(f"{title}: missing {sorted(missing)}")
    x = frame[list(FIELDS)].copy()
    for col in ("accepted_species", "genus"):
        if x[col].isna().any() or x[col].astype(str).str.strip().eq("").any():
            raise ValueError(f"{title}: blank {col}")
        x[col] = x[col].astype(str).str.strip()
    if x.accepted_species.duplicated().any():
        raise ValueError(f"{title}: duplicate accepted species")
    if not (x.accepted_species.str.split().str[0] == x.genus).all():
        raise ValueError(f"{title}: species name and genus disagree")
    for col in TRAITS:
        raw = pd.to_numeric(x[col], errors="coerce")
        if raw.isna().any() or ~raw.isin((0, 1)).all():
            raise ValueError(f"{title}: {col} must be measured binary 0/1")
        x[col] = raw.astype(int)
    return x


def matched_trait_gains(
    donors: pd.DataFrame,
    targets: pd.DataFrame,
    *,
    globally_excluded_species: set[str],
    tau: float = 1.,
    min_test_species: int = 50,
    min_train_species: int = 20,
) -> pd.DataFrame:
    """Use IDENTICAL training and test names for both strictly measured traits.

    The caller must independently verify excluded mainland source entities
    are absent: source IDs do not enter this scoring function. Require full
    all-island-name exclusion, not only target-component species exclusion.
    """
    if min_test_species < 50 or min_train_species < 20:
        raise ValueError("Privacy/support threshold not met")
    if not np.isfinite(tau) or tau <= 0:
        raise ValueError("tau must be finite and positive")
    d = _checked(donors, "donors")
    t = _checked(targets, "targets")
    if len(d) < min_train_species or len(t) < min_test_species:
        raise ValueError("Insufficient observed matched trait carriers")
    excluded = set(globally_excluded_species)
    if not set(t.accepted_species).issubset(excluded):
        raise ValueError("Evaluated species absent from global island exclusion")
    if not set(d.accepted_species).isdisjoint(excluded):
        raise ValueError("Global species-name donor leakage")
    train_n = d.groupby("genus").size()
    target_n = t.genus.map(train_n).fillna(0).to_numpy(dtype=float)
    covered = target_n > 0
    # The most frequent covered genera are selected without using traits.
    covered_size = t.loc[covered, "genus"].value_counts()
    largest = set(covered_size.index[:1])
    top_three = set(covered_size.index[:3])
    records = []
    for trait in TRAITS:
        baseline = float(d[trait].mean())
        train_k = d.groupby("genus")[trait].sum()
        target_k = t.genus.map(train_k).fillna(0).to_numpy(float)
        y = t[trait].to_numpy(float)
        p = (target_k + tau * baseline) / (target_n + tau)
        global_loss = (y - baseline)**2
        genus_loss = (y - p)**2
        gain = global_loss - genus_loss
        by_genus = pd.Series(gain).groupby(t.genus.to_numpy()).mean()
        no_one = ~t.genus.isin(largest).to_numpy()
        no_three = ~t.genus.isin(top_three).to_numpy()
        records.append({
            "trait": trait,
            "n_test_species": int(len(t)),
            "n_train_species": int(len(d)),
            "n_test_genera": int(t.genus.nunique()),
            "n_covered_genera": int(covered_size.size),
            "donor_coverage": float(covered.mean()),
            "top_three_share_of_covered": (
                float(covered_size.iloc[:3].sum()/covered.sum())
                if covered.any() else np.nan),
            "brier_global": float(global_loss.mean()),
            "brier_genus": float(genus_loss.mean()),
            "name_equal_gain": float(gain.mean()),
            "genus_equal_gain": float(by_genus.mean()),
            "covered_only_gain": (
                float(gain[covered].mean()) if covered.any() else np.nan),
            "remove_most_common_covered_genus_gain": (
                float(gain[no_one].mean()) if no_one.any() else np.nan),
            "remove_three_most_common_covered_genera_gain": (
                float(gain[no_three].mean()) if no_three.any() else np.nan),
        })
    return pd.DataFrame(records)
