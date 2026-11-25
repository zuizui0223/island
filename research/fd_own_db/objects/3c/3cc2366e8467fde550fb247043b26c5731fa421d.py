"""Outcome-blind priority for independent SC/SI observations, not imputation.

Locally authorized source-disjoint input only. Returns aggregate evidence
coverage upper bounds, never accepted species/genus/source identities.
Candidate reproductive states must be truly UNKNOWN.
"""
from __future__ import annotations

import pandas as pd

STEPS = (1, 5, 10, 25, 50)
FIELDS = ("accepted_species", "genus")


def _check(frame: pd.DataFrame, label: str) -> pd.DataFrame:
    if not set(FIELDS).issubset(frame.columns):
        raise ValueError(f"{label}: missing fields")
    x = frame[list(FIELDS)].copy()
    for col in FIELDS:
        if x[col].isna().any() or x[col].astype(str).str.strip().eq("").any():
            raise ValueError(f"{label}: blank names")
        x[col] = x[col].astype(str).str.strip()
    if x.accepted_species.duplicated().any():
        raise ValueError(f"{label}: duplicate species")
    if not (x.accepted_species.str.split().str[0] == x.genus).all():
        raise ValueError(f"{label}: genus mismatch")
    return x


def prospective_sc_donor_coverage(
    test_sc: pd.DataFrame,
    donor_sc_direct: pd.DataFrame,
    prospective_sc: pd.DataFrame,
    *,
    all_island_species: set[str],
    minimum_targets: int = 50,
) -> pd.DataFrame:
    """Upper bound after one new measured Direct SC/SI species per selected genus.

    Call only after removing all original island species and evaluated source
    regions from eligible donors. Prospective must NOT contain known Direct
    SC/SI states. Ranking uses target genus counts, never observed SC states.
    """
    if minimum_targets < 50:
        raise ValueError("Privacy minimum is 50 target species")
    test = _check(test_sc, "test_sc")
    donors = _check(donor_sc_direct, "donor_sc_direct")
    cand = _check(prospective_sc, "prospective_sc")
    if len(test) < minimum_targets:
        raise ValueError("Too few tested SC taxa")
    excluded = set(all_island_species)
    if not set(test.accepted_species).issubset(excluded):
        raise ValueError("Target lacks strict all-island-name holdout")
    if not set(donors.accepted_species).isdisjoint(excluded):
        raise ValueError("Existing SC donor leakage")
    if not set(cand.accepted_species).isdisjoint(excluded):
        raise ValueError("Prospective candidate leakage")
    if not set(cand.accepted_species).isdisjoint(set(donors.accepted_species)):
        raise ValueError("Candidate already has Direct SC status")
    counts = test.genus.value_counts()
    existing = set(donors.genus)
    possible = set(cand.genus)
    newly_covered = counts[
        ~counts.index.isin(existing) & counts.index.isin(possible)
    ].sort_values(ascending=False, kind="stable")
    n_test = len(test)
    already_covered = int(counts[counts.index.isin(existing)].sum())
    result = {
        "n_test_sc_species": n_test,
        "n_current_direct_sc_donors": len(donors),
        "n_prospective_species_missing_sc": len(cand),
        "n_prospective_genera_with_uncovered_test_sc": len(newly_covered),
        "current_same_genus_coverage": already_covered/n_test,
        "max_candidate_supported_coverage": (
            already_covered+int(newly_covered.sum()))/n_test,
    }
    for count in STEPS:
        result[f"hypothetical_coverage_after_{count}_new_genera"] = (
            already_covered+int(newly_covered.iloc[:count].sum()))/n_test
    return pd.DataFrame([result])
