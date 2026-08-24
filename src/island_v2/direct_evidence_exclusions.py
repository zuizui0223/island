"""Apply reviewed exclusions to species-direct evidence ledgers."""

from __future__ import annotations

import pandas as pd

DIRECT_EVIDENCE_EXCLUSION_KEY = [
    "accepted_species",
    "trait_name",
    "normalized_value",
    "source_lineage",
]


def apply_direct_evidence_exclusions(
    evidence: pd.DataFrame,
    exclusions: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Remove reviewed exact records or reviewed source-concept families.

    Exact exclusions retain the species x trait x value x lineage contract.
    A source-concept error may instead use ``*`` for species and value together
    with a lineage-prefix key ending in ``*``. This narrow wildcard form lets
    one reviewed ontology decision remove every affected redistribution row.
    """

    required = set(DIRECT_EVIDENCE_EXCLUSION_KEY).union({"reason", "reviewer", "reviewed_at_utc"})
    missing = required.difference(exclusions.columns)
    if missing:
        raise ValueError(f"direct-evidence exclusions missing columns: {sorted(missing)}")
    work = exclusions.copy().fillna("")
    if work[DIRECT_EVIDENCE_EXCLUSION_KEY].duplicated().any():
        raise ValueError("direct-evidence exclusion keys must be unique")
    if work[list(required)].apply(lambda column: column.astype(str).str.strip().eq("").any()).any():
        raise ValueError("direct-evidence exclusions require complete keys and review provenance")
    wildcard = (
        work[DIRECT_EVIDENCE_EXCLUSION_KEY]
        .apply(lambda column: column.astype(str).str.contains("*", regex=False))
        .any(axis=1)
    )
    valid_wildcard = (
        work["accepted_species"].eq("*")
        & work["normalized_value"].eq("*")
        & work["source_lineage"].str.endswith("*")
        & ~work["trait_name"].str.contains("*", regex=False)
        & work["source_lineage"].str.count(r"\*").eq(1)
    )
    if (wildcard & ~valid_wildcard).any():
        raise ValueError(
            "source-concept exclusions require species='*', value='*', "
            "an exact trait, and one trailing lineage-prefix '*'"
        )
    if evidence.empty or work.empty:
        audit = work.copy()
        audit["matched_rows"] = 0
        return evidence.copy(), audit

    frame = evidence.copy().fillna("")
    rejected = pd.Series(False, index=frame.index)
    matched_rows: list[int] = []
    for item in work.to_dict("records"):
        match = frame["trait_name"].eq(item["trait_name"])
        if item["accepted_species"] != "*":
            match &= frame["accepted_species"].eq(item["accepted_species"])
        if item["normalized_value"] != "*":
            match &= frame["normalized_value"].eq(item["normalized_value"])
        if item["source_lineage"].endswith("*"):
            match &= frame["source_lineage"].str.startswith(item["source_lineage"][:-1])
        else:
            match &= frame["source_lineage"].eq(item["source_lineage"])
        matched_rows.append(int(match.sum()))
        rejected |= match
    audit = work.copy()
    audit["matched_rows"] = matched_rows
    return frame.loc[~rejected].reset_index(drop=True), audit
