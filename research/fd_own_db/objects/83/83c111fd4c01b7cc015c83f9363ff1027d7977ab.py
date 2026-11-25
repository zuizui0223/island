"""Pinned-source GloBI original-provider availability and omission gate.

The currently archived 2,228-row GloBI *genus summary* cannot reconstruct
original study/dataset citations. A provider-independent holdout MUST receive
real source-level records with versioned citation and provider identities.
Never simulate sources from aggregate genus reference counts. Nothing returned
from the public audit contains raw source, genus or reference identities.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

ORIGINAL_GLOBI_SHA256 = "095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17"
CHANNELS = ("bombus", "non_bombus_bees", "lepidoptera",
            "diptera", "flower_visiting_birds")
REQUIRED = ("plant_genus", "visitor_channel", "provider_id",
            "reference_id", "argument_type", "evidence_class")


def _validate(original: pd.DataFrame, *, source_sha256: str) -> pd.DataFrame:
    if source_sha256 != ORIGINAL_GLOBI_SHA256:
        raise ValueError("Unverified GloBI source snapshot SHA256")
    missing = set(REQUIRED)-set(original)
    if missing:
        raise ValueError("Original provider/ref keys absent from summarized GloBI: "+str(sorted(missing)))
    x = original[list(REQUIRED)].copy()
    if x.empty:
        raise ValueError("Empty original-source GloBI record frame")
    for col in REQUIRED:
        x[col] = x[col].astype("string").str.strip()
        if x[col].isna().any() or x[col].eq("").any() or x[col].str.lower().isin(("unknown", "nan", "unresolved")).any():
            raise ValueError("Unresolved original provider/interaction evidence")
    if not x.visitor_channel.isin(CHANNELS).all():
        raise ValueError("Unverified visitor guild")
    if not x.argument_type.eq("SUPPORTS").all():
        raise ValueError("REFUTES/unknown interactions are not supporting visits")
    if not x.evidence_class.isin(("pollination", "flower_visit")).all():
        raise ValueError("Unverified pollination/flower-visit evidence")
    return x.drop_duplicates()


def source_coverage_gate(original: pd.DataFrame, *, source_sha256: str) -> dict:
    """Count independently attributable datasets and crossed reference IDs.

    The same citation appearing in multiple provider datasets is not an
    independent test reference. It must be discarded from both sides.
    """
    x = _validate(original, source_sha256=source_sha256)
    n_providers = x.provider_id.nunique()
    if n_providers < 3:
        raise ValueError("Insufficient independent GloBI datasets")
    collision = x.groupby("reference_id").provider_id.nunique()
    repeated = set(collision[collision > 1].index)
    single = x[~x.reference_id.isin(repeated)]
    units = single[["plant_genus", "visitor_channel", "reference_id", "provider_id"]].drop_duplicates()
    size = units.provider_id.value_counts()
    if units.empty:
        raise ValueError("No independent provider/reference evidence remains")
    return {
        "n_raw_source_rows": int(len(x)),
        "n_observed_plant_genera": int(x.plant_genus.nunique()),
        "n_independent_source_providers": int(n_providers),
        "n_cross_provider_reference_collisions": int(len(repeated)),
        "n_uncollided_reference_ids": int(single.reference_id.nunique()),
        "n_genus_channel_reference_units": int(len(units)),
        "largest_provider_reference_fraction": float(size.iloc[0]/len(units)),
        "provider_holdout_completed": False,
    }


def drop_largest_provider(
    original: pd.DataFrame, *, source_sha256: str,
    min_references_per_genus: int = 3,
) -> tuple[pd.DataFrame, dict]:
    """Return PRIVATE revised genus×channel counts after deleting the main dataset.

    Caller must then fit heldout FAMILY models on the EXACT SAME support as
    its matched source reference, and apply distinct-reference/provider gates.
    No effect estimate is made by this routine.
    """
    if min_references_per_genus < 3:
        raise ValueError("Genus evidence threshold must be >=3 independent references")
    x = _validate(original, source_sha256=source_sha256)
    collisions = x.groupby("reference_id").provider_id.nunique()
    shared = set(collisions[collisions > 1].index)
    z = x[~x.reference_id.isin(shared)][
        ["plant_genus", "visitor_channel", "provider_id", "reference_id"]
    ].drop_duplicates()
    rankings = z.provider_id.value_counts()
    if len(rankings) < 3:
        raise ValueError("Not enough original source providers to hold one out")
    after = z[~z.provider_id.eq(rankings.index[0])]
    counts = after.groupby("plant_genus").reference_id.nunique()
    keep = set(counts[counts >= min_references_per_genus].index)
    after = after[after.plant_genus.isin(keep)]
    by = after.groupby(["plant_genus", "visitor_channel"]).reference_id.nunique().unstack(fill_value=0)
    by = by.reindex(columns=CHANNELS, fill_value=0).astype(int)
    return by, {
        "n_original_providers": int(len(rankings)),
        "n_provider_removed": 1,
        "n_cross_provider_references_discarded": int(len(shared)),
        "n_supported_genera_after_omission": int(len(by)),
        "n_channel_reference_units_after_omission": int(by.to_numpy().sum()),
        "heldout_family_prediction_recomputed": False,
    }
