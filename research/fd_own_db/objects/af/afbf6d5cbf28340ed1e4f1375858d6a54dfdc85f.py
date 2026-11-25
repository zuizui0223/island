"""Outcome-blind source-donor independence eligibility on the FULL fd island cohort.

Requires locally authorized island regime and five-nearest GIFT source IDs.
Emits coarse region aggregates only; never export original source labels.
This is a conservative donor-overlap audit, NOT a historical ancestry model.
"""
from __future__ import annotations
import pandas as pd
import numpy as np
from .source_signature_dependency import _component_labels

OUT = (
 "region","n_total_islands","n_labeled_islands","n_source_signatures",
 "components_any_shared_source","largest_component_islands_any_shared",
 "components_at_least_ten_islands","components_three_shared_sources",
 "largest_component_islands_three_shared",
)

def audit_full_cohort(regimes: pd.DataFrame, sources: pd.DataFrame,
                      *, min_islands: int=40, min_signatures: int=5) -> pd.DataFrame:
    """Audit the original eligible island universe *before trait filtering*.

    Do NOT infer absence of ancestral source dispersal from no shared GIFT
    entity labels. Graph connectivity is a conservative source-reuse stress.
    """
    if min_islands<40 or min_signatures<5:
        raise ValueError("Privacy/support threshold too low")
    if not {"island_id","analysis_regime"}.issubset(regimes):
        raise ValueError("Missing island regimes")
    if not {"island_id","source_mode","source_rank","entity_ID"}.issubset(sources):
        raise ValueError("Missing five-source assignments")
    x=regimes[["island_id","analysis_regime"]].copy()
    for col in x:
        if x[col].isna().any() or x[col].astype(str).str.strip().eq("").any():
            raise ValueError("Missing island regime")
        x[col]=x[col].astype(str).str.strip()
    if x.island_id.duplicated().any():
        raise ValueError("An island must have exactly one regime")
    y=sources[["island_id","source_mode","source_rank","entity_ID"]].copy()
    for col in ("island_id","source_mode","entity_ID"):
        if y[col].isna().any() or y[col].astype(str).str.strip().eq("").any():
            raise ValueError("Missing source identity")
        y[col]=y[col].astype(str).str.strip()
    y.source_rank=pd.to_numeric(y.source_rank,errors="coerce")
    if y.source_rank.isna().any():
        raise ValueError("Invalid source rank")
    y=y[y.source_mode.eq("geo_k5") & y.source_rank.between(1,5)].copy()
    if (y.duplicated(["island_id","source_rank"]).any()
            or y.duplicated(["island_id","entity_ID"]).any()):
        raise ValueError("Duplicate donor or source rank")
    s=y.groupby("island_id").entity_ID.agg(lambda z:tuple(sorted(set(z))))
    if not s.map(len).eq(5).all():
        raise ValueError("Each source signature needs five distinct regions")
    x=x.join(s.rename("signature"),on="island_id")
    if x.signature.isna().any():
        raise ValueError("Unlabeled island in full cohort")
    rows=[]
    for region,g in x.groupby("analysis_regime",sort=True):
        if len(g)<min_islands:continue
        sigs=list(g.signature.unique())
        if len(sigs)<min_signatures:continue
        n1,lab1=_component_labels(sigs,1)
        n3,lab3=_component_labels(sigs,3)
        cmap1=dict(zip(sigs,lab1))
        cmap3=dict(zip(sigs,lab3))
        size1=g.signature.map(cmap1).value_counts()
        size3=g.signature.map(cmap3).value_counts()
        rows.append(dict(
            region=region,n_total_islands=int(len(g)),
            n_labeled_islands=int(g.signature.notna().sum()),
            n_source_signatures=int(len(sigs)),
            components_any_shared_source=int(n1),
            largest_component_islands_any_shared=int(size1.max()),
            components_at_least_ten_islands=int((size1>=10).sum()),
            components_three_shared_sources=int(n3),
            largest_component_islands_three_shared=int(size3.max()),
        ))
    return pd.DataFrame(rows,columns=OUT)
