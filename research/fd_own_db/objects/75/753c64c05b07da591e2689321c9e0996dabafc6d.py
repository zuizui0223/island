"""EuPPollNet original STUDY/YEAR field protocol documentary eligibility.

Uses original protected accepted GIFT same-plant visitor×flower survey source
rows (not merely original GloBI provider labels), and public 2025 EuPPollNet
study/year metadata. A study reporting percentages of animal orders is NOT
proof that all guilds were uniformly surveyed or that unreported edges are
observed animal nonvisits. No taxon/island/study/year/site identifiers exported.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from .euppollnet_gift_exact_site_bridge import original_guild

MUST={
 "Study_id","Year","Sampling_method","Sampling_days",
 "Sampling_period","Taxa_recorded","Hymenoptera_percent",
 "Lepidoptera_percent",
}
SOURCE={
 "Study_id","Network_id","Date","island_id","plant",
 "Pollinator_order","Pollinator_family","Pollinator_accepted_name",
}


def study_year_original_protocol_gate(
    same_source_species_flower_visits:pd.DataFrame,
    original_public_study_metadata:pd.DataFrame,
    *,n_original_global_GIFT_floras:int=4505,
)->dict:
    """Study/year methodological opportunity (not animal choice effect)."""
    if n_original_global_GIFT_floras!=4505:
        raise ValueError("GIFT 4505 original flora parent changed")
    if not SOURCE<=set(same_source_species_flower_visits):
        raise ValueError("Original same-plant verified visitor source fields missing")
    if not MUST<=set(original_public_study_metadata):
        raise ValueError("Original European method/protocol metadata columns missing")
    q=same_source_species_flower_visits.copy()
    m=original_public_study_metadata.copy()
    source_year=pd.to_datetime(q.Date,errors="coerce").dt.year
    if source_year.isna().any() or q.Study_id.isna().any():
        raise ValueError("Original source study/dates ambiguous")
    q["_source_year"]=source_year.astype(int)
    original_studies=q[["Study_id","_source_year"]].drop_duplicates()
    m["_source_year"]=pd.to_numeric(m.Year,errors="coerce")
    m=m.loc[m.Study_id.notna() & m._source_year.notna()].copy()
    m["_source_year"]=m._source_year.astype(int)
    fields=("Sampling_method","Sampling_days","Sampling_period",
            "Taxa_recorded","Hymenoptera_percent","Lepidoptera_percent")
    # Source original study×year metadata may be repeated. Count then
    # quarantine non-unique rows rather than silently choose first.
    duplicates=m.groupby(["Study_id","_source_year"]).size()
    ambiguity=set(duplicates.index[duplicates.gt(1)])
    m=m.loc[~pd.MultiIndex.from_frame(
        m[["Study_id","_source_year"]]).isin(ambiguity)]
    match=original_studies.merge(
        m[["Study_id","_source_year",*fields]].copy(),
        on=["Study_id","_source_year"],how="left",validate="one_to_one",
        indicator=True)
    qualified=match.loc[match["_merge"].eq("both")].copy()
    sampling_days=pd.to_numeric(
        qualified.Sampling_days,errors="coerce")
    sampling_period=pd.to_numeric(
        qualified.Sampling_period,errors="coerce")
    hym=pd.to_numeric(
        qualified.Hymenoptera_percent,errors="coerce")
    lep=pd.to_numeric(
        qualified.Lepidoptera_percent,errors="coerce")
    methods=qualified.Sampling_method.fillna(
        "UNREPORTED").astype(str).str.strip()
    known_methods=methods.ne("")&~methods.str.lower().isin(
        ("na","none","unreported","unknown","nan"))
    both_reported=hym.between(0,1)&lep.between(0,1)
    q["_source_guild"]=[original_guild(o,f,a) for o,f,a in zip(
        q.Pollinator_order,q.Pollinator_family,
        q.Pollinator_accepted_name)]
    group=q.groupby(["Study_id","_source_year"],sort=False)
    two_guild_ids={k for k,g in group if len(
        set(g._source_guild)-{"unresolved"})>=2}
    observed_source_pairs=set(map(tuple,
        original_studies[["Study_id","_source_year"]].itertuples(index=False,name=None)))
    if len(observed_source_pairs)!=len(original_studies):
        raise AssertionError("Study/year original positive source identity duplicated")
    return dict(
        status="ORIGINAL_GIFT_EUPPOLLNET_STUDY_YEAR_METHOD_EVIDENCE_AUDIT",
        n_original_parent_island_floras=4505,
        n_original_positive_same_plant_source_visit_rows=int(len(q)),
        n_original_distinct_source_study_ids=int(q.Study_id.nunique()),
        n_original_distinct_source_study_years=int(len(original_studies)),
        n_study_years_with_exact_unique_public_metadata=int(len(qualified)),
        n_study_years_without_unique_original_metadata=int(
            len(original_studies)-len(qualified)),
        n_original_source_metadata_study_year_keys_ambiguous=int(len(ambiguity)),
        n_matched_study_years_with_reported_sampling_method=int(
            known_methods.sum()),
        n_matched_study_years_with_positive_sampling_days=int(
            sampling_days.gt(0).sum()),
        n_matched_study_years_with_positive_sampling_period=int(
            sampling_period.gt(0).sum()),
        n_matched_study_years_with_both_Hymenoptera_and_Lepidoptera_reported_percentages=int(
            both_reported.sum()),
        n_source_study_years_with_positive_links_from_two_animal_guild_classes=int(
            sum(k in two_guild_ids for k in observed_source_pairs)),
        original_method_labels_frequency={
            k:int(v) for k,v in methods[known_methods].value_counts().to_dict().items()},
        original_animal_order_percentage_columns_are_observed_not_sampling_completeness=True,
        original_study_method_is_not_complete_guild_nondetection_protocol=True,
        source_guild_cooccurrence_is_not_standardized_preference=True,
        original_species_direct_plant_traits_not_local_measured=True,
        no_inferred_missing_visits_or_pollen_delivery=True,
        original_private_study_year_plant_island_ids_not_serialized=True,
    )
