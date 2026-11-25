"""Audit whether EuPPollNet FLOWER COUNTS and GIFT visitor plants are identical.

Published Flower_data_merger, study/network/date and positive numeric source
flower count are NECESSARY, NOT sufficient: the original source flower survey
Plant_species must be the SAME unambiguous binomial as the accepted GIFT
flower-visitor plant on the source island. A wrong plant denominator must
never enter visitor preference analyses. This is an outcome-blind source
validation, not a model of visitation choice, flower availability or fitness.

Return only anonymous counts. Source identities, names, raw locations and
linked rows stay inside original SHA-locked CI.
"""
from __future__ import annotations

from collections import Counter
import pandas as pd

from .euppollnet_gift_exact_site_bridge import (
    strict_original_accepted_plant, original_guild,
)
from .fullframe_globi_local_support import genuine_direct_motifs

VISITS={
    "Study_id","Network_id","Date","Flower_data_merger",
    "island_id","plant","origin_status","endemic_status",
    "Pollinator_order","Pollinator_family","Pollinator_accepted_name",
}
FLOWERS={
    "Study_id","Network_id","Day","Month","Year",
    "Flower_data_merger","Plant_species","Flower_count","Units",
}
MOTIFS=("actinomorphic_deep","zygomorphic_deep",
        "tube_bell_deep","white_tube_bell","red_pink_tube_bell")
GUILDS=("bombus","other_bees","butterfly_family",
        "sphingidae_moth","flies_order","other_lepidoptera")


def _clean_nonempty_key(x:pd.Series)->pd.Series:
    s=x.astype("string").str.strip()
    bad=s.str.lower().isin(("","na","nan","none","unknown","null"))
    return s.mask(bad.fillna(True))


def qualified_rows_and_identity_receipt(
    source_local_visits:pd.DataFrame,
    flower_counts:pd.DataFrame,
    *, original_parent_n:int=4505,
) -> tuple[pd.DataFrame,dict]:
    """Private table of exact positive records plus public source waterfall."""
    if original_parent_n!=4505:
        raise ValueError("Original GIFT parent floras changed")
    if not VISITS<=set(source_local_visits) or not FLOWERS<=set(flower_counts):
        raise ValueError("Original source full accepted visit/flower fields missing")
    if source_local_visits[["Study_id","Network_id","island_id","plant"]].isna().any().any():
        raise ValueError("Original accepted botanical visitor/site ID missing")
    visits=source_local_visits.copy()
    fl=flower_counts.copy()
    visits["Flower_data_merger"]=_clean_nonempty_key(visits.Flower_data_merger)
    fl["Flower_data_merger"]=_clean_nonempty_key(fl.Flower_data_merger)
    visits=visits.loc[visits.Flower_data_merger.notna()].copy()
    fl=fl.loc[fl.Flower_data_merger.notna()].copy()
    f_freq=fl.Flower_data_merger.value_counts()
    qualified_keys=set(f_freq.index[f_freq.eq(1)])
    unique=fl.loc[fl.Flower_data_merger.isin(qualified_keys)]
    linked=visits.merge(
        unique,on="Flower_data_merger",how="inner",
        suffixes=("","_flower"),validate="many_to_one")
    after_key=len(linked)

    same_study=(linked.Study_id.eq(linked.Study_id_flower)&
        linked.Network_id.eq(linked.Network_id_flower))
    linked=linked.loc[same_study].copy()
    after_study_network=len(linked)

    # The source day may appear as zero-padded text in visits but numeric
    # tokens in the flower table. Compare actual date components.
    date=pd.to_datetime(linked.Date,errors="coerce")
    yr=pd.to_numeric(linked.Year,errors="coerce")
    mo=pd.to_numeric(linked.Month,errors="coerce")
    dy=pd.to_numeric(linked.Day,errors="coerce")
    same_day=(date.notna()&yr.eq(date.dt.year)&
        mo.eq(date.dt.month)&dy.eq(date.dt.day))
    linked=linked.loc[same_day].copy()
    after_day=len(linked)

    count=pd.to_numeric(linked.Flower_count,errors="coerce")
    unit=_clean_nonempty_key(linked.Units)
    ok=count.notna()&count.gt(0)&unit.notna()
    linked=linked.loc[ok].copy()
    after_positive=len(linked)

    # The core new discrimination: raw flower-survey plant is not silently
    # equated with the visitor plant even when all other keys match.
    original_names=set(linked.Plant_species.dropna().astype(str).unique())
    normalized={x:strict_original_accepted_plant(x) for x in original_names}
    linked["_flower_binomial"]=linked.Plant_species.map(normalized)
    unresolved=int(linked._flower_binomial.isna().sum())
    mismatch=int((linked._flower_binomial.notna()&
                  linked._flower_binomial.ne(linked.plant)).sum())
    matched=linked.loc[
        linked._flower_binomial.notna()&
        linked._flower_binomial.eq(linked.plant)].copy()
    # Reject any original raw plant/visitor that is undecodable, rather than
    # treating replacement-character taxonomy as a valid exact match.
    if len(matched):
        if matched.plant.astype(str).str.contains("\ufffd",regex=False).any():
            raise ValueError("Original visitor plant has undecodable identity")

    receipt=dict(
        original_GIFT_full_botanical_population=4505,
        n_original_same_island_positive_visitor_rows=int(len(source_local_visits)),
        n_positive_rows_with_nonempty_flower_merger=int(len(visits)),
        n_original_flower_count_keys_repeated=int(f_freq.gt(1).sum()),
        n_positive_rows_matching_unique_original_flower_key=int(after_key),
        n_positive_rows_same_published_study_and_network=int(after_study_network),
        n_positive_rows_same_published_study_network_calendar_day=int(after_day),
        n_positive_rows_with_source_nonzero_flower_count_and_unit=int(after_positive),
        n_positive_rows_with_unresolvable_original_flower_plant=int(unresolved),
        n_positive_rows_with_different_source_flower_plant=int(mismatch),
        n_positive_rows_with_exact_same_flower_and_visitor_plant=int(len(matched)),
        n_original_source_studies_with_exact_plant_source_flower_supply=int(
            matched.Study_id.nunique()),
        n_original_GIFT_islands_with_exact_plant_source_flower_supply=int(
            matched.island_id.nunique()),
        n_original_botanical_plants_with_exact_source_flower_supply=int(
            matched.plant.nunique()),
        n_distinct_exact_flower_survey_keys=int(
            matched.Flower_data_merger.nunique()),
        n_original_source_flower_unit_strings_in_exact_plant_positive_rows=int(
            matched.Units.nunique()),
        original_source_species_binomial_not_synonym_fuzzy_reconciled=True,
        original_source_species_mismatch_is_unknown_not_ecological_no_pollination=True,
        source_nonvisitor_flower_availability_NOT_verified=True,
    )
    if (len(matched)+unresolved+mismatch!=after_positive):
        raise AssertionError("Plant matching waterfall does not reconcile")
    return matched,receipt


def _grouped_source_contexts(matched:pd.DataFrame, direct_traits:dict)->dict:
    """Anonymous same-unit source comparability for actual two-axis flowers."""
    q=matched.copy()
    q["_guild"]=[original_guild(o,f,a) for o,f,a in zip(
        q.Pollinator_order,q.Pollinator_family,
        q.Pollinator_accepted_name)]
    motifs=genuine_direct_motifs(direct_traits)
    if set(motifs)!=set(MOTIFS):
        raise ValueError("Source-original same-plant motif ontology drifted")
    cohort=["Study_id","Network_id","Date","island_id","Units"]
    contexts=q.groupby(cohort,dropna=False,sort=False)
    total_contexts=0
    with_two_source_plants=0
    with_two_guilds=0
    with_two_genera=0
    with_two_source_flower_count_units=0
    paired=[]
    for _,group in contexts:
        total_contexts+=1
        plants=set(group.plant)
        guilds=set(group._guild)-{"unresolved"}
        genera={name.split()[0] for name in plants if len(name.split())>=2}
        with_two_source_plants+=int(len(plants)>=2)
        with_two_guilds+=int(len(guilds)>=2)
        with_two_genera+=int(len(genera)>=2)
        with_two_source_flower_count_units+=int(
            group.Flower_data_merger.nunique()>=2)
        for motif in MOTIFS:
            truth=motifs[motif]
            scored=plants&set(truth)
            yes={s for s in scored if truth[s]}
            no=scored-yes
            paired.append((motif,len(scored)>=2,len(yes)>=1,
                len(no)>=1,len(genera)>=2,len(guilds)>=2,
                group.Study_id.iloc[0],group.island_id.iloc[0]))
    counts=[]
    for motif in MOTIFS:
        relevant=[x for x in paired if x[0]==motif]
        comparable=[x for x in relevant if
            x[1] and x[2] and x[3]]
        both_guild=[x for x in comparable if x[5]]
        both_genus=[x for x in both_guild if x[4]]
        counts.append(dict(
            real_same_species_direct_motif=motif,
            n_source_contexts_with_true_and_false_scored_flower_species=len(comparable),
            n_source_contexts_with_flower_states_and_two_documented_guilds=len(both_guild),
            n_source_contexts_with_flower_states_two_guilds_and_two_genera=len(both_genus),
            n_source_studies_with_comparable_flower_states=len({x[6] for x in comparable}),
            n_original_GIFT_islands_with_comparable_flower_states=len({x[7] for x in comparable}),
            source_positive_documentation_not_verified_nondetection=True,
        ))
    return dict(
        n_same_source_study_network_date_island_unit_contexts=total_contexts,
        n_contexts_with_two_botanical_flower_species=with_two_source_plants,
        n_contexts_with_two_documented_visitor_guilds=with_two_guilds,
        n_contexts_with_two_botanical_genera=with_two_genera,
        n_contexts_with_two_distinct_flower_count_survey_keys=with_two_source_flower_count_units,
        source_motif_comparison_opportunities=counts,
        source_repeat_visit_rows_not_independent_same_day_replicates=True,
    )


def original_same_plant_exposure_audit(
    original_visits:pd.DataFrame, flowers:pd.DataFrame,traits:dict,
    *,original_parent_n:int=4505
)->dict:
    matched,waterfall=qualified_rows_and_identity_receipt(
        original_visits,flowers,original_parent_n=original_parent_n)
    contexts=_grouped_source_contexts(matched,traits)
    return dict(
        status="ORIGINAL_GIFT_4505_EUPPOLLNET_SAME_PLANT_EXPOSURE_AND_MOTIF_SUPPORT",
        original_GIFT_parent_botanical_floras=4505,
        exact_source_plant_identity_waterfall=waterfall,
        source_same_study_site_date_and_unit_comparability=contexts,
        source_flower_units_not_commensurate_without_physical_unit_conversion=True,
        original_source_published_nonvisited_botanical_flower_exposure_unverified=True,
        no_visitor_preference_effect_or_pollen_fitness_coefficient=True,
        all_source_original_identifiers_retained_protected_only=True,
    )
