"""Local-source floral assembly vs TRUE pollinator-shift feasibility audit.

Only original high/medium species-direct single-state traits are accepted.
GIFT five-nearest native mainland lists are POTENTIAL sources, not ancestry.
GBIF animal occurrences are not verified visits to a specific flower.
GloBI global genus reference proportions are not local pollinator abundance.
The same accepted species has ONE original floral state across all locations;
therefore this data structure cannot test within-species island trait evolution.
Only broad anonymous aggregate classes are returned.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

TRAITS=("flower_primary_color","floral_form","floral_symmetry")
GUILDS=("bombus","non_bombus_bees","lepidoptera","flower_visiting_birds","diptera")


def _require(frame, cols, name):
    missing = set(cols)-set(frame)
    if missing:
        raise ValueError(f"{name}: required source fields missing {sorted(missing)}")
    return frame[list(cols)].copy()


def source_pool_audit(flora:pd.DataFrame, direct_ledger:pd.DataFrame,
        source_assignments:pd.DataFrame, native_mainland:pd.DataFrame,
        pollinator_evidence:pd.DataFrame, island_geo:pd.DataFrame,
        globi_reference_genera:pd.DataFrame, *,
        min_direct_species:int=30, min_group_islands:int=10,
        min_group_spatial_blocks:int=5, min_globi_refs:int=3):
    """Return (region/group safe aggregates, anonymous cohort receipt).

    Strictly no original island_id, species, genus, provider, study or
    spatial_block identifiers may be released as a result column.
    """
    if (min_direct_species<30 or min_group_islands<10 or
        min_group_spatial_blocks<5 or min_globi_refs<3):
        raise ValueError("Privacy / evidence support gate is too low")
    f=_require(flora,("island_id","accepted_species","origin_status"),"flora")
    if f.isna().any().any() or f.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Duplicate or missing island x plant species or status")
    t=_require(direct_ledger,("accepted_species","trait_name","normalized_value",
                              "evidence_scope","evidence_quality"),"trait ledger")
    t=t[t.trait_name.isin(TRAITS)&t.evidence_scope.eq("species_direct")&
        t.evidence_quality.isin(("medium","high"))].copy()
    if t.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Duplicate original species direct trait")
    t=t[t.normalized_value.notna() &
        ~t.normalized_value.astype(str).str.contains("|",regex=False)]
    x=t.pivot(index="accepted_species",columns="trait_name",values="normalized_value")
    if not set(TRAITS).issubset(x):
        raise ValueError("All three original direct floral traits are required")
    x=x.dropna(subset=list(TRAITS))
    x["open_form"]=x.floral_form.eq("open_radial").astype(float)
    x["actinomorphy"]=x.floral_symmetry.eq("actinomorphic").astype(float)
    x["inconspicuous_colour"]=x.flower_primary_color.eq(
        "green_brown_inconspicuous").astype(float)

    island_f=f[f.accepted_species.isin(x.index)].copy()
    known_n=island_f.groupby("island_id").accepted_species.nunique()
    qualified=set(known_n[known_n>=min_direct_species].index)
    if not qualified:
        return pd.DataFrame(), {
            "n_islands_atleast_30_direct_complete_species":0,
            "n_bombus_nondetected_bird_recorded":0,
            "species_level_pollinator_shift_confirmed":False,
        }
    island_f=island_f[island_f.island_id.isin(qualified)].copy()
    island_f["genus"]=island_f.accepted_species.astype(str).str.split().str[0]
    geo=_require(island_geo,("island_id","analysis_regime","spatial_block"),"island geography")
    if geo.island_id.duplicated().any() or geo.isna().any().any():
        raise ValueError("Invalid original unique island geography")
    island_f=island_f.merge(geo,on="island_id",validate="many_to_one")
    if island_f.island_id.nunique()!=len(qualified):
        raise ValueError("Missing original island spatial classification")

    s=_require(source_assignments,("island_id","source_mode","source_rank","entity_ID"),"GIFT geo_k5")
    s=s[s.island_id.isin(qualified)&s.source_mode.eq("geo_k5")&
        s.source_rank.between(1,5)].copy()
    if s.duplicated(["island_id","source_rank"]).any() or s.duplicated(["island_id","entity_ID"]).any():
        raise ValueError("Duplicate ranked GIFT source region")
    if s.groupby("island_id").entity_ID.nunique().reindex(list(qualified)).fillna(0).ne(5).any():
        raise ValueError("Five distinct geographic source regions required for every island")
    gift=_require(native_mainland,("entity_ID","work_species"),"native mainland GIFT")
    if gift.isna().any().any():
        raise ValueError("Missing GIFT source region or plant name")
    gift=gift[gift.entity_ID.isin(s.entity_ID)].drop_duplicates()
    # Critical memory safeguard: do NOT expand source island x ALL mainland
    # flora names (can be >10 million intermediate rows).
    gift_sets=gift.groupby("entity_ID").work_species.agg(set).to_dict()
    source_sets=s.groupby("island_id").entity_ID.agg(list).to_dict()

    ch=_require(pollinator_evidence,("island_id","channel_id","channel_state","source_state"),"GBIF guild states")
    ch=ch[ch.island_id.isin(qualified)].copy()
    if ch.isna().any().any() or ch.duplicated(["island_id","channel_id"]).any():
        raise ValueError("Duplicated or missing island x pollinator channel")
    if not ch.channel_id.isin(GUILDS).all() or not ch.channel_state.isin(
        ("retained","disrupted","unresolved")).all():
        raise ValueError("Invalid pollinator channel/status")
    if (ch.channel_state.eq("disrupted")&ch.source_state.ne("available")).any():
        raise ValueError("Adequate nondetection requires a verified source opportunity")
    if len(ch)!=len(qualified)*len(GUILDS):
        raise ValueError("All five guild channel evidence states required")
    states=ch.pivot(index="island_id",columns="channel_id",values="channel_state")
    if set(states.index)!=qualified or set(states.columns)!=set(GUILDS) or states.isna().any().any():
        raise ValueError("Unmatched five-guild evidence cohort")
    d=_require(globi_reference_genera,
        ("genus","n_independent_references",*("refs__"+c for c in GUILDS)),
        "global GloBI genus source documentation")
    if d.genus.isna().any() or d.genus.duplicated().any():
        raise ValueError("Duplicate GloBI documented source genus")
    metric=["n_independent_references",*("refs__"+c for c in GUILDS)]
    d[metric]=d[metric].apply(pd.to_numeric,errors="coerce")
    if not np.isfinite(d[metric].to_numpy(float)).all() or (d[metric]<0).any().any():
        raise ValueError("Invalid source reference counts")
    d=d[d.n_independent_references>=min_globi_refs]
    refcols=["refs__"+c for c in GUILDS]
    d=d[d[refcols].sum(axis=1)>0]
    refs=d.set_index("genus")[refcols]
    refs=refs.div(refs.sum(axis=1),axis=0)
    refs.columns=GUILDS

    per_island={}
    candidates=[]
    for island,plant in island_f.groupby("island_id",sort=True):
        source_names=set()
        for entity in source_sets[island]:
            source_names.update(gift_sets.get(entity,set()))
        qualified_source=source_names.intersection(x.index)
        if len(qualified_source)<10:
            continue
        source_form=x.loc[list(qualified_source)]
        island_form=x.loc[plant.accepted_species]
        bg=states.loc[island,"bombus"]
        bird=states.loc[island,"flower_visiting_birds"]
        label=("bombus_undetected_birds_recorded" if
               (bg=="disrupted" and bird=="retained") else
               "both_bombus_and_birds_recorded" if
               (bg=="retained" and bird=="retained") else
               "other_or_unresolved")
        doc_island=refs.reindex(plant.genus).dropna(how="all")
        doc_source=refs.reindex([name.split()[0] for name in qualified_source]).dropna(how="all")
        record={
            "analysis_regime":plant.analysis_regime.iloc[0],
            "spatial_block":plant.spatial_block.iloc[0],
            "group":label,
            "n_direct_species":len(plant),
            "n_direct_source_pool_species":len(source_form),
            "n_same_source_species":int(plant.accepted_species.isin(source_names).sum()),
            "n_native_direct_species":int(plant.origin_status.eq("native").sum()),
            "n_native_same_source_species":int((plant.origin_status.eq("native")&
                                                plant.accepted_species.isin(source_names)).sum()),
        }
        for name in ("open_form","actinomorphy","inconspicuous_colour"):
            record["source_relative_"+name]=float(
                island_form[name].mean()-source_form[name].mean())
        for guild in GUILDS:
            record["documented_source_relative_"+guild] = (
                float(doc_island[guild].mean()-doc_source[guild].mean())
                if len(doc_island)>=5 and len(doc_source)>=10 else np.nan
            )
        per_island[island]=record
        if label=="bombus_undetected_birds_recorded":
            candidates.append((plant,source_names))
    by=pd.DataFrame.from_dict(per_island,orient="index")
    rows=[]
    for (reg,label),sub in by.groupby(["analysis_regime","group"],sort=True):
        if len(sub)<min_group_islands or sub.spatial_block.nunique()<min_group_spatial_blocks:
            continue
        row=dict(region=reg,evidence_group=label,n_islands=int(len(sub)),
                 n_spatial_blocks=int(sub.spatial_block.nunique()),
                 median_n_complete_direct_species=float(sub.n_direct_species.median()),
                 median_n_complete_GIFT_source_species=float(sub.n_direct_source_pool_species.median()),
                 median_n_same_named_source_species=float(sub.n_same_source_species.median()),
                 median_n_native_direct_species=float(sub.n_native_direct_species.median()),
                 share_islands_with_native_same_source_plant=float((sub.n_native_same_source_species>0).mean()))
        for col in [v for v in by if v.startswith(("source_relative_","documented_source_relative_"))]:
            nonmissing=sub[col].dropna()
            row["median_"+col]=float(nonmissing.median()) if len(nonmissing) else np.nan
            row["fraction_positive_"+col]=float((nonmissing>0).mean()) if len(nonmissing) else np.nan
            row["n_complete_"+col]=int(len(nonmissing))
        rows.append(row)
    out=pd.DataFrame.from_records(rows)
    all_cand=pd.concat([plant for plant,_ in candidates],ignore_index=True) if candidates else island_f.iloc[:0]
    stat={
       "n_islands_atleast_30_direct_complete_species":len(qualified),
       "n_bombus_nondetected_bird_recorded":len(candidates),
       "n_candidate_island_species_records":int(len(all_cand)),
       "n_distinct_candidate_species":int(all_cand.accepted_species.nunique()),
       "n_distinct_candidate_genera":int(all_cand.genus.nunique()),
       "n_recorded_native_island_species":int(all_cand.origin_status.eq("native").sum()),
       "n_unresolved_island_origin_records":int(all_cand.origin_status.eq("unresolved").sum()),
       "n_same_named_source_species_records":int(sum(
           plant.accepted_species.isin(source_names).sum() for plant,source_names in candidates)),
       "n_native_same_named_source_species_records":int(sum(
           ((plant.origin_status=="native")&plant.accepted_species.isin(source_names)).sum()
           for plant,source_names in candidates)),
       "source_pool_is_confirmed_historical_ancestor":False,
       "direct_species_traits_have_independent_island_populations":False,
       "local_animal_guild_detection_is_specific_flower_visitation":False,
       "global_GloBI_genus_documentation_is_local_pollinator_abundance":False,
       "confirmed_same_species_evolutionary_pollinator_shift_in_data":False,
    }
    return out,stat
