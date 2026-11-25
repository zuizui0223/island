"""Original global GIFT floral Direct documentation × GBIF background effort.

GBIF channel states alone are not a comparable pollinator community survey:
independent original NEE records additionally carry background occurrence
counts, independent spatial/temporal background units and GBIF datasets.
This audits those documentation opportunities across ALL recorded GIFT
island floras without requiring all-three floral traits or recorded guild
states for every island.

No original island/plant/GIFT/GBIF species/observer/source identities export.
No biological guild absence inferred from a missing background field.
This is NOT estimation of real floral visitation, pollen transfer,
causal guild loss or ecological selection. Correlations are purely descriptive.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .global_partial_trait_flora_fd import strict_single_trait_original
from .global_partial_pollinator_guild_bridge import (
    normalize_guild_evidence,CHANNELS
)

TRAITS=(
  "flower_primary_color","floral_form","floral_symmetry",
  "autonomous_selfing_capacity","inflorescence_display",
  "flower_primary_color__floral_form",
)
EFFORT=(
  "background_record_count","background_spatial_units",
  "background_temporal_units","distinct_dataset_count",
)
RECORD_BIN_EDGES=(0,10,50,200,1000,np.inf)
REGIMES=("GLOBAL_ALL","tropical","northern_midlatitude",
         "northern_high_latitude","southern_extratropical")


def _safe_summary(x:pd.Series)->dict:
    v=pd.to_numeric(x,errors="coerce").dropna().to_numpy(float)
    v=v[np.isfinite(v)]
    if not len(v):return dict(n_nonmissing=0)
    return dict(n_nonmissing=int(len(v)),
                median=float(np.median(v)),
                q10=float(np.quantile(v,.1)),
                q90=float(np.quantile(v,.9)))


def _rank_corr(x:pd.Series,y:pd.Series)->dict|None:
    values=pd.DataFrame({"x":pd.to_numeric(x,errors="coerce"),
                         "y":pd.to_numeric(y,errors="coerce")}).dropna()
    values=values[np.isfinite(values.x)&np.isfinite(values.y)]
    if len(values)<30 or values.x.nunique()<5 or values.y.nunique()<5:
        return None
    rank_x=values.x.rank().to_numpy(float)
    rank_y=values.y.rank().to_numpy(float)
    rho=float(np.corrcoef(rank_x,rank_y)[0,1])
    if not np.isfinite(rho):return None
    return dict(n_original_islands_with_both_effort_and_trait_coverage=len(values),
                spearman_descriptive_rho=rho,
                original_BioGBIF_recording_effort_is_not_standardized_sampling=True)


def _bins(x:pd.Series)->pd.Series:
    return pd.cut(pd.to_numeric(x,errors="coerce"),
       bins=RECORD_BIN_EDGES,labels=False,include_lowest=True,right=False)


def original_global_flora_guild_effort_audit(
    flora:pd.DataFrame,ledger:pd.DataFrame,
    geography:pd.DataFrame,channel:pd.DataFrame
)->dict:
    fields={"island_id","accepted_species","origin_status"}
    if not fields.issubset(flora) or flora.duplicated(
        ["island_id","accepted_species"]).any():
        raise ValueError("Original source GIFT island membership invalid")
    if not flora.origin_status.isin(("native","introduced","unresolved")).all():
        raise ValueError("Original native plant provenance category changed")
    n=flora.island_id.nunique()
    if n<4000:
        raise ValueError("Source scope accidentally collapsed to 928 full triple flora")
    if not {"island_id","analysis_regime","geological_origin","spatial_block"}.issubset(geography):
        raise ValueError("Independent original island geography not present")
    if geography.island_id.duplicated().any():
        raise ValueError("Original independent island geography duplicated")
    req={"island_id","channel_id","channel_state","source_state",*EFFORT}
    if not req.issubset(channel):
        raise ValueError("Original NEE background evidence fields unavailable")
    original=normalize_guild_evidence(channel, set(flora.island_id))
    # Previous provenance gate intentionally keeps ONLY state fields.
    # Reattach independently audited effort columns without inventing
    # background counts for original islands not covered by a source row.
    effort_fields=channel[["island_id","channel_id",*EFFORT]].copy()
    original=original.merge(effort_fields,on=["island_id","channel_id"],
                            how="left",validate="one_to_one")
    for e in EFFORT:
        original[e]=pd.to_numeric(original[e],errors="coerce")
        if (original[e].dropna()<0).any():
            raise ValueError("Original GBIF effort includes negative values")
    # Do NOT replace absent background counts with 0, ever.
    all_names=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    native_names=flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    direct=strict_single_trait_original(ledger)
    bytrait={t:set(direct[t]) for t in direct}
    bytrait["flower_primary_color__floral_form"]=(
        bytrait["flower_primary_color"]&bytrait["floral_form"])
    records=[]
    for island in sorted(all_names):
        row=dict(island_id=island)
        full=all_names[island]
        ind=native_names.get(island,set())
        for trait in TRAITS:
            state=bytrait[trait]
            row["all__"+trait]=len(full&state)/len(full)
            row["native__"+trait]=(
                len(ind&state)/len(ind) if len(ind) else float("nan"))
        records.append(row)
    meta=pd.DataFrame({"island_id":sorted(all_names)}).merge(
        geography[["island_id","analysis_regime","geological_origin"]],
        on="island_id",how="left",validate="one_to_one")
    meta["analysis_regime"]=meta.analysis_regime.fillna("geography_unresolved")
    meta["geological_origin"]=meta.geological_origin.fillna("unresolved")
    obs=original.merge(meta,on="island_id",validate="many_to_one")
    obs=obs.merge(pd.DataFrame(records),on="island_id",validate="many_to_one")
    if len(obs)!=n*len(CHANNELS):
        raise ValueError("Original GIFT flora × GBIF five channel grid truncated")
    panels=[]
    for region in REGIMES:
        allg=obs if region=="GLOBAL_ALL" else obs.loc[
            obs.analysis_regime.eq(region)]
        if allg.empty:continue
        for geological in ("all","oceanic","continental"):
            g=allg if geological=="all" else allg.loc[
                allg.geological_origin.eq(geological)]
            if g.empty:continue
            for guild,one in g.groupby("channel_id"):
                detected=one.loc[one.channel_state.eq("retained")]
                notdetected=one.loc[one.channel_state.eq("disrupted")]
                unresolved=one.loc[one.channel_state.eq("unresolved")]
                bins_d=set(_bins(detected.background_record_count).dropna().astype(int))
                bins_n=set(_bins(notdetected.background_record_count).dropna().astype(int))
                a=dict(
                    region=region,geological_origin=geological,
                    guild=guild,n_all_original_botanical_islands=int(len(one)),
                    n_GBIF_guild_recorded=int(len(detected)),
                    n_guild_source_qualified_NONdetection=int(len(notdetected)),
                    n_guild_unresolved=int(len(unresolved)),
                    n_guild_state_source_missing_not_a_nondetection=int(
                        (~one.source_channel_row_present).sum()),
                    n_background_record_effort_bins_with_both_qualified_states=len(bins_d&bins_n),
                    n_background_record_effort_bins_detected=len(bins_d),
                    n_background_record_effort_bins_qualified_nondetected=len(bins_n),
                )
                for col in EFFORT:
                    a["recorded_"+col]=_safe_summary(detected[col])
                    a["qualified_nondetection_"+col]=_safe_summary(notdetected[col])
                    a["all_evaluable_"+col]=_safe_summary(
                        one.loc[one.channel_state.isin(
                            ("retained","disrupted")),col])
                a["direct_trait_coverage_correlations_with_log1p_GBIF_background_records"]=[]
                effort=one.loc[
                    one.source_state.eq("available")&
                    one.channel_state.isin(("retained","disrupted"))].copy()
                effort["effort_log"]=np.log1p(effort.background_record_count)
                for scope in ("all","native"):
                    for trait in TRAITS:
                        cov=scope+"__"+trait
                        r=_rank_corr(effort["effort_log"],effort[cov])
                        if r:
                            a["direct_trait_coverage_correlations_with_log1p_GBIF_background_records"].append(
                                dict(botanical_flora_status=scope,trait=trait,**r))
                panels.append(a)
    global_bombus=next((x for x in panels if
        x["region"]=="GLOBAL_ALL" and x["geological_origin"]=="all" and
        x["guild"]=="bombus"),None)
    if global_bombus is None or global_bombus["n_all_original_botanical_islands"]!=n:
        raise ValueError("GBIF background denominator no longer global source flora")
    return dict(
        status="ORIGINAL_GLOBAL_4505_DIRECT_TRAIT_AND_GBIF_BACKGROUND_EFFORT_AUDIT",
        n_original_botanical_islands=n,
        n_original_GBIF_island_guild_cells=len(obs),
        n_observation_effort_fields=len(EFFORT),
        original_plant_flora_and_GBIF_background_effort_by_group=panels,
        scientific_limits=dict(
            globally_fixed_original_botanical_flora_NOT_triple_complete=True,
            adequate_nondetection_NOT_historical_pollinator_loss=True,
            GBIF_background_record_count_NOT_equal_standardized_survey_effort=True,
            missing_background_value_not_inferred_zero=True,
            species_direct_trait_not_local_flower_photo=True,
            correlations_descriptive_not_causal_or_multiplicity_controlled=True,
            independent_geological_origin_only_when_verified=True,
            no_botanical_plant_island_animal_or_GIFT_source_names_exported=True,
        ),
    )
