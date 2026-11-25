"""GIFT global 4,505 island flora × independently sourced 5 GBIF guild channels.

Unlike legacy 928-island all-three-floral-traits analysis, this consumes
original island×ONE-trait or TWO-trait DIRECT phenotypes (>=5 accepted
species) and retains all original island floras as coverage denominators.
Missing GBIF source-channel rows are unresolved; NEVER non-detection.

'retained' = qualified island GBIF recorded pollinator channel.
'disrupted' = source-available and effort-qualified non-detection ONLY,
not observed historical decline, lost pollinator or guild replacement.

Original botanically accepted species, original island IDs, GBIF source
locations, GIFT source entities and detailed spatial units remain
private inside original island CI. The return value has NO such IDs.
"""
from __future__ import annotations

from hashlib import sha256

import numpy as np
import pandas as pd

from .global_partial_trait_flora_fd import TRAITS,PAIR_TRAITS,N_RAO

CHANNELS=("bombus","non_bombus_bees","lepidoptera",
          "flower_visiting_birds","diptera")
VALID_STATES=("retained","disrupted","unresolved","structurally_absent")
SOURCE_STATES=("available","unresolved","structurally_absent")
PRIMARY=(
    "flower_primary_color",
    "floral_form",
    "floral_symmetry",
    "autonomous_selfing_capacity",
    "flower_primary_color__floral_form",
    "flower_primary_color__floral_symmetry",
    "floral_form__floral_symmetry",
)
# A priori descriptive tests; no p/q or causal claims without independent
# correction, phylogenetic and detection-effort modeling.
TARGETS=("observed_finite5","observed_genus_equal_finite5",
         "source_relative_finite5","source_relative_genus_equal_finite5",
         "observed_two_trait_pairing","source_relative_pairing")
GATES=dict(min_state_islands=10,min_state_blocks=5,min_state_sources=5)


def _seed(*parts:object)->int:
    key="|".join(map(str,parts))
    return int.from_bytes(sha256(key.encode()).digest()[:8],"big")


def normalize_guild_evidence(table:pd.DataFrame,all_islands:set[str])->pd.DataFrame:
    expected={"island_id","channel_id","channel_state","source_state"}
    if not expected.issubset(table):
        raise ValueError("Original independently qualified GBIF guild source schema missing")
    q=table[list(expected)].copy()
    if q.isna().any().any() or q.duplicated(["island_id","channel_id"]).any():
        raise ValueError("Original GBIF source channel cannot have missing/duplicate identity")
    if not q.channel_id.isin(CHANNELS).all():
        raise ValueError("Unexpected pollinator guild in independently qualified source")
    if not q.channel_state.isin(VALID_STATES).all():
        raise ValueError("Invalid source GBIF observation state")
    if not q.source_state.isin(SOURCE_STATES).all():
        raise ValueError("Invalid GIFT candidate source guild availability state")
    if ((q.channel_state.eq("disrupted") &
         q.source_state.ne("available"))).any():
        raise ValueError("Claimed adequate non-detection without actual available guild source")
    if (q.channel_state.eq("structurally_absent") &
        q.source_state.ne("structurally_absent")).any():
        raise ValueError("No independently proved structurally absent source guild")
    # Keep unknown source rows as unknown, not nonexistent guilds.
    source_ids=set(q.island_id)
    if len(source_ids&all_islands)<1:
        raise ValueError("Original source guild evidence contains NO GIFT flora islands")
    missing=sorted(all_islands-source_ids)
    # Canonical complete five-channel island×guild grid, not five source
    # observations fabricated as absent when a source scan never ran.
    island_frame=pd.DataFrame({"island_id":sorted(all_islands)})
    channel_frame=pd.DataFrame({"channel_id":list(CHANNELS)})
    full=island_frame.merge(channel_frame,how="cross")
    z=full.merge(q,on=["island_id","channel_id"],how="left",
                 validate="one_to_one",indicator=True)
    z["source_channel_row_present"]=z["_merge"].eq("both")
    z["channel_state"]=z.channel_state.fillna("unresolved")
    z["source_state"]=z.source_state.fillna("unresolved")
    z=z.drop(columns="_merge")
    if len(z)!=len(all_islands)*len(CHANNELS):
        raise AssertionError("Original full 5-guild support grid cannot be created")
    if z.duplicated(["island_id","channel_id"]).any():
        raise AssertionError("GBIF guild evidence duplicated after exact original island join")
    return z


def _contrast(q:pd.DataFrame,metric:str,group:str,draws:int,seed:int)->dict|None:
    obs=q.loc[q.channel_state.eq("retained")&q[metric].notna()].groupby(
        group)[metric].mean().dropna().to_numpy(float)
    non=q.loc[q.channel_state.eq("disrupted")&q[metric].notna()].groupby(
        group)[metric].mean().dropna().to_numpy(float)
    if min(len(obs),len(non))<5:
        return None
    rng=np.random.default_rng(seed)
    mean=float(obs.mean()-non.mean())
    a=obs[rng.integers(len(obs),size=(draws,len(obs)))].mean(axis=1)
    b=non[rng.integers(len(non),size=(draws,len(non)))].mean(axis=1)
    return dict(
        observed_GBIF_channel_minus_qualified_nondetection=mean,
        ci95=[float(z) for z in np.quantile(a-b,[.025,.975])],
        n_observed_GBIF_channel_clusters=len(obs),
        n_qualified_nondetection_clusters=len(non),
    )


def source_guild_partial_traits(
    private_island_scores:pd.DataFrame,
    channel_source:pd.DataFrame,*,
    seed:int=20261114,draws:int=999
)->dict:
    if draws<999:
        raise ValueError("Original block/source uncertainty draws fewer than 999")
    needed={"island_id","regime","geological_origin",
            "original_scope","outcome","outcome_type",
            "spatial_block","source_signature",
            "n_original_flora_species","n_recorded_direct_species",
            "observed_direct_coverage",*TARGETS}
    if not needed.issubset(private_island_scores):
        raise ValueError("Global 4,505 island partial-trait original FD private schema missing")
    data=private_island_scores.copy()
    n_original=data.island_id.nunique()
    if n_original<4000:
        raise ValueError("Original global botanical island universe was narrowed to legacy 928")
    if data.duplicated(["island_id","original_scope","outcome"]).any():
        raise ValueError("Global Direct floral island×scope×trait source duplicated")
    endpoints=set(TRAITS)|{"__".join(x) for x in PAIR_TRAITS}
    if set(data.outcome)!=endpoints:
        raise ValueError("One of ten Direct traits or three observed pairs missing")
    if set(data.original_scope)!={"all_recorded","original_native"}:
        raise ValueError("Original GIFT botanically native vs all-recorded scopes conflated")
    if len(data)!=n_original*2*len(endpoints):
        raise ValueError("Original 4,505 island roster must stay intact for every trait")
    if data.n_recorded_direct_species.lt(0).any() or (
        data.n_original_flora_species.lt(0).any()):
        raise ValueError("Original Direct evidence cannot be negative")
    g=normalize_guild_evidence(channel_source,set(data.island_id))
    joined=data.merge(
        g,on="island_id",how="left",validate="many_to_many")
    if len(joined)!=len(data)*len(CHANNELS):
        raise ValueError("Complete original guild qualification grid drifted")
    # Classify true source-context independently of phenotype, and retain
    # unclassified island geology as global botanical study observations.
    joined["analysis_zone"]="global"
    island_keys=data[["island_id","regime","geological_origin"]].drop_duplicates()
    if island_keys.island_id.duplicated().any():
        raise ValueError("Island geographic metadata changed by trait/scope")
    meta=joined[["island_id","regime","geological_origin"]].drop_duplicates()
    n_geology=int(meta.geological_origin.isin(("oceanic","continental")).sum())
    n_original_channel_source=int(g.loc[
        g.source_channel_row_present,"island_id"].nunique())
    panels=[]
    for region in ["GLOBAL_ALL",*sorted(
        str(x) for x in meta.regime.dropna().unique()
        if x not in ("", "geography_unresolved"))]:
        regionframe=(joined if region=="GLOBAL_ALL"
            else joined.loc[joined.regime.eq(region)])
        if regionframe.empty:
            continue
        for geology in ("all","oceanic","continental"):
            if geology=="all":
                frame=regionframe
            else:
                frame=regionframe.loc[
                    regionframe.geological_origin.eq(geology)]
            if frame.empty:
                continue
            for (scope,trait,channel),sub in frame.groupby(
                ["original_scope","outcome","channel_id"],sort=True
            ):
                if len(sub)<10:
                    continue
                sufficient=sub.n_recorded_direct_species.ge(N_RAO)
                observed=sub.channel_state.eq("retained")
                nondetected=sub.channel_state.eq("disrupted")
                sourceknown=sub.source_state.eq("available")
                unresolved=sub.channel_state.eq("unresolved")
                missing=~sub.source_channel_row_present
                row=dict(
                    region=region,geological_origin=geology,
                    original_flora_status=scope,trait_or_pair=trait,
                    pollinator_guild=channel,
                    n_ALL_islands_in_this_original_scope=int(len(sub)),
                    n_original_islands_with_any_direct_trait=int(
                        sub.n_recorded_direct_species.gt(0).sum()),
                    n_islands_with_5plus_direct_species=int(sufficient.sum()),
                    n_islands_original_pollinator_guild_row_present=int(
                        sub.source_channel_row_present.sum()),
                    n_islands_GBIF_guild_row_MISSING_not_absent=int(missing.sum()),
                    n_guild_source_available=int(sourceknown.sum()),
                    n_guild_GBIF_observation_retained=int(observed.sum()),
                    n_guild_qualified_adequate_NONdetection=int(nondetected.sum()),
                    n_guild_unresolved=int(unresolved.sum()),
                    n_source_structurally_absent=int(
                        sub.source_state.eq("structurally_absent").sum()),
                    n_GBIF_observed_and_5plus_direct=int((observed&sufficient).sum()),
                    n_GBIF_adequately_nondetected_and_5plus_direct=int(
                        (nondetected&sufficient).sum()),
                    # From the same original botanical island cohort,
                    # coverage always uses all original island plant members.
                    mean_GBIF_observed_direct_record_coverage=float(
                        sub.loc[observed,"observed_direct_coverage"].mean())
                        if sub.loc[observed,"observed_direct_coverage"].notna().any() else None,
                    mean_GBIF_qualified_nondetection_direct_record_coverage=float(
                        sub.loc[nondetected,"observed_direct_coverage"].mean())
                        if sub.loc[nondetected,"observed_direct_coverage"].notna().any() else None,
                    source_count_not_true_pollinator_richness=True,
                )
                if trait in PRIMARY:
                    for metric in TARGETS:
                        # Some endpoints are inapplicable by design; do not
                        # turn missing joint/trait into zero ecological effect.
                        o=sub.loc[observed&sub[metric].notna()]
                        d=sub.loc[nondetected&sub[metric].notna()]
                        blocks=(o.spatial_block.nunique(),
                                d.spatial_block.nunique())
                        sources=(o.source_signature.nunique(),
                                 d.source_signature.nunique())
                        allowed=(len(o)>=GATES["min_state_islands"] and
                                 len(d)>=GATES["min_state_islands"] and
                                 min(blocks)>=GATES["min_state_blocks"] and
                                 min(sources)>=GATES["min_state_sources"])
                        if not allowed:
                            continue
                        record=dict(
                            metric=metric,
                            n_GBIF_observed_islands=int(len(o)),
                            n_GBIF_qualified_nondetection_islands=int(len(d)),
                            n_spatial_blocks_by_state=list(map(int,blocks)),
                            n_nearest_GIFT_sources_by_state=list(map(int,sources)),
                            block=_contrast(
                                sub,metric,"spatial_block",draws,
                                _seed(seed,region,geology,scope,trait,channel,metric,"block")),
                            first_GIFT_source=_contrast(
                                sub,metric,"source_signature",draws,
                                _seed(seed,region,geology,scope,trait,channel,metric,"source")),
                        )
                        if record["block"] and record["first_GIFT_source"]:
                            row.setdefault("exploratory_guild_floral_trait_associations",[]
                                           ).append(record)
                panels.append(row)
    n_tests=sum(len(p.get("exploratory_guild_floral_trait_associations",[])) for p in panels)
    return dict(
        status="GLOBAL_4505_PARTIAL_DIRECT_TRAIT_FD_GBIF_GUILD_EVIDENCE_BRIDGE",
        n_original_global_botanical_flora_islands=n_original,
        n_original_botanical_islands_with_independent_geological_type=n_geology,
        n_original_guild_evidence_source_islands_present=n_original_channel_source,
        n_original_botanical_islands_without_ANY_guild_source_row=(
            n_original-n_original_channel_source),
        n_original_plant_scope_trait_endpoints=len(data),
        n_original_plant_scope_trait_GBIF_channel_rows=len(joined),
        n_guilds=len(CHANNELS),
        n_original_individual_trait_or_pair_endpoints=len(endpoints),
        n_anonymous_region_geology_scope_trait_guild_panels=len(panels),
        n_exploratory_GBIF_channel_vs_qualified_nondetection_endpoint_contrasts=n_tests,
        regional_and_geological_pollinator_guild_opportunity=panels,
        inference_limits=dict(
            source_guild_retained_means_GBIF_detection_not_visit_or_pollen_transfer=True,
            source_disrupted_means_ADEQUATE_nondetection_not_historical_loss=True,
            source_unresolved_and_no_source_rows_NEVER_treated_as_absence=True,
            original_4505_botanical_flora_membership_not_narrowed_to_928=True,
            ten_individual_traits_and_three_pairs_not_triple_required=True,
            n5_per_trait_Direct_only_and_no_missing_value_imputation=True,
            separate_original_native_vs_all_recorded_flora=True,
            organism_detection_effort_source_and_trait_documentation_confounding=True,
            same_accepted_plant_traits_fixed_across_original_islands=True,
            descriptive_block_and_GIFT_cluster_intervals_not_FDR_corrected=True,
            GIFT_mainland_regions_are_not_ancestral_pollinator_source=True,
            no_proven_pollen_service_pollinator_shift_or_adaptation=True,
            no_source_plant_island_guild_GBIF_identity_ids_exported=True,
        ),
    )
