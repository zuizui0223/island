"""Island endemic vs independent candidate MAINLAND source CONGENER photo support.

A GIFT native mainland source checklist provides potential botanical
source opportunity, NOT direct propagule arrival, the true ancestor or
an iNaturalist mainland-local photograph. Photo-provider errors are never
treated as negative evidence. Original status, source name, accepted
species, island IDs and observations stay inside source-protected CI.

Select before querying photo availability or consulting focal/donor
flower primary colour VALUES. The accepted island species is entirely
excluded from its donor candidate pool. Each selected island endemic
is paired with exactly one source-listed different accepted congener,
with fixed near-source-rank and deterministic SHA tie-breaking.
"""
from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
from typing import Callable

import pandas as pd

from .oceanic_endemic_inat_photo_pilot import (
    LICENSES,MAX_PAIRED_GENERA,original_photo_pair_pilot,
    select_original_endemic_colour_pairs,query_exact_source_taxon,
    official_api_json
)

SEED="fd-island-native-endemic-versus-GIFT-mainland-source-genus-20261114-v1"


def _genus(name:str)->str:
    parts=str(name).split()
    return parts[0] if len(parts)>=2 else ""


def _hash(*tokens:object)->str:
    return sha256((SEED+"|"+ "|".join(map(str,tokens))).encode()).hexdigest()


def source_independent_mainland_controls(
    selected:list[dict],ledger:pd.DataFrame,
    flora:pd.DataFrame,assignments:pd.DataFrame,mainland:pd.DataFrame
)->tuple[list[dict],dict]:
    if len(selected)!=2*MAX_PAIRED_GENERA or len({q["plant"] for q in selected})!=len(selected):
        raise ValueError("Source original 13×2 accepted endemic target names not frozen")
    source_valid=ledger.loc[
        ledger.trait_name.eq("flower_primary_color")&
        ledger.evidence_scope.eq("species_direct")&
        ledger.evidence_quality.isin(("high","medium"))&
        ledger.normalized_value.notna()&
        ~ledger.normalized_value.astype(str).str.contains("|",regex=False)
    ]
    if source_valid.accepted_species.duplicated().any():
        raise ValueError("Direct colour reference has duplicate source species")
    trait_names=set(source_valid.accepted_species)
    if not {"island_id","source_mode","source_rank","entity_ID"}.issubset(assignments):
        raise ValueError("Original ranked GIFT sources missing")
    near=assignments.loc[assignments.source_mode.eq("geo_k5")].copy()
    near["source_rank"]=pd.to_numeric(near.source_rank,errors="raise")
    near=near.loc[near.source_rank.between(1,5)]
    if near.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original source rank collision")
    mapping=near.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    by_source=mainland[["entity_ID","work_species"]].drop_duplicates().groupby(
        "entity_ID").work_species.agg(set).to_dict()
    if flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original island plant records duplicate")
    island_roster=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    pairs=[]
    for q in selected:
        island=str(q["island_id"])
        genus=_genus(q["plant"])
        if q["plant"] not in island_roster.get(island,set()):
            raise ValueError("Original endemic target absent from accepted island flora")
        src=mapping.get(island,())
        if len(src)!=5 or len(set(src))!=5:
            raise ValueError("Original source candidate five GIFT regions changed")
        taken=None
        for rank,entity in enumerate(src,1):
            pool=(by_source.get(entity,set()) & trait_names)-island_roster[island]
            possible=sorted(name for name in pool if _genus(name)==genus and
                            name!=q["plant"])
            if not possible:continue
            chosen=min(possible,key=lambda name:(_hash(q["plant"],rank,name),name))
            taken=dict(
                **q,
                gift_rank=rank,
                botanical_source_congener=chosen,
                candidate_species_count_in_chosen_GIFT_source=len(possible),
            )
            break
        if taken is None:
            raise ValueError("Original source-qualified endemic target lost every GIFT donor")
        pairs.append(taken)
    if len(pairs)!=26 or len({x["plant"] for x in pairs})!=26:
        raise ValueError("Original complete 13×2 island endemic and candidate source pairs lost")
    stats={}
    for origin in ("oceanic","continental"):
        rows=[x for x in pairs if x["origin"]==origin]
        stats[origin]=dict(
            n_source_endemic_accepted_plants=len(rows),
            n_distinct_source_botanical_genus=len({_genus(x["plant"]) for x in rows}),
            n_distinct_candidate_mainland_accepted_congeners=len(
                {x["botanical_source_congener"] for x in rows}),
            n_original_GIFT_first_region_as_source=sum(x["gift_rank"]==1 for x in rows),
            median_GIFT_native_candidate_species_options=float(pd.Series([
                x["candidate_species_count_in_chosen_GIFT_source"] for x in rows]).median()),
        )
    return pairs,dict(
        original_selection_source_only_pre_image_pre_colour=True,
        all_source_controls_are_independently_accepted_name_not_focal_island_species=True,
        n_original_same_genus_pairings=len(pairs),
        n_distinct_mainland_candidate_species=len({x["botanical_source_congener"] for x in pairs}),
        original_candidate_source_rank_statistics_by_geology=stats,
        original_GIFT_native_candidate_source_not_historical_ancestor=True,
    )


def original_mainland_congener_public_photo_gate(
    pairs:list[dict],selected_audit:dict,*,
    client:Callable[[str,dict],dict]=official_api_json
)->dict:
    if len(pairs)!=26 or len({x["plant"] for x in pairs})!=26:
        raise ValueError("Source original island endemic-paired geography sample changed")
    # Query each different GIFT accepted mainland congener once, even if it
    # serves several different endemic source plant comparisons.
    source_names=sorted({x["botanical_source_congener"] for x in pairs})
    source_results={name:query_exact_source_taxon(name,client) for name in source_names}
    # This original frame was frozen before source/photo success was queried.
    # Recheck endemic side for a same-snapshot availability gate.
    target_names=sorted({x["plant"] for x in pairs})
    target_results={name:query_exact_source_taxon(name,client) for name in target_names}
    stats={}
    for origin in ("oceanic","continental"):
        rows=[x for x in pairs if x["origin"]==origin]
        derived=[]
        for row in rows:
            photo=source_results[row["botanical_source_congener"]]
            endemic=target_results[row["plant"]]
            derived.append(dict(
                genus=row["genus"],
                original_source_rank=row["gift_rank"],
                source_taxon_exact=photo["exact_id"] is not None,
                source_public_any_photo=photo.get(
                    "n_qualified_licensed_open_geophoto_on_first_page",0)>0,
                endemic_public_any_photo=endemic.get(
                    "n_qualified_licensed_open_geophoto_on_first_page",0)>0,
                source_provider_error=photo["provider_status"]=="provider_error",
                endemic_provider_error=endemic["provider_status"]=="provider_error",
            ))
        stats[origin]=dict(
            n_frozen_source_endemic_congener_pairs=len(rows),
            n_mainland_GIFT_source_congeners_exact_iNat_taxon=sum(
                x["source_taxon_exact"] for x in derived),
            n_mainland_GIFT_source_congeners_with_global_public_photo=sum(
                x["source_public_any_photo"] for x in derived),
            n_source_genus_pairs_with_endemic_and_mainland_global_photo=sum(
                x["source_public_any_photo"] and x["endemic_public_any_photo"]
                for x in derived),
            n_source_genus_pairs_with_endemic_global_photo=sum(
                x["endemic_public_any_photo"] for x in derived),
            n_mainland_source_provider_errors=sum(
                x["source_provider_error"] for x in derived),
            n_endemic_provider_errors=sum(
                x["endemic_provider_error"] for x in derived),
            n_genera_with_both_global_photo_and_mainland_source_rank_1=sum(
                x["source_public_any_photo"] and x["endemic_public_any_photo"]
                and x["original_source_rank"]==1 for x in derived),
        )
    if sum(v["n_mainland_source_provider_errors"]+
           v["n_endemic_provider_errors"] for v in stats.values())>3:
        raise RuntimeError("Provider error rate invalidates source photo opportunity")
    return dict(
        status="ORIGINAL_13_GENUS_MAINLAND_CONGENER_AND_ISLAND_ENDEMIC_PHOTO_PREFLIGHT",
        frozen_original_source_tips=selected_audit,
        public_photo_opportunity_by_geological_island_type=stats,
        n_unique_source_mainland_accepted_congener_names_queried=len(source_names),
        n_unique_original_endemic_accepted_names_queried=len(target_names),
        source_scope=dict(
            original_island_and_source_plant_species_Direct_colour_qualified=True,
            GIFT_nearest_mainland_list_is_not_ancestral_population=True,
            source_photo_is_GLOBAL_not_confirmed_mainland_locality=True,
            endemic_photo_is_GLOBAL_not_confirmed_source_island_locality=True,
            source_photo_may_not_show_flower=True,
            licensed_geoprivacy_open_research_grade_photo_metadata_only=True,
            no_original_photo_pixels_fetched=True,
            no_original_flower_colour_ROI_measured=True,
            no_botanical_endemic_source_or_public_taxon_identifiers_exported=True,
            no_pollinator_mediated_adaptation_inferred=True,
        )
    )
