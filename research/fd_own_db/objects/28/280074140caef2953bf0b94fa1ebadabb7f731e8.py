"""Bounded iNaturalist ORIGINAL photo-opportunity pilot for island endemic congeners.

The target frame is source-verified native+endemic plants with original
HIGH/MEDIUM species-direct flower colour, a different accepted GIFT
mainland source congener with source-direct colour, and an endemic
congener on both geologically oceanic and continental tropical islands.

Original plant source IDs, geolocation, API taxon/photo/observation IDs
never leave the GitHub protected CI. This does not measure flower pixels,
whether photos display a flower, pollen transfer, within-species island
phenotype or historical plant adaptation.

Taxon choice outcome-blind: fixed SHA sorting over original source
geographic IDs, genus and plant names prior to any iNaturalist API query.
Provider query failures are NEVER interpreted as zero photographs.
"""
from __future__ import annotations
from collections import Counter,defaultdict
from datetime import date
from hashlib import sha256
from time import sleep
from typing import Callable
from urllib.parse import urlencode
from urllib.request import Request,urlopen
import json
import math

import pandas as pd


SALT="fd-genuine-oc-endemic-photos-first-8-paired-genera-20261114"
LICENSES={"cc0","cc-by","cc-by-sa","cc-by-nc","cc-by-nc-sa"}
API_ROOT="https://api.inaturalist.org/v1"
QUALIFIED_MAX_ACCURACY=5000
MAX_PAIRED_GENERA=13


def _genus(taxon:str)->str:
    x=str(taxon).strip().split()
    return x[0] if len(x)>=2 else ""


def _seed(*parts:object)->str:
    return sha256((SALT+"|"+ "|".join(str(x) for x in parts)).encode()).hexdigest()


def select_original_endemic_colour_pairs(
    ledger:pd.DataFrame,flora:pd.DataFrame,geology:pd.DataFrame,
    assignments:pd.DataFrame,mainland:pd.DataFrame,*,
    max_pairs:int=MAX_PAIRED_GENERA
)->tuple[list[dict],dict]:
    if not 1<=max_pairs<=MAX_PAIRED_GENERA:
        raise ValueError("Original bounded first wave exceeds frozen 13 original shared genera")
    needed={"accepted_species","trait_name","normalized_value",
            "evidence_scope","evidence_quality"}
    if not needed.issubset(ledger):
        raise ValueError("Original species-direct source missing")
    phen=ledger.loc[
        ledger.trait_name.eq("flower_primary_color") &
        ledger.evidence_scope.eq("species_direct") &
        ledger.evidence_quality.isin(("high","medium")) &
        ledger.normalized_value.notna() &
        ~ledger.normalized_value.astype(str).str.contains("|",regex=False)
    ]
    if phen.accepted_species.duplicated().any():
        raise ValueError("Original direct colour labels ambiguous")
    trait_plant_names=set(phen.accepted_species)
    if not {"island_id","accepted_species","origin_status",
            "endemic_status","floristic_status"}.issubset(flora):
        raise ValueError("Original independent island endemism missing")
    if flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original island accepted source species duplicate")
    if ((flora.floristic_status.eq("endemic")&
         (~flora.origin_status.eq("native")|
          ~flora.endemic_status.eq("endemic"))).any()):
        raise ValueError("Original endemic status unsupported by native island evidence")
    if not {"island_id","geological_origin","analysis_regime",
            "spatial_block","weigelt_archipelago"}.issubset(geology):
        raise ValueError("Independent geological type and archipelago missing")
    if geology.island_id.duplicated().any():
        raise ValueError("Original island/geology ID repeated")
    eligible=geology.loc[
        geology.geological_origin.isin(("oceanic","continental"))&
        geology.analysis_regime.eq("tropical")]
    targets=flora.loc[
        flora.island_id.isin(eligible.island_id)&
        flora.floristic_status.eq("endemic")&
        flora.accepted_species.isin(trait_plant_names)
    ].merge(eligible[["island_id","geological_origin","analysis_regime",
                      "spatial_block","weigelt_archipelago"]],
        on="island_id",validate="many_to_one",how="inner")
    near=assignments.loc[assignments.source_mode.eq("geo_k5")].copy()
    near["source_rank"]=pd.to_numeric(near.source_rank,errors="raise")
    near=near.loc[near.source_rank.between(1,5)]
    if near.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original GIFT source rank duplicate")
    near_map=near.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    by_source=mainland[["entity_ID","work_species"]].drop_duplicates().groupby(
        "entity_ID").work_species.agg(set).to_dict()
    observed=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    candidates=defaultdict(lambda:defaultdict(list))
    for row in targets.itertuples(index=False):
        refs=near_map.get(str(row.island_id),())
        if len(refs)!=5 or len(set(refs))!=5:continue
        sourced=(set().union(*(by_source.get(x,set()) for x in refs))-
                 observed[str(row.island_id)])
        genus=_genus(row.accepted_species)
        matching={x for x in sourced&trait_plant_names if _genus(x)==genus}
        if not matching or not genus:continue
        candidates[genus][row.geological_origin].append({
            "island_id":str(row.island_id),"plant":str(row.accepted_species),
            "genus":genus,"origin":str(row.geological_origin),
            "spatial_block":str(row.spatial_block),
            "archipelago":str(row.weigelt_archipelago or ""),
        })
    # No FCP/pollination/photo data or floral COLOUR STATE is consulted here.
    shared=sorted(g for g,origins in candidates.items()
                  if "oceanic" in origins and "continental" in origins)
    available=[]
    for genus in shared:
        by=[]
        for origin in ("oceanic","continental"):
            sorted_options=sorted(candidates[genus][origin],key=lambda x:(
                x["archipelago"] in ("","nan","unknown","None"),
                _seed(genus,origin,x["plant"],x["island_id"])))
            by.append(sorted_options[0])
        if by[0]["plant"]==by[1]["plant"]:
            continue
        available.append((genus,by))
    # Historical SALT from first eight remains unchanged: selecting 13
    # keeps exactly the same first eight source pairs and adds five more.
    # Freeze selection by original source opportunity, not provider-photo
    # availability or the actual direct flower colour category.
    available.sort(key=lambda x:_seed("genus",x[0]))
    pairs=available[:max_pairs]
    selected=[]
    for genus,own in pairs:
        for item in own:
            selected.append(dict(item))
    return selected,dict(
        n_original_tropical_shared_genera_with_direct_colour_and_GIFT_congener=len(shared),
        n_original_distinct_source_pair_eligible_genera=len(available),
        n_selected_original_distinct_genera=len(pairs),
        n_selected_original_accepted_endemic_plant_names=len({x["plant"] for x in selected}),
        n_selected_geography_blocks_oceanic=len({x["spatial_block"]
            for x in selected if x["origin"]=="oceanic"}),
        n_selected_geography_blocks_continental=len({x["spatial_block"]
            for x in selected if x["origin"]=="continental"}),
        n_verified_archipelagos_used=len({x["archipelago"] for x in selected
            if x["archipelago"] not in ("","nan","unknown","None")}),
        selection_before_provider_query_does_not_use_colour_value=True,
        accepted_endemic_photo_name_not_in_fcp_4475_assumed=False,
    )


def official_api_json(endpoint:str,params:dict,*,delay:float=1.12)->dict:
    url=API_ROOT+"/"+endpoint.lstrip("/")+"?"+urlencode(params)
    req=Request(url,headers={
        "User-Agent":"zuizui0223-fd-oceanic-endemic-photo-source/0.1 "
                     "(github.com/zuizui0223/fd; research provenance)",
        "Accept":"application/json"})
    try:
        with urlopen(req,timeout=40) as stream:
            data=json.load(stream)
        if not isinstance(data,dict) or not isinstance(data.get("results"),list):
            raise RuntimeError("Official iNaturalist JSON schema changed")
        return data
    finally:
        sleep(max(delay,1.05))


def photo_observation_qualifies(obs:dict,exact_taxon_id:int,accepted_name:str)->bool:
    t=obs.get("taxon") or {}
    if (str(t.get("id"))!=str(exact_taxon_id) or
        str(t.get("name") or "").casefold()!=accepted_name.casefold() or
        str(t.get("rank") or "")!="species" or
        obs.get("quality_grade")!="research" or obs.get("obscured") or
        obs.get("geoprivacy") not in (None,"","open") or
        obs.get("taxon_geoprivacy") not in (None,"","open")):
        return False
    coordinates=(obs.get("geojson") or {}).get("coordinates")
    if (not isinstance(coordinates,list) or len(coordinates)!=2):
        return False
    try:
        longitude,latitude=map(float,coordinates)
        accuracy=float(obs.get("positional_accuracy"))
        recorded=str(obs["observed_on"])
        date.fromisoformat(recorded)
        if not (math.isfinite(longitude) and math.isfinite(latitude)
                and -180<=longitude<=180 and -90<=latitude<=90 and
                0<=accuracy<=QUALIFIED_MAX_ACCURACY):
            return False
    except (TypeError,ValueError,KeyError):
        return False
    return any(str(p.get("license_code") or "").lower() in LICENSES and
        str(p.get("id") or "").isdigit() and
        str(p.get("url") or "").startswith("https://")
        for p in (obs.get("photos") or []))


def query_exact_source_taxon(
    accepted_name:str, client:Callable[[str,dict],dict]
)->dict:
    try:
        taxa=client("taxa",{"q":accepted_name,"rank":"species","per_page":30,
                            "is_active":"true"})
    except Exception:
        return dict(provider_status="provider_error",exact_id=None,
                    n_qualified_licensed_open_geophoto_on_first_page=0)
    names=[r for r in taxa.get("results",[]) if
           str(r.get("name") or "").casefold()==accepted_name.casefold()
           and r.get("rank")=="species" and r.get("is_active") is not False]
    if len(names)!=1 or not str(names[0].get("id") or "").isdigit():
        return dict(provider_status="no_unambiguous_exact_public_taxon",
                    exact_id=None,n_qualified_licensed_open_geophoto_on_first_page=0)
    taxon_id=int(names[0]["id"])
    params={
        "taxon_id":taxon_id,"quality_grade":"research",
        "photos":"true","geo":"true","per_page":30,"page":1,
        "order_by":"id","order":"asc",
        "photo_license":",".join(sorted(LICENSES)),
    }
    try:
        observations=client("observations",params)
    except Exception:
        return dict(provider_status="provider_error",exact_id=taxon_id,
                    n_qualified_licensed_open_geophoto_on_first_page=0)
    first=observations.get("results",[])
    good=sum(photo_observation_qualifies(o,taxon_id,accepted_name) for o in first)
    total=observations.get("total_results")
    return dict(provider_status="queried_exact_original_public_taxon",
                exact_id=taxon_id,
                global_provider_photographed_georeferenced_research_observations=(
                   int(total) if isinstance(total,int) and total>=0 else None),
                n_qualified_licensed_open_geophoto_on_first_page=int(good),
                n_observation_candidates_first_page=len(first))


def original_photo_pair_pilot(
    selected:list[dict],selection_audit:dict,*,
    client:Callable[[str,dict],dict]=official_api_json
)->dict:
    if len(selected)<2 or len(selected)>2*MAX_PAIRED_GENERA:
        raise ValueError("Original bounded 13-genus/two-geology pilot source size invalid")
    if len({x["plant"] for x in selected})!=len(selected):
        raise ValueError("Selected botanical accepted plant names duplicated")
    private=[]
    for src in selected:
        rec=query_exact_source_taxon(src["plant"],client)
        private.append({**src,**rec})
    stats={}
    for origin in ("oceanic","continental"):
        x=[z for z in private if z["origin"]==origin]
        stats[origin]=dict(
            source_selected_endemic_accepted_taxa=len(x),
            n_original_genus_groups=len(set(z["genus"] for z in x)),
            n_unambiguous_public_inat_exact_species=sum(
                z["exact_id"] is not None for z in x),
            n_inat_taxa_queried_for_photos=sum(
                z["provider_status"]=="queried_exact_original_public_taxon" for z in x),
            n_globally_qualified_photo_taxa_from_first_page=sum(
                z.get("n_qualified_licensed_open_geophoto_on_first_page",0)>0
                for z in x),
            total_first_page_qualified_photo_observation_records=sum(
                z.get("n_qualified_licensed_open_geophoto_on_first_page",0) for z in x),
            n_taxa_without_unambiguous_exact_inat_species=sum(
                z["provider_status"]=="no_unambiguous_exact_public_taxon" for z in x),
            n_provider_error_taxa=sum(
                z["provider_status"]=="provider_error" for z in x),
        )
    if sum(stats[x]["n_provider_error_taxa"] for x in stats)>2:
        raise RuntimeError(
           "Public API unstable; do not claim source lack/zero photos from provider errors")
    bygenus=defaultdict(dict)
    for row in private:
        bygenus[row["genus"]][row["origin"]]=row
    n_ready=sum(
        all(bygenus[g].get(o,{}).get(
            "n_qualified_licensed_open_geophoto_on_first_page",0)>0
            for o in ("oceanic","continental"))
        for g in bygenus)
    return dict(
        status="GENUINE_PUBLIC_INAT_EXACT_SOURCE_TAXON_GEOPHOTO_PAIRED_PILOT",
        original_selection_before_photo_query=selection_audit,
        provider_query_counts_by_geology=stats,
        n_shared_genera_with_eligible_public_exact_photo_on_both_geologies=n_ready,
        source_query_protocol=dict(
            official_iNaturalist_v1_taxa_exact_active_species=True,
            official_iNaturalist_v1_observations_first_page_research_geophoto=True,
            client_1request_per_second_pacing_required=True,
            licenses_cc0_by_by_sa_by_nc_by_nc_sa_only=True,
            precise_open_geolocation_5000m_limit=True,
            provider_taxon_lookup_ambiguous_is_NOT_species_absence=True,
            no_original_photo_pixels_downloaded=True,
            no_genuine_flower_photo_ROI_measured=True,
            photo_may_show_leaves_not_flowers=True,
            broad_global_public_inat_photo_presence_NOT_island_locality=True,
            historical_photo_provider_snapshot_NOT_immutable=True,
            restricted_taxon_and_observation_IDs_NOT_exported=True,
            no_ecological_or_evolutionary_adaptation_effect_measured=True,
        ),
    )
