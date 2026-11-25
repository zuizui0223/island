"""Public iNaturalist photo records NEAR independently typed endemic source islands.

This only checks public, non-obscured, currently served photographic
metadata near the frozen GSHHG island centroid. A circle is NOT the
original island polygon and may contain other islands and mainland.
Photo location is not plant native status, observed flowering stage,
a local adaptive phenotype, or measured pollen transfer.

Original island IDs, accepted endemic names, exact coordinates,
taxon/observation/photo IDs never leave the source CI.
"""
from __future__ import annotations
from collections import defaultdict
from math import asin,cos,isfinite,pi,radians,sin,sqrt
import numpy as np

from .oceanic_endemic_inat_photo_pilot import (
    LICENSES,original_photo_pair_pilot,query_exact_source_taxon,
    official_api_json,photo_observation_qualifies
)


def _distance_km(lat1:float,lon1:float,lat2:float,lon2:float)->float:
    a,b,c,d=map(radians,(lat1,lon1,lat2,lon2))
    x=sin((c-a)/2)**2+cos(a)*cos(c)*sin((d-b)/2)**2
    return float(2*6371.0088*asin(min(1,sqrt(max(0,x)))))


def source_island_centroid_radius(area_km2:float)->float:
    """Conservative first-pass circle, explicitly NOT a land polygon."""
    a=float(area_km2)
    if not isfinite(a) or a<=0:
        raise ValueError("GSHHG island area not known")
    return float(min(75.,max(5.,1.5*sqrt(a/pi))))


def public_near_observation_source_gate(
    obs:dict,taxon_id:int,name:str,*,latitude:float,longitude:float,
    radius_km:float
)->bool:
    if not photo_observation_qualifies(obs,taxon_id,name):
        return False
    if obs.get("taxon_geoprivacy") not in (None,"","open"):
        return False
    coords=(obs.get("geojson") or {}).get("coordinates")
    try:
        lat=float(coords[1]);lon=float(coords[0])
    except (TypeError,ValueError,IndexError,KeyError):
        return False
    return _distance_km(latitude,longitude,lat,lon)<=radius_km


def original_near_island_photo_provider_pilot(
    selected:list[dict],selection_audit:dict,*,
    client=official_api_json
)->dict:
    if len(selected)!=26 or len({x["plant"] for x in selected})!=26:
        raise ValueError("Original full 13-genus source pair support must remain fixed")
    if len({x["genus"] for x in selected})!=13:
        raise ValueError("Original 13 paired endemic botanical genera not observed")
    private=[]
    for row in selected:
        name=str(row["plant"])
        t=query_exact_source_taxon(name,client)
        local=dict(
            **row,
            source_taxon_status=t["provider_status"],
            exact_public_taxon_found=t.get("exact_id") is not None,
            global_first_page_licensed_geo_photo=t.get(
                "n_qualified_licensed_open_geophoto_on_first_page",0)>0,
            near_source_centroid_open_photo=False,
            n_near_source_centroid_open_photo_observations_first_page=0,
            local_provider_error=False,
        )
        if t.get("exact_id") is not None and t["provider_status"]==(
            "queried_exact_original_public_taxon"):
            radius=source_island_centroid_radius(row["island_area_km2"])
            params={
                "taxon_id":t["exact_id"],"photos":"true","geo":"true",
                "quality_grade":"research","per_page":30,"page":1,
                "order_by":"id","order":"asc","lat":float(row["island_latitude"]),
                "lng":float(row["island_longitude"]),"radius":round(radius,3),
                "photo_license":",".join(sorted(LICENSES)),
            }
            try:
                resp=client("observations",params)
                n=sum(public_near_observation_source_gate(
                    o,t["exact_id"],name,latitude=float(row["island_latitude"]),
                    longitude=float(row["island_longitude"]),radius_km=radius)
                    for o in resp.get("results",[]))
                local["near_source_centroid_open_photo"]=bool(n)
                local["n_near_source_centroid_open_photo_observations_first_page"]=int(n)
            except Exception:
                local["local_provider_error"]=True
        private.append(local)
    errors=sum(
        r["source_taxon_status"]=="provider_error" or r["local_provider_error"]
        for r in private)
    if errors>2:
        raise RuntimeError("iNaturalist provider failed too often; do not infer photo unavailability")
    bytype={}
    for typ in ("oceanic","continental"):
        one=[x for x in private if x["origin"]==typ]
        bytype[typ]=dict(
            n_original_endemic_accepted_plants_queried=len(one),
            n_source_taxa_unambiguously_matched_to_public_iNat=sum(
                x["exact_public_taxon_found"] for x in one),
            n_global_first_page_licensed_photo_taxa=sum(
                x["global_first_page_licensed_geo_photo"] for x in one),
            n_public_near_island_centroid_qualified_photo_taxa=sum(
                x["near_source_centroid_open_photo"] for x in one),
            n_public_near_island_centroid_qualified_firstpage_observation_records=sum(
                x["n_near_source_centroid_open_photo_observations_first_page"]
                for x in one),
            n_photo_provider_errors=sum(
                x["source_taxon_status"]=="provider_error" or x["local_provider_error"]
                for x in one),
            n_original_distinct_10degree_island_blocks=int(
                len({x["spatial_block"] for x in one})),
        )
    grouped=defaultdict(dict)
    for r in private:
        grouped[r["genus"]][r["origin"]]=r
    n_both=sum(
        all(grouped[genus].get(typ,{}).get("near_source_centroid_open_photo")
            for typ in ("oceanic","continental"))
        for genus in grouped
    )
    return dict(
        status="PUBLIC_INAT_ENDEMIC_CENTROID_RADIUS_LOCALITY_PROXY_PREFLIGHT",
        original_target_selection=selection_audit,
        original_public_geo_photo_provider_by_island_type=bytype,
        n_matched_genera_with_centroid_near_public_photo_in_both_geologies=n_both,
        provider_query_safeguards=dict(
            raw_GSHHG_island_polygon_NOT_USED=True,
            centroid_circle_not_evidence_of_island_presence=True,
            public_open_only_no_obscured_geoprivacy_taxon_or_observation=True,
            approximate_island_radius_min5_max75_km=True,
            maximum_open_public_accuracy_m=5000,
            photo_license_checked_for_CC0_and_selected_CC_terms=True,
            botanical_name_exact_iNaturalist_species_ID_required=True,
            source_taxon_ancestor_not_identified=True,
            geographic_photo_phenotype_not_measured=True,
            no_original_photo_pixels_downloaded=True,
            no_taxon_island_photo_ID_or_precise_coordinates_exported=True,
            provider_taxon_or_photo_absence_does_not_imply_biological_absence=True,
        ),
    )
