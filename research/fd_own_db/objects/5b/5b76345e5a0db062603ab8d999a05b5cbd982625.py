"""Conservative paired island-locality and GIFT-mainland candidate photo support.

Two independent metadata observations are required per pair:
1) an exact endemic species photograph, public/unobscured, within a
   source-island CENTROID radius (NOT original polygon verified), and
2) an exact different source-GIFT native-congener with a licensed public
   photograph SOMEWHERE globally (NOT mainland-locality verified).

DO NOT mistake this as two original floral ROI measurements, adaptation
or an ancestral population comparison. Query species regardless of
whether an earlier worldwide first page happened to include a photo.
Original taxon/island/coordinates/private photo IDs never leave CI.
"""
from __future__ import annotations
from collections import defaultdict
from typing import Callable

from .oceanic_endemic_inat_photo_pilot import (
    LICENSES,official_api_json,query_exact_source_taxon
)
from .oceanic_endemic_near_island_photo import (
    source_island_centroid_radius,public_near_observation_source_gate
)


def original_endemic_centroid_vs_GIFT_mainland_global_photo(
    pairs:list[dict],source_audit:dict,*,
    client:Callable[[str,dict],dict]=official_api_json
)->dict:
    if len(pairs)!=26 or len({x["plant"] for x in pairs})!=26:
        raise ValueError("Not original 26 taxon source pair list")
    expected={x["genus"] for x in pairs}
    if len(expected)!=13 or any({x["origin"] for x in pairs
       if x["genus"]==g}!={"oceanic","continental"} for g in expected):
        raise ValueError("Original paired geology genus identity drift")
    for row in pairs:
        if not {"island_latitude","island_longitude","island_area_km2",
                "spatial_block","archipelago"}.issubset(row):
            raise ValueError("Original verified island geographic data missing")
    source_names=sorted({x["botanical_source_congener"] for x in pairs})
    # Distinct original botanical source plant names queried exactly once.
    source={n:query_exact_source_taxon(n,client) for n in source_names}
    private=[]
    for p in pairs:
        target=query_exact_source_taxon(p["plant"],client)
        supplier=source[p["botanical_source_congener"]]
        radius=source_island_centroid_radius(p["island_area_km2"])
        local_n=0
        local_error=False
        if target["provider_status"]=="queried_exact_original_public_taxon":
            params=dict(
                taxon_id=target["exact_id"],photos="true",geo="true",
                quality_grade="research",per_page=30,page=1,
                order_by="id",order="asc",
                lat=float(p["island_latitude"]),lng=float(p["island_longitude"]),
                radius=round(radius,3),
                photo_license=",".join(sorted(LICENSES)),
            )
            try:
                res=client("observations",params)
                if not isinstance(res.get("results"),list):
                    raise ValueError("Source iNaturalist georeferenced response schema invalid")
                local_n=sum(public_near_observation_source_gate(
                    obs,target["exact_id"],p["plant"],
                    latitude=float(p["island_latitude"]),
                    longitude=float(p["island_longitude"]),
                    radius_km=radius,
                ) for obs in res["results"])
            except Exception:
                local_error=True
        donor_global=(supplier["provider_status"]==
                      "queried_exact_original_public_taxon" and
                      supplier.get("n_qualified_licensed_open_geophoto_on_first_page",0)>0)
        local_ok=(local_n>0)
        private.append(dict(
            genus=p["genus"],origin=p["origin"],spatial_block=p["spatial_block"],
            archipelago=p["archipelago"],island_id=p["island_id"],
            donor_name=p["botanical_source_congener"],target_name=p["plant"],
            donor_photo_globally=bool(donor_global),
            endemic_photo_near_centroid=bool(local_ok),
            source_first_GIFT_region_rank=int(p["gift_rank"]),
            n_source_island_centroid_radius_photo_obs=local_n,
            target_exact=target["exact_id"] is not None,
            supplier_exact=supplier["exact_id"] is not None,
            target_error=(target["provider_status"]=="provider_error" or local_error),
            supplier_error=supplier["provider_status"]=="provider_error",
        ))
    errors=sum(x["target_error"] or x["supplier_error"] for x in private)
    if errors>3:
        raise RuntimeError("Public photo API errors too frequent for reliable source coverage")
    bygeology={}
    for geology in ("oceanic","continental"):
        rows=[x for x in private if x["origin"]==geology]
        valid=[x for x in rows if x["endemic_photo_near_centroid"] and
               x["donor_photo_globally"]]
        bygeology[geology]=dict(
            original_paired_island_endemic_congener_taxa=len(rows),
            n_exact_endemic_iNaturalist_accepted_species=sum(x["target_exact"] for x in rows),
            n_exact_GIFT_mainland_source_congener_species=sum(x["supplier_exact"] for x in rows),
            n_endemic_taxa_with_open_licensed_photo_near_original_island_centroid=sum(
                x["endemic_photo_near_centroid"] for x in rows),
            n_GIFT_native_mainland_candidate_congener_species_with_global_open_photo=sum(
                x["donor_photo_globally"] for x in rows),
            n_pairs_with_centroid_near_endemic_photo_AND_global_GIFT_congener_photo=len(valid),
            n_distinct_10degree_blocks_with_two_sided_photo_opportunity=len({
                x["spatial_block"] for x in valid}),
            n_distinct_verified_archipelagos_with_two_sided_photo_opportunity=len({
                x["archipelago"] for x in valid
                if x["archipelago"] not in ("","nan","unknown","None")}),
            n_pairs_with_both_photo_and_nearest_first_GIFT_source_region=sum(
                x["source_first_GIFT_region_rank"]==1 for x in valid),
            n_local_centroid_near_qualified_endemic_observation_records_on_two_sided_pairs=sum(
                x["n_source_island_centroid_radius_photo_obs"] for x in valid),
            n_provider_error_on_endemic_or_GIFT_candidate=sum(
                x["target_error"] or x["supplier_error"] for x in rows),
        )
    both_genera=sum(all(
        any(x["genus"]==genus and x["origin"]==geo and
            x["endemic_photo_near_centroid"] and x["donor_photo_globally"]
            for x in private)
        for geo in ("oceanic","continental")) for genus in expected)
    return dict(
        status="ORIGINAL_ENDEMIC_CENTROID_NEAR_PLUS_SOURCE_CONGENER_GLOBAL_PHOTO_FEASIBILITY",
        original_source_matching_metadata=source_audit,
        original_protected_13_genus_26_paired_plants=26,
        n_original_different_accepted_GIFT_source_congener_species_queried=len(source),
        support_by_geological_island_class=bygeology,
        n_shared_genera_with_centroid_near_photo_and_source_congener_global_photo_in_BOTH_classes=both_genera,
        qualification=dict(
            source_island_centroid_radius_NOT_polygon_presence=True,
            source_mainland_GIFT_checklist_not_actual_ancestral_source=True,
            source_congener_photo_GLOBALLY_not_proven_mainland_local=True,
            all_public_taxa_and_images_exact_species_only=True,
            original_accepted_endemic_and_GIFT_donor_names_not_exported=True,
            precise_geolocation_and_observer_and_photo_identifiers_not_exported=True,
            no_flower_pixels_or_image_colour_ROI_measured=True,
            no_historical_colonization_or_adaptation_or_pollinator_effect_identified=True,
        ),
    )
