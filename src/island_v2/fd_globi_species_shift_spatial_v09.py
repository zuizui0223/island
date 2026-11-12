"""Spatially verify GloBI v0.9 geocoded exact-plant-species records on
the SAME 21 historically frozen GSHHG candidate island polygon identities.

Original georeferenced GloBI relations are NOT verified pollen transfer
or same-site historical guild shifts; original v0.9 has no event time.
Never infer island locality from global genus documentation, centroid
radius, a substitute shoreline dataset, or a mismatched island polygon.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd


GUILDS=("bombus","non_bombus_bees","lepidoptera",
        "flower_visiting_birds","diptera")
PIN="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"


def file_sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1048576),b""):
            h.update(block)
    return h.hexdigest()


def original_candidate_island_locality_gate(
    recorded_georeferenced_relations:set[tuple],
    bibliographic_providers_by_reference:dict,
    native_plant_island_pairs:pd.DataFrame,
    verified_geometry_path:Path,
    geometry_receipt_path:Path,
    *, dominant_source_archive_label:str|None=None,
)->dict:
    """One-to-one original GSHHG polygon membership preflight.

    Return anonymous counts of cited species×visitor relationships whose
    original georeference lies in the SAME exact 13 candidate islands
    where the accepted plant species is recorded native and source-matched.
    """
    import geopandas as gpd
    from shapely.geometry import Point

    if not {"island_id","accepted_species"}.issubset(native_plant_island_pairs):
        raise ValueError("Missing protected original species×island source-match key")
    p=native_plant_island_pairs[["island_id","accepted_species"]].copy()
    if p.isna().any().any() or p.duplicated().any():
        raise ValueError("Duplicate or unresolved original candidate native plant pairs")
    if (p.accepted_species.nunique()!=54 or
        p.island_id.nunique()!=13 or len(p)!=197):
        raise ValueError("Original 54-species / 13-island / 197-pair cohort drift")
    if not verified_geometry_path.is_file() or not geometry_receipt_path.is_file():
        raise ValueError("Missing exact original geometry artifact and receipt")
    receipt=json.loads(geometry_receipt_path.read_text())
    if (receipt.get("GSHHG_archive_SHA256")!=PIN
        or receipt.get("n_original_focal_islands")!=21
        or receipt.get("n_focal_islands_with_exact_historical_polygon")!=21
        or receipt.get("status") not in (
            "FOCAL_21_POLYGONS_IDENTICAL_BUT_FULL_UNIVERSE_DRIFT",
            "FULL_8265_GLOBAL_GEOMETRY_IDENTICAL")
        or receipt.get("geometry_file_SHA256")!=file_sha256(verified_geometry_path)):
        raise ValueError("Exact GSHHG original-ID/source-hash verification failed")
    polygons=gpd.read_file(verified_geometry_path,layer="islands")
    if (polygons.crs is None or polygons.crs.to_epsg()!=4326 or
        polygons.island_id.isna().any() or polygons.island_id.duplicated().any()):
        raise ValueError("Source original geometry not unique EPSG:4326")
    focal_ids=set(p.island_id)
    if not focal_ids.issubset(set(polygons.island_id)):
        raise ValueError("One or more of the original 13 source-matched focal islands missing geometry")
    # Only exact polygons of the focal original native island records.
    focal=polygons[polygons.island_id.isin(focal_ids)][["island_id","geometry"]]
    if len(focal)!=13 or focal.geometry.isna().any() or focal.geometry.is_empty.any():
        raise ValueError("Original focal island geometry is incomplete")

    georefs=pd.DataFrame.from_records(
        list(recorded_georeferenced_relations),
        columns=("accepted_species","guild","reference_id","provider_id",
                 "evidence_class","latitude","longitude"))
    output={
        "status":"EXACT_FOCAL_ISLAND_POLYGON_OVERLAP_AUDITED_NOT_VERIFIED_POLLINATOR_SHIFT",
        "n_frozen_candidate_species":54,
        "n_strict_native_source_matched_plant_island_rows":197,
        "n_qualified_focal_islands_with_same_native_source_names":13,
        "n_approved_original_candidate_island_polygon_IDs":21,
        "n_original_GloBI_georeferenced_species_provider_reference_records":
            int(len(georefs)),
        "n_source_georeferenced_rows_inside_any_original_13_island_polygon":0,
        "n_source_georeferenced_rows_matching_native_species_on_exact_same_island":0,
        "n_distinct_focal_islands_with_record_matching_native_species":0,
        "n_distinct_original_native_species_with_spatially_matched_guild_record":0,
        "n_same_island_plant_species_with_2plus_documented_guilds":0,
        "n_same_island_plant_species_with_3plus_independent_cited_references":0,
        "n_original_explicit_pollination_relations_inside_island_with_native_species":0,
        "n_distinct_original_georeferenced_site_coordinates_matching_native_plant":0,
        "n_distinct_original_bibliographic_reference_keys_matching_native_plant":0,
        "n_distinct_original_archive_first_labels_matching_native_plant":0,
        "n_original_same_island_plant_pairs_with_2plus_archive_labels":0,
        "n_original_same_island_plant_pairs_with_2plus_distinct_coordinate_sites":0,
        "n_exact_same_island_native_plant_relation_rows_after_globally_dominant_archive_removed":
            0 if dominant_source_archive_label else None,
        "n_native_plant_species_with_record_after_globally_dominant_archive_removed":
            0 if dominant_source_archive_label else None,
        "n_islands_with_record_after_globally_dominant_archive_removed":
            0 if dominant_source_archive_label else None,
        "dominant_archive_identifier_never_exported":True,
        "source_label_removal_is_not_leave_one_independent_study_out":True,
        "n_effective_pollen_deposition_or_seedset_measurements_in_input":0,
        "n_confirmed_same_species_mainland_to_island_pollinator_shifts":0,
        "n_original_mainland_GIFT_source_population_visitor_assessments":0,
        "original_GSHHG_archive_SHA256":PIN,
        "geometry_artifact_SHA256":receipt["geometry_file_SHA256"],
        "original_gshhg_FULL_8265_universe_verified":
            bool(receipt.get("global_full_universe_identity_verified")),
        "locality_analyses_scoped_to_only_original_21_focal_IDs":True,
        "GloBI_original_event_dates_available":False,
        "original_species_island_provider_reference_or_geo_coordinates_exported":False,
        "geospatial_relation_does_not_prove_pollination_effectiveness":True,
        "guilds":[],
    }
    if georefs.empty:
        for g in GUILDS:
            output["guilds"].append({
                "guild":g,
                "n_source_native_focal_plant_species_with_georeferenced_relation":0,
                "n_source_native_focal_islands_with_georeferenced_relation":0,
            })
        return output
    if (georefs.isna().any().any() or
        not georefs.guild.isin(GUILDS).all() or
        (~georefs.evidence_class.isin(("flower_visit","pollination"))).any()
        or not np.isfinite(georefs[["longitude","latitude"]].to_numpy(float)).all()
        or georefs.longitude.abs().gt(180).any()
        or georefs.latitude.abs().gt(90).any()):
        raise ValueError("Invalid original cited species-visit georeference state")
    multi={ref for ref,prov in bibliographic_providers_by_reference.items()
           if len(prov)>1}
    georefs=georefs[~georefs.reference_id.isin(multi)].copy()
    points=gpd.GeoDataFrame(georefs,
            geometry=gpd.points_from_xy(georefs.longitude,georefs.latitude),
            crs="EPSG:4326")
    overlap=gpd.sjoin(points,focal,how="inner",predicate="within")
    # A coordinate may fall in overlapping original landmask polygons.
    # Such ambiguous site evidence is unqualified, never counted twice.
    n_per_site=overlap.groupby(overlap.index).island_id.nunique()
    ambiguous=set(n_per_site[n_per_site>1].index)
    overlap=overlap[~overlap.index.isin(ambiguous)].copy()
    output["n_source_georeferenced_rows_inside_any_original_13_island_polygon"]=int(len(overlap))
    if overlap.empty:
        for g in GUILDS:
            output["guilds"].append({
                "guild":g,
                "n_source_native_focal_plant_species_with_georeferenced_relation":0,
                "n_source_native_focal_islands_with_georeferenced_relation":0,
            })
        return output
    matched=overlap.merge(p,on=["island_id","accepted_species"],
                          how="inner",validate="many_to_one")
    output["n_source_georeferenced_rows_matching_native_species_on_exact_same_island"]=len(matched)
    output["n_distinct_focal_islands_with_record_matching_native_species"]=int(matched.island_id.nunique())
    output["n_distinct_original_native_species_with_spatially_matched_guild_record"]=int(matched.accepted_species.nunique())
    if not matched.empty:
        plant_island=matched.groupby(["island_id","accepted_species"])
        output["n_distinct_original_georeferenced_site_coordinates_matching_native_plant"]=int(
            matched[["island_id","longitude","latitude"]].drop_duplicates().shape[0])
        output["n_distinct_original_bibliographic_reference_keys_matching_native_plant"]=int(
            matched.reference_id.nunique())
        output["n_distinct_original_archive_first_labels_matching_native_plant"]=int(
            matched.provider_id.nunique())
        output["n_original_same_island_plant_pairs_with_2plus_archive_labels"]=int(
            plant_island.provider_id.nunique().ge(2).sum())
        output["n_original_same_island_plant_pairs_with_2plus_distinct_coordinate_sites"]=int(
            plant_island.apply(lambda frame:frame[["latitude","longitude"]].drop_duplicates().shape[0],
                               include_groups=False).ge(2).sum())
        output["n_same_island_plant_species_with_2plus_documented_guilds"]=int(
            plant_island.guild.nunique().ge(2).sum())
        output["n_same_island_plant_species_with_3plus_independent_cited_references"]=int(
            plant_island.reference_id.nunique().ge(3).sum())
        output["n_original_explicit_pollination_relations_inside_island_with_native_species"]=int(
            matched.evidence_class.eq("pollination").sum())
    if dominant_source_archive_label is not None:
        rest=matched.loc[matched.provider_id.ne(dominant_source_archive_label)]
        output["n_exact_same_island_native_plant_relation_rows_after_globally_dominant_archive_removed"]=int(len(rest))
        output["n_native_plant_species_with_record_after_globally_dominant_archive_removed"]=int(rest.accepted_species.nunique())
        output["n_islands_with_record_after_globally_dominant_archive_removed"]=int(rest.island_id.nunique())
    for g in GUILDS:
        part=matched[matched.guild.eq(g)]
        rest_g=part[part.provider_id.ne(dominant_source_archive_label)] if dominant_source_archive_label is not None else None
        output["guilds"].append({
            "guild":g,
            "n_source_native_focal_plant_species_with_georeferenced_relation":
                int(part.accepted_species.nunique()),
            "n_source_native_focal_islands_with_georeferenced_relation":
                int(part.island_id.nunique()),
            "n_source_native_focal_plant_species_after_largest_archive_removed":
                int(rest_g.accepted_species.nunique()) if rest_g is not None else None,
        })
    return output
