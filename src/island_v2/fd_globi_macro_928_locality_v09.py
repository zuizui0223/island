"""Macro-scale exact-plant-species local GloBI documentation coverage.

A rigorously restricted **coverage preflight** over the frozen 928
Direct-colour/form/symmetry-complete island floras. Global GloBI source
relationships must have original bibliography/provider, coordinates within
the historical exact GSHHG polygon, and the SAME accepted plant name in that
island flora. This DOES NOT demonstrate pollen transfer, a local pollination
syndrome, historical change, or actual island-mainland pollinator shift.

The GSHHG 8,265-ID universe has 9 expired-source hash mismatches. Every
nonreproducible polygon is marked UNRESOLVED; it is not silently treated
as a biological absence. Original identities are never returned.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import json

import pandas as pd
import numpy as np

from island_v2.fd_globi_species_shift_spatial_v09 import file_sha256, PIN

TRAITS=("flower_primary_color","floral_form","floral_symmetry")
GUILDS=("bombus","non_bombus_bees","lepidoptera",
        "flower_visiting_birds","diptera")
ACCEPTED_ORIGINS={"native","introduced","unresolved","endemic"}


def strict_direct_928_flora(
        original_direct:pd.DataFrame, original_flora:pd.DataFrame
)->tuple[pd.DataFrame,dict]:
    required={"accepted_species","trait_name","normalized_value",
              "evidence_scope","evidence_quality"}
    if not required.issubset(original_direct):
        raise ValueError("Missing original Direct floral source contract")
    if not {"island_id","accepted_species","origin_status"}.issubset(original_flora):
        raise ValueError("Original island flora/origin evidence missing")
    d=original_direct.loc[
        original_direct.trait_name.isin(TRAITS) &
        original_direct.evidence_scope.eq("species_direct") &
        original_direct.evidence_quality.isin(("high","medium")) &
        original_direct.normalized_value.notna() &
        ~original_direct.normalized_value.astype(str).str.contains("|",regex=False)
    ]
    if d.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Ambiguous original species-direct floral trait")
    x=d.pivot(index="accepted_species",columns="trait_name",
              values="normalized_value")
    if not set(TRAITS).issubset(x):
        raise ValueError("All 3 original direct trait categories required")
    complete=x.dropna(subset=list(TRAITS))
    flora=original_flora.loc[
        original_flora.accepted_species.isin(complete.index),
        ["island_id","accepted_species","origin_status"]].copy()
    if flora.isna().any().any() or flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Incomplete or duplicate island-by-plant original flora")
    if not flora.origin_status.isin(ACCEPTED_ORIGINS).all():
        raise ValueError("Unknown plant origin category")
    counts=flora.groupby("island_id").accepted_species.nunique()
    eligible=set(counts[counts>=30].index)
    if len(eligible)!=928:
        raise ValueError("Original exactly 928 Direct triple-complete flora islands drifted")
    flora=flora[flora.island_id.isin(eligible)].copy()
    return flora,{
        "n_original_eligible_islands":len(eligible),
        "n_original_complete_plant_species":int(flora.accepted_species.nunique()),
        "n_original_island_plant_species_rows":int(len(flora)),
        "n_strict_native_island_plant_species_rows":
            int(flora.origin_status.eq("native").sum()),
        "n_origin_unresolved_island_plant_species_rows":
            int(flora.origin_status.eq("unresolved").sum()),
        "strict_native_claim_applied_to_unresolved_origins":False,
        "plant_or_island_names_in_receipt":False,
    }


def macro_exact_928_locality_gate(
        georeferenced_plant_visitor_relations:set[tuple],
        bibliographic_providers_by_reference:dict,
        original_928_flora:pd.DataFrame,
        verified_geometry_path:Path,
        verified_geometry_receipt_path:Path,
        *, largest_global_archive_first_label:str|None=None,
)->dict:
    """Source-provenance-preserving polygon × exact plant name match.

    All input original accepted species, citations, raw coordinates, locality,
    original GSHHG polygon IDs and provider labels remain private.
    """
    import geopandas as gpd

    if not {"island_id","accepted_species","origin_status"}.issubset(original_928_flora):
        raise ValueError("Missing original Direct flora join keys")
    p=original_928_flora[["island_id","accepted_species","origin_status"]].copy()
    if (p.isna().any().any() or p.duplicated(["island_id","accepted_species"]).any()
            or p.island_id.nunique()!=928):
        raise ValueError("Original 928 Direct island flora not verified")
    if not verified_geometry_path.is_file() or not verified_geometry_receipt_path.is_file():
        raise ValueError("Original GSHHG geometry and receipt not available")
    receipt=json.loads(verified_geometry_receipt_path.read_text())
    if (receipt.get("GSHHG_archive_SHA256")!=PIN
        or receipt.get("n_original_focal_islands")!=21
        or receipt.get("n_focal_islands_with_exact_historical_polygon")!=21
        or receipt.get("geometry_file_SHA256")!=file_sha256(verified_geometry_path)
        or receipt.get("status") not in (
            "FULL_8265_GLOBAL_GEOMETRY_IDENTICAL",
            "FOCAL_21_POLYGONS_IDENTICAL_BUT_FULL_UNIVERSE_DRIFT")):
        raise ValueError("Geometry is not source-pin and original-island-ID verified")
    geometry=gpd.read_file(verified_geometry_path,layer="islands")
    if (geometry.crs is None or geometry.crs.to_epsg()!=4326
        or geometry.island_id.isna().any() or geometry.island_id.duplicated().any()):
        raise ValueError("Original island polygons are not unique WGS84")
    eligible=set(p.island_id)
    polygons=geometry[geometry.island_id.isin(eligible)][["island_id","geometry"]]
    retained=set(polygons.island_id)
    excluded=eligible-retained
    if polygons.empty:
        raise ValueError("No compatible original 928-island geometry")
    original_pairs=p.loc[p.island_id.isin(retained)].copy()
    raw=pd.DataFrame.from_records(
        list(georeferenced_plant_visitor_relations),
        columns=("accepted_species","guild","reference_id","provider_id",
                 "evidence_class","latitude","longitude"))
    out={
        "status":"MACRO_928_FLORA_EXACT_GLOBI_LOCALITY_EVIDENCE_PREFLIGHT",
        "n_original_Direct_triple_complete_islands":928,
        "n_islands_matching_original_GSHHG_geometric_id":int(len(retained)),
        "n_islands_unmatchable_to_original_geometry":int(len(excluded)),
        "n_original_island_plant_pairs_in_geometrically_evaluable_islands":int(len(original_pairs)),
        "n_original_GloBI_cited_georeferenced_candidate_relations":int(len(raw)),
        "n_georeferenced_relations_inside_evaluable_polygon":0,
        "n_georeferenced_relations_on_exact_same_island_flora_species":0,
        "n_distinct_islands_with_species_locality_match":0,
        "n_distinct_plant_species_with_locality_match":0,
        "n_distinct_island_species_pairs_with_record":0,
        "n_distinct_source_citations_in_matched_relation_rows":0,
        "n_distinct_original_source_archive_labels_in_matched_rows":0,
        "n_island_species_pairs_with_2plus_documented_visitor_guilds":0,
        "n_island_species_pairs_with_3plus_source_bibliography_keys":0,
        "n_island_species_pairs_with_2plus_distinct_archive_first_labels":0,
        "n_island_species_pairs_with_2plus_distinct_coordinate_sites":0,
        "n_native_island_species_pairs_with_2plus_documented_visitor_guilds":0,
        "n_unresolved_origin_island_species_pairs_with_2plus_documented_visitor_guilds":0,
        "n_georeferenced_same_native_plant_relation_units_with_pollination_label":0,
        "n_source_archive_labels_accounting_for_local_cited_units":0,
        "largest_archive_fraction_in_local_cited_units":None,
        "n_island_species_pairs_with_record_after_globally_largest_archive_omitted":
            0 if largest_global_archive_first_label is not None else None,
        "n_islands_with_same_plant_cited_relations_after_globally_largest_archive_omitted":
            0 if largest_global_archive_first_label is not None else None,
        "n_native_island_species_pairs_with_record_after_globally_largest_archive_omitted":
            0 if largest_global_archive_first_label is not None else None,
        "n_native_islands_with_cited_plant_after_globally_largest_archive_omitted":
            0 if largest_global_archive_first_label is not None else None,
        "n_island_species_pairs_with_2plus_guilds_after_globally_largest_archive_omitted":
            0 if largest_global_archive_first_label is not None else None,
        "n_native_island_species_pairs_with_2plus_guilds_after_globally_largest_archive_omitted":
            0 if largest_global_archive_first_label is not None else None,
        "n_same_original_island_plant_pairs_with_bird_and_another_cited_guild":0,
        "n_native_original_island_plant_pairs_with_bird_and_another_cited_guild":0,
        "global_largest_source_is_selected_before_locality_or_trait_results":True,
        "original_archive_first_label_is_not_proven_independent_study":True,
        "n_strict_native_flora_island_species_pairs_with_visitor_relation":0,
        "n_origin_unresolved_flora_island_species_pairs_with_visitor_relation":0,
        "n_islands_with_native_plant_visitation_reference":0,
        "n_biological_pollinator_shifts_demonstrated":0,
        "n_effective_pollen_transfer_or_seedset_records":0,
        "candidate_locality_does_not_demonstrate_effective_pollination":True,
        "unreproducible_original_polygon_IDs_never_called_absence":True,
        "GSHHG_original_file_sha256":PIN,
        "GSHHG_original_polygon_file_sha256":file_sha256(verified_geometry_path),
        "GloBI_source_has_event_time":False,
        "plant_island_geography_provider_reference_IDs_released":False,
        "per_guild":[],
    }
    if raw.empty:
        for guild in GUILDS:
            out["per_guild"].append(dict(
                documented_guild=guild,n_islands_with_same_flora_species_locality=0,
                n_locality_matched_plant_species=0,
                n_strict_native_islands_with_same_species_guild=0))
        return out
    if (raw.isna().any().any()
            or (~raw.guild.isin(GUILDS)).any()
            or (~raw.evidence_class.isin(("flower_visit","pollination"))).any()
            or not np.isfinite(raw[["latitude","longitude"]].to_numpy(float)).all()
            or raw.latitude.abs().gt(90).any()
            or raw.longitude.abs().gt(180).any()):
        raise ValueError("Original GloBI matched species site/reference data invalid")
    duplicates={ref for ref,providers in bibliographic_providers_by_reference.items()
                if len(providers)>1}
    raw=raw.loc[~raw.reference_id.isin(duplicates)]
    points=gpd.GeoDataFrame(raw,
        geometry=gpd.points_from_xy(raw.longitude,raw.latitude),crs="EPSG:4326")
    located=gpd.sjoin(points,polygons,predicate="within",how="inner")
    conflicted=located.groupby(located.index).island_id.nunique()
    located=located.loc[~located.index.isin(set(conflicted[conflicted>1].index))]
    out["n_georeferenced_relations_inside_evaluable_polygon"]=int(len(located))
    linked=located.merge(original_pairs,on=["island_id","accepted_species"],
                        how="inner",validate="many_to_one")
    out["n_georeferenced_relations_on_exact_same_island_flora_species"]=int(len(linked))
    out["n_distinct_islands_with_species_locality_match"]=int(linked.island_id.nunique())
    out["n_distinct_plant_species_with_locality_match"]=int(linked.accepted_species.nunique())
    out["n_distinct_island_species_pairs_with_record"]=int(
        linked[["island_id","accepted_species"]].drop_duplicates().shape[0])
    out["n_distinct_source_citations_in_matched_relation_rows"]=int(linked.reference_id.nunique())
    out["n_distinct_original_source_archive_labels_in_matched_rows"]=int(linked.provider_id.nunique())
    native=linked.loc[linked.origin_status.eq("native")]
    unresolved=linked.loc[linked.origin_status.eq("unresolved")]
    out["n_strict_native_flora_island_species_pairs_with_visitor_relation"]=int(
        native[["island_id","accepted_species"]].drop_duplicates().shape[0])
    out["n_origin_unresolved_flora_island_species_pairs_with_visitor_relation"]=int(
        unresolved[["island_id","accepted_species"]].drop_duplicates().shape[0])
    out["n_islands_with_native_plant_visitation_reference"]=int(native.island_id.nunique())
    if not linked.empty:
        by=linked.groupby(["island_id","accepted_species"])
        out["n_island_species_pairs_with_2plus_documented_visitor_guilds"]=int(
            by.guild.nunique().ge(2).sum())
        out["n_island_species_pairs_with_3plus_source_bibliography_keys"]=int(
            by.reference_id.nunique().ge(3).sum())
        out["n_island_species_pairs_with_2plus_distinct_archive_first_labels"]=int(
            by.provider_id.nunique().ge(2).sum())
        out["n_island_species_pairs_with_2plus_distinct_coordinate_sites"]=int(
            by.apply(lambda q:q[["latitude","longitude"]].drop_duplicates().shape[0],
                     include_groups=False).ge(2).sum())
        out["n_native_island_species_pairs_with_2plus_documented_visitor_guilds"]=int(
            native.groupby(["island_id","accepted_species"]).guild.nunique().ge(2).sum())
        out["n_unresolved_origin_island_species_pairs_with_2plus_documented_visitor_guilds"]=int(
            unresolved.groupby(["island_id","accepted_species"]).guild.nunique().ge(2).sum())
        out["n_georeferenced_same_native_plant_relation_units_with_pollination_label"]=int(
            native.evidence_class.eq("pollination").sum())
        local_providers=linked.provider_id.value_counts()
        out["n_source_archive_labels_accounting_for_local_cited_units"]=int(len(local_providers))
        out["largest_archive_fraction_in_local_cited_units"]=float(local_providers.iloc[0]/len(linked))
        bird="flower_visiting_birds"
        bridge=by.guild.agg(set)
        pairs_with_birds_and_other={k for k,v in bridge.items() if bird in v and len(v)>1}
        out["n_same_original_island_plant_pairs_with_bird_and_another_cited_guild"]=len(pairs_with_birds_and_other)
        origin=p.set_index(["island_id","accepted_species"]).origin_status
        out["n_native_original_island_plant_pairs_with_bird_and_another_cited_guild"]=int(
            sum(origin.get(k)=="native" for k in pairs_with_birds_and_other))
    if largest_global_archive_first_label is not None:
        after=linked.loc[linked.provider_id.ne(largest_global_archive_first_label)]
        after_native=after.loc[after.origin_status.eq("native")]
        out["n_island_species_pairs_with_record_after_globally_largest_archive_omitted"]=int(
            after[["island_id","accepted_species"]].drop_duplicates().shape[0])
        out["n_islands_with_same_plant_cited_relations_after_globally_largest_archive_omitted"]=int(
            after.island_id.nunique())
        out["n_native_island_species_pairs_with_record_after_globally_largest_archive_omitted"]=int(
            after_native[["island_id","accepted_species"]].drop_duplicates().shape[0])
        out["n_native_islands_with_cited_plant_after_globally_largest_archive_omitted"]=int(
            after_native.island_id.nunique())
        out["n_island_species_pairs_with_2plus_guilds_after_globally_largest_archive_omitted"]=int(
            after.groupby(["island_id","accepted_species"]).guild.nunique().ge(2).sum())
        out["n_native_island_species_pairs_with_2plus_guilds_after_globally_largest_archive_omitted"]=int(
            after_native.groupby(["island_id","accepted_species"]).guild.nunique().ge(2).sum())
    for guild in GUILDS:
        sub=linked.loc[linked.guild.eq(guild)]
        nsub=native.loc[native.guild.eq(guild)]
        after_guild=(linked.loc[
            linked.guild.eq(guild) &
            linked.provider_id.ne(largest_global_archive_first_label)]
            if largest_global_archive_first_label is not None else None)
        after_native=(after_guild.loc[after_guild.origin_status.eq("native")]
            if after_guild is not None else None)
        out["per_guild"].append({
            "documented_guild":guild,
            "n_original_island_plant_pairs_with_documented_guild":int(
                sub[["island_id","accepted_species"]].drop_duplicates().shape[0]),
            "n_native_island_plant_pairs_with_documented_guild":int(
                nsub[["island_id","accepted_species"]].drop_duplicates().shape[0]),
            "n_original_island_plant_pairs_with_guild_after_globally_largest_source_removed":
                int(after_guild[["island_id","accepted_species"]].drop_duplicates().shape[0])
                if after_guild is not None else None,
            "n_native_island_plant_pairs_with_guild_after_globally_largest_source_removed":
                int(after_native[["island_id","accepted_species"]].drop_duplicates().shape[0])
                if after_native is not None else None,
            "n_islands_with_same_flora_species_locality":
                int(sub.island_id.nunique()),
            "n_locality_matched_plant_species":
                int(sub.accepted_species.nunique()),
            "n_strict_native_islands_with_same_species_guild":
                int(nsub.island_id.nunique()),
        })
    return out
