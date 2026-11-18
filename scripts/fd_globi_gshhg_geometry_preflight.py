"""Regenerate ONLY the pinned original 8,265 GSHHG island polygons.

Historical GSHHG artifact 8066083419 expired in October 2026.
Never substitute Natural Earth, approximate centroids or republish an
unverified geometry as an exact geographic visitor observation.
This preflight runs independently of GloBI flower interaction labels.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import geopandas as gpd
import pandas as pd
import yaml

from island_v2.gshhg_source import build_from_gshhg, checksum

EXPECTED_GSHHG="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
SOURCE_URL="https://github.com/GenericMappingTools/gshhg-gmt/releases/download/2.3.7/gshhg-shp-2.3.7.zip"
EXPECTED_N_ISLANDS=8265

def original_focal_island_ids(ledger_path,flora_path,channel_path):
    """Build 21 candidate IDs BEFORE any geometry or GloBI point inspection."""
    d=pd.read_csv(ledger_path,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"])
    t=d[d.trait_name.isin(("flower_primary_color","floral_form","floral_symmetry"))
        & d.evidence_scope.eq("species_direct")
        & d.evidence_quality.isin(("high","medium"))
        & d.normalized_value.notna()
        & ~d.normalized_value.astype(str).str.contains("|",regex=False)]
    if t.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Original Direct floral trait duplicated")
    wide=t.pivot(index="accepted_species",columns="trait_name",
                 values="normalized_value")
    known=set(wide.dropna(subset=[
        "flower_primary_color","floral_form","floral_symmetry"]).index)
    flora=pd.read_csv(flora_path,usecols=["island_id","accepted_species"])
    sample=flora[flora.accepted_species.isin(known)]
    counts=sample.groupby("island_id").accepted_species.nunique()
    eligible=set(counts[counts>=30].index)
    if len(eligible)!=928:
        raise ValueError("Original 928 Direct flora islands drifted")
    q=pd.read_csv(channel_path,usecols=[
        "island_id","channel_id","channel_state","source_state"])
    q=q[q.island_id.isin(eligible)]
    if q.duplicated(["island_id","channel_id"]).any():
        raise ValueError("Original qualified guild state duplicated")
    if not q.loc[q.channel_state.eq("disrupted"),"source_state"].eq("available").all():
        raise ValueError("Adequate nondetection without source opportunity")
    states=q.pivot(index="island_id",columns="channel_id",values="channel_state")
    if not {"bombus","flower_visiting_birds"}.issubset(states.columns):
        raise ValueError("Missing original Bombus/bird evidence")
    focal=set(states.index[states.bombus.eq("disrupted") &
                           states.flower_visiting_birds.eq("retained")])
    if len(focal)!=21:
        raise ValueError("Original candidate island count not exactly 21")
    return focal


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-source-assignments",required=True,type=Path)
    p.add_argument("--original-direct-ledger",required=True,type=Path)
    p.add_argument("--original-island-flora",required=True,type=Path)
    p.add_argument("--original-qualified-channel-states",required=True,type=Path)
    p.add_argument("--output-dir",required=True,type=Path)
    p.add_argument("--config",default="config/island_source_gshhg.yml",type=Path)
    a=p.parse_args()
    config=yaml.safe_load(a.config.read_text())
    if (config["source_version"]!="2.3.7" or config["resolution"]!="h"
            or float(config["min_area_km2"])!=5.0
            or float(config["mainland_area_threshold_km2"])!=7000000.0):
        raise ValueError("Frozen GSHHG 2.3.7 h area contract mismatch")
    a.output_dir.mkdir(parents=True,exist_ok=True)
    units, original_source=build_from_gshhg(
        SOURCE_URL,config,a.output_dir)
    digest=checksum(original_source)
    if digest!=EXPECTED_GSHHG:
        raise ValueError("Regenerated GSHHG archive SHA256 differs from source pin")
    if not units.island_id.is_unique:
        raise ValueError("Regenerated GSHHG island polygon IDs not unique")
    # Compare source geometry to the EXACT original island ID set before
    # deciding whether a rebuilt geometry can be used. A different count
    # is not an invitation to relax the historical identity contract.
    if len(units)!=EXPECTED_N_ISLANDS:
        print(json.dumps({
            "geometry_rebuild_diagnostic":"COUNT_NOT_EQUAL_ORIGINAL_GIFT_SOURCE_UNIVERSE",
            "n_rebuilt_island_polygons":int(len(units)),
            "n_expected_islands":EXPECTED_N_ISLANDS,
            "GSHHG_archive_SHA256":digest,
            "source_polygon_ids_released":False,
        },sort_keys=True),flush=True)
    if units.crs is None or units.crs.to_epsg()!=4326:
        raise ValueError("Original geometry CRS must be EPSG:4326")
    source=pd.read_csv(a.original_source_assignments,
                       compression="infer",
                       usecols=["island_id","source_mode","source_rank"])
    candidates=source[source.source_mode.eq("geo_k5")]
    if candidates.empty or candidates.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original GIFT geo_k5 source region assignment invalid")
    original_ids=set(candidates.island_id)
    observed_ids=set(units.island_id)
    focal=original_focal_island_ids(
        a.original_direct_ledger,a.original_island_flora,
        a.original_qualified_channel_states)
    missing=original_ids-observed_ids
    extra=observed_ids-original_ids
    print(json.dumps({
        "diagnostic":"COMPLETE_8265_GLOBAL_ID_GATE" if not(missing or extra)
                     else "HISTORIC_GLOBAL_ID_DRIFT_SUBSET_GATE",
        "n_original_GIFT_islands":len(original_ids),
        "n_rebuilt_islands":len(observed_ids),
        "n_original_only_ids":len(missing),
        "n_rebuild_only_ids":len(extra),
        "n_original_focal_candidates":len(focal),
        "n_original_focal_candidates_missing_polygon":len(focal-observed_ids),
        "focal_island_IDs_exposed":False,
    },sort_keys=True),flush=True)
    if not focal.issubset(observed_ids):
        raise ValueError("Some original focal candidate polygons missing: reject all alleged locality inference")
    if candidates.groupby("island_id").source_rank.nunique().ne(5).any():
        raise ValueError("Frozen source regions must have five distinct ranks")
    # Only polygons bearing EXACT original WKB-derived island IDs survive.
    # Non-focal global hash drift does not authorize a full-universe replay.
    use=units[units.island_id.isin(original_ids)].copy()
    gpkg=a.output_dir/"original_gshhg_v237_h_candidate21_exact_shared_ids.gpkg"
    use.to_file(gpkg,layer="islands",driver="GPKG")
    saved=gpd.read_file(gpkg,layer="islands")
    if set(saved.island_id)!=set(use.island_id) or len(saved)!=len(use):
        raise ValueError("GeoPackage round-trip changed original island hashes")
    if not focal.issubset(set(saved.island_id)):
        raise ValueError("Original focal island lost on GeoPackage round-trip")
    status=("FULL_8265_GLOBAL_GEOMETRY_IDENTICAL" if not(missing or extra)
            else "FOCAL_21_POLYGONS_IDENTICAL_BUT_FULL_UNIVERSE_DRIFT")

    receipt={
        "status":status,
        "GSHHG_version":"2.3.7",
        "GSHHG_archive_SHA256":digest,
        "GSHHG_resolution":"h",
        "source_backend":"GSHHG_NO_FALLBACK",
        "original_frozen_island_universe_artifact":"8066083419_EXPIRED",
        "original_frozen_island_count":EXPECTED_N_ISLANDS,
        "n_reconstructed_polygon_units":int(len(saved)),
        "n_original_GIFT_geo_k5_island_ids":int(len(original_ids)),
        "n_IDs_missing_from_geometry":int(len(missing)),
        "n_extra_ids_in_unfiltered_rebuild":int(len(extra)),
        "n_original_focal_islands":int(len(focal)),
        "n_focal_islands_with_exact_historical_polygon":int(len(focal)),
        "global_full_universe_identity_verified":not bool(missing or extra),
        "spatial_inference_allowed_only_for_verified_21_candidate_polygon_IDs":True,
        "n_original_GIFT_source_rows":int(len(candidates)),
        "n_source_species_pollinator_records_tested":0,
        "GloBI_locality_or_pollination_shift_confirmed":False,
        "geometry_file_SHA256":checksum(gpkg),
        "island_or_botanical_species_names_in_receipt":False,
    }
    (a.output_dir/"original_gshhg_geometry_receipt.json").write_text(
        json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))


if __name__=="__main__":
    main()
