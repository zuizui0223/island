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


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-source-assignments",required=True,type=Path)
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
    if original_ids!=observed_ids:
        print(json.dumps({
            "geometry_rebuild_diagnostic":"ID_SET_NOT_EQUAL_ORIGINAL_GIFT_SOURCE_UNIVERSE",
            "n_rebuilt_polygons":int(len(observed_ids)),
            "n_original_GIFT_island_IDs":int(len(original_ids)),
            "n_source_IDs_not_in_rebuild":int(len(original_ids-observed_ids)),
            "n_rebuilt_IDs_not_in_original_source":int(len(observed_ids-original_ids)),
            "GSHHG_archive_SHA256":digest,
            "source_polygon_ids_released":False,
        },sort_keys=True),flush=True)
        raise ValueError(
            "Exact historical island IDs differ; STOP all purported original-island GloBI geolocation inference")
    if candidates.groupby("island_id").source_rank.nunique().ne(5).any():
        raise ValueError("Frozen source regions must have five ranks per island")
    gpkg=a.output_dir/"original_gshhg_v237_h_8265_verified.gpkg"
    units.to_file(gpkg,layer="islands",driver="GPKG")
    saved=gpd.read_file(gpkg,layer="islands")
    if set(saved.island_id)!=original_ids or len(saved)!=EXPECTED_N_ISLANDS:
        raise ValueError("GeoPackage round-trip altered island identities")
    receipt={
        "status":"ORIGINAL_GSHHG_GEOMETRY_RECONSTRUCTED_SAME_8265_IDS",
        "GSHHG_version":"2.3.7",
        "GSHHG_archive_SHA256":digest,
        "GSHHG_resolution":"h",
        "source_backend":"GSHHG_NO_FALLBACK",
        "original_frozen_island_universe_artifact":"8066083419_EXPIRED",
        "original_frozen_island_count":EXPECTED_N_ISLANDS,
        "n_reconstructed_polygon_units":int(len(saved)),
        "n_original_GIFT_geo_k5_island_ids":int(len(original_ids)),
        "n_IDs_missing_from_geometry":0,
        "n_extra_ids_in_geometry":0,
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
