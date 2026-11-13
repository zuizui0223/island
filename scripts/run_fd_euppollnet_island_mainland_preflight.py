"""Original EuPPollNet SAME-PLANT island vs verified CONTINENT feasibility.

Unlike absence from the 8,265 island polygons, continental membership is
validated POSITIVELY by the source-verified GSHHG 2.3.7 high-res L1 mainland
polygons above the frozen 7,000,000 km2 landmass threshold. Small or
unverified islands stay 'unresolved', never synthetic mainland controls.
Source metadata and all row-level study, plant, site and island IDs stay
inside source-rights-safe CI. No flower visitor change is inferred.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from island_v2.gshhg_source import (
    download, extract_l1, dissolve_ids, checksum
)
from fd.pollinator_functional_diversity import direct_triple
from fd.euppollnet_island_overlap import validate_original_visits
from fd.euppollnet_island_mainland import (
    original_same_plant_island_mainland_opportunity
)

SOURCE_URL="https://github.com/GenericMappingTools/gshhg-gmt/releases/download/2.3.7/gshhg-shp-2.3.7.zip"
GSHHG_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
EU_INTERACTION_SHA="cde17692dc66946002e19026433d9db429357d28873a369d884941f71131db61"
EU_METADATA_SHA="adb396f31e7e7e3257425872b92d56c1051ca05849c6228a39675a71364d4833"
SITE_KEYS=["Study_id","Network_id","Latitude","Longitude"]


def original_mainland_polygons(workdir:Path)->gpd.GeoDataFrame:
    archive=workdir/"source"/"gshhg-shp-2.3.7.zip"
    download(SOURCE_URL,archive)
    if checksum(archive)!=GSHHG_SHA:
        raise ValueError("Original GSHHG 2.3.7 ZIP checksum mismatch")
    shape=extract_l1(archive,workdir/"original-L1-h","h")
    land=gpd.read_file(shape)
    if land.crs is None:
        land=land.set_crs(4326)
    pieces=dissolve_ids(land).to_crs(4326)
    area_km2=pieces.to_crs(6933).area/1e6
    continental=pieces.loc[area_km2.gt(7_000_000.0),["geometry"]].copy()
    if not (3<=len(continental)<=9):
        raise ValueError("Original source continental landmass universe unexpected")
    return continental


def _sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def original_site_classes(visits:pd.DataFrame,geometries:gpd.GeoDataFrame,
                          original_geography:pd.DataFrame, continents:gpd.GeoDataFrame
                          )->pd.DataFrame:
    sites=visits[SITE_KEYS].drop_duplicates().copy()
    if sites.duplicated(SITE_KEYS).any():
        raise ValueError("Same original study-network-site key duplicated")
    points=gpd.GeoDataFrame(
       sites,geometry=gpd.points_from_xy(sites.Longitude,sites.Latitude),
       crs="EPSG:4326")
    islands=geometries[["island_id","geometry"]].copy()
    island_join=gpd.sjoin(points,islands,predicate="within",how="inner")
    if not island_join.empty:
        ambiguous=island_join.groupby(SITE_KEYS).island_id.nunique()
        bad=ambiguous.loc[ambiguous.gt(1)].index
        if len(bad):
            raise ValueError("Source location assigned to multiple original island polygons")
    is_map=island_join[SITE_KEYS+["island_id"]].drop_duplicates()
    eligible_geo=original_geography[["island_id","analysis_regime","spatial_block"]].copy()
    if eligible_geo.island_id.duplicated().any() or eligible_geo.isna().any().any():
        raise ValueError("Original protected FD island geography not unique")
    is_map=is_map.merge(eligible_geo,on="island_id",how="inner",validate="many_to_one")
    major=gpd.sjoin(points,continents,predicate="within",how="inner")
    mainland_keys=major[SITE_KEYS].drop_duplicates()
    types=sites.merge(is_map,on=SITE_KEYS,how="left",validate="one_to_one")
    types["habitat_class"]="unresolved"
    is_island=types.island_id.notna()
    continent_hit=(pd.MultiIndex.from_frame(types[SITE_KEYS]).isin(
       pd.MultiIndex.from_frame(mainland_keys[SITE_KEYS])))
    if (is_island & continent_hit).any():
        raise ValueError("Overlap between original island and continental polygons")
    types.loc[is_island,"habitat_class"]="fd_island"
    types.loc[~is_island & continent_hit,"habitat_class"]="mainland"
    # No continental polygon == unresolved, NEVER an island-inferred mainland.
    return types[SITE_KEYS+["habitat_class","island_id",
                            "analysis_regime","spatial_block"]]


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--euppollnet-interactions",type=Path,required=True)
    p.add_argument("--euppollnet-metadata",type=Path,required=True)
    p.add_argument("--original-progressive-root",type=Path,required=True)
    p.add_argument("--verified-gshhg-gpkg",type=Path,required=True)
    p.add_argument("--verified-gshhg-receipt",type=Path,required=True)
    p.add_argument("--source-cache",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    if _sha256(a.euppollnet_interactions)!=EU_INTERACTION_SHA:
        raise ValueError("Source-corrected EuPPollNet 623k interactions changed")
    if _sha256(a.euppollnet_metadata)!=EU_METADATA_SHA:
        raise ValueError("Original source protocol metadata changed")
    receipt=json.loads(a.verified_gshhg_receipt.read_text())
    if (receipt.get("GSHHG_archive_SHA256")!=GSHHG_SHA
        or receipt.get("original_frozen_island_count")!=8265
        or receipt.get("n_original_focal_islands")!=21
        or receipt.get("n_focal_islands_with_exact_historical_polygon")!=21
        or receipt.get("status") not in (
           "FULL_8265_GLOBAL_GEOMETRY_IDENTICAL",
           "FOCAL_21_POLYGONS_IDENTICAL_BUT_FULL_UNIVERSE_DRIFT")
        or _sha256(a.verified_gshhg_gpkg)!=receipt.get("geometry_file_SHA256")):
        raise ValueError("Original island GSHHG history/polygon identity invalid")
    islands=gpd.read_file(a.verified_gshhg_gpkg,layer="islands")
    if islands.crs is None or islands.crs.to_epsg()!=4326 or islands.island_id.duplicated().any():
        raise ValueError("Original island WKB geography invalid")
    root=a.original_progressive_root
    direct=pd.read_csv(root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        usecols=["accepted_species","trait_name","normalized_value",
                 "evidence_scope","evidence_quality"])
    known=set(direct_triple(direct).index)
    original_flora=pd.read_csv(root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        usecols=["island_id","accepted_species","origin_status"])
    flora=original_flora.loc[original_flora.accepted_species.isin(known)].copy()
    if flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original licensed flora contains duplicate botanical species")
    n=flora.groupby("island_id").accepted_species.nunique()
    eligible=set(n.loc[n>=30].index)
    if len(eligible)!=928:
        raise ValueError("Original fixed fd Direct three-trait 928 island cohort changed")
    flora=flora.loc[flora.island_id.isin(eligible)]
    geo=pd.read_csv(root/"fixed/isolation/results/purpose_shortest_island_data.csv",
        usecols=["island_id","analysis_regime","spatial_block"])
    geo=geo.loc[geo.island_id.isin(eligible)].copy()
    continents=original_mainland_polygons(a.source_cache)
    original=pd.read_csv(a.euppollnet_interactions,compression="gzip",
             dtype=str,encoding="utf-8",encoding_errors="replace",low_memory=False)
    if len(original)!=623476:
        raise ValueError("Original v1.3.0 de-duplicated flower-visit count changed")
    text_fields=original.select_dtypes(include=["object"]).columns
    encoding_bad=pd.Series(False,index=original.index)
    for f in text_fields:
        encoding_bad |= original[f].fillna("").str.contains("\ufffd",regex=False)
    if encoding_bad.mean()>.01:
        raise ValueError("Original public source decoding error >1%")
    original=original.loc[~encoding_bad].copy()
    visits=validate_original_visits(original)
    metadata=pd.read_csv(a.euppollnet_metadata,dtype=str,encoding="utf-8")
    classes=original_site_classes(visits,islands,geo,continents)
    receipt=original_same_plant_island_mainland_opportunity(
        original,classes,flora,metadata)
    receipt.update(
        original_full_interaction_rows=623476,
        original_encoding_excluded_rows=int(encoding_bad.sum()),
        original_source_visit_studies=int(metadata.Study_id.nunique()),
        original_GSHHG_2_3_7_h_SHA256=GSHHG_SHA,
        original_928_fd_islands=928,
        original_flora_islands_exact_GSHHG_geometry=int(len(eligible & set(islands.island_id))),
        n_source_network_sites_with_positive_continent_membership=int(
          classes.habitat_class.eq("mainland").sum()),
        n_source_network_sites_with_exact_fd_island_membership=int(
          classes.habitat_class.eq("fd_island").sum()),
        n_source_network_sites_geographically_unresolved=int(
          classes.habitat_class.eq("unresolved").sum()),
        no_fake_mainland_from_island_record_absence=True,
        original_plant_study_island_identity_never_released=True)
    serialized=json.dumps(receipt,sort_keys=True)
    if any('"'+key+'"' in serialized for key in (
       "Study_id","Network_id","Latitude","Longitude","Plant_accepted_name",
       "Pollinator_accepted_name","island_id","accepted_species","study_id","plant")):
        raise ValueError("Original source/site identity present in aggregate")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print(json.dumps({
        "status":receipt["status"],
        "n_fd_island_visit_sites":receipt["n_source_network_sites_with_exact_fd_island_membership"],
        "n_verified_continent_visit_sites":receipt["n_source_network_sites_with_positive_continent_membership"],
        "n_exact_plant_study_year_pairs":receipt["n_matching_same_study_year_species_island_mainland"] if "n_matching_same_study_year_species_island_mainland" in receipt
            else receipt["n_same_study_year_plant_candidate_groups_any_distance"],
        "n_nearby_native_pairs":receipt.get("n_eligible_strict_native_pair_groups_500km",0),
    },sort_keys=True))
if __name__=="__main__":
    main()
