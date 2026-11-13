"""Outcome-blind GloPL experimental pollen-limitation × fd floral FD gate.

NEVER read PL_Effect_Size or treatment outcomes in this preflight. Source
GloPL bytes, exact original GSHHG island polygons and original Direct flora
must match their frozen provenance. Only anonymized region counts are public.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import geopandas as gpd
import pandas as pd

from island_v2.chapter1_h5_glopl_island_overlap_preflight import (
    load_config, read_glopl_metadata, match_glopl_to_islands
)
from fd.pollinator_functional_diversity import direct_triple

GL_OPL_SHA="264ca8c2237f126b3c209fcba039e19b291324ddf6e7a96992dc8235a55984de"
GSHHG_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"


def digest(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--glopl",type=Path,required=True)
    p.add_argument("--glopl-config",type=Path,required=True)
    p.add_argument("--progressive-root",type=Path,required=True)
    p.add_argument("--original-gshhg-gpkg",type=Path,required=True)
    p.add_argument("--original-gshhg-receipt",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()

    if digest(a.glopl)!=GL_OPL_SHA:
        raise ValueError("GloPL original source SHA mismatch")
    receipt=json.loads(a.original_gshhg_receipt.read_text())
    if (receipt.get("GSHHG_archive_SHA256")!=GSHHG_SHA
        or receipt.get("original_frozen_island_count")!=8265
        or receipt.get("n_original_focal_islands")!=21
        or receipt.get("n_focal_islands_with_exact_historical_polygon")!=21
        or receipt.get("status") not in (
          "FULL_8265_GLOBAL_GEOMETRY_IDENTICAL",
          "FOCAL_21_POLYGONS_IDENTICAL_BUT_FULL_UNIVERSE_DRIFT")
        or digest(a.original_gshhg_gpkg)!=receipt.get("geometry_file_SHA256")):
        raise ValueError("GSHHG historical source geometry provenance mismatch")

    cfg=load_config(a.glopl_config)
    # Only geography/study/taxon columns are read: no pollen limitation.
    glopl=read_glopl_metadata(a.glopl,cfg)
    geom=gpd.read_file(a.original_gshhg_gpkg,layer="islands")
    matched,audit=match_glopl_to_islands(glopl,geom,cfg)
    root=a.progressive_root
    ledger=pd.read_csv(root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
       usecols=["accepted_species","trait_name","normalized_value",
                "evidence_scope","evidence_quality"])
    complete=direct_triple(ledger)
    names=set(complete.index)
    flora=pd.read_csv(root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
                      usecols=["island_id","accepted_species","origin_status"])
    flora=flora.loc[flora.accepted_species.isin(names)].copy()
    rich=flora.groupby("island_id").accepted_species.nunique()
    elig=set(rich[rich>=30].index)
    if len(elig)!=928:
        raise ValueError("Original strict Direct 928 island-flora universe drift")
    flora=flora.loc[flora.island_id.isin(elig)]
    geo=pd.read_csv(root/"fixed/isolation/results/purpose_shortest_island_data.csv",
       usecols=["island_id","analysis_regime","spatial_block"])
    if geo.island_id.duplicated().any():
        raise ValueError("Original geographic block id duplicates")
    original=matched.loc[matched.island_id.isin(elig)].copy()
    if original.empty:
        classified=pd.DataFrame(columns=[
          "island_id","study_key","site_key","Species_accepted_names",
          "analysis_regime","spatial_block","origin_status"])
    else:
        original=original.merge(geo,
            on="island_id",how="left",validate="many_to_one")
        classified=original.merge(flora.rename(
           columns={"accepted_species":"Species_accepted_names"}),
           on=["island_id","Species_accepted_names"],
           how="left",validate="many_to_one")
    region=[]
    if not classified.empty:
        for name,part in classified.groupby("analysis_regime"):
            strict=part.loc[part.origin_status.eq("native")]
            exact=part.loc[part.origin_status.notna()]
            islands=part.island_id.nunique()
            studies=part.study_key.nunique()
            blocks=part.spatial_block.nunique()
            region.append(dict(regime=str(name),
               n_flora_islands_with_experimental_geolocations=int(islands),
               n_spatial_blocks=int(blocks),
               n_original_GloPL_publication_keys=int(studies),
               n_GloPL_records_in_original_flora_islands=int(len(part)),
               n_direct_species_by_island_exact_name_records=int(len(exact)),
               n_direct_species_by_island_strict_native_records=int(len(strict)),
               n_islands_with_exact_direct_plant_experimental_record=int(exact.island_id.nunique()),
               n_islands_with_strict_native_direct_plant_experimental_record=int(strict.island_id.nunique()),
               enough_support_to_consider_plant_side_association=bool(
                    islands>=10 and blocks>=8 and studies>=5),
               effects_or_service_mechanism_tested=False))
    result=dict(
      status="GLOPL_EXPERIMENTAL_SERVICE_FD_COLOCATION_ONLY",
      original_GloPL_source_sha256=GL_OPL_SHA,
      GloPL_source_original_rows=int(len(glopl)),
      n_original_GloPL_rows_matching_any_verified_GSHHG_island=int(audit["n_rows_within_any_frozen_island"]),
      n_original_GloPL_sites_matching_any_verified_GSHHG_island=int(audit["n_unique_sites_within_any_frozen_island"]),
      n_original_GloPL_islands_matching_any_verified_GSHHG=int(matched.island_id.nunique()),
      n_original_direct_eligible_islands=928,
      n_original_fd_islands_with_GloPL=int(classified.island_id.nunique()),
      n_original_fd_plant_exact_source_GloPL_records=int(classified.origin_status.notna().sum()) if not classified.empty else 0,
      regions=region,
      n_regions_meeting_pinned_GloPL_to_flora_support=sum(
          r["enough_support_to_consider_plant_side_association"] for r in region),
      GloPL_pollen_effect_column_never_read_in_preflight=True,
      pollen_service_effect_on_fd_or_evolution_not_identified=True,
      all_original_taxon_island_study_geolocations_kept_private=True,
    )
    if any('"'+k+'"' in json.dumps(result) for k in (
      "island_id","Species_accepted_names","study_key","site_key",
      "Latitude","Longitude","PL_Effect_Size","accepted_species")):
        raise ValueError("Original source private identity/outcome leaking")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
      "n_GloPL_original_rows":result["GloPL_source_original_rows"],
      "n_GloPL_islands_original_GSHHG":result["n_original_GloPL_islands_matching_any_verified_GSHHG"],
      "n_GloPL_islands_exact_fd_flora":result["n_original_fd_islands_with_GloPL"],
      "n_GloPL_records_exact_same_direct_plant":result["n_original_fd_plant_exact_source_GloPL_records"],
      "n_regions_identifiable":result["n_regions_meeting_pinned_GloPL_to_flora_support"]
    },sort_keys=True))
if __name__=="__main__":
    main()
