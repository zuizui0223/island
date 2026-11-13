"""Reproduce FD real plant–pollinator interaction overlap from pinned EuPPollNet.

Original v1.3.0 623,476 interaction rows have corrected duplicates. Read the
public source privately, require original GSHHG polygons and the exact
rights-protected island Direct flora. Upload ONLY anonymous group counts.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import geopandas as gpd
import pandas as pd

from fd.pollinator_functional_diversity import direct_triple
from fd.euppollnet_island_overlap import (
    SOURCE_COLUMNS, validate_original_visits, audit_original_island_euppollnet
)
from fd.euppollnet_guild_loss_simulation import original_island_network_guild_ablation

EXPECTED_GEOMETRY_ORIGINAL_ZIP_SHA = "8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"


def sha256(p: Path) -> str:
    h=hashlib.sha256()
    with p.open("rb") as f:
        for z in iter(lambda:f.read(1024*1024),b""):
            h.update(z)
    return h.hexdigest()


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--euppollnet-interactions",type=Path,required=True)
    p.add_argument("--euppollnet-metadata",type=Path,required=True)
    p.add_argument("--original-progressive-root",type=Path,required=True)
    p.add_argument("--original-gshhg-gpkg",type=Path,required=True)
    p.add_argument("--original-gshhg-receipt",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    meta=json.loads(a.original_gshhg_receipt.read_text())
    if (meta.get("GSHHG_archive_SHA256")!=EXPECTED_GEOMETRY_ORIGINAL_ZIP_SHA
        or meta.get("original_frozen_island_count")!=8265
        or meta.get("n_original_focal_islands")!=21
        or meta.get("n_focal_islands_with_exact_historical_polygon")!=21
        or meta.get("status") not in (
          "FULL_8265_GLOBAL_GEOMETRY_IDENTICAL",
          "FOCAL_21_POLYGONS_IDENTICAL_BUT_FULL_UNIVERSE_DRIFT")
        or sha256(a.original_gshhg_gpkg)!=meta.get("geometry_file_SHA256")):
        raise ValueError("Original GSHHG source identity or exact polygon SHA failed")
    g=gpd.read_file(a.original_gshhg_gpkg,layer="islands")
    if (g.crs is None or g.crs.to_epsg()!=4326
        or g.island_id.isna().any() or g.island_id.duplicated().any()):
        raise ValueError("Original GSHHG original-island polygons invalid")

    root=a.original_progressive_root
    direct=pd.read_csv(root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        usecols=["accepted_species","trait_name","normalized_value",
                 "evidence_scope","evidence_quality"])
    complete=set(direct_triple(direct).index)
    original_flora=pd.read_csv(
        root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        usecols=["island_id","accepted_species","origin_status"])
    original_flora=original_flora.loc[original_flora.accepted_species.isin(complete)]
    eligibility=original_flora.groupby("island_id").accepted_species.nunique()
    eligible=set(eligibility[eligibility>=30].index)
    if len(eligible)!=928:
        raise ValueError("Original frozen 928 direct triple-morphology eligible islands drifted")
    floral=original_flora.loc[original_flora.island_id.isin(eligible)]
    geo=pd.read_csv(root/"fixed/isolation/results/purpose_shortest_island_data.csv",
                    usecols=["island_id","analysis_regime","spatial_block"])
    # SOURCE original data checks are strict: user-facing older 1,162,109
    # record description predates the actual deduplicated v1.3.0 source.
    source=pd.read_csv(a.euppollnet_interactions,
                       compression="gzip",dtype=str,low_memory=False,
                       encoding="utf-8",encoding_errors="replace")
    n_original_source_rows=int(len(source))
    if not (600000<=n_original_source_rows<=650000):
        raise ValueError("EuPPollNet original v1.3.0 no longer has corrected-row support")
    # A few accented original fields can include legacy invalid UTF8 bytes.
    # Count such rows, but never guess or use corrupted botanical IDs,
    # coordinates, study IDs, dates, taxa or visitor names in the join.
    text_fields=source.select_dtypes(include=["object"]).columns
    bad_encoding=pd.Series(False,index=source.index)
    for k in text_fields:
        bad_encoding |= source[k].fillna("").str.contains("\ufffd",regex=False)
    n_bad_encoding=int(bad_encoding.sum())
    if n_bad_encoding/n_original_source_rows > 0.01:
        raise ValueError("More than 1% EuPPollNet rows have undecodable UTF8: stop")
    source=source.loc[~bad_encoding].copy()
    source_meta=pd.read_csv(a.euppollnet_metadata,dtype=str,
                            encoding="utf-8",encoding_errors="strict")
    if not {"Study_id","Sampling_method"}.issubset(source_meta):
        raise ValueError("EuPPollNet v1.3.0 metadata schema drift")
    visits=validate_original_visits(source)
    # Exact site coordinates for each published study and network, no
    # centroids or radius-based assignment.
    keys=["Study_id","Network_id","Latitude","Longitude"]
    positions=visits[keys].drop_duplicates().copy()
    pts=gpd.GeoDataFrame(positions,
        geometry=gpd.points_from_xy(positions.Longitude,positions.Latitude),
        crs="EPSG:4326")
    loc=gpd.sjoin(pts,g[["island_id","geometry"]],predicate="within",how="inner")
    if not loc.empty:
        ambiguous=loc.groupby(keys).island_id.nunique()
        ambiguous=set(ambiguous[ambiguous>1].index)
        if ambiguous:
            index=loc.set_index(keys).index
            loc=loc.loc[~index.isin(ambiguous)]
    site_to_island=pd.DataFrame(
        loc[[*keys,"island_id"]],columns=[*keys,"island_id"]).drop_duplicates()
    if site_to_island.duplicated(keys).any():
        raise ValueError("Conflicting site across exact original island polygons")
    result=audit_original_island_euppollnet(
       source,site_to_island,floral,geo)
    # Model WHAT COULD HAPPEN to observed flower-visit links if a sampled
    # guild disappeared, not actual decline, true substitution, or fitness.
    simulated=original_island_network_guild_ablation(
        source,site_to_island,floral,direct_triple(direct),geo,source_meta)
    result["hypothetical_within_network_visitor_guild_removal"]=simulated
    result.update(
        source_public_license="CC_BY_4_0",
        source_zenodo_record="10.5281/zenodo.15183272",
        source_public_repository="joseblanuza/EuPPollNet",
        source_git_ref="v1.3.0",
        original_source_interaction_csv_SHA256=sha256(a.euppollnet_interactions),
        original_source_metadata_csv_SHA256=sha256(a.euppollnet_metadata),
        original_input_interaction_rows=n_original_source_rows,
        original_source_rows_with_invalid_UTF8_excluded=n_bad_encoding,
        source_rows_qualified_after_encoding_quality=int(len(source)),
        original_source_studies=int(source_meta.Study_id.nunique()),
        original_source_longitude_latitude_sites=int(len(positions)),
        n_source_visiting_sites_on_original_GSHHG_island_polygons=int(len(site_to_island)),
        n_original_GSHHG_island_ids_with_visiting_sites=int(site_to_island.island_id.nunique()),
        n_original_928_flora_ids_without_verified_geometry=int(
            len(eligible-set(g.island_id))),
        original_GSHHG_source_zip_SHA256=EXPECTED_GEOMETRY_ORIGINAL_ZIP_SHA,
        exact_island_georeferencing_subject_to_upstream_site_precision=True,
        absolute_pollinator_loss_not_demonstrated=True,
        dated_flower_visitor_events_not_pollination_effectiveness=True,
    )
    raw=json.dumps(result,sort_keys=True)
    forbidden=['"island_id"','"Study_id"','"Network_id"',
       '"Plant_accepted_name"','"Pollinator_accepted_name"',
       '"Latitude"','"Longitude"','"accepted_species"']
    if any(k in raw for k in forbidden):
        raise ValueError("Raw original identifier fields leaked into public report")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
       "status":result["status"],
       "source_rows":result["original_input_interaction_rows"],
       "source_studies":result["original_source_studies"],
       "original_gshhg_islands_with_sites":result["n_original_GSHHG_island_ids_with_visiting_sites"],
       "fd_islands_with_sites":result["n_original_928_flora_islands_with_exact_network_site"],
       "fd_direct_plant_visits":result.get("n_direct_species_joined_interaction_rows",0),
       "n_source_repeated_site_years":result.get("n_islands_with_two_years_site",0),
       "source_disjoint_regions_eligible":result.get("n_regions_meeting_10_islands_5blocks_3studies",0),
       "simulated_network_scenarios":simulated.get("n_qualified_network_simulation_scenarios",0)
    },sort_keys=True))
if __name__=="__main__":
    main()
