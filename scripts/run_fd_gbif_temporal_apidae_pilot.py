"""Live GBIF Search *documentation* temporal feasibility pilot for original fd floras.

Source-pinned island species/trait eligibility and GSHHG polygon identity.
Strictly bounded API calls. Reads original raw GBIF species/dataset/year records
only in process memory and publishes aggregate feasibility by region. This
cannot establish actual pollinator decline, visit rate or floristic evolution.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import time

import geopandas as gpd
import httpx
import numpy as np
import pandas as pd
from shapely.geometry import Point

from fd.pollinator_functional_diversity import direct_triple
from fd.gbif_temporal_opportunity import (
    PERIODS, select_outcome_blind_islands,
    aggregate_repeated_dataset_year_evidence,
)
from island_v2.chapter1_nee_island_search_observation import exact_query_wkts
from island_v2.gbif_blocks import resolve_query_taxon

EXPECTED_GSHHG = "8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
MAX_RAW = 900
API = "https://api.gbif.org/v1/occurrence/search"
PAGE = 300


def _text(x):
    return "" if x is None else str(x).strip()


def _retry_search(client,params,retries=3):
    error=None
    for attempt in range(retries):
        try:
            response=client.get(API,params=params)
            if response.status_code in (429,500,502,503,504):
                error=RuntimeError("GBIF rate-limited or temporary upstream outage")
                time.sleep(1.5*(attempt+1))
                continue
            response.raise_for_status()
            data=response.json()
            if not isinstance(data,dict) or not isinstance(data.get("results"),list):
                raise ValueError("Invalid GBIF Search JSON payload")
            return data
        except (httpx.RequestError,httpx.HTTPStatusError,ValueError) as exc:
            error=exc
            time.sleep(1.5*(attempt+1))
    raise RuntimeError("GBIF Search API request failed after retries") from error


def _source_record(rec,geom,window,island,seen):
    """Quality filter; no background row becomes a survey or pseudo-absence."""
    if _text(rec.get("occurrenceStatus")).upper()=="ABSENT":
        return None
    if _text(rec.get("basisOfRecord")).upper()=="LIVING_SPECIMEN":
        return None
    establishment=_text(rec.get("establishmentMeans")).upper()
    if any(w in establishment for w in ("CAPTIVE","CULTIVATED","INTRODUCED")):
        return None
    gbifid=_text(rec.get("key") or rec.get("gbifID"))
    dataset=_text(rec.get("datasetKey"))
    genus=_text(rec.get("genus"))
    if not gbifid or not dataset or not genus or gbifid in seen:
        return None
    try:
        year=int(rec.get("year"))
        lat=float(rec.get("decimalLatitude"))
        lon=float(rec.get("decimalLongitude"))
    except (ValueError,TypeError):
        return None
    if not (PERIODS[window][0]<=year<=PERIODS[window][1]):
        return None
    if not(np.isfinite(lat) and np.isfinite(lon) and -90<=lat<=90 and -180<=lon<=180):
        return None
    if not geom.contains(Point(lon,lat)):
        return None
    # The source taxonKey is the whole Apidae background, not a taxon
    # composition learned from trait data or a claim of visitor effectiveness.
    seen.add(gbifid)
    return dict(island_id=island,period=window,dataset_key=dataset,
                year=year,gbif_id=gbifid,genus=genus)


def scan_window(client,*,geometry,parts,taxon_key,island,window):
    """Exhaustive only if ALL exact polygon components reach endOfRecords."""
    rows=[]
    examined=0
    complete=True
    truncated=False
    error=""
    seen=set()
    for polygon_wkt in parts:
        offset=0
        while True:
            remaining=MAX_RAW-examined
            if remaining<=0:
                complete=False
                truncated=True
                break
            limit=min(PAGE,remaining)
            params={
                "taxonKey":str(taxon_key),
                "hasCoordinate":"true",
                "hasGeospatialIssue":"false",
                "geometry":polygon_wkt,
                "year":f"{PERIODS[window][0]},{PERIODS[window][1]}",
                "limit":limit,"offset":offset,
            }
            try:
                page=_retry_search(client,params=params)
                recs=page["results"]
                if len(recs)>limit:
                    raise ValueError("GBIF returned records beyond requested page limit")
                examined+=len(recs)
                for rec in recs:
                    val=_source_record(rec,geometry,window,island,seen)
                    if val is not None:
                        rows.append(val)
                ended=bool(page.get("endOfRecords")) or len(recs)<limit
                if ended:
                    break
                if not recs:
                    raise ValueError("GBIF Search empty page without endOfRecords")
                offset+=len(recs)
            except Exception as exc:
                # Never leak raw API error strings with identifiers into output.
                error=type(exc).__name__
                complete=False
                break
        if not complete:
            break
    return rows,dict(island_id=island,period=window,
                     complete=bool(complete),truncated=bool(truncated),
                     api_error=error,n_raw_examined=examined,
                     n_usable_records=len(rows))


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--original-progressive-root",type=Path,required=True)
    a.add_argument("--verified-gpkg",type=Path,required=True)
    a.add_argument("--geometry-receipt",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    z=a.parse_args()
    proof=json.loads(z.geometry_receipt.read_text())
    if (proof.get("GSHHG_archive_SHA256")!=EXPECTED_GSHHG
        or proof.get("original_frozen_island_count")!=8265
        or proof.get("n_original_focal_islands")!=21
        or proof.get("n_focal_islands_with_exact_historical_polygon")!=21
        or proof.get("status") not in ("FULL_8265_GLOBAL_GEOMETRY_IDENTICAL",
                                      "FOCAL_21_POLYGONS_IDENTICAL_BUT_FULL_UNIVERSE_DRIFT")):
        raise ValueError("Original source GSHHG polygon source-hash contract failed")
    from island_v2.fd_globi_species_shift_spatial_v09 import file_sha256
    if file_sha256(z.verified_gpkg)!=proof.get("geometry_file_SHA256"):
        raise ValueError("Original exact GSHHG GeoPackage checksum mismatch")
    geometry=gpd.read_file(z.verified_gpkg,layer="islands")
    if (geometry.crs is None or geometry.crs.to_epsg()!=4326
        or geometry.island_id.isna().any() or geometry.island_id.duplicated().any()):
        raise ValueError("Original exact GSHHG polygon geometry not WGS84 unique")
    root=z.original_progressive_root
    direct=pd.read_csv(root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        usecols=["accepted_species","trait_name","normalized_value",
                 "evidence_scope","evidence_quality"])
    names=set(direct_triple(direct).index)
    flora=pd.read_csv(root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        usecols=["island_id","accepted_species"])
    flora=flora.loc[flora.accepted_species.isin(names)]
    n=flora.groupby("island_id").accepted_species.nunique()
    elig=set(n[n>=30].index)
    if len(elig)!=928:
        raise ValueError("Locked fd 928 Direct triple-evaluable island sample drifted")
    geo=pd.read_csv(root/"fixed/isolation/results/purpose_shortest_island_data.csv",
        usecols=["island_id","analysis_regime","spatial_block"])
    selected=select_outcome_blind_islands(geo,elig,set(geometry.island_id))
    private=geo.loc[geo.island_id.isin(selected)].copy()
    polygons=geometry.loc[geometry.island_id.isin(selected)].set_index("island_id")
    if len(polygons)!=80:
        raise ValueError("Some outcome-blind selected original polygons not recovered")
    out_rows=[]
    scans=[]
    with httpx.Client(timeout=75,follow_redirects=True,headers={
       "User-Agent":"fd-island-floral-temporal-feasibility/20261113 (GBIF research pilot)"
    }) as client:
        resolved=resolve_query_taxon(client,"Apidae","FAMILY",None,None)
        taxon_key=resolved["usage_key"]
        if not taxon_key:raise ValueError("GBIF current accepted Apidae taxon could not be resolved")
        for island in selected:
            geom=polygons.loc[island,"geometry"]
            try:
                parts=exact_query_wkts(geom)
            except Exception:
                parts=[]
            for window in PERIODS:
                if not parts:
                    scans.append(dict(island_id=island,period=window,
                         complete=False,truncated=False,api_error="OriginalGeometryFailure",
                         n_raw_examined=0,n_usable_records=0))
                    continue
                rows,meta=scan_window(client,geometry=geom,parts=parts,
                   taxon_key=taxon_key,island=island,window=window)
                out_rows.extend(rows)
                scans.append(meta)
    raw=pd.DataFrame.from_records(out_rows,columns=(
        "island_id","period","dataset_key","year","gbif_id","genus"))
    meta=pd.DataFrame.from_records(scans)
    result=aggregate_repeated_dataset_year_evidence(raw,meta,private)
    result.update(
        original_928_flora_preselected_before_GBIF_outcomes=True,
        source_provenance="verified_original_GSHHG_2.3.7_h_same_WKB",
        GSHHG_zip_sha256=EXPECTED_GSHHG,
        source_trait_checkpoint="island_original_progressive_artifact_10058653212",
        GBIF_v1_search_live_snapshot_not_DOI_pinned=True,
        n_original_flora_islands_unmatchable_to_GSHHG=len(elig-set(geometry.island_id)),
        not_actual_pollinator_abundance_fitness_or_ecological_decline=True,
        time_matched_plant_phenotype_evolution_unavailable=True,
    )
    if any(k in json.dumps(result) for k in (
       '"island_id"','"gbif_id"','"dataset_key"','"accepted_species"',
       '"longitude"','"latitude"')):
        raise ValueError("Private identity field found in aggregate output")
    z.output.parent.mkdir(parents=True,exist_ok=True)
    z.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":result["status"],
        "n_islands":result["fixed_islands"],
        "n_qualifying_repeated_dataset_islands":
            result["n_islands_with_persistent_dataset_year_evidence"],
        "n_regions_reaching_support":
            result["n_regimes_reaching_predeclared_20_island_8_block_gate"],
        "by_region":[{k:v for k,v in row.items() if k.startswith("n_") or k in ("regime","status")}
                     for row in result["by_regime"]],
    },sort_keys=True))
if __name__=="__main__":
    main()
