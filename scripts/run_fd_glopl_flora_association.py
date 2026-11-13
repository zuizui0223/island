"""Exploratory *island-level* GloPL hand-pollen response vs floral-source FD.

This runner is only executed with original source-pinned verified inputs.
It first freezes independent island and potential GIFT mainland plant-morph
outcomes and checks sufficient study/island support *before reading GloPL
effect sizes*. The GloPL experiment species are NOT matched to those
direct flower trait names (original co-location preflight found zero).
No rights-restricted original taxa/island/site/study IDs leave CI.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import geopandas as gpd
import numpy as np
import pandas as pd

from island_v2.chapter1_h5_glopl_island_overlap_preflight import (
    load_config, read_glopl_metadata, match_glopl_to_islands
)
from fd.pollinator_functional_diversity import direct_triple, sample_morphology
from fd.glopl_floral_fd_association import audit_crosssectional_pollen_fd

GLOPL_PIN="264ca8c2237f126b3c209fcba039e19b291324ddf6e7a96992dc8235a55984de"
GSHHG_PIN="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
PRIMARY=("marginal_rao","joint_rao")
ALSO=("open_radial_share",)
NRAREF=30


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for blk in iter(lambda:f.read(1<<20),b""):
            h.update(blk)
    return h.hexdigest()


def _seed(name):
    return int.from_bytes(hashlib.sha256(("fd-PL-v1:"+str(name)).encode()).digest()[:8],"big")


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--glopl",type=Path,required=True)
    a.add_argument("--glopl-config",type=Path,required=True)
    a.add_argument("--progressive-root",type=Path,required=True)
    a.add_argument("--original-gshhg-gpkg",type=Path,required=True)
    a.add_argument("--original-gshhg-receipt",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    p=a.parse_args()
    if sha256(p.glopl)!=GLOPL_PIN:
        raise ValueError("Original GloPL 2018 data file changed: reject")
    receipt=json.loads(p.original_gshhg_receipt.read_text())
    if (receipt.get("GSHHG_archive_SHA256")!=GSHHG_PIN
        or receipt.get("original_frozen_island_count")!=8265
        or receipt.get("n_original_focal_islands")!=21
        or receipt.get("n_focal_islands_with_exact_historical_polygon")!=21
        or receipt.get("status") not in (
            "FULL_8265_GLOBAL_GEOMETRY_IDENTICAL",
            "FOCAL_21_POLYGONS_IDENTICAL_BUT_FULL_UNIVERSE_DRIFT")
        or sha256(p.original_gshhg_gpkg)!=receipt.get("geometry_file_SHA256")):
        raise ValueError("Original GSHHG source/polygon SHA does not match")
    config=load_config(p.glopl_config)
    metadata=read_glopl_metadata(p.glopl,config) # no effect or treatment values
    geography=gpd.read_file(p.original_gshhg_gpkg,layer="islands")
    if geography.crs is None or geography.crs.to_epsg()!=4326:
        raise ValueError("Wrong original WGS84 geometry")
    matched,overlap=match_glopl_to_islands(metadata,geography,config)
    root=p.progressive_root
    ledger=pd.read_csv(root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
       usecols=["accepted_species","trait_name","normalized_value",
                "evidence_scope","evidence_quality"])
    traits=direct_triple(ledger)
    florida=pd.read_csv(root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
       usecols=["island_id","accepted_species","origin_status"])
    florida=florida.loc[florida.accepted_species.isin(traits.index)]
    if florida.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original flora has duplicate island species")
    coverage=florida.groupby("island_id").accepted_species.nunique()
    eligible=set(coverage[coverage>=30].index)
    if len(eligible)!=928:
        raise ValueError("Original source Direct morphology 928 island cohort drifted")
    g=pd.read_csv(root/"fixed/isolation/results/purpose_shortest_island_data.csv",
       usecols=["island_id","analysis_regime","spatial_block",
               "log_distance_to_continent_km","log_island_area_km2"])
    if g.island_id.duplicated().any():
        raise ValueError("Duplicate frozen geography island")
    gmap=g.set_index("island_id")
    source=pd.read_csv(root/"fixed/source/result/island_source_assignments.csv.gz",
       usecols=["island_id","source_mode","source_rank","entity_ID"])
    source=source.loc[source.island_id.isin(eligible)&source.source_mode.eq("geo_k5")&
       source.source_rank.between(1,5)].copy()
    if source.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original GIFT source rank duplicated")
    sources=source.sort_values(["island_id","source_rank"]).groupby("island_id").entity_ID.agg(tuple).to_dict()
    if any(len(sources.get(i,()))!=5 or len(set(sources[i]))!=5 for i in eligible):
        raise ValueError("Frozen five GIFT geographic source regions incomplete")
    pool=pd.read_csv(root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
       usecols=["entity_ID","work_species"])
    pool=pool.loc[pool.work_species.isin(traits.index)].drop_duplicates()
    gift=pool.groupby("entity_ID").work_species.agg(set).to_dict()
    glopl_eligible=matched.loc[
       matched.island_id.isin(eligible)&matched.study_key.fillna("").astype(str).str.len().gt(0)
    ]
    islands=sorted(set(glopl_eligible.island_id))
    floras={i:part for i,part in florida.loc[
       florida.island_id.isin(islands)].groupby("island_id")}
    records=[]
    for island in islands:
        if island not in gmap.index:
            continue
        name=sorted(set(floras[island].accepted_species))
        if len(name)<30:
            continue
        candidate=set().union(*(gift.get(j,set()) for j in sources[island]))
        candidate&=set(traits.index)
        if len(candidate)<30:
            continue
        candidate=sorted(candidate)
        local=sample_morphology(traits.loc[name].to_numpy(dtype=str),
            repeats=NRAREF,rng=np.random.default_rng(_seed(island+":island")))
        mainland=sample_morphology(traits.loc[candidate].to_numpy(dtype=str),
            repeats=NRAREF,rng=np.random.default_rng(_seed(island+":mainland")))
        record=dict(island_id=island,
          analysis_regime=str(gmap.loc[island,"analysis_regime"]),
          spatial_block=str(gmap.loc[island,"spatial_block"]),
          log_distance_to_continent_km=float(
              gmap.loc[island,"log_distance_to_continent_km"]),
          log_island_area_km2=float(gmap.loc[island,"log_island_area_km2"]))
        for t in (*PRIMARY,*ALSO):
            record["source_delta_"+t]=float(local[t]-mainland[t])
        records.append(record)
    fd=pd.DataFrame.from_records(records)
    observed=glopl_eligible.loc[
       glopl_eligible.island_id.isin(fd.island_id)].copy()
    nblocks=observed.merge(g[["island_id","spatial_block"]],
        on="island_id",how="left").spatial_block.nunique()
    counts=observed.groupby("analysis_regime").island_id.nunique().to_dict() if "analysis_regime" in observed else (
        observed.merge(g[["island_id","analysis_regime"]],on="island_id").groupby(
          "analysis_regime").island_id.nunique().to_dict()
    )
    # Crucial source/experimental geography support gate before inspecting
    # ANY pollen-supplementation effect sizes.
    gate=(len(fd)>=25 and nblocks>=15 and
          observed.study_key.nunique()>=20 and
          counts.get("northern_midlatitude",0)>=10 and
          counts.get("tropical",0)>=10)
    payload=dict(
       status="OUTCOME_BLIND_ISLAND_FLORAL_FD_EXPERIMENTAL_COVERAGE",
       n_source_supported_islands_with_GloPL_preflight=int(len(fd)),
       n_original_spatial_blocks_before_effects=int(nblocks),
       n_published_GloPL_studies_before_effects=int(observed.study_key.nunique()),
       region_counts_before_effects={k:int(v) for k,v in counts.items()},
       original_GloPL_source_sha256=GLOPL_PIN,
       original_GSHHG_source_zip_sha256=GSHHG_PIN,
       direct_source_taxon_flora_population=928,
       original_hand_pollen_measurements_on_same_flower_species_count=0,
       genuine_pollinator_decline_or_shift_measured=False,
       pre_effect_outcome_blind_gate_passed=bool(gate),
       n_species_rarefaction_per_island_and_mainland=30,
       n_rarefaction_draws=NRAREF,
       exact_study_island_species_geography_not_exported=True)
    if not gate:
        payload["status"]="INSUFFICIENT_GLOPL_ISLAND_FD_SOURCE_SUPPORT_NO_EFFECTS_READ"
        payload["effect_sizes_read"]=False
    else:
        effects=pd.read_csv(p.glopl,encoding="latin-1",
            usecols=["PL_Effect_Size","PL_Effect_Size_Type2"],dtype=str).reset_index(drop=True)
        if len(effects)!=len(matched):
            raise ValueError("GloPL original metadata/effect rows not identical")
        effects["row_id"]=np.arange(len(effects),dtype=int)
        effect=observed[["row_id","island_id","study_key"]].merge(
            effects,on="row_id",how="left",validate="one_to_one")
        effect["PL_Effect_Size"]=pd.to_numeric(effect.PL_Effect_Size,errors="coerce")
        effect=effect.loc[np.isfinite(effect.PL_Effect_Size)].copy()
        if effect.empty:
            payload["status"]="NO_FINITE_GLOPL_EXPERIMENTAL_EFFECTS"
            payload["effect_sizes_read"]=True
        else:
            cells=effect.groupby(["island_id","study_key"],as_index=False).agg(
              PL_Effect_Size_mean=("PL_Effect_Size","mean"))
            payload["all_measurement_types"]=audit_crosssectional_pollen_fd(fd,cells)
            # Prior original H5 reported narrower supplemental-only sensitivity.
            # This subset is NOT allowed to rescue a failed main claim.
            supplemental=effect.loc[effect.PL_Effect_Size_Type2.str.strip().eq("Sup")]
            if len(supplemental):
                supp=supplemental.groupby(["island_id","study_key"],as_index=False).agg(
                    PL_Effect_Size_mean=("PL_Effect_Size","mean"))
                payload["supplemental_only"]=audit_crosssectional_pollen_fd(fd,supp)
            else:
                payload["supplemental_only"]={"status":"NO_SUPPLEMENTAL_RECORDS"}
            payload["effect_sizes_read"]=True
            payload["status"]="EXPLORATORY_CROSSSECTIONAL_GLOPL_ISLAND_FUNCTIONAL_ASSOCIATION"
    raw=json.dumps(payload,sort_keys=True)
    forbidden=['"island_id"','"study_key"','"accepted_species"','"entity_ID"',
               '"Species_accepted_names"','"Longitude"','"Latitude"']
    if any(a in raw for a in forbidden):
        raise ValueError("Rights-protected original identities leaked into GloPL public result")
    p.output.parent.mkdir(parents=True,exist_ok=True)
    p.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],"source_islands":payload["n_source_supported_islands_with_GloPL_preflight"],
        "blocks":payload["n_original_spatial_blocks_before_effects"],
        "read_effects":payload.get("effect_sizes_read",False),
        "outcomes":payload.get("all_measurement_types",{}).get("status","not_tested"),
        "supplemental":payload.get("supplemental_only",{}).get("status","not_tested")
    },sort_keys=True))
if __name__=="__main__":
    main()
