"""Source-pinned TRUE dated tree phylogenetic beta vs floral categorical beta.

The original 928-island FD Direct trait flora and protected GIFT native source
rosters never leave original island CI. Published V.PhyloMaker2 3 actual
backbones are read from pinned original source but not redistributed.
Only externally safe anonymous geographic-region/source-cluster aggregate
relationships are published, not source taxon/plant/island coordinates.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from fd.colxr_family_opportunity import strict_candidate_family_map
from fd.true_backbone_phylo_source_match import _strict_source_taxon_groups
from fd.phylo_beta_functional_beta_macro import (
    OriginalTreeSubtree,macro_original_phylo_functional_beta
)
from run_fd_gbotb_direct_phylogeny_preflight import direct_triplets

PINNED_GBOTB={
    "TPL":("9ac19e04d6016c9e7ee8758524996f40a089271b2fd34371a5dbd607b5dbf2a3",1187),
    "LCVP":("30fd7e9e67bdfff427b336d3b35983f63e284ef118e91f76a7f51a884c9b4e5b",1162),
    "WP":("7fc97970fcd9e69f6506b0c860b1296f4a797635193a8ce4bcb052a83fd1f965",1160)
}
TAXON_CROSSWALK_SHA="079d9c27b0ca3aa51f0c7310f43ad4901b07a319d6063eab5a7d62b291bdd9e3"


def _sha(path:Path)->str:
    dig=hashlib.sha256()
    with path.open("rb") as src:
        for batch in iter(lambda:src.read(1<<20),b""):
            dig.update(batch)
    return dig.hexdigest()


def _safe_json(v):
    if isinstance(v,dict):
        return {str(k):_safe_json(x) for k,x in v.items()}
    if isinstance(v,list):
        return [_safe_json(x) for x in v]
    if isinstance(v,(float,np.float64)) and not math.isfinite(v):
        return None
    if isinstance(v,np.integer):
        return int(v)
    return v


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",required=True,type=Path)
    p.add_argument("--colxr-crosswalk",required=True,type=Path)
    p.add_argument("--published-backbone-checkout",required=True,type=Path)
    p.add_argument("--version",choices=("TPL","LCVP","WP"),required=True)
    p.add_argument("--workdir",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--original-island-repeats",type=int,default=20)
    p.add_argument("--profile-k-null",type=int,default=199)
    p.add_argument("--bootstrap",type=int,default=999)
    a=p.parse_args()
    sha,expected=PINNED_GBOTB[a.version]
    if _sha(a.published_backbone_checkout/f"data/GBOTB.extended.{a.version}.rda")!=sha:
        raise ValueError("Published original dated GBOTB source tree SHA mismatch")
    if _sha(a.colxr_crosswalk)!=TAXON_CROSSWALK_SHA:
        raise ValueError("Original CoL XR strict unreviewed source family candidacy changed")
    orig=a.original_root
    direct_frame=pd.read_csv(orig/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        usecols=["accepted_species","trait_name","normalized_value",
                 "evidence_scope","evidence_quality"],dtype=str)
    original_flora=pd.read_csv(orig/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        usecols=["island_id","accepted_species","origin_status"],dtype=str)
    flora,gate=strict_direct_928_flora(direct_frame,original_flora)
    if gate["n_original_eligible_islands"]!=928:
        raise ValueError("Original high/medium single-state Direct 928 islands changed")
    triple=direct_triplets(direct_frame)
    family,quality=strict_candidate_family_map(
         pd.read_csv(a.colxr_crosswalk,dtype=str),set(triple.index))
    if len(triple)!=2234 or len(family)!=2190:
        raise ValueError("Original Direct 2234 and strict candidate family 2190 changed")
    private=a.workdir/"private_no_export"
    private.mkdir(parents=True,exist_ok=True)
    input_names=private/"original_fd_direct_taxa_with_family.tsv"
    pd.DataFrame({
        "accepted_species":list(triple.index),
        "source_family":[family.get(x,"") for x in triple.index],
        "source_genus":[x.split()[0] for x in triple.index],
        "rarity_tier":["audit_only"]*len(triple)
    }).to_csv(input_names,index=False,sep="\t")
    placements=private/"original_published_tree_tip_names.tsv"
    subprocess.run(["Rscript","scripts/fd_gbotb_original_tip_match.R",
         str(a.published_backbone_checkout/"data"),
         str(input_names),str(placements)],check=True,timeout=120)
    ptable=pd.read_csv(placements,dtype=str,sep="\t")
    names=set(ptable.loc[
        ptable.tree_version.eq(a.version)&
        ptable.supported_species_placement.eq("exact_backbone_tip"),
        "accepted_species"])
    if len(names)!=expected:
        raise ValueError("Actual dated GBOTB original tip acceptance drift")
    edges_file=private/f"{a.version}_source_rooted_edges.tsv"
    tips_file=private/f"{a.version}_source_exact_tip_nodes.tsv"
    subprocess.run(["Rscript","scripts/fd_gbotb_export_private_edges.R",
         str(a.published_backbone_checkout/"data"),
         str(input_names),a.version,str(edges_file),str(tips_file)],
         check=True,timeout=120)
    edges=pd.read_csv(edges_file,sep="\t")
    tips=pd.read_csv(tips_file,sep="\t",dtype={"accepted_species":str,"tip_node_id":np.int64})
    tips=tips.loc[tips.accepted_species.isin(names)].copy()
    if len(tips)!=expected or tips.accepted_species.duplicated().any():
        raise ValueError("Original actual tree tip-to-accepted-plant mapping ambiguous")
    tipmap=dict(zip(tips.accepted_species,
                    tips.tip_node_id.astype(int),strict=True))
    tree=OriginalTreeSubtree(
        edges[["parent_node_id","child_node_id"]].to_numpy(dtype=int),
        edges.original_branch_length.to_numpy(float),
        tips.tip_node_id.to_numpy(dtype=int)
    )
    geo=pd.read_csv(orig/"fixed/isolation/results/purpose_shortest_island_data.csv",
        usecols=["island_id","analysis_regime","spatial_block",
                 "log_island_area_km2","log_distance_to_continent_km"],dtype=str)
    for c in ("log_island_area_km2","log_distance_to_continent_km"):
        geo[c]=pd.to_numeric(geo[c],errors="raise")
    donor=pd.read_csv(orig/"fixed/source/result/island_source_assignments.csv.gz",
        usecols=["island_id","source_mode","source_rank","entity_ID"],dtype=str)
    mainland=pd.read_csv(orig/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        usecols=["entity_ID","work_species"],dtype=str)
    result=macro_original_phylo_functional_beta(
        triple,flora,geo,donor,mainland,tree,
        sorted(names),tipmap,family,tree_version=a.version,
        repeats=a.original_island_repeats,
        reference_null_draws=a.profile_k_null,bootstrap=a.bootstrap
    )
    result["original_published_tree_source_sha256"]=sha
    result["original_COLXR_crosswalk_sha256"]=TAXON_CROSSWALK_SHA
    result["original_tree_source_version_SHA_verified"]=True
    result["source_accepted_tip_name_family_genus_ambiguities_rejected"]=True
    result["no_original_true_phylogeny_backbone_or_taxon_record_redistributed"]=True
    result=_safe_json(result)
    content=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    forbidden=('island_id','accepted_species','source_genus',
        'source_family','candidate_family','entity_ID','source_signature',
        'spatial_block','work_species','family','genus','tip_node_id')
    if any('"'+term+'"' in content for term in forbidden):
        raise ValueError("Private original source/taxon/island identity leaked into public summary")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(content+"\n")
    print(json.dumps({
        "status":result["status"],"tree":a.version,
        "original_eligible_flora":result["n_original_928_direct_island_population"],
        "original_true_tree_tips":result["n_original_true_phylogenetic_tip_species"],
        "island_scopes_scored":result["n_original_private_island_scope_matched_rows"],
        "regions":[{
           "scope":r["scope"],"region":r["regime"],
           "n_original_islands":r["n_original_eligible_islands"],
           "n_original_blocks":r["n_physical_spatial_blocks"],
           "n_source_groups":r["n_first_nearest_GIFT_source_signatures"],
           "supported":r["passed_20islands_10blocks_10source_support"],
           "median_phylogenetic_beta":r["median_original_phylo_beta"],
           "median_functional_beta":r["median_original_functional_beta"]
        }for r in result["region_support_and_categorical_beta_relationships"]]
    },sort_keys=True))
if __name__=="__main__":
    main()
