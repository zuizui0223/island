"""Original SHA-pinned observed dated LCVP phylogeny matching island floral rarity.

Only DIRECT source taxa with exact published tree tips, accepted species name
and independently verified original family and genus agreement are admitted.
Never impute missing tips from taxonomy; do not treat source GIFT as ancestors.
Entire original plant/tip/tree-distance matrix is kept private in this CI job.
"""
from __future__ import annotations
import argparse,hashlib,json,math,subprocess
from pathlib import Path
import numpy as np
import pandas as pd

from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from fd.colxr_family_opportunity import strict_candidate_family_map
from fd.true_backbone_phylo_distance import true_tip_patristic_distance
from fd.true_backbone_phylo_source_match import original_true_phylogeny_match

from run_fd_gbotb_direct_phylogeny_preflight import direct_triplets

PINNED_TREE_SHA={
    "TPL":"9ac19e04d6016c9e7ee8758524996f40a089271b2fd34371a5dbd607b5dbf2a3",
    "LCVP":"30fd7e9e67bdfff427b336d3b35983f63e284ef118e91f76a7f51a884c9b4e5b",
    "WP":"7fc97970fcd9e69f6506b0c860b1296f4a797635193a8ce4bcb052a83fd1f965"
}
EXPECTED_TRUE_TIPS={"TPL":1187,"LCVP":1162,"WP":1160}
EXPECTED_FULL_TIPS={"TPL":74529,"LCVP":73420,"WP":72570}
PINNED_COLXR_SHA="079d9c27b0ca3aa51f0c7310f43ad4901b07a319d6063eab5a7d62b291bdd9e3"


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda:stream.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def cleanup(result):
    if isinstance(result,dict):return {str(k):cleanup(v) for k,v in result.items()}
    if isinstance(result,list):return [cleanup(x) for x in result]
    if isinstance(result,float) and not math.isfinite(result):return None
    return result


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--original-root",type=Path,required=True)
    a.add_argument("--crosswalk",type=Path,required=True)
    a.add_argument("--backbone",type=Path,required=True)
    a.add_argument("--workdir",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    a.add_argument("--tree-version",choices=("TPL","LCVP","WP"),default="LCVP")
    a.add_argument("--rarefactions",type=int,default=20)
    a.add_argument("--reference-null",type=int,default=199)
    a.add_argument("--bootstrap",type=int,default=999)
    args=a.parse_args()
    if sha256(args.crosswalk)!=PINNED_COLXR_SHA:
        raise ValueError("Frozen original CoL XR source taxonomy changed")
    version=args.tree_version
    tree_sha=PINNED_TREE_SHA[version]
    if sha256(args.backbone/f"data/GBOTB.extended.{version}.rda")!=tree_sha:
        raise ValueError("True dated original source tree hash altered")
    root=args.original_root
    ledger=pd.read_csv(root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        usecols=["accepted_species","trait_name","normalized_value",
                 "evidence_scope","evidence_quality"],dtype=str)
    flora=pd.read_csv(root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        usecols=["island_id","accepted_species","origin_status"],dtype=str)
    eligible,receipt=strict_direct_928_flora(ledger,flora)
    if receipt["n_original_eligible_islands"]!=928:
        raise ValueError("Original Direct source 928 island population drift")
    direct=direct_triplets(ledger)
    family,audit=strict_candidate_family_map(
        pd.read_csv(args.crosswalk,dtype=str),set(direct.index))
    if len(family)!=2190:
        raise ValueError("Original accepted species family source is not frozen")
    private=args.workdir/"private"
    private.mkdir(parents=True,exist_ok=True)
    info=private/"original_direct_species_names_and_accepted_families.tsv"
    pd.DataFrame({
        "accepted_species":list(direct.index),
        "source_family":[family.get(x,"") for x in direct.index],
        "source_genus":[x.split()[0] for x in direct.index],
        "rarity_tier":["audit_only"]*len(direct),
    }).to_csv(info,index=False,sep="\t")
    placements=private/"original_true_tip_placement.tsv"
    subprocess.run(["Rscript","scripts/fd_gbotb_original_tip_match.R",
       str(args.backbone/"data"),str(info),str(placements)],
       check=True,timeout=120)
    p=pd.read_csv(placements,sep="\t",dtype=str)
    true=set(p.loc[p.tree_version.eq(version)&
            p.supported_species_placement.eq("exact_backbone_tip"),
            "accepted_species"])
    if len(true)!=EXPECTED_TRUE_TIPS[version]:
        raise ValueError("Original published true phylogenetic tip count changed")
    edge_file=private/f"dated_{version}_original_edges.tsv"
    tip_file=private/f"{version}_source_tip_accepted_species_node.tsv"
    subprocess.run(["Rscript","scripts/fd_gbotb_export_private_edges.R",
       str(args.backbone/"data"),str(info),version,
       str(edge_file),str(tip_file)],
       check=True,timeout=120)
    edges=pd.read_csv(edge_file,sep="\t")
    tips=pd.read_csv(tip_file,sep="\t")
    if set(edges.columns)!={"parent_node_id","child_node_id",
                             "original_branch_length"}:
        raise ValueError("Published phylogeny graph source schema invalid")
    tips=tips.loc[tips.accepted_species.isin(true)].copy()
    if len(tips)!=EXPECTED_TRUE_TIPS[version] or tips.accepted_species.duplicated().any():
        raise ValueError("One original accepted taxon is not EXACT observed tree tip")
    ordered=tips.sort_values("accepted_species")
    tree_tip_nodes=ordered.tip_node_id.to_numpy(dtype=np.int32)
    name_to_tip=dict(zip(ordered.accepted_species,
                         ordered.tip_node_id,strict=True))
    D=true_tip_patristic_distance(
       edges[["parent_node_id","child_node_id"]].to_numpy(dtype=np.int64),
       edges.original_branch_length.to_numpy(dtype=float),
       tree_tip_nodes,n_original_tree_tips=EXPECTED_FULL_TIPS[version])
    if np.any(D<0) or not np.isfinite(D).all():
        raise ValueError("Original true phylogenetic distance matrix invalid")
    geo=pd.read_csv(root/"fixed/isolation/results/purpose_shortest_island_data.csv",
       usecols=["island_id","analysis_regime","spatial_block"],dtype=str)
    donors=pd.read_csv(root/"fixed/source/result/island_source_assignments.csv.gz",
       usecols=["island_id","source_mode","source_rank","entity_ID"],dtype=str)
    gift=pd.read_csv(root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
       usecols=["entity_ID","work_species"],dtype=str)
    result=original_true_phylogeny_match(
       direct,eligible,geo,donors,gift,
       name_to_tip,ordered.accepted_species.tolist(),D,
       tree_version=version,original_accepted_family=family,
       rarefactions=args.rarefactions,
       reference_draws=args.reference_null,bootstrap=args.bootstrap)
    result["original_taxonomy_source_family_candidates"]=len(family)
    result["original_GBOTB_tree_SHA256"]=tree_sha
    result["original_accepted_species_exact_tip_names_agreed_to_CoL"]=len(true)
    result["original_phylogeny_tree_tips_with_distances"]=len(ordered)
    result["external_plant_phylogeny_GPL2_source_not_redistributed"]=True
    result=cleanup(result)
    txt=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    if any('"'+word+'"' in txt for word in (
        "accepted_species","source_genus","source_family",
        "island_id","work_species","family","genus","entity_ID",
        "spatial_block","source_signature","tip_label")):
        raise ValueError("Original source taxa/phylo names have leaked")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(txt+"\n")
    print(json.dumps({
       "status":result["status"],
       "n_original_true_tips":len(ordered),
       "n_original_eligible_flora":928,
       "tree_version":version,
       "by_region":[{
           "region":x["regime"],"scope":x["scope"],
           "n_true_tree_tips":x["n_dated_true_tree_tips_per_side"],
           "n_real_source_comparable_islands":x["n_original_tree_tip_eligible_islands"],
           "n_blocks":x["n_physical_blocks"],
           "n_nearest_source_regions":x["n_first_GIFT_source_signatures"],
           "reportable":x["regional_effect_reportable_20island_10block_10sources"],
       } for x in result["region_phylogenetic_tip_matching_support"]],
    },sort_keys=True))
if __name__=="__main__":
    main()
