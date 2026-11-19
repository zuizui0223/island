"""Original protected 928-island Direct/GIFT/WCVP true-tip floral pairing.

Runs independent published TPL/LCVP/WorldPlants GBOTB exact-tree-tip
source qualification. All original island taxon names, original tree edges,
dated source distances, WCVP/CoL identities and GIFT candidate source
names stay in CI. Only anonymous regional conditionally bootstrapped
results can leave the source runner.

NO historical pollinator decline, pollen effectiveness, true ancestral
population or within-species island floral adaptation is measured.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd

from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from fd.colxr_family_opportunity import strict_candidate_family_map
from fd.regional_native_provenance_volume import regional_native_candidates
from fd.regional_native_flora_fd import TRAITS
from fd.true_backbone_phylo_distance import true_tip_patristic_distance
from fd.true_tip_native_pairing import original_true_tip_native_pairing

TREE_SHA={
    "TPL":"9ac19e04d6016c9e7ee8758524996f40a089271b2fd34371a5dbd607b5dbf2a3",
    "LCVP":"30fd7e9e67bdfff427b336d3b35983f63e284ef118e91f76a7f51a884c9b4e5b",
    "WP":"7fc97970fcd9e69f6506b0c860b1296f4a797635193a8ce4bcb052a83fd1f965"
}
N_TRUE_TIPS={"TPL":1187,"LCVP":1162,"WP":1160}
N_FULL_TIPS={"TPL":74529,"LCVP":73420,"WP":72570}
COLXR_SHA="079d9c27b0ca3aa51f0c7310f43ad4901b07a319d6063eab5a7d62b291bdd9e3"


def verify(path:Path,expected:str):
    if not path.is_file():
        raise ValueError("Original frozen input missing")
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    if h.hexdigest()!=expected:
        raise ValueError("Original rights-restricted input hash mismatch")


def original_direct_triples(ledger:pd.DataFrame)->pd.DataFrame:
    x=ledger.loc[
        ledger.trait_name.isin(TRAITS)&
        ledger.evidence_scope.eq("species_direct")&
        ledger.evidence_quality.isin(("high","medium"))&
        ledger.normalized_value.notna()&
        ~ledger.normalized_value.astype(str).str.contains("|",regex=False)
    ]
    if x.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Original accepted species Direct trait conflict")
    d=x.pivot(index="accepted_species",columns="trait_name",
              values="normalized_value").dropna(subset=list(TRAITS))[list(TRAITS)]
    if len(d)!=2234:
        raise ValueError("Original Direct colour-form-symmetry species changed")
    return d


def run(root:Path,crosswalk:Path,backbone:Path,wcvp_status:Path,
        island_tdwg:Path,private:Path,*,tree_version:str,seed:int)->dict:
    if tree_version not in TREE_SHA:
        raise ValueError("Unsupported original published source backbone")
    verify(crosswalk,COLXR_SHA)
    verify(backbone/"data"/f"GBOTB.extended.{tree_version}.rda",
           TREE_SHA[tree_version])
    ledger=pd.read_csv(
        root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        dtype=str,usecols=["accepted_species","trait_name","normalized_value",
                            "evidence_scope","evidence_quality"])
    original_flora=pd.read_csv(
        root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        dtype=str,usecols=["island_id","accepted_species","origin_status"])
    eligible,gate=strict_direct_928_flora(ledger,original_flora)
    if gate["n_original_eligible_islands"]!=928:
        raise ValueError("Original original-island Direct roster altered")
    direct=original_direct_triples(ledger)
    col=pd.read_csv(crosswalk,dtype=str)
    family,fa=strict_candidate_family_map(col,set(direct.index))
    if len(family)!=2190 or fa["n_candidate_families"]<50:
        raise ValueError("Source-checked CoL XR family candidates changed")
    native=pd.read_csv(
        wcvp_status,dtype=str,
        usecols=["accepted_species","native_l3_codes"]).drop_duplicates()
    mappings=pd.read_csv(
        island_tdwg,dtype=str,
        usecols=["island_id","tdwg_l3_code","tdwg_match_status"])
    if len(mappings)!=8265 or native.accepted_species.nunique()!=55137:
        raise ValueError("Original pinned WCVP native/TDWG mapping changed")
    classified=regional_native_candidates(eligible,mappings,native)
    n=classified.groupby("island_id")[
        ["strict_source_native","regional_native_compatible"]].sum()
    if (n.strict_source_native>=30).sum()!=87 or (
        (n.regional_native_compatible>=30).sum()!=357):
        raise ValueError("Original native/WCVP candidate cohort eligibility drift")
    private.mkdir(parents=True,exist_ok=True)
    src=private/"source_original_2234_direct_taxon_family.tsv"
    pd.DataFrame({
        "accepted_species":list(direct.index),
        "source_family":[family.get(i,"") for i in direct.index],
        "source_genus":[str(i).split()[0] for i in direct.index],
        "rarity_tier":["audit_only"]*len(direct),
    }).to_csv(src,sep="\t",index=False)
    placements=private/"exact_observed_source_tree_tip_support.tsv"
    subprocess.run([
        "Rscript","scripts/fd_gbotb_original_tip_match.R",
        str(backbone/"data"),str(src),str(placements)
    ],check=True,timeout=150)
    p=pd.read_csv(placements,sep="\t",dtype=str)
    valid=set(p.loc[p.tree_version.eq(tree_version)&
                  p.supported_species_placement.eq("exact_backbone_tip"),
                  "accepted_species"])
    if len(valid)!=N_TRUE_TIPS[tree_version]:
        raise ValueError("Original accepted species true tip count changed")
    edgesfile=private/f"{tree_version}_original_true_tree_edges.tsv"
    tipsfile=private/f"{tree_version}_original_exact_tip_names.tsv"
    subprocess.run([
        "Rscript","scripts/fd_gbotb_export_private_edges.R",
        str(backbone/"data"),str(src),tree_version,
        str(edgesfile),str(tipsfile)
    ],check=True,timeout=150)
    e=pd.read_csv(edgesfile,sep="\t")
    t=pd.read_csv(tipsfile,sep="\t",dtype={"accepted_species":str,"tip_node_id":int})
    t=t.loc[t.accepted_species.isin(valid)].copy()
    if len(t)!=N_TRUE_TIPS[tree_version] or t.accepted_species.duplicated().any():
        raise ValueError("Original published tree exact tip source mismatch")
    ordered=t.sort_values("accepted_species")
    original_tree_nodes=ordered.tip_node_id.to_numpy(dtype=np.int32)
    D=true_tip_patristic_distance(
        e[["parent_node_id","child_node_id"]].to_numpy(dtype=np.int64),
        e.original_branch_length.to_numpy(dtype=float),
        original_tree_nodes,n_original_tree_tips=N_FULL_TIPS[tree_version]
    )
    tip_order=ordered.accepted_species.tolist()
    name_to_tip={name:i for i,name in enumerate(tip_order)}
    geo=pd.read_csv(
        root/"fixed/isolation/results/purpose_shortest_island_data.csv",
        dtype=str,usecols=["island_id","analysis_regime","spatial_block"])
    geo=geo.loc[geo.island_id.isin(classified.island_id)].drop_duplicates("island_id")
    if len(geo)!=928 or geo.spatial_block.isna().any():
        raise ValueError("Original geography failed independent island support")
    assign=pd.read_csv(
        root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    assign=assign.loc[assign.island_id.isin(classified.island_id)&
                      assign.source_mode.eq("geo_k5")].copy()
    assign["source_rank"]=pd.to_numeric(assign.source_rank,errors="raise")
    assign=assign.loc[assign.source_rank.between(1,5)]
    if assign.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original candidate GIFT nearest source ranks duplicate")
    near=assign.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    if any(len(near.get(i,()))!=5 or len(set(near[i]))!=5
           for i in classified.island_id.unique()):
        raise ValueError("Original source donor region opportunity missing")
    mainland=pd.read_csv(
        root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    mainland=mainland.loc[
        mainland.work_species.isin(direct.index)].drop_duplicates()
    donors=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    result=original_true_tip_native_pairing(
        direct,classified,geo,near,donors,
        family,name_to_tip,tip_order,D,
        tree_version=tree_version,repeats=24,seed=seed,cluster_draws=999)
    if result["original_eligible_islands"]!=928 or (
        result["n_original_true_tip_source_species"]!=N_TRUE_TIPS[tree_version]):
        raise ValueError("Original tree source result integrity failed")
    result["original_tree_SHA256"]=TREE_SHA[tree_version]
    result["original_colxr_SHA256"]=COLXR_SHA
    result["original_family_source_candidate_n"]=len(family)
    result["tree_and_flora_source_are_rights_restricted_no_identifiers_exported"]=True
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-root",type=Path,required=True)
    ap.add_argument("--crosswalk",type=Path,required=True)
    ap.add_argument("--backbone",type=Path,required=True)
    ap.add_argument("--wcvp-status",type=Path,required=True)
    ap.add_argument("--island-tdwg",type=Path,required=True)
    ap.add_argument("--workdir",type=Path,required=True)
    ap.add_argument("--tree-version",choices=("TPL","LCVP","WP"),required=True)
    ap.add_argument("--seed",type=int,default=20261114)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    result=run(
        args.original_root,args.crosswalk,args.backbone,args.wcvp_status,
        args.island_tdwg,args.workdir/"private",
        tree_version=args.tree_version,seed=args.seed)
    text=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    forbidden=("accepted_species","island_id","source_family","source_genus",
      "work_species","entity_ID","source_signature","spatial_block",
      "tip_label","tip_node_id","native_l3_codes","candidate_family",
      "tdwg_l3_code","family_map","longitude","latitude")
    if any('"'+k+'"' in text for k in forbidden):
        raise ValueError("Original restricted tree/island/source/taxon identity escaped")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(text+"\n")
    panels=result["regional_conditional_phylo_pairing"]
    print(json.dumps({
        "status":result["status"],
        "tree_version":args.tree_version,
        "seed":args.seed,
        "n_original_eligible_islands":928,
        "n_published_true_tip_Direct_species":result["n_original_true_tip_source_species"],
        "n_region_condition_panels":len(panels),
        "n_exploratory_supported_panels":sum(
            x["exploratory_support_10island_5block_5source"] for x in panels),
        "n_strict_supported_panels":sum(
            x["strict_support_20island_10block_10source"] for x in panels),
    },sort_keys=True))


if __name__=="__main__":main()
