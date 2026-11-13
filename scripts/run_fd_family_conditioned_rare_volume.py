"""Protected original-source CoL XR candidate-family floral FD experiment.

Uses two frozen source archives only:
 - original 928-island Direct high/medium floral plant trait checkpoint,
 - CoL XR 2026 alpha2 crosswalk: unpromoted but exact accepted plant
   family candidate that agrees with independent original reported family.

Outputs only regional family turnover and family-shared floral category
reference volumes, no plant/island/genus/family names or coordinates.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
import pandas as pd

from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from fd.discrete_functional_volume import TRAITS
from fd.colxr_family_opportunity import strict_candidate_family_map
from fd.family_conditioned_rare_volume import original_family_rare_support

SHA_COLXR="079d9c27b0ca3aa51f0c7310f43ad4901b07a319d6063eab5a7d62b291bdd9e3"


def sha(path:Path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for data in iter(lambda:f.read(1<<20),b""):
            h.update(data)
    return h.hexdigest()


def _direct(ledger):
    x=ledger.loc[
        ledger.trait_name.isin(TRAITS)&
        ledger.evidence_scope.eq("species_direct")&
        ledger.evidence_quality.isin(("high","medium"))&
        ledger.normalized_value.notna()&
        ~ledger.normalized_value.astype(str).str.contains("|",regex=False)
    ]
    if x.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Original strict species trait evidence conflict")
    r=x.pivot(index="accepted_species",columns="trait_name",
        values="normalized_value").dropna(subset=list(TRAITS))[list(TRAITS)]
    if len(r)!=2234:
        raise ValueError("Original global 2234 source species Direct trio changed")
    return r


def _json(x):
    if isinstance(x,dict):
        return {str(k):_json(v) for k,v in x.items()}
    if isinstance(x,list):
        return [_json(v) for v in x]
    if isinstance(x,float) and not math.isfinite(x):
        return None # Unsupported group -> missing, never pseudo-zero.
    return x


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-progressive-root",type=Path,required=True)
    p.add_argument("--colxr-crosswalk",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--rarefactions",type=int,default=40)
    p.add_argument("--reference-null",type=int,default=399)
    z=p.parse_args()
    if sha(z.colxr_crosswalk)!=SHA_COLXR:
        raise ValueError("Frozen 2026 COL XR source accepted classification candidate altered")
    base=z.original_progressive_root
    direct=pd.read_csv(base/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        usecols=["accepted_species","trait_name","normalized_value",
                 "evidence_scope","evidence_quality"])
    floras=pd.read_csv(base/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        usecols=["island_id","accepted_species","origin_status"])
    original,proof=strict_direct_928_flora(direct,floras)
    if proof["n_original_eligible_islands"]!=928:
        raise ValueError("Original 928 independent original polygon IDs drifted")
    triple=_direct(direct)
    tax=pd.read_csv(z.colxr_crosswalk,dtype=str)
    family,checks=strict_candidate_family_map(tax,set(triple.index))
    if len(family)<1900:
        raise ValueError("Insufficient original source-verified COL family candidate coverage")
    geo=pd.read_csv(base/"fixed/isolation/results/purpose_shortest_island_data.csv",
        usecols=["island_id","analysis_regime","spatial_block"])
    donor=pd.read_csv(base/"fixed/source/result/island_source_assignments.csv.gz",
        usecols=["island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(base/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        usecols=["entity_ID","work_species"])
    result=original_family_rare_support(triple,original,geo,donor,
        mainland,family,rarefactions=z.rarefactions,
        reference_draws=z.reference_null,bootstraps=1999,seed=20261113)
    result["strict_CoLXR_candidate_family_mapping_coverage"]={
        "n_original_2234_plant_species":2234,
        "n_original_COLXR_strict_family_candidate_species":len(family),
        "n_candidate_families":checks["n_candidate_families"],
        "candidate_taxonomy_reviewed_by_human":False,
        "candidate_taxonomy_not_genetically_resolved":True,
        "colxr_source_crosswalk_sha256":SHA_COLXR}
    assert result["n_original_eligible_island_flora"]==928
    assert result["inference"]["family_labels_are_unreviewed_exact_accepted_COL_XR_candidates"]
    assert result["inference"]["different_genus_condition_restricts_to_different_genus_in_same_family"]
    assert result["no_original_species_family_source_island_or_coordinates_exported"]
    assert result["inference"]["not_causal_selection_or_phylogenetic_effect"]
    result=_json(result)
    txt=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    forbidden=("island_id","accepted_species","family","genus","family_map",
          "work_species","entity_ID","source_signature","spatial_block",
          "longitude","latitude","candidate_family","submitted_name")
    if any('"'+s+'"' in txt for s in forbidden):
        raise ValueError("Original taxon/island/family/source ID leaked in aggregate")
    z.output.parent.mkdir(parents=True,exist_ok=True)
    z.output.write_text(txt+"\n")
    print(json.dumps({
        "status":result["status"],
        "original_928_islands":result["n_original_eligible_island_flora"],
        "candidate_taxonomy_matched_species":len(family),
        "shapley_panels":len(result["family_shapley_rare_category_macro_panels"]),
        "matched_panels":[{
            "scope":r["scope"],"regime":r["regime"],
            "same_family_count":r["n_exact_shared_families_per_side"],
            "require_different_genus":r["condition"]=="same_family_different_genus",
            "n_islands":r["n_family_matched_islands"],
            "n_blocks":r["n_blocks"],
            "n_source":r["n_first_GIFT_source_signatures"],
            "supported":r["reportable_under_10island_5block_5source_gate"],
        }for r in result["same_family_island_mainland_matched_panels"]],
    },sort_keys=True))
if __name__=="__main__":main()
