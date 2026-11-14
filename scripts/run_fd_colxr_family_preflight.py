"""Original source-pinned CoL XR exact plant family candidate availability.

This is *not* a fitted family/phylogenetic floral selection effect.
The CoL XR alpha2 source is reviewed as candidate, not independently
accepted taxonomy. It can support a rigorous labeled sensitivity only.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import pandas as pd
from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from fd.colxr_family_opportunity import (
    strict_candidate_family_map,family_source_support
)
from fd.discrete_functional_volume import TRAITS

PINNED_CROSSWALK_SHA="079d9c27b0ca3aa51f0c7310f43ad4901b07a319d6063eab5a7d62b291bdd9e3"


def _sha256(path:Path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for s in iter(lambda:f.read(1<<20),b""):
            h.update(s)
    return h.hexdigest()


def _original_direct(ledger):
    t=ledger.loc[
       ledger.trait_name.isin(TRAITS)&
       ledger.evidence_scope.eq("species_direct")&
       ledger.evidence_quality.isin(("high","medium"))&
       ledger.normalized_value.notna()&
       ~ledger.normalized_value.astype(str).str.contains("|",regex=False)
    ]
    if t.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Conflicting directly measured source evidence")
    x=t.pivot(index="accepted_species",columns="trait_name",
        values="normalized_value").dropna(subset=list(TRAITS))[list(TRAITS)]
    if len(x)!=2234:
        raise ValueError("Original 2234 species direct flower phenotype cohort drifted")
    return x


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--original-root",type=Path,required=True)
    a.add_argument("--colxr-crosswalk",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    z=a.parse_args()
    if _sha256(z.colxr_crosswalk)!=PINNED_CROSSWALK_SHA:
        raise ValueError("Frozen CoL XR 2026 alpha2 crosswalk checksum mismatch")
    root=z.original_root
    ledger=pd.read_csv(root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        usecols=["accepted_species","trait_name","normalized_value",
                 "evidence_scope","evidence_quality"])
    flora=pd.read_csv(root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        usecols=["island_id","accepted_species","origin_status"])
    eligible,gate=strict_direct_928_flora(ledger,flora)
    if gate["n_original_eligible_islands"]!=928:
        raise ValueError("Original 928 island locked roster missing")
    direct=_original_direct(ledger)
    colxr=pd.read_csv(z.colxr_crosswalk,dtype=str)
    mapping,quality=strict_candidate_family_map(colxr,set(direct.index))
    geography=pd.read_csv(root/"fixed/isolation/results/purpose_shortest_island_data.csv",
        usecols=["island_id","analysis_regime","spatial_block"])
    source=pd.read_csv(root/"fixed/source/result/island_source_assignments.csv.gz",
        usecols=["island_id","source_mode","source_rank","entity_ID"])
    gift=pd.read_csv(root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        usecols=["entity_ID","work_species"])
    score=family_source_support(direct,eligible,geography,source,gift,mapping)
    result={
        "status":"COLXR_2026_ACCEPTED_EXACT_CANDIDATE_FAMILY_928_OPPORTUNITY_ONLY",
        "original_taxonomy":"COL_XR_2026-2.0.0-alpha2-auto-candidate-not-promoted",
        "original_taxonomy_crosswalk_sha256":PINNED_CROSSWALK_SHA,
        "mapping_quality":quality,
        "original_source_support":score,
        "no_original_family_taxon_island_or_source_identity_exported":True,
        "no_fitted_family_effect_or_evolutionary_tree":True,
    }
    text=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    if any('"'+key+'"' in text for key in ("island_id","accepted_species",
        "family","genus","entity_ID","work_species","candidate_family",
        "submitted_name","spatial_block","source_signature")):
        raise ValueError("Rights-protected name or original source label leaked")
    z.output.parent.mkdir(parents=True,exist_ok=True)
    z.output.write_text(text+"\n")
    print(json.dumps({
        "status":result["status"],
        "n_exact_COLXR_family_qualified_species":len(mapping),
        "n_colxr_original_species_name_matches":quality["n_source_colxr_name_matched"],
        "n_family_source_panels":len(score["region_family_support"]),
        "panels":[{
            "scope":x["scope"],"region":x["regime"],
            "eligible":x["n_original_eligible_flora_islands"],
            "family5_shared":x["n_family_5_shared_identity"],
            "family5_diff_genus":x["n_family_5_different_genus_identity"],
            "family15_equal_richness":x["n_family_15_equal_richness"],
            "n_blocks_family5":x["n_distinct_physical_blocks_with_family5"]
        }for x in score["region_family_support"]],
    },sort_keys=True))
if __name__=="__main__":main()
