"""Source-rights protected original 928-island equal genus-RICHNESS trial.

Unlike original genus-identity matched cohort, one Direct flower-complete
plant species per distinct genus gives the same genus count on island and
GIFT source WITHOUT requiring genus identities to coincide. No original
plant names, genus names, GIFT entities, locations or islands are exported.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import pandas as pd

from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from fd.discrete_functional_volume import TRAITS
from fd.genus_richness_rare_volume import original_genus_richness_normalization


def _complete_direct(ledger:pd.DataFrame)->pd.DataFrame:
    req={"accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"}
    if not req.issubset(ledger):
        raise ValueError("Missing original Direct trait evidence contract")
    x=ledger.loc[
        ledger.trait_name.isin(TRAITS)&
        ledger.evidence_scope.eq("species_direct")&
        ledger.evidence_quality.isin(("high","medium"))&
        ledger.normalized_value.notna()&
        ~ledger.normalized_value.astype(str).str.contains("|",regex=False)
    ]
    if x.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Original source evidence contains duplicate/conflicting direct trait")
    t=x.pivot(index="accepted_species",columns="trait_name",
         values="normalized_value").dropna(subset=list(TRAITS))[list(TRAITS)]
    if len(t)!=2234: raise ValueError("Original global Direct 2234 species changed")
    return t


def _json(v):
    if isinstance(v,dict):return {str(k):_json(x) for k,x in v.items()}
    if isinstance(v,list):return [_json(x) for x in v]
    if isinstance(v,float) and not math.isfinite(v):return None
    return v


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--progressive-root",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--repeats",type=int,default=40)
    p.add_argument("--reference-null",type=int,default=399)
    a=p.parse_args()
    base=a.progressive_root
    direct=pd.read_csv(base/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
         usecols=["accepted_species","trait_name","normalized_value",
                  "evidence_scope","evidence_quality"])
    flora=pd.read_csv(base/"fixed/canonical/input/chapter1_status_flora.csv.gz",
         usecols=["island_id","accepted_species","origin_status"])
    original,gate=strict_direct_928_flora(direct,flora)
    if gate["n_original_eligible_islands"]!=928:
        raise ValueError("Original protected 928-island universe changed")
    triple=_complete_direct(direct)
    geo=pd.read_csv(base/"fixed/isolation/results/purpose_shortest_island_data.csv",
         usecols=["island_id","analysis_regime","spatial_block"])
    source=pd.read_csv(base/"fixed/source/result/island_source_assignments.csv.gz",
         usecols=["island_id","source_mode","source_rank","entity_ID"])
    gift=pd.read_csv(base/"fixed/source/source/gift_native_mainland_flora.csv.gz",
         usecols=["entity_ID","work_species"])
    r=original_genus_richness_normalization(
         triple,original,geo,source,gift,
         repeats=a.repeats,reference_draws=a.reference_null,
         bootstrap=1999,seed=20261113)
    assert r["n_original_eligible_islands"]==928
    assert r["original_complete_Direct_species"]==2234
    assert r["original_global_profile_count"]==100
    assert r["design"]["one_different_accepted_plant_species_per_genus"]
    assert r["design"]["same_genus_number_but_genus_IDENTITY_NOT_matched"]
    assert r["no_original_source_genus_taxon_island_or_GIFT_name_exported"]
    assert r["no_actual_pollinator_loss_or_fitness_inferred"]
    r=_json(r)
    txt=json.dumps(r,indent=2,sort_keys=True,allow_nan=False)
    forbidden=('island_id','accepted_species','genus','family',
              'entity_ID','work_species','spatial_block','source_signature')
    if any('"'+x+'"' in txt for x in forbidden):
        raise ValueError("Original protected plant/genera/source IDs leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(txt+"\n")
    print(json.dumps({
        "status":r["status"],
        "panels":[{
            "scope":x["scope"],"region":x["regime"],
            "n_distinct_genera":x["n_distinct_genera_per_side"],
            "n_islands":x["n_supported_islands"],
            "n_blocks":x["n_spatial_blocks"],
            "n_nearest_GIFT_sources":x["n_first_nearest_GIFT_source_signatures"],
            "reported":x["regional_comparison_reportable"],
        } for x in r["equal_genus_richness_panels"]]
    },sort_keys=True))
if __name__=="__main__":main()
