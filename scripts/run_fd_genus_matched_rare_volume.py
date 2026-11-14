"""Protected original 928-island source-genus controlled rare coverage cohort.

This runner measures a NEW, separate 15-species exact-genus sample, not
an after-the-fact n=30 genus-adjustment. The ability to match source genera
is decided only by accepted plant names, not visitor or morphology outcomes.
No original accepted taxon names, GIFT entities, island IDs, genera, or
coordinates appear in the published region-level JSON.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import pandas as pd

from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from fd.discrete_functional_volume import TRAITS
from fd.genus_matched_rare_volume import original_genus_matched_rare_reanalysis


def _direct_complete(direct:pd.DataFrame)->pd.DataFrame:
    x=direct.loc[
        direct.trait_name.isin(TRAITS)&
        direct.evidence_scope.eq("species_direct")&
        direct.evidence_quality.isin(("high","medium"))&
        direct.normalized_value.notna()&
        ~direct.normalized_value.astype(str).str.contains("|",regex=False)
    ]
    if x.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Original Direct species x trait conflicted")
    result=x.pivot(index="accepted_species",columns="trait_name",
                    values="normalized_value").dropna(subset=list(TRAITS))[list(TRAITS)]
    if len(result)!=2234:
        raise ValueError("Original 2234 direct species full reference changed")
    return result


def _json(value):
    if isinstance(value,dict):
        return {str(k):_json(v) for k,v in value.items()}
    if isinstance(value,list):
        return [_json(v) for v in value]
    if isinstance(value,float) and not math.isfinite(value):
        return None
    return value


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
        raise ValueError("Original 928 eligible island macro population mismatch")
    triple=_direct_complete(direct)
    geo=pd.read_csv(base/"fixed/isolation/results/purpose_shortest_island_data.csv",
        usecols=["island_id","analysis_regime","spatial_block"])
    sources=pd.read_csv(base/"fixed/source/result/island_source_assignments.csv.gz",
        usecols=["island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(base/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        usecols=["entity_ID","work_species"])
    output=original_genus_matched_rare_reanalysis(
        triple,original,geo,sources,mainland,
        repeats=a.repeats,reference_draws=a.reference_null,
        bootstrap=1999,seed=20261113)
    assert output["n_original_928_eligible_islands"]==928
    assert output["n_original_complete_Direct_reference_species"]==2234
    assert output["design"]["identical_genus_counts_per_island_source_draw"]
    assert output["design"]["original_exact_local_plant_names_removed_from_all_candidate_sources"]
    assert output["protected_taxon_species_island_genus_or_GIFT_entity_never_exported"]
    assert output["actual_pollinator_loss_or_evolution_not_measured"]
    output=_json(output)
    text=json.dumps(output,indent=2,sort_keys=True,allow_nan=False)
    if any('"'+x+'"' in text for x in [
        "accepted_species","island_id","source_signature","spatial_block",
        "work_species","entity_ID","genus","family","species"]):
        raise ValueError("Original source genus species or island identity exported")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(text+"\n")
    print(json.dumps({
      "status":output["status"],
      "region_support":[{
        "scope":r["scope"],"region":r["regime"],
        "eligible_original_islands":r["n_original_direct_30_species_islands"],
        "genus_matched_15_islands":r["n_with_complete_same_genus_match_15_islands"],
        "genus_matched_30_islands":r["n_with_complete_same_genus_match_30_islands"],
        "spatial_blocks":r["n_original_physical_blocks_match15"],
        "first_GIFT_source_signatures":r["n_first_GIFT_source_signatures_match15"],
        "regional_effect_estimable":r["eligible_to_report_conditional_genus15_comparison"],
      }for r in output["same_genus_source_support_by_scope_region"]],
      "n_regions_with_genus_matched_effects":output["n_genus_matched_source_effect_region_panels_reported"],
    },sort_keys=True))
if __name__=="__main__":
    main()
