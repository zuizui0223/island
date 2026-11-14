"""Source-rights-safe original 4505 floral richness raw versus GIFT-reference
same-cohort observational guild comparisons; not 928 flora complete cases.

No three-trait completion threshold, no guessed absence for missing guild
rows, and no false interpretation of qualified non-detection as guild loss.
Source botanical names, exact animal occurrences, island IDs, GIFT source
regions and original channel/pollinator identities stay in original CI.

First compute the full original traitwise and pairwise source phenotype
on 4,505 botanical floras (928 triple-complete sensitivity only), then
join five existing independently source-qualified GBIF channel records
to the SAME island keys. Only anonymous aggregate JSON is exported.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import pandas as pd

from fd.global_partial_trait_flora_fd import original_global_partial_traits
from fd.global_guild_source_same_cohort_audit import source_cohort_and_guild_robustness


def original_sha256(path:Path)->str:
    h=sha256()
    with path.open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):
            h.update(part)
    return h.hexdigest()


def run(root:Path,geology:Path,channel_path:Path,*,seed:int)->dict:
    ledger_path=root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora_path=root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    if original_sha256(ledger_path)!=(
        "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745"):
        raise ValueError("Frozen original botanical Direct trait source changed")
    if original_sha256(flora_path)!=(
        "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a"):
        raise ValueError("Original 4,505 source floras changed")
    ledger=pd.read_csv(ledger_path,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"])
    flora=pd.read_csv(flora_path,dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status"])
    if flora.island_id.nunique()!=4505:
        raise ValueError("Original botanical flora population 4505 drifted")
    geo=pd.read_csv(geology,dtype=str,usecols=[
        "island_id","geological_origin","analysis_regime","spatial_block"])
    if (geo.island_id.nunique()<8200 or
        geo.geological_origin.eq("continental").sum()!=176 or
        geo.geological_origin.eq("oceanic").sum()!=258):
        raise ValueError("Original independently published geological source drift")
    assignments=pd.read_csv(
        root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(
        root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    channels=pd.read_csv(channel_path,dtype=str,
        usecols=["island_id","channel_id","channel_state","source_state"])
    original_main,private_rows=original_global_partial_traits(
        ledger,flora,geo,assignments,mainland,
        seed=seed,bootstrap_draws=999,return_private_island_rows=True)
    if (original_main["n_original_global_distinct_island_floras"]!=4505 or
        original_main["n_original_three_trait_complete30_islands_retained_for_SENSITIVITY"]!=928 or
        len(private_rows)!=117130):
        raise ValueError("4505 / 928 historic botanical source frame changed")
    result=source_cohort_and_guild_robustness(
        private_rows,channels,seed=seed,draws=999)
    if result["n_original_individual_botanical_islands"]!=4505:
        raise AssertionError("GBIF source join lost global island plant floras")
    if len(private_rows)!=117130:
        raise AssertionError("Original 10+3 traits lost original island×native status resolution")
    result["original_numeric_seed"]=seed
    result["original_evidence_SHA256"]=dict(
        GBIF_guild_source_qualifier_sha256=original_sha256(channel_path),
        original_species_direct_ledger_sha256=original_sha256(ledger_path),
        original_GIFT_flora_sha256=original_sha256(flora_path))
    result["global_10_trait_coverage_preflight"]=dict(
        n_global_flora_islands=4505,
        n_historical_triple_complete30_sensitivity_only=928,
        original_direct_species_by_trait=original_main["n_original_direct_plant_species_per_trait"],
        original_botanical_geology_intersection_islands=original_main[
            "n_original_distinct_islands_with_geological_origin_typed"],
    )
    result["all_original_source_taxa_animal_island_and_GIFT_IDs_stayed_inside_CI"]=True
    return result


def main()->None:
    a=argparse.ArgumentParser()
    a.add_argument("--original-root",type=Path,required=True)
    a.add_argument("--geology-crosswalk",type=Path,required=True)
    a.add_argument("--channel-evidence-csv",type=Path,required=True)
    a.add_argument("--seed",type=int,default=20261114)
    a.add_argument("--output",type=Path,required=True)
    args=a.parse_args()
    summary=run(args.original_root,args.geology_crosswalk,
                args.channel_evidence_csv,seed=args.seed)
    body=json.dumps(summary,ensure_ascii=False,sort_keys=True,indent=2,
                    allow_nan=False)
    forbidden=(
        "island_id","accepted_species","entity_ID","work_species",
        "spatial_block","source_signature","source_genus",
        "source_family","source_region_id","longitude","latitude",
        "taxon_id","GBIF_id","source_observation_id","original_photo_id")
    if any('"'+name+'"' in body for name in forbidden):
        raise ValueError("GBIF plant/pollinator restricted original evidence leaked")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+"\n",encoding="utf-8")
    print(json.dumps(dict(
        status=summary["status"],
        n_global_botanical_islands=summary["n_original_individual_botanical_islands"],
        n_guild_source_matched_analytical_panels=summary["n_source_cohort_contrast_panels"],
        n_source_same_island_panels_with_independent_support=summary[
            "n_source_same_island_panels_with_support"],
        n_eight_slot_source_block_mixing_hypotheses=summary[
            "n_original_primary_eight_slot_testable_with_mixed_blocks"],
        historical_guild_decline_inferred=False,
    ),sort_keys=True))


if __name__=="__main__":
    main()
