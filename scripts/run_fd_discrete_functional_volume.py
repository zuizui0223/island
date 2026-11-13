"""Original-source protected 928-island categorical finite-volume replay.

Only original Direct high/medium flower colour, form and symmetry are used.
GIFT original five-nearest candidate native mainland-source lists are
comparators, NOT proven ancestral populations. Upload only aggregate JSON.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from fd.discrete_functional_volume import (
    TRAITS,original_categorical_flora_volume
)


def _source_direct_triples(direct:pd.DataFrame)->pd.DataFrame:
    needed={"accepted_species","trait_name","normalized_value",
            "evidence_scope","evidence_quality"}
    if not needed.issubset(direct):
        raise ValueError("Original source Direct ledger changed")
    d=direct.loc[
        direct.trait_name.isin(TRAITS) &
        direct.evidence_scope.eq("species_direct") &
        direct.evidence_quality.isin(("high","medium")) &
        direct.normalized_value.notna() &
        ~direct.normalized_value.astype(str).str.contains("|",regex=False)
    ].copy()
    if d.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Original Direct species×trait conflict")
    x=d.pivot(index="accepted_species",columns="trait_name",values="normalized_value")
    if not set(TRAITS).issubset(x):
        raise ValueError("Original Direct morphology 3 trait columns missing")
    x=x.dropna(subset=list(TRAITS))[list(TRAITS)]
    if not 2000<=len(x)<=3000:
        raise ValueError("Original source 2,234 complete Direct phenotype species cohort drifted")
    return x


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--progressive-root",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--repeats",type=int,default=25)
    p.add_argument("--null-shuffles",type=int,default=19)
    p.add_argument("--random-profiles",type=int,default=199)
    a=p.parse_args()
    root=a.progressive_root
    direct=pd.read_csv(root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        usecols=["accepted_species","trait_name","normalized_value",
                 "evidence_scope","evidence_quality"])
    original_flora=pd.read_csv(root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        usecols=["island_id","accepted_species","origin_status"])
    matched_flora, original_gate=strict_direct_928_flora(direct,original_flora)
    if original_gate["n_original_eligible_islands"]!=928:
        raise ValueError("Locked original 928-island Direct macro-flora check failed")
    triple=_source_direct_triples(direct)
    if not len(triple)>=original_gate["n_original_complete_plant_species"]:
        raise ValueError("Original Direct species reference smaller than eligible-island cohort")
    if len(triple)!=2234:
        raise ValueError("Global triple-complete source population must be original 2,234 species")
    geo=pd.read_csv(root/"fixed/isolation/results/purpose_shortest_island_data.csv",
        usecols=["island_id","analysis_regime","spatial_block"])
    source=pd.read_csv(root/"fixed/source/result/island_source_assignments.csv.gz",
        usecols=["island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        usecols=["entity_ID","work_species"])
    result=original_categorical_flora_volume(
        triple,matched_flora,geo,source,mainland,
        n=30,repeats=a.repeats,null_shuffles=a.null_shuffles,
        random_profiles=a.random_profiles,bootstrap=1999)
    assert result["n_928_direct_islands"]==928
    assert result["n_all_recorded_islands_scored"]==928
    assert result["measurement"]["radius0_is_only_unique_observed_profile_richness_over_K"]
    assert result["measurement"]["radius1_interpreted_as_one_trait_neighbourhood_not_continuous_Euclidean_volume"]
    assert result["measurement"]["null_preserves_observed_colour_margins_and_form_symmetry_pairs"]
    assert result["no_implicit_source_species_or_island_taxa_export"]
    assert result["not_a_species_independent_or_phylogenetically_adjusted_null"]
    assert not result["demonstrated_causal_pollinator_loss_or_flower_evolution"]
    text=json.dumps(result,indent=2,sort_keys=True)
    forbidden=('"accepted_species"','"island_id"','"spatial_block"','"entity_ID"',
               '"work_species"','"floral_form"','"floral_symmetry"')
    if any(k in text for k in forbidden):
        raise ValueError("Rights-protected original source identity leaked into aggregate")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(text+"\n")
    print(json.dumps({
        "status":result["status"],
        "n_islands":result["n_all_recorded_islands_scored"],
        "n_reference_species":result["n_original_reference_direct_plant_species"],
        "n_reference_profiles":result["n_original_distinct_triple_profiles"],
        "n_native_islands":result["n_original_strict_native_islands_scored"],
        "by_region":[{
            "region":x["regime"],
            "n":x["n_islands"],
            "source_all":x["n_source_all_supported"],
            "source_disjoint":x["n_source_disjoint_supported"],
            "median_r0":x["median_occupancy_r0"],
            "median_r1":x["median_metric_volume_r1"],
            "median_r2":x["median_metric_volume_r2"],
            "r1_saturated_frac":x["fraction_islands_r1_near_saturation"],
        } for x in result["floral_metric_volume_by_geographic_regime"]
          if x["scope"]=="all_recorded"]
    },sort_keys=True))
if __name__=="__main__":
    main()
