"""Focused ORIGINAL GIFT 4505 botanical flower form/joint architecture source TRAINING block-out.

This is the exact SAME strict all13 original source TRAINING exclusion code
and 4,505 island denominator but focused pre-declared three endpoints to
secure the main floral regionality test before full all13 work finishes:
floral_form, true SAME-species colour×symmetry, true SAME-species colour×form.
Both provenance scopes and all FOUR original GSHHG geographic regimes;
three donor-training deletion modes remain. Source SHA-pinned and anonymous.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_global_geo_source_training_holdout import (
    analyse,sha,MODES
)

FOCUS={"floral_form","flower_primary_color__floral_symmetry",
       "flower_primary_color__floral_form"}
PINS=dict(
    traits="9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
    flora="6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
    assign="8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44")


def main():
    ap=argparse.ArgumentParser()
    for k in ("original-root","geography","output"):
        ap.add_argument("--"+k,type=Path,required=True)
    q=ap.parse_args()
    files=dict(
        traits=q.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        flora=q.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        assign=q.original_root/"fixed/source/result/island_source_assignments.csv.gz",
    )
    for k,f in files.items():
        if sha(f)!=PINS[k]:
            raise ValueError("Frozen original full 4505 botanical source changed: "+k)
    states=strict_traits(pd.read_csv(files["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(files["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    assign=pd.read_csv(files["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    geo=pd.read_csv(q.geography,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    if geo.island_id.nunique()!=8264:raise ValueError("Original GSHHG geography changed")
    x=analyse(states,flora,geo,assign,nboot=499,
              endpoint_filter=FOCUS)
    assert x["n_original_botanical_GIFT_island_floras"]==4505
    assert x["n_actual_trait_endpoint_families_tested"]==3
    assert len(x["original_source_training_eligibility_and_scoring"])==72
    output=json.dumps(x,sort_keys=True,indent=2,allow_nan=False)
    for field in ('"island_id"','"accepted_species"','"entity_ID"',
                  '"source_cluster"','"spatial_block"','"island_latitude"',
                  '"island_longitude"','"source_genus"'):
        if field in output:
            raise ValueError("Private original botanical source name leaked")
    q.output.parent.mkdir(parents=True,exist_ok=True)
    q.output.write_text(output+"\n",encoding="utf-8")
    print("ORIGINAL_4505_FOCUSED_SOURCE_TRAINING_EXCLUSION="+json.dumps({
        "status":x["status"],"global_botanical_islands":4505,
        "source_modes":list(MODES),"trait_endpoints":sorted(FOCUS),
        "n_source_reported_cells":len(x["original_source_training_eligibility_and_scoring"]),
        "n_qualified_both_group_modes":sum(
            r.get("status")=="ORIGINAL_STRICT_TRAINING_EXCLUSION_BOTH_GROUPS_QUALIFIED"
            for r in x["original_source_training_scored_results"]),
        "all_source_training_same_group_botanical_names_excluded":True,
    },sort_keys=True))
    for z in x["original_source_training_scored_results"]:
        if z["original_scope"]!="all_recorded" or z[
            "original_source_region"] not in (
                "northern_midlatitude","tropical","southern_extratropical"):
            continue
        print("ORIGINAL_4505_FOCUSED_SPATIAL_HOLDOUT="+json.dumps({
            k:z.get(k) for k in (
                "original_source_region","original_trait_endpoint",
                "source_training_exclusion_mode","status",
                "n_original_true_direct_n5_eligible_focal_islands",
                "n_original_original_GIFT_islands_with_strict_common_focal_scores",
                "island_equal_strict_local_minus_external_Brier",
                "island_equal_same_focal_baseline_local_minus_external_Brier",
                "physical_10degree_block","first_candidate_source_group")
        },sort_keys=True))
if __name__=="__main__":
    main()
