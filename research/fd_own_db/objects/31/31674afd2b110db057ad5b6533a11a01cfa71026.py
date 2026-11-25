"""Verify complete anonymous original GIFT equal-N local/external trait-source sensitivity.

Input ONLY a source-verified GitHub Actions original anonymized ZIP/JSON, never
the protected 4,505-island source botanical names. Recompute independently the
very strict ALL 3 outcome-blind hash-seed × BOTH proper Brier/log losses ×
BOTH physical-block/first-source-group signed interval gate. All 78
13 genuine Direct endpoint × 3 macroregion × 2 botanical-status cells
must remain present, including unsupported cells and negative evidence.

This is source training sample-size sensitivity; the three hash samples may
overlap and are NOT independent ecological replication. It does NOT isolate
pollinator effects, island adaptation, or verified mainland ancestry.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import zipfile

REGIONS = ("northern_midlatitude", "tropical", "southern_extratropical")
SCOPES = ("all_recorded", "original_native")
GROUPS = ("physical_block", "first_source_group")
SEEDS = (0, 1, 2)
PARENT = 4505
DIRECT_ENDPOINTS = frozenset((
    "flower_primary_color", "floral_form", "floral_symmetry",
    "flower_size_class", "tube_depth_class", "inflorescence_display",
    "self_incompatibility", "autonomous_selfing_capacity",
    "mating_system", "cleistogamy",
    "flower_primary_color__floral_form",
    "flower_primary_color__floral_symmetry",
    "floral_form__floral_symmetry",
))
N_TRAITS = len(DIRECT_ENDPOINTS)
N_ROWS = N_TRAITS * len(REGIONS) * len(SCOPES)
STUDY_STATUS = "ORIGINAL_GIFT_STAGE2_BALANCED_SOURCE_LOCAL_EXTERNAL_BLIND_SENSITIVITY"
SENSITIVE_KEYS = frozenset({
    "island_id", "accepted_species", "entity_ID", "genus_name",
    "spatial_block", "first_source", "source_genus",
    "latitude", "longitude", "island_latitude", "island_longitude",
})


def ensure_anon(obj):
    """Fail closed if actual source species/island/source-group identifiers leak."""
    if isinstance(obj, dict):
        if SENSITIVE_KEYS.intersection(obj):
            raise ValueError("Protected GIFT plant/island/source identifying key")
        for v in obj.values():
            ensure_anon(v)
    elif isinstance(obj, list):
        for v in obj:
            ensure_anon(v)


def read_source(path: Path, expected_digest: str | None = None) -> dict:
    if expected_digest is not None and hashlib.sha256(path.read_bytes()).hexdigest() != expected_digest:
        raise ValueError("Original source artifact digest mismatch")
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as archive:
            names=[x for x in archive.namelist() if x.lower().endswith(".json")]
            if len(names)!=1:
                raise ValueError("Expected exactly one anonymous original GIFT result")
            data=json.loads(archive.read(names[0]))
    else:
        data=json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data,dict):
        raise ValueError("Expected anonymous source dict root")
    ensure_anon(data)
    return data


def score_gate(row: dict, *, beneficiary: str) -> str:
    """Both proper scores × both geographical clusters × all three hash seeds."""
    if beneficiary not in ("external","local"):
        raise ValueError("Unregistered source reference direction")
    if not row.get("support_both_independent_source_groups"):
        return "unsupported"
    ci=[]
    for seed in SEEDS:
        for group in GROUPS:
            frame=row.get(group,{})
            if frame.get("supported_for_inference") is not True:
                raise ValueError("Source geographical group qualification disagreement")
            for loss in ("Brier","logloss"):
                if beneficiary=="external":
                    key=f"balanced_external_seed{seed}_gain_vs_balanced_local_seed{seed}_{loss}"
                else:
                    key=f"balanced_local_seed{seed}_gain_vs_balanced_external_seed{seed}_{loss}"
                interval=frame.get(key,{}).get("bootstrap_95")
                if (not isinstance(interval,list) or len(interval)!=2
                        or not all(isinstance(v,(float,int)) and math.isfinite(v) for v in interval)
                        or interval[0]>interval[1]):
                    raise ValueError("Original equal-N source interval missing or malformed: "+key)
                ci.append(interval)
    if all(z[0]>0 for z in ci):
        return "positive"
    if all(z[1]<0 for z in ci):
        return "negative"
    return "ambiguous"


def audit(data: dict) -> dict:
    if (data.get("status")!=STUDY_STATUS
            or data.get("original_parent_GIFT_botanical_floras")!=PARENT
            or data.get("original_distinct_Direct_endpoints_total")!=N_TRAITS):
        raise ValueError("Original source identity, parent or endpoint family changed")
    if (data.get("matched_source_direct_name_sample_seeds")!=list(SEEDS)
            or data.get("matched_source_direct_name_per_side_cap")!=2000
            or data.get("original_focal_trait_cohorts_and_nonoracle_predictions_unchanged") is not True
            or data.get("direct_local_vs_external_blind_losses_bootstrapped_on_identical_source_carriers") is not True):
        raise ValueError("Original equal-N species donor sampling protocol drift")
    endpoints=data.get("frozen_pilot_endpoints")
    rows=data.get("regional_source_heldout_results")
    if (not isinstance(endpoints,list) or len(endpoints)!=N_TRAITS
            or len(set(endpoints))!=N_TRAITS or set(endpoints)!=DIRECT_ENDPOINTS
            or not isinstance(rows,list) or len(rows)!=N_ROWS):
        raise ValueError("Missing required complete original 13×3×2 source panel table")
    seen=set()
    result=[]
    count=Counter()
    for row in rows:
        scope,region,endpoint=row["scope"],row["region"],row["trait_endpoint"]
        key=scope,region,endpoint
        if (scope not in SCOPES or region not in REGIONS or endpoint not in endpoints
                or key in seen):
            raise ValueError("Original source scope/region/trait key absent or duplicated")
        seen.add(key)
        if (row.get("original_parent_botanical_floras")!=PARENT
                or row.get("original_full_accepted_name_group_ban_on_both_models") is not True
                or row.get("original_focal_genus_oracle_NOT_used_in_predictor") is not True):
            raise ValueError("Botanical parent or no-oracle name donor embargo lost")
        n=row.get("n_original_target_islands_scored")
        if not isinstance(n,int) or n<0:
            raise ValueError("Original target botanical island counts invalid")
        ext=score_gate(row,beneficiary="external")
        local=score_gate(row,beneficiary="local")
        if (row.get("balanced_external_better_all_three_source_samples_both_losses_groups")
                is not (ext=="positive")
                or row.get("balanced_local_better_all_three_source_samples_both_losses_groups")
                is not (local=="positive")):
            raise ValueError("Stored equal-N source verdict disagrees with independent 12-interval gate")
        if ext=="positive" and local=="positive":
            raise ValueError("Both signed local and external reference improvements supported")
        donor_n=row.get("balanced_direct_source_names_per_side_mean")
        donor_ratio=row.get("external_local_original_source_n_ratio_mean")
        if row.get("support_both_independent_source_groups"):
            if (not isinstance(donor_n,(float,int)) or not math.isfinite(donor_n)
                    or donor_n<20 or donor_n>2000
                    or not isinstance(donor_ratio,(float,int)) or not math.isfinite(donor_ratio)
                    or donor_ratio<=0):
                raise ValueError("Original equal source-N/ratio support metadata missing")
        doc=dict(region=region,scope=scope,trait_endpoint=endpoint,
                 n_original_scored_islands=n,
                 both_geographical_source_groups_eligible=bool(
                     row.get("support_both_independent_source_groups")),
                 mean_balanced_distinct_source_names_each_side=donor_n,
                 mean_external_to_local_source_names_ratio=donor_ratio,
                 balanced_external_source_advantage=ext,
                 balanced_local_source_advantage=local)
        result.append(doc)
        count[(scope,region,"external_advantage",ext)]+=1
        count[(scope,region,"local_advantage",local)]+=1
    if len(seen)!=N_ROWS:
        raise ValueError("Incomplete original botanical source grid")
    return {
        "source_status":STUDY_STATUS,
        "original_parent_GIFT_islands":PARENT,
        "n_endpoint_trait_family":N_TRAITS,
        "n_original_region_status_panels":N_ROWS,
        "n_source_hash_samples":len(SEEDS),
        "n_geographical_groupings":len(GROUPS),
        "source_verdict_counts":{"|".join(map(str,k)):v for k,v in sorted(count.items())},
        "panels":sorted(result,key=lambda x:(x["scope"],x["region"],x["trait_endpoint"])),
        "nonclaims":"Not three biological replications, not observed ancestral mainland trait evolution, island adaptation or confirmed local pollinator selection; exploratory nonmultiplicity-adjusted intervals.",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("source_zip_or_json",type=Path)
    p.add_argument("--sha256")
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    r=audit(read_source(a.source_zip_or_json,a.sha256))
    text=json.dumps(r,indent=2,sort_keys=True,ensure_ascii=False)+"\n"
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(text,encoding="utf-8")
    print(text)


if __name__=="__main__":
    main()
