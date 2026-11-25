"""Independent source-safe GIFT target species/genus-unequal-repetition robustness audit.

Reads *only* anonymous original source artifact JSON/ZIP, never protected
accepted species names, source entities or coordinates.

The original data reuse the same taxon-level Direct trait state on multiple
islands. Reweighting their paired Brier/log losses by distinct target accepted
species, then by target genus, is essential SENSITIVITY to pseudoreplication,
not equivalent to original 10°-block-equal or candidate-source-equal inference.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import zipfile

PARENT=4505
ENDPOINTS=frozenset((
    "flower_primary_color","floral_form","floral_symmetry",
    "flower_size_class","tube_depth_class","inflorescence_display",
    "self_incompatibility","autonomous_selfing_capacity","mating_system",
    "cleistogamy","flower_primary_color__floral_form",
    "flower_primary_color__floral_symmetry","floral_form__floral_symmetry"
))
REGIONS=("northern_midlatitude","tropical","southern_extratropical")
SCOPES=("all_recorded","original_native")
BASELINES=("regional_blind","external_blind","regional_genus_prior")
TARGET_UNITS=("species_equal","genus_equal")
STATUS="ORIGINAL_GIFT_STAGE2_TARGET_ACCEPTED_SPECIES_AND_GENUS_CLUSTER_SENSITIVITY"
PRIVATE_KEYS=frozenset(("island_id","accepted_species","entity_ID",
                        "genus_name","source_genus","spatial_block","first_source",
                        "latitude","longitude","island_latitude","island_longitude"))


def _no_private_ids(value):
    if isinstance(value,dict):
        if PRIVATE_KEYS.intersection(value):
            raise ValueError("Protected original botanical plant/source identity key")
        for v in value.values():
            _no_private_ids(v)
    elif isinstance(value,list):
        for v in value:
            _no_private_ids(v)


def read_original(path:Path, *, sha256:str|None=None) -> dict:
    b=path.read_bytes()
    if sha256 is not None and hashlib.sha256(b).hexdigest()!=sha256:
        raise ValueError("SHA256 of original ZIP/JSON differs")
    if path.suffix.lower()==".zip":
        with zipfile.ZipFile(path) as archive:
            files=[x for x in archive.namelist() if x.endswith(".json")]
            if len(files)!=1:
                raise ValueError("Expected exactly one anonymous GIFT JSON artifact")
            obj=json.loads(archive.read(files[0]))
    else:
        obj=json.loads(b)
    if not isinstance(obj,dict):
        raise ValueError("Expected original GIFT JSON object")
    _no_private_ids(obj)
    return obj


def signed_pair(row:dict, unit:str, baseline:str) -> str:
    if unit not in TARGET_UNITS or baseline not in BASELINES:
        raise ValueError("Unknown source target unit or prior")
    result=row["target_species_genus_reweighting"][unit]
    if not result.get("supported_for_inference"):
        return "unsupported"
    cis=[]
    for loss in ("Brier","logloss"):
        key=f"environment_source_nearest_gain_vs_{baseline}_{loss}"
        v=result.get(key,{}).get("bootstrap_95")
        if (not isinstance(v,list) or len(v)!=2
                or not all(isinstance(x,(float,int)) and math.isfinite(x) for x in v)
                or v[0]>v[1]):
            raise ValueError("Missing original taxon-cluster proper-score signed interval")
        cis.append(v)
    if all(x[0]>0 for x in cis):
        return "positive"
    if all(x[1]<0 for x in cis):
        return "negative"
    return "ambiguous"


def audit(original:dict, prior:dict|None=None) -> dict:
    if (original.get("original_parent_GIFT_botanical_floras")!=PARENT
            or original.get("status")!=STATUS
            or original.get("original_distinct_Direct_endpoints_total")!=13
            or set(original.get("frozen_pilot_endpoints",[]))!=ENDPOINTS
            or len(original.get("frozen_pilot_endpoints",[]))!=len(ENDPOINTS)
            or original.get("unique_target_accepted_species_repetitions_not_independent") is not True
            or original.get("species_equal_and_genus_equal_are_new_weighting_estimands") is not True
            or original.get("source_island_based_original_proper_scores_unchanged") is not True):
        raise ValueError("Original 4505-island all13 source study provenance changed")
    rows=original.get("regional_source_heldout_results")
    if not isinstance(rows,list) or len(rows)!=78:
        raise ValueError("Incomplete 13×3×2 original source panel output")
    pairs={}
    if prior is not None:
        if prior.get("original_parent_GIFT_botanical_floras")!=PARENT:
            raise ValueError("Previous GIFT population is not original 4505 islands")
        prior_rows=prior.get("regional_source_heldout_results",[])
        if len(prior_rows)!=78:
            raise ValueError("Historical prior original source must contain 78 panels")
        pairs={(r["scope"],r["region"],r["trait_endpoint"]):r for r in prior_rows}
        if len(pairs)!=78:
            raise ValueError("Historical original region/trait statuses duplicated")
    seen=set()
    found=[]
    counts=Counter()
    for r in rows:
        key=r["scope"],r["region"],r["trait_endpoint"]
        if (key in seen or key[0] not in SCOPES or key[1] not in REGIONS
                or key[2] not in ENDPOINTS):
            raise ValueError("Invalid or duplicate original source region×status×trait")
        seen.add(key)
        if (r.get("original_parent_botanical_floras")!=PARENT
                or r.get("original_full_accepted_name_group_ban_on_both_models") is not True
                or r.get("original_focal_genus_oracle_NOT_used_in_predictor") is not True):
            raise ValueError("Original botanical parent/donor embargo or no-genus-oracle lost")
        v=r.get("target_species_genus_reweighting")
        if not isinstance(v,dict) or v.get("target_accepted_names_absent_from_export") is not True:
            raise ValueError("Taxon name original source-private counter missing")
        sn,gn,occ=(v.get("n_distinct_target_species"),
                    v.get("n_distinct_target_genera"),
                    v.get("n_species_island_occurrences"))
        # Original supported=0 cohort exports no occurrence ratio because
        # there were no eligible Direct carriers; this is a genuine missing
        # analysis, NOT a zero ecological effect.
        if (sn,gn)==(0,0) and r.get("n_original_target_islands_scored")==0:
            if occ is not None and occ!=0:
                raise ValueError("Unsupported original cohort has invented occurrences")
            occ=0
        if (not all(isinstance(z,int) and z>=0 for z in (sn,gn,occ))
                or not 0<=gn<=sn<=occ):
            raise ValueError("Impossible protected target species/genus incidence counts")
        if sn and not math.isclose(v["mean_island_occurrences_per_distinct_species"],
                                   occ/sn,rel_tol=0,abs_tol=1e-10):
            raise ValueError("Species occurrence duplication ratio inconsistent")
        cohort=r["n_original_target_islands_scored"]
        if (not isinstance(cohort,int) or cohort<0
                or (cohort>0 and (sn<5 or occ<cohort*5))
                or (cohort==0 and (sn!=0 or occ!=0))):
            raise ValueError("Target Direct original cohort count or incidence inconsistency")
        doc=dict(scope=key[0],region=key[1],trait_endpoint=key[2],
                 n_original_target_islands=cohort,
                 n_unique_direct_target_species=sn,
                 n_unique_direct_target_genera=gn,
                 n_direct_species_island_occurrences=occ,
                 mean_appearances_per_unique_target_species=(
                     v.get("mean_island_occurrences_per_distinct_species")),
                 original_geographic_group_support=r["support_both_independent_source_groups"])
        for unit in TARGET_UNITS:
            frame=v.get(unit,{})
            if frame.get("supported_for_inference") is True:
                if (frame.get("n_original_group_units") != (sn if unit=="species_equal" else gn)
                        or frame.get("resampling_unit") !=
                           ("distinct_accepted_species" if unit=="species_equal"
                            else "distinct_accepted_genus")):
                    raise ValueError("Target species/genus independent source grouping count mismatch")
            for baseline in BASELINES:
                outcome=signed_pair(r,unit,baseline)
                stored_flag=r.get(
                    f"{unit}_positive_local_external_pooled",{}).get(baseline)
                if stored_flag is not (outcome=="positive"):
                    raise ValueError("Source stored unique taxon positive gate inconsistent")
                doc[f"{unit}_gain_over_{baseline}"]=outcome
                counts[(key[0],key[1],unit,baseline,outcome)]+=1
        if prior is not None:
            b=pairs.get(key)
            if b is None or b["n_original_target_islands_scored"]!=cohort:
                raise ValueError("Original island source cohort denominator drift")
            for grp in ("physical_block","first_source_group"):
                for baseline in BASELINES:
                    for loss in ("Brier","logloss"):
                        label=f"environment_source_nearest_gain_vs_{baseline}_{loss}"
                        old=b.get(grp,{}).get(label)
                        now=r.get(grp,{}).get(label)
                        if (old is None)!=(now is None):
                            raise ValueError("Original proper-score report missing")
                        if old is not None:
                            if not math.isclose(old["group_equal_mean"],
                                                now["group_equal_mean"],
                                                rel_tol=0,abs_tol=1e-10):
                                raise ValueError("Original source proper scores changed")
        found.append(doc)
    if len(seen)!=78:
        raise ValueError("Original 78 panels not fully covered")
    return dict(original_parent=PARENT,original_panels=len(found),
                n_source_status_direct_traits=len(ENDPOINTS),
                original_island_equal_loss_forecasts_unchanged=prior is not None,
                source_verdict_counts={"|".join(k):v for k,v in sorted(counts.items())},
                all_anonymous_panels=sorted(
                    found,key=lambda x:(x["scope"],x["region"],x["trait_endpoint"])),
                limitations="Unique species and genus weighting are new estimands; no source identities, real within-population evolution, full phylogenetic independence, effective pollinators or independent mainland ancestors.")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("source_original_json_or_zip",type=Path)
    ap.add_argument("--source-sha256")
    ap.add_argument("--historical-original-category",type=Path)
    ap.add_argument("--output",type=Path)
    args=ap.parse_args()
    out=audit(read_original(args.source_original_json_or_zip,sha256=args.source_sha256),
              read_original(args.historical_original_category)
              if args.historical_original_category else None)
    raw=json.dumps(out,sort_keys=True,indent=2)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(raw,encoding="utf-8")
    print(raw)


if __name__=="__main__":
    main()
