"""Anonymous source audit: taxonomy, Direct documentation, congener transfer.

Given a successful protected-source GIFT 4,505 botanical parent original
source JSON/ZIP, independently verify all 78 13-Direct-endpoint×3
macroregion×2 nested botanical-status panels. Its signed residual-vector
pieces must sum to the full predicted-minus-empirical categorical error;
signed projections are algebraic, NOT additive causal variance shares.

Require original isotope-free source cohort and older original no-oracle
floristic category prediction scores unchanged; no accepted botanical name,
genus identifier, source group ID or exact coordinates may leave original CI.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import zipfile

EXPECTED_PARENT=4505
ENDPOINTS=frozenset((
    "flower_primary_color","floral_form","floral_symmetry",
    "flower_size_class","tube_depth_class","inflorescence_display",
    "self_incompatibility","autonomous_selfing_capacity",
    "mating_system","cleistogamy",
    "flower_primary_color__floral_form",
    "flower_primary_color__floral_symmetry",
    "floral_form__floral_symmetry",
))
REGIONS=("northern_midlatitude","tropical","southern_extratropical")
SCOPES=("all_recorded","original_native")
GROUPS=("physical_block","first_source_group")
METRICS=(
    "taxonomic_genus_forecast_residual",
    "direct_trait_documentation_residual",
    "within_genus_documented_trait_residual",
)
STATUS="ORIGINAL_GIFT_STAGE2_TAXONOMY_DIRECT_ANNOTATION_TRAIT_TRANSFER_RESIDUAL"
SENSITIVE=frozenset(("island_id","accepted_species","entity_ID","spatial_block",
                     "first_source","genus_name","source_genus",
                     "latitude","longitude","island_latitude","island_longitude"))


def read_original(path:Path, *, sha256:str|None=None)->dict:
    raw=path.read_bytes()
    if sha256 is not None and hashlib.sha256(raw).hexdigest()!=sha256:
        raise ValueError("Original source ZIP SHA256 digest differs")
    if path.suffix.lower()==".zip":
        with zipfile.ZipFile(path) as z:
            names=[n for n in z.namelist() if n.endswith(".json")]
            if len(names)!=1:
                raise ValueError("Expected one source-anonymous JSON in ZIP")
            d=json.loads(z.read(names[0]))
    else:
        d=json.loads(raw)
    if not isinstance(d,dict):
        raise ValueError("Invalid source-anonymous root")
    _reject_sensitive(d)
    return d


def _reject_sensitive(obj):
    if isinstance(obj,dict):
        if set(obj).intersection(SENSITIVE):
            raise ValueError("Protected botanical source species/genus/island identity leaked")
        for x in obj.values():
            _reject_sensitive(x)
    elif isinstance(obj,list):
        for x in obj:
            _reject_sensitive(x)


def _index(data:dict)->dict:
    if (data.get("original_parent_GIFT_botanical_floras")!=EXPECTED_PARENT
            or set(data.get("frozen_pilot_endpoints",[]))!=ENDPOINTS
            or len(data.get("frozen_pilot_endpoints",[]))!=13):
        raise ValueError("Original 4,505 GIFT parent/13 genuine Direct traits changed")
    rows=data.get("regional_source_heldout_results",[])
    if len(rows)!=78:
        raise ValueError("Original 13×3×2 protected source panels incomplete")
    result={}
    for row in rows:
        key=row.get("scope"),row.get("region"),row.get("trait_endpoint")
        if (key in result or key[0] not in SCOPES or key[1] not in REGIONS
                or key[2] not in ENDPOINTS
                or row.get("original_parent_botanical_floras")!=EXPECTED_PARENT
                or row.get("original_full_accepted_name_group_ban_on_both_models") is not True
                or row.get("original_focal_genus_oracle_NOT_used_in_predictor") is not True):
            raise ValueError("Duplicate, unregistered or unprotected original GIFT source cell")
        result[key]=row
    return result


def verify(data:dict, old:dict|None=None, *, tol:float=1e-10)->dict:
    if (data.get("status")!=STATUS
            or data.get("direct_residual_vector_identity_verified_per_original_target_island") is not True
            or data.get("residual_projection_signed_not_causal_explained_variance") is not True
            or data.get("observed_true_focal_genera_used_only_in_post_prediction_audit") is not True):
        raise ValueError("Original descriptive residual-identity provenance missing")
    rows=_index(data)
    old_rows=_index(old) if old else None
    results=[]
    counts=Counter()
    for key in sorted(rows):
        r=rows[key]
        original=old_rows.get(key) if old_rows is not None else None
        if (r.get("source_residual_projection_is_not_a_causal_selection_partition") is not True
                or r.get("whole_flora_and_direct_carrier_genus_used_only_after_prediction") is not True):
            raise ValueError("Residual is not protected as post-prediction noncausal diagnostic")
        n=r.get("n_original_target_islands_scored")
        if not isinstance(n,int) or n<0:
            raise ValueError("Invalid original accepted target island count")
        if old_rows is not None:
            if original is None or original["n_original_target_islands_scored"]!=n:
                raise ValueError("Original Direct target cohort changed")
            if r["support_both_independent_source_groups"]!=original["support_both_independent_source_groups"]:
                raise ValueError("Original geographic donor support changed")
        out=dict(scope=key[0],region=key[1],endpoint=key[2],
                 n_scored_original_target_islands=n,
                 original_source_geographical_support=bool(r["support_both_independent_source_groups"]))
        for g in GROUPS:
            block=r.get(g,{})
            if old_rows is not None:
                for base in ("regional_blind","external_blind","regional_genus_prior"):
                    for loss in ("Brier","logloss"):
                        label=f"environment_source_nearest_gain_vs_{base}_{loss}"
                        a=block.get(label)
                        b=original.get(g,{}).get(label)
                        if (a is None)!=(b is None):
                            raise ValueError("Previous GIFT category proper-score presence changed")
                        if a is not None and not math.isclose(
                                a["group_equal_mean"],b["group_equal_mean"],
                                rel_tol=0,abs_tol=tol):
                            raise ValueError("Previous GIFT source category proper scores changed")
            if block.get("supported_for_inference"):
                measures={}
                for stem in METRICS:
                    for suffix in ("squared_norm","signed_projection_on_total_error"):
                        label=f"{stem}_{suffix}"
                        m=block.get(label,{}).get("group_equal_mean")
                        if not isinstance(m,(float,int)) or not math.isfinite(m):
                            raise ValueError("Missing original source residual component")
                        measures[label]=float(m)
                    if measures[f"{stem}_squared_norm"] < -tol:
                        raise ValueError("Negative true squared residual magnitude")
                for label in ("observed_minus_predicted_direct_trait_squared_norm",
                              "observed_direct_trait_plugin_Q"):
                    v=block.get(label,{}).get("group_equal_mean")
                    if not isinstance(v,(float,int)) or not math.isfinite(v):
                        raise ValueError("Missing original observed Direct trait Q/total error")
                    measures[label]=float(v)
                total=measures["observed_minus_predicted_direct_trait_squared_norm"]
                projection=sum(measures[f"{stem}_signed_projection_on_total_error"]
                               for stem in METRICS)
                if not math.isclose(total,projection,rel_tol=0,abs_tol=tol):
                    raise ValueError("Original residual vector signed projections do not add")
                if (total < -tol or not 0-tol<=measures["observed_direct_trait_plugin_Q"]<=1+tol):
                    raise ValueError("Original category Q or source residual outside valid interval")
                norms={name:measures[f"{name}_squared_norm"] for name in METRICS}
                largest=max(norms,key=norms.get)
                out[g+"_largest_norm_component"]=largest
                out[g+"_squared_error"]=total
                out[g+"_observed_plugin_Q"]=measures["observed_direct_trait_plugin_Q"]
                for stem in METRICS:
                    out[g+"_"+stem+"_squared_norm"]=norms[stem]
                    out[g+"_"+stem+"_signed_projection"]=measures[
                        f"{stem}_signed_projection_on_total_error"]
                counts[(key[0],key[1],g,"largest_component",largest)]+=1
            else:
                out[g+"_largest_norm_component"]="unsupported"
                counts[(key[0],key[1],g,"support","unsupported")]+=1
        results.append(out)
    return dict(original_GIFT_botanical_parent=EXPECTED_PARENT,
                n_full_Direct_endpoints=13,
                n_anonymous_source_panels=len(results),
                old_original_proper_category_scores_unchanged=old is not None,
                counts={"|".join(k):v for k,v in sorted(counts.items())},
                source_anonymous_panels=results,
                interpretation="Algebraic noncausal components. Squared magnitudes do not add; signed projections do. True focal botanical species/genus used only AFTER no-oracle prediction. Not adaptation, pollinator effect or ancestral mainland data.")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("source_json_or_zip",type=Path)
    p.add_argument("--sha256")
    p.add_argument("--prior_original_json_or_zip",type=Path)
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    source=read_original(args.source_json_or_zip,sha256=args.sha256)
    prior=read_original(args.prior_original_json_or_zip) if args.prior_original_json_or_zip else None
    r=verify(source,prior)
    content=json.dumps(r,sort_keys=True,indent=2,ensure_ascii=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(content,encoding="utf-8")
    print(content)


if __name__=="__main__":
    main()
