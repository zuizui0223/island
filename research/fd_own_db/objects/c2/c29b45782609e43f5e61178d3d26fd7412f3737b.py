"""Independently audit original 4,505-GIFT high/medium Direct trait grading.

Consumes only the protected original source-anonymous JSON/ZIP artifact,
never source species names, botanical genus identities, island IDs or any
geographical coordinates.

The source study compares DISTINCT accepted plant species having high vs
medium species-Direct evidence categories, genus-stratifying permutation
of the evidence grade with within-genus trait-state labels fixed. It DOES
NOT compare repeat measurements of the SAME species and cannot identify
actual biological trait heterogeneity vs documentation bias.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import zipfile

ENDPOINTS=frozenset((
    "flower_primary_color","floral_form","floral_symmetry",
    "flower_size_class","tube_depth_class","inflorescence_display",
    "self_incompatibility","autonomous_selfing_capacity","mating_system",
    "cleistogamy","flower_primary_color__floral_form",
    "flower_primary_color__floral_symmetry","floral_form__floral_symmetry",
))
REGIONS=("northern_midlatitude","tropical","southern_extratropical")
SCOPES=("all_recorded","original_native")
SOURCE_STATUS="ORIGINAL_GIFT4505_SPECIES_DIRECT_QUALITY_WITHIN_GENUS_NULL"
PRIVATE=frozenset((
    "island_id","accepted_species","entity_ID","source_genus",
    "genus_name","genus","spatial_block","first_source",
    "latitude","longitude","island_latitude","island_longitude"
))


def _no_private(obj):
    if isinstance(obj,dict):
        if PRIVATE.intersection(obj):
            raise ValueError("Protected original plant/genus/island identity exported")
        for v in obj.values():
            _no_private(v)
    elif isinstance(obj,list):
        for v in obj:
            _no_private(v)


def read_anon(path:Path,expected_sha256:str|None=None)->dict:
    buf=path.read_bytes()
    if expected_sha256 and hashlib.sha256(buf).hexdigest()!=expected_sha256:
        raise ValueError("Original GIFT source artifact SHA256 digest changed")
    if path.suffix.lower()==".zip":
        with zipfile.ZipFile(path) as archive:
            items=[p for p in archive.namelist() if p.endswith(".json")]
            if len(items)!=1:
                raise ValueError("Exactly one original anonymous JSON expected")
            doc=json.loads(archive.read(items[0]))
    else:
        doc=json.loads(buf)
    if not isinstance(doc,dict):
        raise ValueError("Expected original anonymous JSON root object")
    _no_private(doc)
    return doc


def audit(data:dict)->dict:
    if (data.get("status")!=SOURCE_STATUS
        or data.get("original_parent_GIFT_botanical_island_floras")!=4505
        or data.get("original_parent_geography_islands")!=8264
        or data.get("n_genuine_Direct_endpoints")!=len(ENDPOINTS)
        or data.get("n_region_scope_endpoint_panels")!=78
        or data.get("quality_ontology")!="species_direct_high_vs_medium_once_per_species"
        or data.get("joint_endpoint_high_grade_requires_both_same_species_traits_high") is not True
        or data.get("original_species_name_and_genus_IDs_never_exported") is not True
        or data.get("no_independent_within_species_or_local_island_phenotype_replicates") is not True
        or data.get("BH_q_exploratory_across_all_qualified_78_nested_tests") is not True):
        raise ValueError("Original botanical evidence quality and noncausal source provenance drift")
    raw=data.get("panels")
    if not isinstance(raw,list) or len(raw)!=78:
        raise ValueError("Expected complete 13×3×2 GIFT evidence quality source panels")
    seen=set()
    out=[]
    counts=Counter()
    pvals=[]
    for row in raw:
        key=row.get("scope"),row.get("region"),row.get("trait_endpoint")
        if (key in seen or key[0] not in SCOPES or key[1] not in REGIONS
                or key[2] not in ENDPOINTS
                or row.get("original_GIFT_botanical_parent_island_floras")!=4505
                or row.get("high_vs_medium_are_different_plant_species_not_repeated_measurements") is not True
                or row.get("biological_within_genus_variation_vs_annotation_not_identifiable") is not True):
            raise ValueError("Duplicate/unregistered/noncausal original region/trait/status panel")
        seen.add(key)
        fields=(
            "regional_original_parent_botanical_floras",
            "source_region_original_distinct_botanical_species",
            "n_original_distinct_accepted_species",
            "n_true_species_Direct_coverage",
            "n_direct_quality_high",
            "n_direct_quality_medium",
            "n_source_distinct_genera_with_Direct",
            "n_direct_names_without_certified_genus",
            "n_within_genus_quality_matched_genera",
            "n_different_category_states",
        )
        if not all(isinstance(row.get(name),int) and row[name]>=0 for name in fields):
            raise ValueError("Original botanical evidence source cohort integers malformed")
        n=row["n_original_distinct_accepted_species"]
        scored=row["n_true_species_Direct_coverage"]
        hi=row["n_direct_quality_high"]
        med=row["n_direct_quality_medium"]
        gen=row["n_source_distinct_genera_with_Direct"]
        un=row["n_direct_names_without_certified_genus"]
        matched=row["n_within_genus_quality_matched_genera"]
        if (row["source_region_original_distinct_botanical_species"]!=n
                or not 0<=scored<=n
                or hi+med!=scored
                or not 0<=matched<=gen<=scored
                or not 0<=un<=scored):
            raise ValueError("Original Direct coverage, unique-species and genus denominators inconsistent")
        frac=row.get("high_grade_fraction_among_Direct")
        cov=row.get("original_distinct_species_Direct_coverage_fraction")
        if scored:
            if (frac is None or cov is None or
                abs(frac-hi/scored)>1e-10 or abs(cov-scored/n)>1e-10):
                raise ValueError("Direct high/medium quality coverage proportions malformed")
        elif frac is not None or cov is not None and n==0:
            raise ValueError("Nonexistent original source Direct trait coverage fabricated")
        qualified=row.get("within_genus_permutation_qualified")
        if not isinstance(qualified,bool):
            raise ValueError("Original quality permutation support missing")
        tv=row.get("pooled_high_medium_category_total_variation")
        within=row.get("within_genus_equal_mean_high_medium_total_variation")
        for v in (tv,within):
            if v is not None and (not isinstance(v,(int,float))
                                  or not math.isfinite(v) or v < -1e-10 or v>1+1e-10):
                raise ValueError("Observed quality category TV invalid")
        p=row.get("within_genus_quality_label_permutation_upper_tail_p_exploratory")
        q=row.get("within_genus_quality_label_permutation_BH_q_exploratory")
        obs=row.get("within_genus_observed_minus_permutation_mean_TV")
        mean=row.get("within_genus_null_mean_total_variation")
        if qualified:
            if (hi<30 or med<30 or matched<10 or
                not all(isinstance(z,(int,float)) and math.isfinite(z)
                        for z in (p,q,obs,mean)) or
                not 0<=p<=1 or not 0<=q<=1 or
                abs(obs-(within-mean))>1e-10):
                raise ValueError("Incomplete matched-genus quality permutation original result")
            pvals.append((key,float(p),float(q)))
        elif any(z is not None for z in (p,q,obs,mean)):
            raise ValueError("Original unsupported quality test has invented inference")
        doc={
            "scope":key[0],"region":key[1],"trait_endpoint":key[2],
            "n_original_distinct_species":n,
            "n_direct_covered_unique_species":scored,
            "n_high":hi,"n_medium":med,
            "n_matched_genera_high_medium":matched,
            "pooled_high_medium_TV":tv,
            "within_genus_quality_TV":within,
            "within_genus_permutation_qualified":qualified,
            "genus_conditioned_null_mean_TV":mean,
            "within_genus_observed_minus_null_TV":obs,
            "one_sided_permutation_p":p,
            "nested_family_exploratory_BH_q":q,
        }
        out.append(doc)
        counts[(key[0],key[1],"qualified" if qualified else "unsupported")]+=1
    if len(seen)!=78:
        raise ValueError("Incomplete original source endpoints")
    if data.get("n_permutation_qualified_panels")!=len(pvals):
        raise ValueError("Original 78-cell source permutation qualification count changed")

    # Rederive BH across ALL qualified source panels, never only the
    # subset of nominally positive findings. Underlying nested botanical
    # category/status endpoints share species: q remains exploratory.
    order=sorted(pvals,key=lambda x:x[1])
    m=len(order)
    corrected={}
    bound=1.0
    for i in range(m-1,-1,-1):
        key,p,q=order[i]
        bound=min(bound,p*m/(i+1))
        corrected[key]=bound
    for key,p,q in pvals:
        if not math.isclose(q,corrected[key],rel_tol=0,abs_tol=1e-10):
            raise ValueError("Original all-qualified-family BH q mismatch")
        if q<.05:
            counts[(key[0],key[1],"exploratory_BH_q_below_0_05")]+=1
    return {
        "source_status":SOURCE_STATUS,
        "GIFT_botanical_parent":4505,
        "geography_parent":8264,
        "source_endpoints":13,
        "n_full_region_status_panels":len(out),
        "n_qualified_within_genus_quality_tests":len(pvals),
        "scope_region_source_support_counts":{
            "|".join(t):v for t,v in sorted(counts.items())},
        "anonymous_panels":sorted(out,key=lambda r:(
            r["scope"],r["region"],r["trait_endpoint"])),
        "nonclaims":(
            "high/medium are evidence grades on DIFFERENT species, not "
            "repeated measures within species; even conditional on genus "
            "association may reflect biological sorting, coding, sampling "
            "or source-quality bias. BH across nested traits is exploratory. "
            "No within-island trait evolution, guild shift, phylogenetic "
            "independence, true mainland ancestor or causal detection effect."
        )
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("original_source_zip_or_json",type=Path)
    p.add_argument("--sha256")
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    result=audit(read_anon(args.original_source_zip_or_json,args.sha256))
    raw=json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(raw,encoding="utf-8")
    print(raw)


if __name__=="__main__":
    main()
