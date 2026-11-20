"""GIFT 4505 source-Direct evidence-quality / category mismatch sensitivity.

This is NOT duplicate within-species phenotype measurement. The protected
ledger has ONE accepted species × trait normalized state and quality grade.
High/medium label comparisons among different species, including within
the same genus, can reveal SOURCE DOCUMENTATION ASSOCIATION but cannot
separate true biological within-genus variation from taxonomy/coding/effort.

Restore original source SHA-pinned botanical flora, GSHHG region crosswalk and
species-Direct evidence ledger in protected island GitHub Actions. Export ONLY
anonymized 13 genuine endpoints × 3 regions × 2 nested statuses = 78 panels.
Never export original names, taxa, species, island IDs, source entities,
coordinates or botanical category names.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from fd_island_functional_beta import TRAITS, PAIRS

REGIONS=("northern_midlatitude","tropical","southern_extratropical")
SCOPES=("all_recorded","original_native")
QUALITY=("high","medium")
ENDPOINTS=TRAITS+tuple(a+"__"+b for a,b in PAIRS)
SEED=20261121
PERMUTATIONS=199
MIN_HIGH_MEDIUM_PER_GENUS=2
MIN_MATCHED_GENERA=10
MIN_QUALIFIED_SPECIES=30
PIN={
    "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
    "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
    "assign":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44",
}


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(2**20),b""):
            h.update(block)
    return h.hexdigest()


def direct_grade_states(table:pd.DataFrame)->dict:
    required={"accepted_species","trait_name","normalized_value",
              "evidence_scope","evidence_quality"}
    if not required<=set(table):
        raise ValueError("Missing original source Direct evidence columns")
    q=table.loc[
        table.trait_name.isin(TRAITS)
        & table.evidence_scope.eq("species_direct")
        & table.evidence_quality.isin(QUALITY)
        & table.normalized_value.notna()
    ].copy()
    q["normalized_value"]=q.normalized_value.astype(str).str.strip()
    q=q.loc[
        q.normalized_value.ne("")
        & ~q.normalized_value.str.contains("|",regex=False)
        & ~q.normalized_value.str.lower().isin(
            ("none","nan","null","unknown","unresolved"))
    ]
    if (q.accepted_species.isna().any()
            or q.duplicated(["accepted_species","trait_name"]).any()):
        raise ValueError("Original evidence has duplicated or unidentified species×trait")
    state={}
    grade={}
    for name in TRAITS:
        x=q.loc[q.trait_name.eq(name)]
        state[name]=dict(zip(x.accepted_species,x.normalized_value))
        grade[name]=dict(zip(x.accepted_species,x.evidence_quality))
    for a,b in PAIRS:
        label=a+"__"+b
        shared=sorted(set(state[a])&set(state[b]))
        state[label]={name:(state[a][name],state[b][name]) for name in shared}
        # A true two-axis Direct pair has high quality ONLY if BOTH
        # independent Direct trait labels are high on the SAME plant name.
        grade[label]={
            name:"high" if grade[a][name]=="high" and grade[b][name]=="high"
            else "medium" for name in shared
        }
    if set(state)!=set(ENDPOINTS) or len(state)!=13:
        raise AssertionError("Expected 10 singles + 3 true same-species pairs")
    return {name:{"state":state[name],"quality":grade[name]}
            for name in ENDPOINTS}


def botanical_genus(species:str)->str:
    chunks=species.split()
    if len(chunks)<2 or not chunks[0] or not chunks[0][0].isupper():
        raise ValueError("Invalid protected source accepted botanical binomial")
    return chunks[0]


def variation_panel(
    species:set[str], endpoint:dict,
    *, universe_count:int, seed:int, nperm:int=PERMUTATIONS,
    min_genus:int=MIN_MATCHED_GENERA) -> dict:
    """Quality-conditioned category TV, genus-stratified permutation null.

    Original accepted species are distinct before counting: the same name
    across 100 islands does NOT count as 100 independent phenotypes.
    No actual source name/genus/category label is emitted.
    """
    if not isinstance(species,set) or universe_count<0 or not species:
        # A true empty group yields explicitly no support.
        if species:
            raise ValueError("Unexpected nonempty original species group")
    values=endpoint["state"]
    grades=endpoint["quality"]
    if set(values)!=set(grades):
        raise ValueError("Original category/quality mappings misaligned")
    n_original=len(species)
    if n_original>universe_count or nperm<50:
        raise ValueError("Original source cohort size or permutation count invalid")
    scored=sorted(species&set(values))
    categories=sorted(set(values.values()),key=repr)
    codes={category:i for i,category in enumerate(categories)}
    K=len(codes)
    high=[name for name in scored if grades[name]=="high"]
    med=[name for name in scored if grades[name]=="medium"]
    if len(high)+len(med)!=len(scored) or len(scored)>n_original:
        raise AssertionError("Direct evidence category/quality original denominator drift")
    group=defaultdict(list)
    for name in scored:
        group[botanical_genus(name)].append((codes[values[name]],grades[name]=="high"))
    matched=[]
    for genus,items in sorted(group.items()):
        idx=np.array([p[0] for p in items],dtype=int)
        flag=np.array([p[1] for p in items],dtype=bool)
        nh=int(flag.sum())
        nm=len(items)-nh
        if nh>=MIN_HIGH_MEDIUM_PER_GENUS and nm>=MIN_HIGH_MEDIUM_PER_GENUS:
            matched.append((idx,nh))
    def tv(a,b):
        if len(a)==0 or len(b)==0:
            return None
        return float(.5*np.abs(
            np.bincount(a,minlength=K)/len(a)
            -np.bincount(b,minlength=K)/len(b)).sum())
    pooled=None
    if high and med:
        pooled=tv(np.array([codes[values[z]] for z in high]),
                  np.array([codes[values[z]] for z in med]))
    observed=None
    null_mean=null95=pvalue=centred=None
    supported=(
        len(matched)>=min_genus and len(scored)>=MIN_QUALIFIED_SPECIES
        and len(high)>=MIN_QUALIFIED_SPECIES
        and len(med)>=MIN_QUALIFIED_SPECIES
    )
    if matched:
        obs=[]
        for idx,nh in matched:
            flags=np.array([a==0 for a in range(len(idx))],dtype=bool)
            # Actual label positions must be preserved separately; repeated
            # sorted source names are private and never serialized.
        matched_real=[]
        for genus,items in sorted(group.items()):
            idx=np.array([p[0] for p in items],dtype=int)
            flags=np.array([p[1] for p in items],dtype=bool)
            if flags.sum()>=MIN_HIGH_MEDIUM_PER_GENUS and (~flags).sum()>=MIN_HIGH_MEDIUM_PER_GENUS:
                matched_real.append((idx,flags))
                obs.append(tv(idx[flags],idx[~flags]))
        observed=float(np.mean(obs))
        if supported:
            rng=np.random.default_rng(seed)
            null=np.empty(nperm,float)
            for i in range(nperm):
                vals=[]
                for idx,flags in matched_real:
                    nh=int(flags.sum())
                    chosen=rng.choice(len(idx),size=nh,replace=False)
                    mask=np.zeros(len(idx),dtype=bool)
                    mask[chosen]=True
                    vals.append(tv(idx[mask],idx[~mask]))
                null[i]=float(np.mean(vals))
            null_mean=float(null.mean())
            null95=float(np.quantile(null,.95))
            centred=observed-null_mean
            pvalue=float((1+sum(null>=observed-1e-12))/(1+nperm))
    return {
        "n_original_distinct_accepted_species":n_original,
        "n_true_species_Direct_coverage":len(scored),
        "original_distinct_species_Direct_coverage_fraction":(
            len(scored)/n_original if n_original else None),
        "n_direct_quality_high":len(high),
        "n_direct_quality_medium":len(med),
        "high_grade_fraction_among_Direct":(
            len(high)/len(scored) if scored else None),
        "n_source_distinct_genera_with_Direct":len(group),
        "n_within_genus_quality_matched_genera":len(matched),
        "n_different_category_states":K,
        "pooled_high_medium_category_total_variation":pooled,
        "within_genus_equal_mean_high_medium_total_variation":observed,
        "within_genus_permutation_qualified":bool(supported),
        "within_genus_null_mean_total_variation":null_mean,
        "within_genus_null_95th_total_variation":null95,
        "within_genus_observed_minus_permutation_mean_TV":centred,
        "within_genus_quality_label_permutation_upper_tail_p_exploratory":pvalue,
        "high_vs_medium_are_different_plant_species_not_repeated_measurements":True,
        "biological_within_genus_variation_vs_annotation_not_identifiable":True,
    }


def bh_qvalues(rows:list[dict]) -> None:
    valid=[i for i,r in enumerate(rows)
           if r["within_genus_permutation_qualified"]]
    if not valid:
        return
    order=sorted(valid,key=lambda k:rows[k][
        "within_genus_quality_label_permutation_upper_tail_p_exploratory"])
    m=len(valid)
    last=1.0
    for pos in range(m-1,-1,-1):
        i=order[pos]
        p=rows[i]["within_genus_quality_label_permutation_upper_tail_p_exploratory"]
        last=min(last,p*m/(pos+1))
        rows[i]["within_genus_quality_label_permutation_BH_q_exploratory"]=float(last)


def run_audit(flora:pd.DataFrame,geo:pd.DataFrame,traits:pd.DataFrame,
              *,nperm:int=PERMUTATIONS)->dict:
    if (not {"island_id","accepted_species","origin_status"}<=set(flora)
            or not {"island_id","analysis_regime"}<=set(geo)):
        raise ValueError("Original botanical/geo schema missing")
    if (flora.island_id.nunique()!=4505
            or flora.duplicated(["island_id","accepted_species"]).any()
            or geo.island_id.nunique()!=8264
            or geo.duplicated("island_id").any()):
        raise ValueError("Original GIFT 4505 floras / GSHHG 8264 geography changed")
    f=flora.merge(geo[["island_id","analysis_regime"]],
                  on="island_id",how="left",validate="many_to_one")
    if len(f)!=len(flora):
        raise AssertionError("Original botanical flora row count drift on geo join")
    states=direct_grade_states(traits)
    outputs=[]
    for region in REGIONS:
        subset=f.loc[f.analysis_regime.eq(region)]
        regional_parent=int(subset.island_id.nunique())
        for scope in SCOPES:
            scoped=(subset if scope=="all_recorded"
                    else subset.loc[subset.origin_status.eq("native")])
            species=set(scoped.accepted_species.dropna().astype(str))
            for name in ENDPOINTS:
                seed=int.from_bytes(hashlib.sha256(
                    (region+scope+name+str(SEED)).encode()).digest()[:8],"big")
                diagnostic=variation_panel(
                    species,states[name],universe_count=len(species),
                    seed=seed,nperm=nperm)
                outputs.append({
                    "region":region,
                    "scope":scope,
                    "trait_endpoint":name,
                    "original_GIFT_botanical_parent_island_floras":4505,
                    "regional_original_parent_botanical_floras":regional_parent,
                    "source_region_original_distinct_botanical_species":len(species),
                    **diagnostic,
                })
    if len(outputs)!=78 or len(set((
        r["region"],r["scope"],r["trait_endpoint"]) for r in outputs))!=78:
        raise AssertionError("Original complete 13x3x2 evidence panel family drift")
    bh_qvalues(outputs)
    return {
        "status":"ORIGINAL_GIFT4505_SPECIES_DIRECT_QUALITY_WITHIN_GENUS_NULL",
        "original_parent_GIFT_botanical_island_floras":4505,
        "original_parent_geography_islands":8264,
        "n_genuine_Direct_endpoints":13,
        "n_region_scope_endpoint_panels":78,
        "n_within_genus_permutations":nperm,
        "quality_ontology":"species_direct_high_vs_medium_once_per_species",
        "joint_endpoint_high_grade_requires_both_same_species_traits_high":True,
        "original_species_name_and_genus_IDs_never_exported":True,
        "no_independent_within_species_or_local_island_phenotype_replicates":True,
        "BH_q_exploratory_across_all_qualified_78_nested_tests":True,
        "n_permutation_qualified_panels":sum(
            r["within_genus_permutation_qualified"] for r in outputs),
        "panels":outputs,
        "noncausal_ceiling":"Quality labels may track biological category differences among species within genera, source provenance, annotation choices and sampling. This test does not identify ecological pollinator selection, actual island morphological evolution or independent mainland ancestry.",
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-root",type=Path,required=True)
    ap.add_argument("--geography",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--nperm",type=int,default=PERMUTATIONS)
    args=ap.parse_args()
    paths={
        "flora":args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "traits":args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "assign":args.original_root/"fixed/source/result/island_source_assignments.csv.gz",
    }
    for k,path in paths.items():
        if sha(path)!=PIN[k]:
            raise ValueError("Original source botanical SHA changed: "+k)
    flora=pd.read_csv(paths["flora"],dtype=str,
                      usecols=["island_id","accepted_species","origin_status"])
    geo=pd.read_csv(args.geography,dtype=str,
                    usecols=["island_id","analysis_regime"])
    traits=pd.read_csv(paths["traits"],dtype=str,
                       usecols=["accepted_species","trait_name",
                                "normalized_value","evidence_scope","evidence_quality"])
    result=run_audit(flora,geo,traits,nperm=args.nperm)
    raw=json.dumps(result,sort_keys=True,indent=2,allow_nan=False)
    for private in ('"island_id":','"accepted_species":','"entity_ID":',
                    '"spatial_block":','"first_source":','"genus_name":',
                    '"latitude":','"longitude":','"source_genus":'):
        if private in raw:
            raise AssertionError("Protected original botanical identity leaked")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(raw+"\n",encoding="utf-8")
    for r in result["panels"]:
        print("GIFT_ORIGINAL_DIRECT_QUALITY_ANONYMOUS="+json.dumps(r,sort_keys=True))
    print("GIFT_ORIGINAL_DIRECT_QUALITY_SUMMARY="+json.dumps({
        "parent":4505,"panels":78,
        "qualified":result["n_permutation_qualified_panels"],
        "BH_q_below_0_05":sum(
            r.get("within_genus_quality_label_permutation_BH_q_exploratory",1.0)<.05
            for r in result["panels"]),
    },sort_keys=True))


if __name__=="__main__":
    main()
