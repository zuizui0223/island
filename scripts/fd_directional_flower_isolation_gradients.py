"""Original island directional flower-phenotype share gradients vs continent isolation.

Source: 4505 botanical island floras, one original high/medium Direct species-wide
trait codebook, and independently joined GSHHG geographic distance/area/blocks.
No three-trait complete-case gate, no KDE volume, no pollinator identity claim.

The fixed 6 primary slopes test original structural sector share as a function
of log2(1+distance_to_continent_km), adjusting log(area), documentation ratio,
number of observed Direct species and absolute latitude. Every source species
in a given island is ONE recorded carrier. Physical 10degree block-cluster
robust inference and source-cluster separate sensitivity. The in-region
within-genus reference is IN-SAMPLE, not causal or independent history.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from collections import Counter,defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from fd_island_functional_beta import strict_traits

MIN_DIRECT=5
MIN_ISLANDS=30
MIN_BLOCKS=10
N_BOOT=999
MAIN=[
    ("northern_midlatitude","floral_symmetry","bilateral_zygomorphic","decrease"),
    ("northern_midlatitude","floral_form","tube_or_bell","decrease"),
    ("northern_midlatitude","floral_form","open_radial","nondecrease"),
    ("tropical","floral_symmetry","bilateral_zygomorphic","decrease"),
    ("tropical","floral_form","tube_or_bell","nondecrease"),
    ("northern_midlatitude","tube_depth_class","deepest","decrease")]
EXTRA=[
    ("tropical","tube_depth_class","deepest"),
    ("southern_extratropical","floral_symmetry","bilateral_zygomorphic"),
    ("southern_extratropical","floral_form","tube_or_bell"),
    ("southern_extratropical","floral_form","open_radial"),
    ("southern_extratropical","tube_depth_class","deepest"),
    ("highlatitude_ge60_nested","floral_symmetry","bilateral_zygomorphic"),
    ("highlatitude_ge60_nested","floral_form","tube_or_bell"),
    ("highlatitude_ge60_nested","floral_form","open_radial"),
    ("highlatitude_ge60_nested","tube_depth_class","deepest")]
EXPECTED={
    "floral_symmetry":{"actinomorphic","asymmetric","zygomorphic"},
    "floral_form":{"bell_campanulate","bilabiate","brush_puff",
                   "composite_head","funnel_trumpet","open_radial",
                   "other_described","papilionaceous","salverform",
                   "spurred","tubular","urn_urceolate"},
    "tube_depth_class":{"absent_or_open","shallow","intermediate","deep"}}


def sector(trait,state,sector_name):
    if trait=="floral_symmetry" and sector_name=="bilateral_zygomorphic":
        return state=="zygomorphic"
    if trait=="floral_form" and sector_name=="tube_or_bell":
        return state in {"tubular","bell_campanulate"}
    if trait=="floral_form" and sector_name=="open_radial":
        return state=="open_radial"
    if trait=="tube_depth_class" and sector_name=="deepest":
        return state=="deep"
    raise ValueError("Undeclared original floral directional sector")


def holm(p):
    """All six predeclared ORIGINAL tests, not observed positive subset."""
    x=np.asarray(p,float)
    if x.ndim!=1 or not np.isfinite(x).all() or (x<0).any() or (x>1).any():
        raise ValueError("Predeclared p-values must all be calculable")
    order=np.argsort(x,kind="stable")
    adjusted=np.empty(len(x))
    past=0.
    for rank,ix in enumerate(order):
        past=max(past,min(1.,(len(x)-rank)*x[ix]))
        adjusted[ix]=past
    return adjusted


def weighted_cluster_fit(frame:pd.DataFrame,*,cluster="spatial_block",
                         n_boot=N_BOOT,seed=20261115):
    """Analytic block-equal OLS and CR1 group-sandwich t uncertainty.

    Each original geographic group contributes ONE equal-weight regression
    vote regardless of whether it has many floras. Bivariate trait fractions
    are NOT interpreted as independent species outcomes.
    """
    use=["response","log2_isolation_plus1","log1p_area",
         "direct_coverage","log_direct_n","absolute_latitude",cluster]
    if any(k not in frame for k in use):
        raise ValueError("Original data frame lacks required response/covariates")
    s=frame.dropna(subset=use).copy()
    G=s[cluster].nunique()
    n=len(s)
    p=6
    if (n<MIN_ISLANDS or G<MIN_BLOCKS or
            s["log2_isolation_plus1"].nunique()<8):
        return dict(status="INSUFFICIENT_ORIGINAL_GEOGRAPHIC_SUPPORT",
                    n_islands=n,n_clusters=G)
    X=np.column_stack((
        np.ones(n),
        s["log2_isolation_plus1"].to_numpy(float),
        s["log1p_area"].to_numpy(float),
        s["direct_coverage"].to_numpy(float),
        s["log_direct_n"].to_numpy(float),
        s["absolute_latitude"].to_numpy(float)))
    y=s.response.to_numpy(float)
    if np.any(~np.isfinite(X)) or np.any(~np.isfinite(y)):
        raise ValueError("Invalid finite source data")
    # Aggregate one X'WX contribution per cluster, weighting each island 1/n_g.
    mats=[];rhs=[];n_by=[]
    for _,ids in s.groupby(cluster,sort=True).indices.items():
        Z=X[ids];q=y[ids]
        mats.append(Z.T@Z/len(ids))
        rhs.append(Z.T@q/len(ids))
        n_by.append(len(ids))
    mats=np.stack(mats)
    rhs=np.stack(rhs)
    normal=mats.sum(axis=0)
    if np.linalg.matrix_rank(normal, tol=1e-9)!=p:
        return dict(status="RANK_DEFICIENT_ORIGINAL_GEOGRAPHIC_SUPPORT",
                    n_islands=n,n_clusters=G)
    try:
        beta=np.linalg.solve(normal,rhs.sum(axis=0))
    except np.linalg.LinAlgError:
        return dict(status="RANK_DEFICIENT_ORIGINAL_GEOGRAPHIC_SUPPORT",
                    n_islands=n,n_clusters=G)
    H=np.linalg.inv(normal)
    score=rhs-np.einsum("ijk,k->ij",mats,beta)
    meat=score.T@score
    V=H@meat@H.T * G/(G-1)* (n-1)/(n-p)
    se=float(np.sqrt(max(0,V[1,1])))
    if not np.isfinite(se) or se<1e-11:
        return dict(status="NONIDENTIFIABLE_CLUSTER_SLOPE_STANDARD_ERROR",
                    n_islands=n,n_clusters=G,
                    estimate_per_doubling=float(beta[1]))
    tval=float(beta[1]/se)
    pval=float(2*student_t.sf(abs(tval),df=G-1))
    rng=np.random.default_rng(seed)
    b=[]
    for _ in range(n_boot):
        idx=rng.integers(G,size=G)
        normal_star=mats[idx].sum(0)
        if np.linalg.matrix_rank(normal_star,tol=1e-9)!=p:continue
        coeff=np.linalg.solve(normal_star,rhs[idx].sum(0))
        b.append(float(coeff[1]))
    if len(b)<max(50,int(.8*n_boot)):
        return dict(status="INSUFFICIENT_BLOCK_BOOTSTRAP_IDENTIFICATION",
                    n_islands=n,n_clusters=G,n_success_boot=len(b))
    boot=np.asarray(b,float)
    loo=[]
    for i in range(G):
        lhs=normal-mats[i]
        if np.linalg.matrix_rank(lhs,tol=1e-9)==p:
            loo.append(float(np.linalg.solve(lhs,rhs.sum(0)-rhs[i])[1]))
    return dict(
        status="ORIGINAL_ISLAND_BLOCK_EQUAL_DESCRIPTIVE_REGRESSION",
        n_islands=n,n_clusters=G,
        n_islands_per_cluster_min=int(min(n_by)),
        n_islands_per_cluster_median=float(np.median(n_by)),
        n_islands_per_cluster_max=int(max(n_by)),
        estimate_per_doubling=float(beta[1]),
        block_robust_cluster_se=se,
        block_robust_two_sided_t_p_exploratory=pval,
        block_resample_quantile_95=[
            float(np.quantile(boot,.025)),float(np.quantile(boot,.975))],
        fraction_block_bootstrap_slopes_negative=float(np.mean(boot<0)),
        n_valid_block_resamples=len(boot),
        delete_one_block_min=float(min(loo)) if loo else None,
        delete_one_block_max=float(max(loo)) if loo else None,
        baseline_mean_actual_sector_share=float(y.mean()),
        mean_original_doc_coverage=float(s.direct_coverage.mean()),
        mean_direct_carriers=float(np.exp(s.log_direct_n).mean()))


def expected_same_region_genus_sector(sets:list[set],state:dict,all_region:list[set],
                                      *,trait:str,sector_name:str):
    """Exact expected sector prevalence under in-region whole genus label shuffle.

    Conditional on observed genus carrier counts and ORIGINAL regional
    measured accepted species pool, not an independent source or selection null.
    """
    names=set().union(*all_region)&set(state)
    group=defaultdict(list)
    for name in names:
        genus=name.split()[0] if len(name.split())>=2 else "fixed:"+name
        group[genus].append(int(sector(trait,state[name],sector_name)))
    p={g:float(np.mean(v)) for g,v in group.items()}
    out=[]
    for roster in sets:
        values=[p[name.split()[0] if len(name.split())>=2 else "fixed:"+name]
                for name in sorted(roster)]
        out.append(float(np.mean(values)) if values else float("nan"))
    return np.asarray(out)


def build_island_rows(states,flora,geo,assign):
    if not {"island_id","accepted_species","origin_status"}<=set(flora):
        raise ValueError("Original name/status columns missing")
    if flora.duplicated(["island_id","accepted_species"]).any() or (
        not flora.origin_status.isin(("native","introduced","unresolved")).all()):
        raise ValueError("Original species occurrence source invalid")
    need={"island_id","analysis_regime","spatial_block",
          "distance_to_continent_km","area_km2","island_latitude"}
    if not need<=set(geo) or geo.island_id.duplicated().any():
        raise ValueError("Missing independent GSHHG source isolation/area")
    for t,allowed in EXPECTED.items():
        if set(states[t].values())!=allowed:
            raise ValueError("Unrecognized original Direct trait ontology: "+t)
    names=sorted(flora.island_id.unique())
    allmap=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    native=flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    g=geo.set_index("island_id").reindex(names)
    source_rank=assign.loc[assign.source_mode.eq("geo_k5")].copy()
    source_rank.source_rank=pd.to_numeric(source_rank.source_rank,errors="raise")
    if source_rank.duplicated(["island_id","source_mode","source_rank"]).any():
        raise ValueError("Original geography native source ranks duplicated")
    first=source_rank.loc[source_rank.source_rank.eq(1)].set_index(
        "island_id").entity_ID.to_dict()
    reports={}
    for scope,lookup in (("all_recorded",allmap),("original_native",native)):
        for region,trait,label,_ in [(x[0],x[1],x[2],"") for x in MAIN]+[
                                (r,t,l,"") for r,t,l in EXTRA]:
            region_ids=[
                i for i in names if (
                    g.at[i,"analysis_regime"]==region
                    if region!="highlatitude_ge60_nested"
                    else pd.notna(g.at[i,"island_latitude"]) and
                         abs(float(g.at[i,"island_latitude"]))>=60)]
            subset=lookup
            regional_rosters=[lookup.get(i,set()) for i in region_ids]
            targets=[]
            valid=[]
            for i in region_ids:
                orig=lookup.get(i,set())
                known=orig&set(states[trait])
                if len(known)<MIN_DIRECT:continue
                v=g.loc[i]
                try:
                    distance=float(v.distance_to_continent_km)
                    area=float(v.area_km2)
                    lat=float(v.island_latitude)
                except (TypeError,ValueError):
                    continue
                if (not np.isfinite([distance,area,lat]).all() or
                        distance<0 or area<=0 or abs(lat)>90 or
                        pd.isna(v.spatial_block)):
                    continue
                successes=sum(sector(trait,states[trait][sp],label)
                              for sp in known)
                targets.append(known)
                valid.append(dict(
                    response=successes/len(known),
                    n_success=int(successes),n_direct=len(known),
                    log2_isolation_plus1=math.log2(1+distance),
                    log1p_area=math.log1p(area),
                    direct_coverage=len(known)/len(orig),
                    log_direct_n=math.log(len(known)),
                    absolute_latitude=abs(lat),
                    spatial_block=str(v.spatial_block),
                    source_cluster=first.get(i) if first.get(i) else None))
            if not valid:
                reports[(scope,region,trait,label)]=pd.DataFrame()
                continue
            frame=pd.DataFrame(valid)
            # Retain same exact original eligible target island ordering.
            genus_mean=expected_same_region_genus_sector(
                targets,states[trait],regional_rosters,trait=trait,
                sector_name=label)
            frame["genus_relative_response"]=frame.response-genus_mean
            frame["genus_expected_sector_share"]=genus_mean
            reports[(scope,region,trait,label)]=frame
    return reports


def analyze(states,flora,geo,assign,n_boot=N_BOOT):
    original_n=flora.island_id.nunique()
    if original_n!=4505:
        raise ValueError("Original botanical 4505 island universe changed")
    rows=build_island_rows(states,flora,geo,assign)
    data=[]
    for key,frame in rows.items():
        scope,regime,trait,label=key
        is_primary=(scope=="all_recorded" and
                    any((regime,trait,label)==(x[0],x[1],x[2]) for x in MAIN))
        row=dict(scope=scope,region=regime,trait=trait,sector=label,
                 primary_directional_family=is_primary,
                 original_population_4505=original_n,
                 n_floras_eligible_with_verified_sector_isolation=len(frame))
        if len(frame)<MIN_ISLANDS:
            row["status"]="INSUFFICIENT_DIRECT_SECTOR_ISOLATION_OPPORTUNITY"
            data.append(row);continue
        for field,prefix in (("response","observed_share"),
                             ("genus_relative_response","in_region_genus_relative")):
            q=frame.copy()
            q["response"]=frame[field]
            for group,suffix in (("spatial_block","physical_block"),
                                 ("source_cluster","GIFT_first_source")):
                stat=weighted_cluster_fit(
                    q,cluster=group,n_boot=n_boot,
                    seed=int.from_bytes(hashlib.sha256(
                        (scope+regime+trait+label+prefix+suffix).encode()).digest()[:8],"big"))
                for f,val in stat.items():
                    row[prefix+"__"+suffix+"__"+f]=val
        row["status"]=row.get("observed_share__physical_block__status")
        data.append(row)
    # Predeclared Holm family requires all SIX raw, named primary outcomes.
    main=[]
    for region,trait,sector_name,expected_direction in MAIN:
        matching=[r for r in data if
                  (r["scope"],r["region"],r["trait"],r["sector"])==
                  ("all_recorded",region,trait,sector_name)]
        if len(matching)!=1:raise AssertionError("Original primary sector missing")
        q=matching[0]
        main.append((q,expected_direction))
    if all(r["observed_share__physical_block__status"]==
           "ORIGINAL_ISLAND_BLOCK_EQUAL_DESCRIPTIVE_REGRESSION"
           for r,d in main):
        adjusted=holm([r["observed_share__physical_block__block_robust_two_sided_t_p_exploratory"]
                       for r,d in main])
        for (row,d),p in zip(main,adjusted):
            row["original_6_sector_family_Holm_adjusted_exploratory_two_sided_p"]=float(p)
            row["hypothesis_direction_before_isolation_gradient"]=d
    else:
        for r,d in main:
            r["original_6_sector_family_Holm_adjusted_exploratory_two_sided_p"]=None
            r["hypothesis_direction_before_isolation_gradient"]=d
    return dict(
        status="ORIGINAL_DIRECTIONAL_FLORAL_REGION_ISOLATION_GRADIENTS",
        n_original_botanical_island_floras=original_n,
        original_primary_region_sector_hypotheses=[dict(
            region=a,trait=b,sector=c,prediction_direction=d) for a,b,c,d in MAIN],
        n_primary_predeclared_original_tests=len(MAIN),
        n_secondary_region_scope_sector_diagnostics=len(data)-len(MAIN),
        original_island_morphological_sector_slope_results=data,
        claims=dict(
            original_population_not_narrowed_from_4505=True,
            marginal_trait_per_island_Direct_carrier_denominators=True,
            zygomorphic_not_assumed_bee_pollination=True,
            real_original_deep_tube_bin_not_flower_form_proxy=True,
            isolation_is_GSHHG_distance_to_continent_not_GIFT_source=True,
            genus_conditional_reference_in_sample_not_evolutionary_null=True,
            bootstrap_and_source_groups_not_independent_phylogenetic_replications=True,
            regression_association_not_category_loss_or_fitness_causation=True,
            no_original_island_species_genus_GIFT_entity_geolocation_exported=True))


def sha(path):
    d=hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):
            d.update(part)
    return d.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for x in ("original-root","geology","output"):
        p.add_argument("--"+x,type=Path,required=True)
    p.add_argument("--n-boot",type=int,default=N_BOOT)
    a=p.parse_args()
    orig={
        "traits":a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "flora":a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "assign":a.original_root/"fixed/source/result/island_source_assignments.csv.gz"}
    pins={
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        "assign":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44"}
    if a.n_boot<99:raise ValueError("Cannot weaken geographic source bootstrap gate")
    for k,path in orig.items():
        if sha(path)!=pins[k]:
            raise ValueError("Original verified source trait status and assignments changed")
    states=strict_traits(pd.read_csv(orig["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value","evidence_scope",
        "evidence_quality"]))
    flora=pd.read_csv(orig["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    geo=pd.read_csv(a.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block",
        "distance_to_continent_km","area_km2","island_latitude"])
    assignment=pd.read_csv(orig["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    if geo.island_id.nunique()!=8264:
        raise ValueError("Independent original geographic island universe changed")
    result=analyze(states,flora,geo,assignment,n_boot=a.n_boot)
    content=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    for forbidden in ('"island_id"','"accepted_species"','"entity_ID"',
                      '"spatial_block"','"source_signature"',
                      '"island_latitude"','"island_longitude"'):
        if forbidden in content:raise ValueError("Original botanical identifiers leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(content+"\n")
    print(json.dumps(dict(status=result["status"],
                          n_original_floras=4505,
                          n_source_primary_tests=len(MAIN),
                          n_anonymous_sector_result_rows=len(result[
                            "original_island_morphological_sector_slope_results"]),
                          no_pollination_selection_inference=True)))


if __name__=="__main__":main()
