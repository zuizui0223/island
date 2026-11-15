"""Geologically oceanic vs continental islands: source-audited floral sector contrasts.

The FULL original 4505 GIFT island botanical flora is the population frame.
Only separately independently GEOLOGY-TYPED botanical islands (425 expected)
enter origin contrasts. Original GSHHG continent distance is an adjusted
covariate, NEVER used to assign geological origin. Two source-direct trait
axes can be joined only when BOTH states exist on the SAME original accepted
species. Native-known and all-recorded scopes remain nonindependent cohorts.

Observational island proportion, in-region genus-conditional expected share,
and regional genus-conditional exact absences are distinct estimands.
The in-region genus reference contains focal Direct-labelled plants and is
NOT a historical source ancestral trait or causal evolutionary null.
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
from fd_directional_flower_space_missing_profile_null import (
    prepare_regional_reference,null_expectation_one_island,
)
from fd_directional_flower_isolation_gradients import holm

ORIGINS=("continental","oceanic")
REGIONS=("tropical","northern_midlatitude","southern_extratropical",
         "GLOBAL_POOLED_NONINDEPENDENT")
SCOPES=("all_recorded","original_native")
MIN_SOURCE_DIRECT=5
MIN_ISLANDS_PER_ORIGIN=20
MIN_BLOCKS_PER_ORIGIN=5
MIN_BLOCKS_ALL=10
N_BOOT=999

SECTORS=(
    ("bilateral_zygomorphic","floral_symmetry","zygomorphic",True),
    ("tube_or_bell","floral_form","tube_or_bell",True),
    ("open_radial","floral_form","open_radial",True),
    ("deepest_tube","tube_depth_class","deep",True),
    ("strict_self_compatible","self_incompatibility","SC",True),
    ("strict_autonomous_selfing","autonomous_selfing_capacity","autonomous",True),
    ("joint_bilateral_tube_or_bell","floral_form__floral_symmetry","joint_special",True),
    ("joint_open_radial_actinomorphic","floral_form__floral_symmetry","joint_open",False))
STRICT_BINARY={
    "self_incompatibility":{"SC","SI"},
    "autonomous_selfing_capacity":{"autonomous","absent"},
}


def valid_sector_state(endpoint:str,label,value):
    if endpoint=="floral_symmetry" and label=="zygomorphic":
        return value=="zygomorphic"
    if endpoint=="floral_form" and label=="tube_or_bell":
        return value in {"tubular","bell_campanulate"}
    if endpoint=="floral_form" and label=="open_radial":
        return value=="open_radial"
    if endpoint=="tube_depth_class" and label=="deep":
        return value=="deep"
    if endpoint=="self_incompatibility" and label=="SC":
        return value=="SC"
    if endpoint=="autonomous_selfing_capacity" and label=="autonomous":
        return value=="autonomous"
    if endpoint=="floral_form__floral_symmetry" and label=="joint_special":
        return value[1]=="zygomorphic" and value[0] in {
            "tubular","bell_campanulate"}
    if endpoint=="floral_form__floral_symmetry" and label=="joint_open":
        return value==("open_radial","actinomorphic")
    raise ValueError("Unknown predeclared biological functional sector")


def sector_roster(states:dict,endpoint:str):
    if endpoint=="floral_form__floral_symmetry":
        a,b=states["floral_form"],states["floral_symmetry"]
        return {sp:(a[sp],b[sp]) for sp in set(a)&set(b)}
    d=states[endpoint]
    if endpoint in STRICT_BINARY:
        return {sp:v for sp,v in d.items() if v in STRICT_BINARY[endpoint]}
    return d


def genus_ref(expected_labels:dict,region_union:set,label:str,endpoint:str):
    """One full observed binary/joint trait label per source region species."""
    pool=region_union &set(expected_labels)
    groups=defaultdict(list)
    for sp in pool:
        g=sp.split()[0] if len(sp.split())>=2 else "fixed:"+sp
        groups[g].append(float(valid_sector_state(endpoint,label,expected_labels[sp])))
    return {g:float(np.mean(v)) for g,v in groups.items()}


def genus_expected_share(species:set,ref:dict):
    p=[]
    for s in species:
        g=s.split()[0] if len(s.split())>=2 else "fixed:"+s
        if g not in ref:raise AssertionError("Original measured species outside reference")
        p.append(ref[g])
    return float(np.mean(p))


def source_group_cluster_fit(df:pd.DataFrame,response:str,cluster:str,
                             *,n_boot:int=N_BOOT,seed:int=20261115):
    """Origin coefficient under equal original physical 10-degree block weights."""
    required={"origin","spatial_block",
              "x_continent_distance","x_area","x_coverage",
              "x_direct_n","x_latitude",response,cluster}
    if not required<=set(df):raise ValueError("Original model fields incomplete")
    d=df.dropna(subset=list(required)).copy()
    n=len(d)
    groupnames=sorted(d[cluster].unique())
    G=len(groupnames)
    ocean=d.origin.eq("oceanic")
    continent=d.origin.eq("continental")
    nb=dict(
        n_islands=n,n_oceanic=int(ocean.sum()),n_continental=int(continent.sum()),
        n_clusters=G,n_oceanic_clusters=d.loc[ocean,cluster].nunique(),
        n_continental_clusters=d.loc[continent,cluster].nunique(),
        n_physical_blocks_with_both_origin_types=int(sum(
            len(set(frame.origin))==2 for _,frame in d.groupby("spatial_block"))))
    if (n<2*MIN_ISLANDS_PER_ORIGIN or
        nb["n_oceanic"]<MIN_ISLANDS_PER_ORIGIN or
        nb["n_continental"]<MIN_ISLANDS_PER_ORIGIN or
        nb["n_oceanic_clusters"]<MIN_BLOCKS_PER_ORIGIN or
        nb["n_continental_clusters"]<MIN_BLOCKS_PER_ORIGIN or
        G<MIN_BLOCKS_ALL):
        return dict(**nb,status="INSUFFICIENT_GEOLOGY_ORIGIN_SOURCE_SUPPORT")
    covariates=["x_continent_distance","x_area","x_coverage",
                "x_direct_n","x_latitude"]
    X=np.column_stack([np.ones(n),ocean.to_numpy(float),
                       *[d[x].to_numpy(float) for x in covariates]])
    y=d[response].to_numpy(float)
    p=X.shape[1]
    if np.any(~np.isfinite(X)) or np.any(~np.isfinite(y)):
        raise ValueError("Nonfinite protected geology category response")
    block_X=[];block_y=[]
    for block,idx in d.groupby(cluster,sort=True).indices.items():
        Z=X[idx]
        block_X.append(Z.T@Z/len(idx))
        block_y.append(Z.T@y[idx]/len(idx))
    M=np.stack(block_X)
    b=np.stack(block_y)
    gram=M.sum(0)
    if np.linalg.matrix_rank(gram,tol=1e-9)<p:
        return dict(**nb,status="ORIGINAL_COVARIATES_RANK_DEFICIENT")
    beta=np.linalg.solve(gram,b.sum(0))
    v=np.linalg.inv(gram)
    S=b-np.einsum("ijk,k->ij",M,beta)
    vc=v@(S.T@S)@v.T *G/(G-1)*(n-1)/(n-p)
    se=np.sqrt(max(float(vc[1,1]),0.))
    if se<1e-10 or not np.isfinite(se):
        return dict(**nb,status="NONIDENTIFIABLE_ORIGIN_CLUSTER_SE")
    pval=float(2*student_t.sf(abs(beta[1]/se),df=G-1))
    rng=np.random.default_rng(seed)
    estimates=[]
    for _ in range(n_boot):
        ix=rng.integers(0,G,size=G)
        trial=M[ix].sum(0)
        if np.linalg.matrix_rank(trial,tol=1e-9)<p:continue
        estimates.append(float(np.linalg.solve(trial,b[ix].sum(0))[1]))
    if len(estimates)<int(.8*n_boot):
        return dict(**nb,status="BOOTSTRAP_ORIGIN_COEFFICIENT_UNIDENTIFIABLE",
                    n_valid_bootstrap=len(estimates))
    values=np.asarray(estimates)
    loo=[]
    for j in range(G):
        trial=gram-M[j]
        if np.linalg.matrix_rank(trial,tol=1e-9)==p:
            loo.append(float(np.linalg.solve(trial,b.sum(0)-b[j])[1]))
    means=df.groupby("origin")[response].mean()
    group_equal=df.groupby([cluster,"origin"])[response].mean().reset_index()
    group_means=group_equal.groupby("origin")[response].mean()
    cov_matrix=np.corrcoef(X[:,1:].T)[0,1:]
    return dict(**nb,
                status="EXPLORATORY_ADJUSTED_GEOLOGICAL_ORIGIN_CATEGORY_CONTRAST",
                adjusted_oceanic_minus_continental_fraction=float(beta[1]),
                origin_cluster_CR1_robust_SE=float(se),
                origin_cluster_two_sided_t_exploratory_p=pval,
                origin_cluster_999resample_quantile95=[
                    float(np.quantile(values,.025)),float(np.quantile(values,.975))],
                n_valid_cluster_resamples=len(values),
                origin_cluster_resample_positive_fraction=float(np.mean(values>0)),
                delete_one_cluster_min=float(min(loo)) if loo else None,
                delete_one_cluster_max=float(max(loo)) if loo else None,
                raw_island_equal_oceanic_minus_continental=float(
                    means.get("oceanic",float("nan"))-
                    means.get("continental",float("nan"))),
                raw_group_equal_oceanic_minus_continental=float(
                    group_means.get("oceanic",float("nan"))-
                    group_means.get("continental",float("nan"))),
                max_absolute_origin_covariate_correlation=float(
                    np.nanmax(np.abs(cov_matrix))))


def analyze(states:dict,flora:pd.DataFrame,geo:pd.DataFrame,
            assignment:pd.DataFrame,*,n_boot:int=N_BOOT)->dict:
    if not {"island_id","accepted_species","origin_status"}<=set(flora):
        raise ValueError("Missing original island species source")
    if flora.duplicated(["island_id","accepted_species"]).any() or not (
        flora.origin_status.isin(("native","introduced","unresolved")).all()):
        raise ValueError("Original botanical island names/status invalid")
    geo_fields={"island_id","geological_origin","analysis_regime",
                "distance_to_continent_km","area_km2","island_latitude",
                "spatial_block"}
    if not geo_fields<=set(geo) or geo.island_id.duplicated().any():
        raise ValueError("Independent geological source mapping incomplete")
    if not {"island_id","source_mode","source_rank","entity_ID"}<=set(assignment):
        raise ValueError("Independent native source candidates missing")
    full=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    native=flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    if len(full)!=4505:raise ValueError("Original population of 4505 changed")
    idx=geo.set_index("island_id").reindex(sorted(full))
    typed=idx.geological_origin.isin(ORIGINS)
    if int(typed.sum())!=425:
        raise ValueError("Source geological-type botanical flora coverage drift")
    near=assignment.loc[assignment.source_mode.eq("geo_k5")].copy()
    near.source_rank=pd.to_numeric(near.source_rank,errors="raise")
    if near.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original GIFT source ranks duplicated")
    first=near.loc[near.source_rank.eq(1)].set_index("island_id").entity_ID.to_dict()
    codebook={}
    for name,endpoint,label,primary in SECTORS:
        raw=sector_roster(states,endpoint)
        if not raw:raise ValueError("Original Direct phenotype sector empty")
        codebook[name]=(endpoint,label,raw,primary)
    controls=defaultdict(dict)
    for scope,source in (("all_recorded",full),("original_native",native)):
        # In-region genus nulls built using ALL source-available islands
        # within the same broad region, NOT selected typed-only islands.
        for region in set(REGIONS)-{"GLOBAL_POOLED_NONINDEPENDENT"}:
            region_all=[id for id in full if idx.at[id,"analysis_regime"]==region]
            pool=set().union(*(source.get(x,set()) for x in region_all))
            for name,(endpoint,label,states_sector,_) in codebook.items():
                ref=genus_ref(states_sector,pool,label,endpoint)
                controls[(scope,region,name)]=ref
    original_typed=[]
    out=[]
    n_source_by_scope={}
    for scope,source in (("all_recorded",full),("original_native",native)):
        for region in REGIONS:
            included=[id for id in full if typed.loc[id] and (
                 idx.at[id,"analysis_regime"]==region
                 if region!="GLOBAL_POOLED_NONINDEPENDENT"
                 else idx.at[id,"analysis_regime"] in REGIONS[:3])]
            n_geo={origin:sum(idx.at[id,"geological_origin"]==origin
                              for id in included) for origin in ORIGINS}
            original_typed.append(dict(scope=scope,region=region,
                n_typed_original_island_floras=len(included),
                n_original_oceanic=n_geo["oceanic"],
                n_original_continental=n_geo["continental"]))
            for name,(endpoint,label,trait_values,primary) in codebook.items():
                frames=[]
                observed_count={origin:0 for origin in ORIGINS}
                for id in included:
                    ro=source.get(id,set())
                    actual=ro&set(trait_values)
                    if len(actual)<MIN_SOURCE_DIRECT:continue
                    obs=idx.loc[id]
                    try:
                        dist=float(obs.distance_to_continent_km)
                        area=float(obs.area_km2)
                        lat=float(obs.island_latitude)
                    except (TypeError,ValueError):continue
                    if not np.isfinite([dist,area,lat]).all() or (
                        dist<0 or area<=0 or abs(lat)>90 or
                        pd.isna(obs.spatial_block)):
                        continue
                    region_ref=str(obs.analysis_regime)
                    ref=controls.get((scope,region_ref,name))
                    if not ref:continue
                    p=float(np.mean([valid_sector_state(endpoint,label,trait_values[sp])
                                     for sp in actual]))
                    g=genus_expected_share(actual,ref)
                    first_source=first.get(id)
                    frames.append(dict(
                        origin=str(obs.geological_origin),
                        region=region_ref,
                        spatial_block=str(obs.spatial_block),
                        first_source=first_source,
                        response=p,genus_relative_response=p-g,
                        n_source_Direct=len(actual),
                        documentable_fraction=len(actual)/len(ro),
                        x_continent_distance=math.log2(1+dist),
                        x_area=math.log1p(area),
                        x_coverage=len(actual)/len(ro),
                        x_direct_n=math.log(len(actual)),
                        x_latitude=abs(lat)))
                    observed_count[str(obs.geological_origin)]+=1
                frame=pd.DataFrame(frames)
                report=dict(scope=scope,region=region,sector=name,
                    trait_endpoint=endpoint,original_frozen_primary_tropical=(
                        scope=="all_recorded" and region=="tropical" and primary),
                    n_original_geologically_typed_floras=len(included),
                    n_geologically_typed_oceanic=n_geo["oceanic"],
                    n_geologically_typed_continental=n_geo["continental"],
                    n_source_direct_n5_valid_oceanic=observed_count["oceanic"],
                    n_source_direct_n5_valid_continental=observed_count["continental"])
                if frame.empty:
                    report["status"]="INSUFFICIENT_ORIGINAL_SOURCE_DIRECT_SECTOR"
                    out.append(report);continue
                for response,prefix in (
                    ("response","observed"),
                    ("genus_relative_response","in_region_genus_relative")):
                    for cluster,label_group in (
                        ("spatial_block","physical_block"),
                        ("first_source","GIFT_first_mainland_source")):
                        if cluster not in frame:continue
                        estimate=source_group_cluster_fit(
                            frame,response,cluster,n_boot=n_boot,
                            seed=int.from_bytes(hashlib.sha256(
                                (scope+region+name+prefix+label_group).encode()
                                ).digest()[:8],"big"))
                        for key,val in estimate.items():
                            report[prefix+"__"+label_group+"__"+key]=val
                report["status"]=report.get("observed__physical_block__status")
                out.append(report)
    # Explicitly preserve seven hypothesis family even when original
    # typed island counts cannot identify all seven outcomes.
    primary_rows=[z for z in out if z["original_frozen_primary_tropical"]]
    if len(primary_rows)!=7:
        raise AssertionError("Original geology tropical primary family truncated")
    for variant in ("observed","in_region_genus_relative"):
        p=[]
        for record in primary_rows:
            name=variant+"__physical_block__"
            valid=(record.get(name+"status")==
                   "EXPLORATORY_ADJUSTED_GEOLOGICAL_ORIGIN_CATEGORY_CONTRAST")
            p.append(record[name+"origin_cluster_two_sided_t_exploratory_p"]
                     if valid else 1.)
        padj=holm(p)
        for record,x in zip(primary_rows,padj):
            record[variant+"__Holm7_full_family_p"]=float(x)
            record[variant+"__original_main_family_test_identifiable"]=(
                record.get(variant+"__physical_block__status")==
                "EXPLORATORY_ADJUSTED_GEOLOGICAL_ORIGIN_CATEGORY_CONTRAST")
    return dict(
        status="ORIGINAL_GEOLOGICAL_ORIGIN_DIRECTIONAL_FLORAL_SECTORS",
        n_original_botanical_island_floras=4505,
        n_independently_geologically_typed_source_floras=int(typed.sum()),
        n_typed_physical_origin_oceanic=int(sum(
            idx.loc[typed,"geological_origin"]=="oceanic")),
        n_typed_physical_origin_continental=int(sum(
            idx.loc[typed,"geological_origin"]=="continental")),
        n_tropical_frozen_main_sectors=len(primary_rows),
        original_geology_support_by_scope_and_region=original_typed,
        original_sector_oceanic_minus_continental_contrasts=out,
        claims=dict(
            study_universe_all_4505_not_geology_only=True,
            distance_not_substitute_for_true_geological_type=True,
            geology_matching_audited_Roeble_Zizka_not_LGM_unconnected=True,
            original_direct_carrier_fraction_and_documentation_denominators=True,
            floral_structure_not_actual_pollination_visitor_guild=True,
            self_compatible_not_autonomous_selfing=True,
            joint_sector_has_real_same_species_morphological_pairs=True,
            seven_frozen_tropical_holm_hypotheses_include_failed_source_tests=True,
            in_region_genus_source_reference_contains_focal_flora=True,
            type_source_coverage_nonrandom_and_geographically_limited=True,
            no_adaptive_evolution_extinction_or_pollination_service_inference=True,
            no_private_island_plant_genus_source_entity_or_coordinates_exported=True))


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as fh:
        for c in iter(lambda:fh.read(1<<20),b""):h.update(c)
    return h.hexdigest()


def main():
    parser=argparse.ArgumentParser()
    for key in ("original-root","geology","output"):
        parser.add_argument("--"+key,type=Path,required=True)
    p=parser.parse_args()
    source={
        "traits":p.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "flora":p.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "assign":p.original_root/"fixed/source/result/island_source_assignments.csv.gz"}
    pins={
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        "assign":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44"}
    for k,path in source.items():
        if sha(path)!=pins[k]:raise ValueError("Original source geological Direct trait fingerprint drift "+k)
    trait=strict_traits(pd.read_csv(source["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(source["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    geology=pd.read_csv(p.geology,dtype=str,usecols=[
        "island_id","geological_origin","analysis_regime",
        "distance_to_continent_km","area_km2","island_latitude","spatial_block"])
    mainland=pd.read_csv(source["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    if geology.island_id.nunique()!=8264:
        raise ValueError("Independently audited original 8264 geographical islands changed")
    result=analyze(trait,flora,geology,mainland)
    if (result["n_original_botanical_island_floras"]!=4505 or
            result["n_independently_geologically_typed_source_floras"]!=425):
        raise ValueError("Original typed geology observed population drift")
    body=json.dumps(result,sort_keys=True,indent=2,allow_nan=False)
    for key in ('"island_id"','"accepted_species"','"spatial_block"',
                '"entity_ID"','"island_latitude"','"island_longitude"',
                '"source_signature"'):
        if key in body:
            raise ValueError("Original private botanical/geological entity escaped")
    p.output.parent.mkdir(parents=True,exist_ok=True)
    p.output.write_text(body+"\n")
    print(json.dumps(dict(status=result["status"],
                          n_original_source_floras=4505,
                          n_real_geologically_typed_floras=425,
                          n_predeclared_main_tropical_tests=7,
                          source_origin_noncausal=True)))


if __name__=="__main__":main()
