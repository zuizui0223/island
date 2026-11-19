"""Exact 4505-island floral category beta algebra: genus assembly vs residual.

One true original high/medium species Direct endpoint at a time (10 single,
3 REAL same-species pairs, 2 original floristic status scopes), same broad
region, zero overlapping accepted species in FULL GIFT botanical floras, and
different original physically indexed 10deg GSHHG blocks. The original
within-genus category reference INCLUDES focal records; this is a purely
DESCRIPTIVE identity, never independent heldout prediction or adaptation.

||p_i-p_j||² = ||m_i-m_j||² + ||e_i-e_j||² + 2Δm·Δe.
Cross term MAY BE NEGATIVE; do not suppress it or turn genus term into a
causal 'fraction of adaptation explained'. Never treat millions of dyads as
independent island/evolutionary replicates. Export only anonymous aggregates.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from fd_island_functional_beta import (
    strict_traits,TRAITS,PAIRS,membership,prepare_endpoint,
    geographical_distance)

MIN_PAIRS,MIN_ISLANDS,MIN_BLOCKS=100,30,10
BANDS=(("all_distances",0.,float("inf")),
       ("less_than_2000km",0.,2000.),
       ("at_least_2000km",2000.,float("inf")))
REGIONS=("northern_midlatitude","tropical","southern_extratropical")


def source_endpoints(states:dict)->dict:
    if not set(TRAITS)<=set(states):raise ValueError("Original 10 Direct traits missing")
    result={t:states[t] for t in TRAITS}
    for a,b in PAIRS:
        result[a+"__"+b]={s:(states[a][s],states[b][s])
                          for s in set(states[a])&set(states[b])}
    if len(result)!=13:raise AssertionError("Main 10 + 3 endpoint family drift")
    return result


def exact_components(P:np.ndarray,Q:np.ndarray)->tuple[np.ndarray,...]:
    """Full exact categorical squared-distance algebra with signed cross."""
    if P.ndim!=2 or P.shape!=Q.shape or not np.isfinite(P+Q).all():
        raise ValueError("Observed category frequencies and genus mean shape invalid")
    E=P-Q
    d2=cdist(P,P,"sqeuclidean")
    genus=cdist(Q,Q,"sqeuclidean")
    residual=cdist(E,E,"sqeuclidean")
    # Cross term is the only potentially negative component.
    cross=d2-genus-residual
    # Independent direct check for > machine error.
    if len(P)>=2:
        i=np.triu_indices(min(len(P),40),1)
        dp=P[:40][i[0]]-P[:40][i[1]]
        dq=Q[:40][i[0]]-Q[:40][i[1]]
        de=E[:40][i[0]]-E[:40][i[1]]
        check=2*np.einsum("ij,ij->i",dq,de)
        if not np.allclose(check,cross[:40,:40][i],atol=1e-9):
            raise AssertionError("Genus/residual interaction algebra is wrong")
    if np.max(np.abs((genus+residual+cross)-d2))>1e-9:
        raise AssertionError("Category squared distance reconstruction fails")
    return d2,genus,residual,cross


class Aggregate:
    """Pair, island and physical-block-equal summaries of exact signed identity."""
    KEYS=("observed_D2","genus_mean_component_D2",
          "within_genus_residual_component_D2","signed_genus_residual_cross_D2")
    def __init__(self,block):
        b=np.asarray(block,dtype=int)
        if b.ndim!=1 or (b<0).any():raise ValueError("Unknown physical block")
        self.block=b
        self.n=len(b)
        self.ng=int(max(b)+1)
        self.count=0
        self.sums=np.zeros(4)
        self.node_count=np.zeros(self.n,int)
        self.node_sums=np.zeros((self.n,4))
        self.block_count=np.zeros(self.ng,int)
        self.block_sums=np.zeros((self.ng,4))
    def add(self,i,j,m):
        if len(i)==0:return
        if m.shape!=(len(i),4) or not np.isfinite(m).all():
            raise ValueError("Incomplete original pair algebra")
        if np.any(i==j) or np.any(self.block[i]==self.block[j]):
            raise ValueError("Reused island or identical block in cross-block comparison")
        if not np.allclose(m[:,0],m[:,1]+m[:,2]+m[:,3],atol=1e-9):
            raise AssertionError("Original pair additive decomposition fails")
        self.count+=len(i)
        self.sums+=m.sum(0)
        both=np.r_[i,j]
        values=np.concatenate((m,m),axis=0)
        self.node_count+=np.bincount(both,minlength=self.n)
        np.add.at(self.node_sums,both,values)
        b=self.block[both]
        self.block_count+=np.bincount(b,minlength=self.ng)
        np.add.at(self.block_sums,b,values)
    def summary(self,*,min_pairs=MIN_PAIRS,min_islands=MIN_ISLANDS,
                min_blocks=MIN_BLOCKS):
        active=self.node_count>0
        blocks=np.unique(self.block[active])
        if (self.count<min_pairs or active.sum()<min_islands or
                len(blocks)<min_blocks):
            return None
        pair=self.sums/self.count
        island=self.node_sums[active]/self.node_count[active,None]
        block=np.array([island[self.block[active]==b].mean(0)
                        for b in blocks])
        row=dict(n_species_disjoint_cross_block_pairs=self.count,
                 n_distinct_pair_end_islands=int(active.sum()),
                 n_distinct_physical_10degree_blocks=len(blocks))
        for idx,key in enumerate(self.KEYS):
            row["pair_equal_"+key]=float(pair[idx])
            row["island_equal_"+key]=float(island[:,idx].mean())
            row["physical_block_equal_"+key]=float(block[:,idx].mean())
            loo=[]
            for b in blocks:
                # Each accepted pair belongs to exactly two DIFFERENT blocks,
                # so its contribution is deleted once for each endpoint block.
                remaining=self.count-self.block_count[b]
                if remaining:
                    loo.append(float((self.sums[idx]-self.block_sums[b,idx])/remaining))
            row["delete_one_physical_block_pair_mean_"+key+"_min"]=(
                float(min(loo)) if loo else None)
            row["delete_one_physical_block_pair_mean_"+key+"_max"]=(
                float(max(loo)) if loo else None)
        for prefix in ("pair_equal_","island_equal_","physical_block_equal_"):
            terms=[row[prefix+k] for k in self.KEYS]
            if abs(terms[0]-sum(terms[1:]))>1e-8:
                raise AssertionError("Signed genus/within-genus reconciliation fails")
        row["pair_equal_genus_component_to_observed_ratio_NOT_causes"]=(
            float(pair[1]/pair[0]) if pair[0]>1e-10 else None)
        row["pair_equal_cross_term_may_be_negative"]=bool(pair[3]<0)
        return row


def analyse(states,flora,geo,*,batch=96,
            min_pairs=MIN_PAIRS,min_islands=MIN_ISLANDS,min_blocks=MIN_BLOCKS):
    if batch<1:raise ValueError("Invalid distance matrix chunk size")
    if not {"island_id","accepted_species","origin_status"}<=set(flora):
        raise ValueError("Original GIFT species incidence missing")
    if flora.duplicated(["island_id","accepted_species"]).any() or not (
        flora.origin_status.isin(("native","introduced","unresolved")).all()):
        raise ValueError("Original accepted plant/status identity invalid")
    need={"island_id","analysis_regime","spatial_block",
          "island_latitude","island_longitude"}
    if not need<=set(geo) or geo.island_id.duplicated().any():
        raise ValueError("Original physical geography invalid")
    full=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    native=flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    if len(full)!=4505:raise ValueError("The complete original GIFT flora changed")
    geo=geo.set_index("island_id")
    # Retain EVERY original GIFT flora's true GSHHG source regime.
    # Three original broad regimes remain primary; all additional,
    # independently labelled geographies are explicitly SECONDARY.
    region_of={
        i:(str(geo.at[i,"analysis_regime"])
           if i in geo.index and pd.notna(geo.at[i,"analysis_regime"])
           and str(geo.at[i,"analysis_regime"]).strip()
           else "geography_unresolved")
        for i in full}
    from collections import Counter
    region_counts=Counter(region_of.values())
    all_regions=list(REGIONS)+sorted(set(region_counts)-set(REGIONS))
    ep=source_endpoints(states)
    coverage={}
    results=[]
    eligibility=[]
    for scope,source in (("all_recorded",full),("original_native",native)):
        for name,phenotypes in ep.items():
            support=sum(len(source.get(i,set())&set(phenotypes))>=5 for i in full)
            coverage[scope+"__"+name]=int(support)
            for region in all_regions:
                region_ids=[i for i in sorted(full) if region_of[i]==region]
                region_roster=[source.get(i,set()) for i in region_ids]
                source_pool=set().union(*region_roster)&set(phenotypes)
                carrier=[source.get(i,set())&set(phenotypes) for i in region_ids]
                use=[k for k,v in enumerate(carrier) if len(v)>=5]
                entry=dict(scope=scope,trait_endpoint=name,region=region,
                    n_original_botanical_region_floras=len(region_ids),
                    n_original_direct_n5_islands=len(use),
                    n_in_region_DISTINCT_Direct_species=len(source_pool))
                if len(use)<2:
                    eligibility.append(entry)
                    continue
                original_ids=[region_ids[k] for k in use]
                use_roster=[carrier[k] for k in use]
                prepared=prepare_endpoint(use_roster,
                    {s:phenotypes[s] for s in source_pool})
                entry.update(
                    n_source_unique_genus_groups=prepared["n_global_genus_groups"],
                    n_source_variable_genus_groups=prepared["n_variable_genus_groups"])
                full_roster=[full[x] for x in original_ids]
                origin_species=sorted(set().union(*full_roster))
                original_A=membership(full_roster,origin_species)
                reg_geo=geo.reindex(original_ids)
                lat=pd.to_numeric(reg_geo.island_latitude,errors="coerce").to_numpy()
                lon=pd.to_numeric(reg_geo.island_longitude,errors="coerce").to_numpy()
                bcodes,blabels=pd.factorize(np.asarray(
                    reg_geo.spatial_block.to_numpy(),dtype=object))
                n=len(original_ids)
                if np.any(bcodes<0):
                    # Unresolved block islands remain in original eligibility
                    # but not in physical-block-different pair comparisons.
                    block_for_agg=np.where(bcodes<0,len(blabels),bcodes)
                else:
                    block_for_agg=bcodes
                sinks={label:Aggregate(block_for_agg) for label,_,_ in BANDS}
                P=np.asarray(prepared["P"],float)
                Q=np.asarray(prepared["genus_mean"],float)
                E=P-Q
                for start in range(0,n,batch):
                    stop=min(n,start+batch)
                    irow,jcol=np.nonzero(
                        np.arange(n)[None,:]>np.arange(start,stop)[:,None])
                    if not len(irow):continue
                    i=irow+start;j=jcol
                    shared=(original_A[start:stop]@original_A.T).toarray()
                    km=geographical_distance(lat,lon,start,stop)[irow,jcol]
                    mask=(shared[irow,jcol]==0)&(bcodes[i]>=0)&(
                        bcodes[j]>=0)&(bcodes[i]!=bcodes[j])&(
                        np.isfinite(km))&(
                        np.isfinite(lat[i])&np.isfinite(lat[j]))&(
                        np.isfinite(lon[i])&np.isfinite(lon[j]))&(
                        (np.abs(lat[i])<=90)&(np.abs(lat[j])<=90))&(
                        (np.abs(lon[i])<=180)&(np.abs(lon[j])<=180))
                    if not mask.any():continue
                    dp=P[i]-P[j]
                    dq=Q[i]-Q[j]
                    de=E[i]-E[j]
                    observed=np.einsum("ij,ij->i",dp,dp)
                    genus=np.einsum("ij,ij->i",dq,dq)
                    residual=np.einsum("ij,ij->i",de,de)
                    cross=2*np.einsum("ij,ij->i",dq,de)
                    metrics=np.column_stack([observed,genus,residual,cross])
                    if not np.allclose(observed,genus+residual+cross,
                                       atol=1e-9):
                        raise AssertionError("Original source pair identity fails")
                    for band,low,high in BANDS:
                        take=mask&(km>=low)&(km<high)
                        if take.any():
                            sinks[band].add(i[take],j[take],metrics[take])
                for band,_,_ in BANDS:
                    row=sinks[band].summary(
                        min_pairs=min_pairs,min_islands=min_islands,
                        min_blocks=min_blocks)
                    if row is not None:
                        results.append(dict(**entry,geographic_pair_band=band,
                                            **row))
                eligibility.append(entry)
    for coverage_key,expected_n in coverage.items():
        scope,endpoint=coverage_key.split("__",1)
        observed_n=sum(x["n_original_direct_n5_islands"]
                       for x in eligibility if x["scope"]==scope and
                       x["trait_endpoint"]==endpoint)
        if observed_n!=expected_n:
            raise AssertionError("Some original 4505 source floras were silently omitted")
    if (coverage["all_recorded__flower_primary_color__floral_form"]!=2324 or
        coverage["original_native__flower_primary_color__floral_form"]!=327):
        raise ValueError("Frozen 2324/327 original Direct paired endpoint drift")
    return dict(
        status="ORIGINAL_GLOBAL_13_ENDPOINT_GENUS_ASSEMBLY_EXACT_BETA_GEOMETRY",
        n_original_botanical_island_floras=len(full),
        n_actual_single_Direct_traits=10,
        n_actual_Direct_joint_two_trait_endpoints=3,
        original_endpoint_coverage_n5=coverage,
        source_region_genus_opportunity=eligibility,
        original_4505_botanical_GSHHG_region_universe=[
            dict(region=r,n_original_botanical_floras=int(region_counts[r]),
                 primary_preset_region=bool(r in REGIONS))
            for r in all_regions],
        original_exact_signed_category_distance_decompositions=results,
        claims=dict(
            pair_additive_squares_and_signed_cross_always_reconcile=True,
            genus_source_reference_is_in_sample_and_includes_focal_species=True,
            residual_fraction_is_not_percentage_of_adaptation=True,
            same_region_different_10degree_blocks=True,
            all_original_full_flora_accepted_species_disjointness_required=True,
            original_all_4505_islands_retained_despite_geo_unresolved=True,
            no_old_928_island_complete_case_criterion=True,
            original_13_trait_endpoints_no_unmeasured_cross_products=True,
            source_native_and_all_recorded_nested_not_independent=True,
            dyad_counts_not_independent_evolutionary_replications=True,
            no_real_pollinator_shift_within_species_adaptation_or_ancestral_colonization=True,
            no_private_island_species_genus_GIFT_source_geographic_ids_exported=True))


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for k in ("original-root","geology","output"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    source={
        "traits":a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "flora":a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz"}
    pins={
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a"}
    for key,path in source.items():
        if sha(path)!=pins[key]:
            raise ValueError("Original 4505 source species-direct archive changed: "+key)
    states=strict_traits(pd.read_csv(source["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(source["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    geo=pd.read_csv(a.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block",
        "island_latitude","island_longitude"])
    if geo.island_id.nunique()!=8264:
        raise ValueError("Original GSHHG 8264 physical island frame changed")
    report=analyse(states,flora,geo)
    payload=json.dumps(report,indent=2,sort_keys=True,allow_nan=False)
    for word in ('"island_id"','"accepted_species"','"entity_ID"',
                 '"spatial_block"','"source_signature"',
                 '"island_latitude"','"island_longitude"'):
        if word in payload:raise ValueError("Original botanical identity leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(payload+"\n",encoding="utf-8")
    print(json.dumps(dict(status=report["status"],
                          n_original_floras=4505,
                          n_anonymous_decomposition_rows=len(
                            report["original_exact_signed_category_distance_decompositions"]))))


if __name__=="__main__":
    main()
