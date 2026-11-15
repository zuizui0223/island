"""Outcome-blind 10degree-block-disjoint matching for exploratory floral beta scale signal.

The primary question was selected AFTER the 15-Nov geographical reference
diagnostic: tropical, 500–2,000km, zero shared full-flora species, species
Direct colour+form, different 10degree spatial blocks. This is POST HOC.
All 4,505 original botanical islands remain in source universe; no simulated
trait outcomes, no pollination effectiveness or evolved syndrome assignment.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx

from fd_island_functional_beta import (
    strict_traits, prepare_endpoint, compare_block,
    geographical_distance, membership,
)
from fd_island_beta_regional_genus import _regional_donor_audit

N_DRAWS=199
MIN_PAIRS,MIN_ISLANDS,MIN_BLOCKS=100,30,10
STRATA=(("tropical","500_to_2000km",500.,2000.),
        ("northern_midlatitude","500_to_2000km",500.,2000.))
GATES=("n5","both_n10","both_direct_fraction_ge_0p10")


def match_blocks(block_a, block_b, *, seed:int=1)->np.ndarray:
    """Select ONE original dyad per block-disjoint edge, WITHOUT outcome access.

    Maximum cardinality is computed on the spatial-block graph, then one
    eligible original island pair within each selected block pair is sampled.
    All tie weights are seeded independently from the floral response.
    """
    a=np.asarray(block_a)
    b=np.asarray(block_b)
    if len(a)!=len(b):raise ValueError("Unpaired edge endpoints")
    if not len(a):return np.array([],dtype=int)
    if np.any(a==b):raise ValueError("Same physical block not allowed")
    candidates={}
    for i,(x,y) in enumerate(zip(a,b)):
        edge=tuple(sorted((int(x),int(y))))
        candidates.setdefault(edge,[]).append(i)
    rng=np.random.default_rng(seed)
    graph=nx.Graph()
    for edge in sorted(candidates):
        graph.add_edge(*edge,weight=float(rng.uniform(0.001,1.)))
    selected=nx.algorithms.matching.max_weight_matching(
        graph,maxcardinality=True,weight="weight")
    picks=[]
    for edge in sorted(tuple(sorted((int(x),int(y)))) for x,y in selected):
        entries=candidates[edge]
        picks.append(int(entries[int(rng.integers(len(entries)))]))
    picks.sort()
    chosen=np.array(picks,dtype=int)
    endpoints=np.r_[a[chosen],b[chosen]]
    if len(endpoints)!=len(set(endpoints)):
        raise AssertionError("A geographic block reused in a match")
    return chosen


def graph_support(block_a,block_b)->dict:
    g=nx.Graph()
    g.add_edges_from(zip(block_a,block_b))
    comp=sorted((len(x) for x in nx.connected_components(g)),reverse=True)
    return dict(n_candidate_block_pairs=g.number_of_edges(),
                n_candidate_blocks=g.number_of_nodes(),
                n_connected_spatial_block_components=len(comp),
                largest_component_block_count=comp[0] if comp else 0)


def summarize(edges,*,draws=N_DRAWS,min_pairs=MIN_PAIRS,
              min_islands=MIN_ISLANDS,min_blocks=MIN_BLOCKS):
    """Only supports and model-residual summaries are exported, never island IDs."""
    required={"i","j","block_i","block_j","obs","expected_g","expected_r"}
    if not required<=set(edges):raise ValueError("Incomplete exact original dyads")
    if not len(edges["i"]):return None
    if draws<1:raise ValueError("Must assess ≥1 original spatial reselection")
    n_pairs=len(edges["i"])
    active_islands=set(np.r_[edges["i"],edges["j"]])
    block_a=np.asarray(edges["block_i"],int)
    block_b=np.asarray(edges["block_j"],int)
    if len(set(np.r_[block_a,block_b]))<min_blocks or n_pairs<min_pairs or len(active_islands)<min_islands:
        return None
    obs=np.asarray(edges["obs"])
    expected_g=np.asarray(edges["expected_g"])
    expected_r=np.asarray(edges["expected_r"])
    if not all(np.isfinite(x).all() for x in (obs,expected_g,expected_r)):
        raise ValueError("Non-finite original floral D² response")
    residual_g=obs-expected_g
    residual_r=obs-expected_r
    support=graph_support(block_a,block_b)
    samples=[]
    sizes=[]
    for seed in range(draws):
        chosen=match_blocks(block_a,block_b,seed=20261115+seed)
        sizes.append(len(chosen))
        samples.append((float(np.mean(residual_g[chosen])),
                        float(np.mean(residual_r[chosen])),
                        float(np.mean((residual_r-residual_g)[chosen]))))
    samples=np.asarray(samples)
    if np.any(~np.isfinite(samples)):raise AssertionError("Empty independent matching")
    return dict(
        n_candidate_pairs=int(n_pairs),
        n_eligible_distinct_islands=len(active_islands),
        n_eligible_distinct_spatial_blocks=len(set(np.r_[block_a,block_b])),
        **support,
        candidate_pair_mean_global_residual=float(np.mean(residual_g)),
        candidate_pair_mean_regional_residual=float(np.mean(residual_r)),
        candidate_pair_mean_region_minus_global=float(np.mean(residual_r-residual_g)),
        n_outcome_blind_matching_variants=draws,
        matched_block_disjoint_pairs_min=int(min(sizes)),
        matched_block_disjoint_pairs_max=int(max(sizes)),
        matched_block_disjoint_pairs_median=float(np.median(sizes)),
        matched_global_residual_median=float(np.median(samples[:,0])),
        matched_regional_residual_median=float(np.median(samples[:,1])),
        matched_regional_residual_p05=float(np.quantile(samples[:,1],.05)),
        matched_regional_residual_p95=float(np.quantile(samples[:,1],.95)),
        matched_regional_mean_positive_fraction=float(np.mean(samples[:,1]>0)),
        matched_reference_shift_median=float(np.median(samples[:,2])),
        selected_dyad_margins_not_statistical_independent_replications=True,
    )


def analyze(states, flora, geo, *, draws=N_DRAWS,
            min_pairs=MIN_PAIRS,min_islands=MIN_ISLANDS,min_blocks=MIN_BLOCKS):
    if not {'island_id','accepted_species','origin_status'}<=set(flora):
        raise ValueError("Missing original flora")
    if flora.duplicated(['island_id','accepted_species']).any() or (
            not flora.origin_status.isin(('native','introduced','unresolved')).all()):
        raise ValueError("Source botanical roster or origins invalid")
    if not {'island_id','analysis_regime','spatial_block',
            'island_latitude','island_longitude'}<=set(geo) or geo.island_id.duplicated().any():
        raise ValueError("Independent geography invalid")
    all_f=flora.groupby('island_id').accepted_species.agg(set).to_dict()
    native=flora.loc[flora.origin_status.eq('native')].groupby(
        'island_id').accepted_species.agg(set).to_dict()
    islands=sorted(all_f)
    joined=geo.set_index('island_id').reindex(islands)
    colour=states['flower_primary_color']
    form=states['floral_form']
    both=set(colour)&set(form)
    states_pair={sp:(colour[sp],form[sp]) for sp in both}
    out=[]
    coverage={}
    for scope,flora_scope in (('all_recorded',all_f),('original_native',native)):
        eligible=[i for i in islands if len(flora_scope.get(i,set())&both)>=5]
        coverage[scope]=len(eligible)
        for region,band,low,high in STRATA:
            filtered=[i for i in eligible
                      if joined.at[i,'analysis_regime']==region]
            if len(filtered)<2:continue
            small_rosters=[flora_scope.get(i,set())&both for i in filtered]
            used=set().union(*small_rosters)
            direct_reference={sp:states_pair[sp] for sp in used}
            ref=prepare_endpoint(small_rosters,direct_reference)
            source_audit=_regional_donor_audit(states_pair,used)
            n=len(filtered)
            # All-original-flora species-sharing remains the exclusion unit,
            # even for a native-only morphological endpoint.
            source_full=[all_f[i] for i in filtered]
            labels=sorted(set().union(*source_full))
            Full=membership(source_full,labels)
            shares=(Full@Full.T).toarray()
            idx_geo=joined.reindex(filtered)
            latitude=pd.to_numeric(idx_geo.island_latitude,errors='coerce').to_numpy()
            longitude=pd.to_numeric(idx_geo.island_longitude,errors='coerce').to_numpy()
            block_id,unique_block=pd.factorize(idx_geo.spatial_block.to_numpy())
            n_direct=np.array([len(r) for r in small_rosters])
            denom=np.array([len(flora_scope.get(i,set())) for i in filtered])
            frac=n_direct/denom
            mat=compare_block(ref,0,n)
            observed=mat['squared_frequency_distance']
            regional_expected=mat['genus_null_squared_distance']
            # Global reference includes all ORIGINAL Direct trait species, not
            # just focal-region species or the names in selected island pairs.
            gp=prepare_endpoint(small_rosters,states_pair)
            global_expected=compare_block(gp,0,n)['genus_null_squared_distance']
            distance=geographical_distance(latitude,longitude,0,n)
            left,right=np.triu_indices(n,1)
            eligible_graph=(
                (shares[left,right]==0)&
                (block_id[left]>=0)&(block_id[right]>=0)&
                (block_id[left]!=block_id[right])&
                np.isfinite(distance[left,right])&
                (distance[left,right]>=low)&(distance[left,right]<high)&
                np.isfinite(latitude[left])&np.isfinite(latitude[right])&
                np.isfinite(longitude[left])&np.isfinite(longitude[right])&
                (np.abs(latitude[left])<=90)&(np.abs(latitude[right])<=90)&
                (np.abs(longitude[left])<=180)&(np.abs(longitude[right])<=180))
            if not eligible_graph.any():continue
            for gate in GATES:
                keep=eligible_graph.copy()
                if gate=='both_n10':
                    keep&=(n_direct[left]>=10)&(n_direct[right]>=10)
                elif gate=='both_direct_fraction_ge_0p10':
                    keep&=(frac[left]>=.10)&(frac[right]>=.10)
                if not keep.any():continue
                ii,jj=left[keep],right[keep]
                ed=dict(i=ii,j=jj,
                        block_i=block_id[ii],block_j=block_id[jj],
                        obs=observed[ii,jj],expected_g=global_expected[ii,jj],
                        expected_r=regional_expected[ii,jj])
                report=summarize(ed,draws=draws,min_pairs=min_pairs,
                                 min_islands=min_islands,min_blocks=min_blocks)
                if report:
                    out.append(dict(scope=scope,region=region,
                                    distance_class=band,support=gate,
                                    **source_audit,**report))
    return dict(
        status="FLORAL_BETA_BLOCK_DISJOINT_SPATIAL_STRESS",
        n_original_botanical_islands=len(islands),
        eligible_n5_direct_islands_by_scope=coverage,
        comparisons=out,
        claims=dict(
            primary_tropical_intermediate_distance_chosen_posthoc=True,
            comparison_requires_full_flora_species_disjointness=True,
            selected_island_dyads_are_unique_in_physical_blocks=True,
            selection_uses_no_flower_trait_outcomes=True,
            199_matching_reselections_not_independent_geographic_histories=True,
            regional_reference_is_in_sample_genus_pool=True,
            field_pollination_and_local_adaptation_not_identified=True,
            no_private_original_taxon_island_block_ids_exported=True,
        ))


def main():
    p=argparse.ArgumentParser()
    for key in ('original-root','geology','output'):
        p.add_argument('--'+key,type=Path,required=True)
    args=p.parse_args()
    source={
        'traits':args.original_root/'fixed/source/current/evidence/direct_only_trait_ledger.csv.gz',
        'flora':args.original_root/'fixed/canonical/input/chapter1_status_flora.csv.gz'}
    expected={
        'traits':'9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745',
        'flora':'6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a'}
    for k,f in source.items():
        h=hashlib.sha256()
        with f.open('rb') as fp:
            for chunk in iter(lambda:fp.read(1<<20),b''):h.update(chunk)
        if h.hexdigest()!=expected[k]:raise ValueError('Source original SHA mismatch')
    traits=pd.read_csv(source['traits'],dtype=str,usecols=[
        'accepted_species','trait_name','normalized_value',
        'evidence_scope','evidence_quality'])
    flora=pd.read_csv(source['flora'],dtype=str,usecols=[
        'island_id','accepted_species','origin_status'])
    geo=pd.read_csv(args.geology,dtype=str)
    if flora.island_id.nunique()!=4505:raise ValueError('Not full botanical frame')
    result=analyze(strict_traits(traits),flora,geo)
    if result['eligible_n5_direct_islands_by_scope']['all_recorded']!=2324:
        raise ValueError('Source direct colour/form opportunity changed')
    # Pair count is a qualification identity: audit the exact original cohort
    # before allowing any descriptive inference.
    primary=[x for x in result['comparisons'] if x['scope']=='all_recorded'
             and x['region']=='tropical' and x['support']=='n5']
    if len(primary)!=1 or primary[0]['n_candidate_pairs']!=318:
        raise ValueError('Frozen tropical same-region 500–2000 km pair universe changed')
    payload=json.dumps(result,sort_keys=True,indent=2,allow_nan=False)
    for s in ('"island_id"','"accepted_species"','"spatial_block"',
              '"island_latitude"','"island_longitude"','"source_signature"'):
        if s in payload:raise ValueError('Private original entity key in output')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(payload+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],
                          n_anonymous_comparisons=len(result['comparisons']),
                          biological_replication_not_assumed=True)))


if __name__=='__main__':
    main()
