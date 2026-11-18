"""Island floristic composition vs recorded floral-category beta-diversity.

10 independent high/medium Direct traits and 3 observed two-trait endpoints.
A plant's categorical trait is species-wide; none of these comparisons
measure population evolution, adaptive shifts, or effective pollination.
Only anonymous island-pair category aggregates can leave original CI.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import cdist

TRAITS = ('flower_primary_color','floral_form','floral_symmetry',
          'flower_size_class','tube_depth_class','inflorescence_display',
          'self_incompatibility','autonomous_selfing_capacity',
          'mating_system','cleistogamy')
PAIRS = (('flower_primary_color','floral_form'),
         ('flower_primary_color','floral_symmetry'),
         ('floral_form','floral_symmetry'))
MIN_DIRECT = 5

def strict_traits(ledger: pd.DataFrame) -> dict:
    fields = {'accepted_species','trait_name','normalized_value',
              'evidence_scope','evidence_quality'}
    if not fields <= set(ledger):
        raise ValueError('Missing original Direct source fields')
    q = ledger.loc[ledger.trait_name.isin(TRAITS) &
                   ledger.evidence_scope.eq('species_direct') &
                   ledger.evidence_quality.isin(('high','medium')) &
                   ledger.normalized_value.notna()].copy()
    q['normalized_value'] = q.normalized_value.astype(str).str.strip()
    q = q.loc[q.normalized_value.ne('') &
              ~q.normalized_value.str.contains('|',regex=False) &
              ~q.normalized_value.str.lower().isin(('none','nan','null','unknown','unresolved'))]
    if q.accepted_species.isna().any() or q.duplicated(['accepted_species','trait_name']).any():
        raise ValueError('Unusable or duplicated accepted species×trait source row')
    return {trait:q.loc[q.trait_name.eq(trait)].set_index('accepted_species')
            .normalized_value.to_dict() for trait in TRAITS}

def incidence_partition(shared, left, right):
    """Sørensen dissimilarity = Baselga Simpson turnover + nestedness-resultant.

    Categorical inclusion is not an inference of historical functional loss.
    """
    a,l,r = np.broadcast_arrays(np.asarray(shared,float),np.asarray(left,float),
                                np.asarray(right,float))
    if (not np.isfinite(a+l+r).all() or np.any(l<=0) or np.any(r<=0) or
        np.any(a<0) or np.any(a>np.minimum(l,r)+1e-9)):
        raise ValueError('Invalid nonempty original species/category intersections')
    total = np.clip(1-2*a/(l+r),0,1)
    turnover = np.clip(1-a/np.minimum(l,r),0,1)
    nestedness = np.maximum(0,total-turnover)
    return total,turnover,nestedness

def membership(rosters:list[set], species:list[str]) -> sparse.csr_matrix:
    lookup = {name:i for i,name in enumerate(species)}
    rows=[];cols=[]
    for i,roster in enumerate(rosters):
        for name in sorted(roster):
            if name in lookup:
                rows.append(i);cols.append(lookup[name])
    return sparse.csr_matrix((np.ones(len(rows)),(rows,cols)),
                             shape=(len(rosters),len(species)))

def accepted_botanical_genus(name: str) -> str | None:
    """Quarantine unresolved/non-binomial accepted-name records for strict genus exclusion."""
    if not isinstance(name, str):
        return None
    parts = name.strip().split()
    if (len(parts)<2 or not parts[0][0].isupper()
            or not parts[0].replace("-", "").isalpha()
            or parts[1].lower() in {"sp.", "spp.", "cf.", "aff.", "?"}):
        return None
    return parts[0]

def certified_genus_sets(rosters: list[set]) -> tuple[list[set], np.ndarray]:
    """Full-roster certification: unknown genus on either side cannot prove non-sharing."""
    genera, complete = [], []
    for roster in rosters:
        named = [accepted_botanical_genus(s) for s in roster]
        genera.append({g for g in named if g is not None})
        complete.append(bool(roster) and all(g is not None for g in named))
    return genera, np.asarray(complete, dtype=bool)

def prepare_endpoint(rosters:list[set], states:dict) -> dict:
    """Original species masks and original recorded category combinations.

    One observed two-trait category state is one actual observation; no
    imagined cross-product of mutually exclusive/missing trait states.
    """
    if not states or not rosters:raise ValueError('Source states and rosters required')
    source = set(states)
    if any(not isinstance(r,set) or not r or not r<=source for r in rosters):
        raise ValueError('Every included island has real directly recorded species')
    species=sorted(states)
    categories=sorted(set(states.values()),key=repr)
    cat={c:i for i,c in enumerate(categories)}
    code=np.array([cat[states[s]] for s in species],int)
    A=membership(rosters,species)
    n=np.asarray(A.sum(1)).ravel()
    C=sparse.csr_matrix((np.ones(len(species)),(np.arange(len(species)),code)),
                        shape=(len(species),len(categories)))
    counts=(A@C).toarray()
    P=counts/n[:,None]
    totals=np.bincount(code,minlength=len(categories))
    global_h=1-np.square(totals/len(species)).sum()
    genera={};gids=[]
    for sp in species:
        chunks=sp.split()
        key=('genus',chunks[0]) if len(chunks)>=2 else ('fixed_original_species',sp)
        if key not in genera:genera[key]=len(genera)
        gids.append(genera[key])
    gids=np.asarray(gids,int)
    G=sparse.csr_matrix((np.ones(len(species)),(np.arange(len(species)),gids)),
                        shape=(len(species),len(genera)))
    table=(G.T@C).toarray()
    sizes=table.sum(1)
    gp=table/sizes[:,None]
    H=1-np.square(gp).sum(1)
    variance=np.zeros(len(sizes))
    np.divide(H,sizes-1,out=variance,where=sizes>1)
    w=sizes*variance
    B=(A@G).tocsr()
    mean=np.asarray(B@gp)/n[:,None]
    trace=(np.asarray(A@w[gids]).ravel()-
           np.asarray(B.power(2)@variance).ravel())/n**2
    return dict(A=A,B=B,n=n,counts=counts,P=P,categories=categories,
                global_n=len(species),global_h=float(global_h),
                genus_mean=mean,genus_trace=np.maximum(trace,0),
                species_variance_weight=w[gids],
                genus_variance_weight=variance,
                n_global_genus_groups=len(genera),
                n_variable_genus_groups=int(np.sum(np.count_nonzero(table,axis=1)>1)))

def compare_block(prepared:dict,start:int,stop:int) -> dict:
    """Exact expected squared categorical-frequency dissimilarity.

    Null labels use either the observed global categorical frequencies or
    permute entire observed categorical profiles within documented genus.
    Both nulls preserve the ORIGINAL species-sharing mask across islands.
    """
    z=prepared;A=z['A'];B=z['B'];n=z['n'];P=z['P']
    if not 0<=start<stop<=len(n):raise ValueError('Invalid original row block')
    shared=(A[start:stop]@A.T).toarray()
    denom=n[start:stop,None]*n[None,:]
    tax,_,_=incidence_partition(shared,n[start:stop,None],n[None,:])
    presence=(z['counts']>0).astype(float)
    cat_shared=presence[start:stop]@presence.T
    nc=presence.sum(1)
    sor,turn,nest=incidence_partition(cat_shared,nc[start:stop,None],nc[None,:])
    squared=cdist(P[start:stop],P,metric='sqeuclidean')
    expectation=np.zeros_like(squared)
    if z['global_n']>1:
        expectation=z['global_h']*z['global_n']/(z['global_n']-1)*(
            1/n[start:stop,None]+1/n[None,:]-2*shared/denom)
    covariance=((A[start:stop].multiply(z['species_variance_weight'])@A.T).toarray()-
                (B[start:stop].multiply(z['genus_variance_weight'])@B.T).toarray())/denom
    genus=cdist(z['genus_mean'][start:stop],z['genus_mean'],
                metric='sqeuclidean')
    genus+=z['genus_trace'][start:stop,None]+z['genus_trace'][None,:]-2*covariance
    if np.min(genus)<-1e-9:raise ArithmeticError('Negative genus expectation')
    return dict(overlap=np.clip(1-.5*cdist(P[start:stop],P,metric='cityblock'),0,1),
        observed_taxonomic_sorensen=tax,observed_shared_species=shared,
        category_sorensen=sor,category_turnover=turn,
        category_nestedness=nest,squared_frequency_distance=squared,
        global_null_squared_distance=np.maximum(expectation,0),
        genus_null_squared_distance=np.maximum(genus,0))

class PairSummary:
    def __init__(self,n:int):
        self.n_pairs=0;self.sums={};self.counts={}
        self.node_count=np.zeros(n,np.int64)
        self.node_overlap=np.zeros(n);self.node_tax=np.zeros(n)
    def add(self,i,j,metrics):
        if len(i)==0:return
        self.n_pairs+=len(i)
        both=np.r_[i,j]
        self.node_count+=np.bincount(both,minlength=len(self.node_count))
        self.node_overlap+=np.bincount(both,
            weights=np.tile(metrics['functional_overlap'],2),
            minlength=len(self.node_count))
        self.node_tax+=np.bincount(both,
            weights=np.tile(metrics['taxonomic_similarity_full'],2),
            minlength=len(self.node_count))
        for name,values in metrics.items():
            valid=np.isfinite(values)
            self.sums[name]=self.sums.get(name,0.)+float(np.sum(values[valid]))
            self.counts[name]=self.counts.get(name,0)+int(valid.sum())
    def export(self,blocks):
        active=self.node_count>0
        out=dict(n_pairs=int(self.n_pairs),n_unique_islands=int(active.sum()),
            n_unique_geographic_blocks=int(len({str(b) for b in
                np.asarray(blocks,dtype=object)[active] if pd.notna(b)})),
            island_equal_mean_functional_overlap=float(
                np.mean(self.node_overlap[active]/self.node_count[active])),
            island_equal_mean_taxonomic_similarity_full=float(
                np.mean(self.node_tax[active]/self.node_count[active])))
        for name,total in self.sums.items():
            out['mean_'+name]=total/self.counts[name] if self.counts[name] else None
        return out

def geographical_distance(lat,lon,start,stop):
    lat=np.radians(lat);lon=np.radians(lon)
    dlat=lat[start:stop,None]-lat[None,:]
    dlon=lon[start:stop,None]-lon[None,:]
    a=np.sin(dlat/2)**2+np.cos(lat[start:stop,None])*np.cos(
        lat[None,:])*np.sin(dlon/2)**2
    return 2*6371.0088*np.arcsin(np.sqrt(np.clip(a,0,1)))

def analyse(states:dict,flora:pd.DataFrame,geography:pd.DataFrame,*,
            pairs=PAIRS,block_size=96) -> dict:
    need={'island_id','accepted_species','origin_status'}
    if not need<=set(flora) or flora.duplicated(['island_id','accepted_species']).any():
        raise ValueError('Original accepted island flora must be unique')
    if flora[list(need)].isna().any().any() or not flora.origin_status.isin(
        ('native','introduced','unresolved')).all():
        raise ValueError('Original source status missing or inconsistent')
    geo_cols={'island_id','analysis_regime','geological_origin','spatial_block'}
    if not geo_cols<=set(geography) or geography.island_id.duplicated().any():
        raise ValueError('Independent published island source geography invalid')
    if block_size<1:raise ValueError('Block size must be positive')
    islands=sorted(flora.island_id.unique())
    total_islands=len(islands)
    geo=geography.set_index('island_id').reindex(islands)
    missing=int(geo.analysis_regime.isna().sum())
    regions=geo.analysis_regime.fillna('geography_unresolved').astype(str).to_numpy()
    origins=geo.geological_origin.fillna('unresolved').astype(str).to_numpy()
    blocks=geo.spatial_block.to_numpy(dtype=object)
    endpoints=list(states.items())
    for a,b in pairs:
        if a not in states or b not in states:raise ValueError('Unknown pair Direct source')
        common=set(states[a])&set(states[b])
        endpoints.append((a+'__'+b,{
            k:(states[a][k],states[b][k]) for k in common}))
    coverage=[];summaries=[];composition=[]
    allmap=flora.groupby('island_id').accepted_species.agg(set).to_dict()
    native=flora.loc[flora.origin_status.eq('native')].groupby(
        'island_id').accepted_species.agg(set).to_dict()
    for scope,original in (('all_recorded',allmap),('original_native',native)):
        rosters=[original.get(i,set()) for i in islands]
        allspecies=sorted(set().union(*rosters))
        fullA=membership(rosters,allspecies)
        fullN=np.asarray(fullA.sum(1)).ravel()
        # Source accepted species identities only: no phylogenetic or genus imputation.
        full_genus_rosters, full_genus_certified = certified_genus_sets(rosters)
        full_genus_names = sorted(set().union(*full_genus_rosters))
        fullG = membership(full_genus_rosters, full_genus_names)
        for endpoint,st in endpoints:
            names=set(st)
            candidate=[r&names for r in rosters]
            n=np.array([len(r) for r in candidate])
            use=np.flatnonzero(n>=MIN_DIRECT)
            nn=len(use)
            record=dict(scope=scope,endpoint=endpoint,
                n_original_islands=total_islands,
                n_any_direct_islands=int((n>0).sum()),
                n_n5_islands=nn,
                n_n10_islands=int((n>=10).sum()),
                n_n30_islands=int((n>=30).sum()),
                n_empty_original_scope_islands=int((fullN==0).sum()),
                n_unordered_n5_pairs=nn*(nn-1)//2,
                n_global_direct_endpoint_species=len(st))
            coverage.append(record)
            if nn<2:continue
            p=prepare_endpoint([candidate[i] for i in use],st)
            direct_genus_rosters,direct_genus_certified=certified_genus_sets(
                [candidate[i] for i in use])
            directG=membership(direct_genus_rosters,
                sorted(set().union(*direct_genus_rosters)))
            FG=fullG[use].tocsr()
            certified_full=full_genus_certified[use]
            record['n_certified_full_flora_genus_islands']=int(certified_full.sum())
            record['n_certified_direct_genus_islands']=int(
                direct_genus_certified.sum())
            record['n_reference_genera']=p['n_global_genus_groups']
            record['n_reference_genera_with_category_variation']=p['n_variable_genus_groups']
            reg=regions[use];org=origins[use];block=blocks[use];obsN=n[use]
            F=fullA[use];FN=fullN[use]
            for region in ['GLOBAL_ALL',*sorted(set(reg)-{'geography_unresolved'})]:
                rm=np.ones(nn,bool) if region=='GLOBAL_ALL' else reg==region
                for origin in ('all','oceanic','continental'):
                    mask=rm if origin=='all' else rm&(org==origin)
                    if mask.sum()<2:continue
                    props=p['P'][mask].mean(0)
                    composition.append(dict(scope=scope,endpoint=endpoint,
                        region=region,geological_origin=origin,
                        n_islands=int(mask.sum()),
                        categories=[dict(
                            state=list(c) if isinstance(c,tuple) else c,
                            island_equal_mean_recorded_species_fraction=float(z))
                          for c,z in zip(p['categories'],props)]))
            # Streaming upper triangle prevents exporting original pair IDs.
            collectors={}
            for start in range(0,nn,block_size):
                stop=min(nn,start+block_size)
                mat=compare_block(p,start,stop)
                fullshared=(F[start:stop]@F.T).toarray()
                fullgenusshared=(FG[start:stop]@FG.T).toarray()
                directgenusshared=(directG[start:stop]@directG.T).toarray()
                fulltax=incidence_partition(fullshared,FN[start:stop,None],
                                            FN[None,:])[0]
                ii,jj=np.nonzero(np.arange(nn)[None,:] >
                                 np.arange(start,stop)[:,None])
                I=ii+start;J=jj
                if not len(I):continue
                metrics={
                    'functional_overlap':mat['overlap'][ii,jj],
                    'taxonomic_similarity_full':1-fulltax[ii,jj],
                    'taxonomic_similarity_observed_mask':
                         1-mat['observed_taxonomic_sorensen'][ii,jj],
                    'category_sorensen':mat['category_sorensen'][ii,jj],
                    'category_turnover':mat['category_turnover'][ii,jj],
                    'category_nestedness':mat['category_nestedness'][ii,jj],
                    'squared_frequency_distance':mat['squared_frequency_distance'][ii,jj],
                    'global_null_squared_distance':mat['global_null_squared_distance'][ii,jj],
                    'genus_null_squared_distance':mat['genus_null_squared_distance'][ii,jj]
                }
                metrics['pair_recorded_species_fraction']=.5*(
                    obsN[I]/FN[I]+obsN[J]/FN[J])
                metrics['absolute_recording_fraction_difference']=np.abs(
                    obsN[I]/FN[I]-obsN[J]/FN[J])
                metrics['global_null_excess_squared_distance']=(
                    metrics['squared_frequency_distance']-
                    metrics['global_null_squared_distance'])
                metrics['genus_null_excess_squared_distance']=(
                    metrics['squared_frequency_distance']-
                    metrics['genus_null_squared_distance'])
                metrics['genus_null_expected_distance_positive']=(
                    metrics['genus_null_squared_distance']>1e-12).astype(float)
                if {'island_latitude','island_longitude'}<=set(geo):
                    lat=pd.to_numeric(geo.island_latitude,errors='coerce').to_numpy()[use]
                    lon=pd.to_numeric(geo.island_longitude,errors='coerce').to_numpy()[use]
                    metrics['geographic_distance_km']=geographical_distance(
                        lat,lon,start,stop)[ii,jj]
                no_shared=fullshared[ii,jj]==0
                full_genus_same=fullgenusshared[ii,jj]>0
                direct_genus_same=directgenusshared[ii,jj]>0
                certified=(certified_full[I]&certified_full[J])
                certified_direct=(direct_genus_certified[I]&
                                  direct_genus_certified[J])
                no_shared_full_genus=certified&~full_genus_same
                no_shared_direct_genus=certified_direct&~direct_genus_same
                distinct_blocks=pd.notna(block[I])&pd.notna(block[J]) & (
                    block[I]!=block[J])
                sample30=(obsN[I]>=30)&(obsN[J]>=30)
                type_masks={
                    'all':np.ones(len(I),bool),
                    'oceanic_oceanic':(org[I]=='oceanic')&(org[J]=='oceanic'),
                    'continental_continental':(org[I]=='continental')&
                                               (org[J]=='continental'),
                    'oceanic_continental':(
                      ((org[I]=='oceanic')&(org[J]=='continental')) |
                      ((org[I]=='continental')&(org[J]=='oceanic')))
                }
                for region in ['GLOBAL_ALL',*sorted(set(reg)-{'geography_unresolved'})]:
                    rm=np.ones(len(I),bool) if region=='GLOBAL_ALL' else (
                        (reg[I]==region)&(reg[J]==region))
                    for kind,typ in type_masks.items():
                        base=rm&typ
                        if not base.any():continue
                        for subset,extra in (
                            ('all',np.ones(len(I),bool)),
                            ('no_shared_full_flora_species',no_shared),
                            ('no_shared_full_flora_species_genus_shared',
                             no_shared&certified&full_genus_same),
                            ('no_shared_full_flora_genera_strict',
                             no_shared_full_genus),
                            ('no_shared_direct_genera',
                             no_shared_direct_genus),
                            ('genus_disjoint_different_blocks',
                             no_shared_full_genus&distinct_blocks),
                            ('both_n30_direct',sample30),
                            ('different_spatial_blocks',distinct_blocks)):
                            keep=base&extra
                            if not keep.any():continue
                            key=(region,kind,subset)
                            if key not in collectors:
                                collectors[key]=PairSummary(nn)
                            collectors[key].add(I[keep],J[keep],{
                                key:v[keep] for key,v in metrics.items()})
            for (region,kind,subset),agg in sorted(collectors.items()):
                summaries.append(dict(scope=scope,endpoint=endpoint,
                    region=region,geology_pair=kind,subset=subset,
                    **agg.export(block)))
            global_row=next(x for x in summaries
                if x['scope']==scope and x['endpoint']==endpoint and
                   x['region']=='GLOBAL_ALL' and
                   x['geology_pair']=='all' and x['subset']=='all')
            if global_row['n_pairs']!=nn*(nn-1)//2:
                raise AssertionError('Original complete upper triangle lost pairs')
    return dict(status='ISLAND_FUNCTIONAL_COMPOSITION_BETA_DESCRIPTIVE',
        n_original_islands=total_islands,n_without_geography=missing,
        n_geologically_typed_original_flora=int(np.isin(
            origins,('oceanic','continental')).sum()),
        endpoint_coverage=coverage,pair_summaries=summaries,
        regional_recorded_trait_composition=composition,
        claims=dict(no_complete_three_trait_gate=True,
            no_missing_trait_imputation=True,
            no_intraspecific_trait_change_measured=True,
            no_pair_count_treated_as_independent_replications=True,
            overlap_uses_species_frequency_not_individual_abundance=True,
            taxonomic_full_and_same_trait_mask_both_reported=True,
            pairwise_trait_states_are_observed_combinations_only=True,
            genus_null_permutates_whole_observed_profile_for_paired_endpoints=True,
            null_expectations_condition_on_shared_species_and_fixed_recording_mask=True,
            genus_null_is_global_recorded_pool_not_ancestral_source=True,
            genus_disjoint_uses_entire_full_flora_not_only_direct_species=True,
            unresolved_accepted_genus_quarantined_not_assumed_distinct=True,
            genus_disjoint_not_phylogenetic_evolutionary_independence=True,
            no_null_pvalues_or_confidence_intervals_inferred_from_expectation=True,
            geological_groups_descriptive_not_climate_or_area_adjusted=True,
            category_turnover_and_nestedness_not_historical_replacement_or_loss=True,
            no_continuous_hypervolume_calculated=True,
            no_realized_pollinator_identity_or_function_inferred=True,
            no_private_taxon_island_or_source_identifiers_exported=True))

def reconcile_baseline(result,baseline):
    if result['n_original_islands']!=baseline['n_original_global_distinct_island_floras']:
        raise ValueError('Original full botanical flora sample changed')
    for x in result['endpoint_coverage']:
        old=[r for r in baseline['marginal_trait_and_pairwise_coverage']
             if r['source_plant_status_scope']==x['scope'] and
                r['trait_name']==x['endpoint']]
        if sum(r['n_islands_with_at_least_five_direct_species']
               for r in old)!=x['n_n5_islands']:
            raise ValueError('Original endpoint species mask changed')
    return dict(n_reconciled_endpoints=len(result['endpoint_coverage']),
                original_n5_island_counts_unchanged=True)

def main():
    ap=argparse.ArgumentParser()
    for key in ('original-root','geology','baseline','output'):
        ap.add_argument('--'+key,type=Path,required=True)
    ap.add_argument('--block-size',type=int,default=96)
    args=ap.parse_args()
    paths=dict(traits=args.original_root/
        'fixed/source/current/evidence/direct_only_trait_ledger.csv.gz',
        flora=args.original_root/
        'fixed/canonical/input/chapter1_status_flora.csv.gz')
    fingerprints={
        'traits':'9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745',
        'flora':'6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a'}
    sha={}
    for name,path in paths.items():
        h=hashlib.sha256()
        with path.open('rb') as inp:
            for chunk in iter(lambda:inp.read(1<<20),b''):h.update(chunk)
        sha[name]=h.hexdigest()
        if sha[name]!=fingerprints[name]:
            raise ValueError('Original private source checksum changed')
    traits=pd.read_csv(paths['traits'],dtype=str,usecols=[
        'accepted_species','trait_name','normalized_value',
        'evidence_scope','evidence_quality'])
    flora=pd.read_csv(paths['flora'],dtype=str,usecols=[
        'island_id','accepted_species','origin_status'])
    geography=pd.read_csv(args.geology,dtype=str)
    if flora.island_id.nunique()!=4505:raise ValueError('Original 4505 botanical floras changed')
    if geography.geological_origin.isin(('oceanic','continental')).sum()!=434:
        raise ValueError('Independently verified geology sample changed')
    result=analyse(strict_traits(traits),flora,geography,
                   block_size=args.block_size)
    baseline=json.loads(args.baseline.read_text())
    result['source_baseline_reconciliation']=reconcile_baseline(result,baseline)
    result['source_sha256']=sha
    if result['n_geologically_typed_original_flora']!=425:
        raise ValueError('Original geology×botanical flora subset changed')
    body=json.dumps(result,sort_keys=True,indent=2,allow_nan=False)
    for term in ('island_id','accepted_species','entity_ID','work_species',
                 'spatial_block','source_signature'):
        if '"'+term+'"' in body:
            raise ValueError('Private original taxon/island/source key in result')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+'\n',encoding='utf-8')
    # Only anonymous pre-fixed pooled/regional summaries reach logs; never island/species keys.
    endpoint='flower_primary_color__floral_form'
    selected={'all','no_shared_full_flora_species',
        'no_shared_full_flora_species_genus_shared',
        'no_shared_full_flora_genera_strict','no_shared_direct_genera',
        'genus_disjoint_different_blocks'}
    highlights=[{k:row[k] for k in ('scope','region','subset','n_pairs',
        'n_unique_islands','n_unique_geographic_blocks',
        'mean_functional_overlap','mean_squared_frequency_distance',
        'mean_genus_null_squared_distance',
        'mean_genus_null_excess_squared_distance')
        if k in row} for row in result['pair_summaries']
        if row['endpoint']==endpoint and row['region'] in
            ('GLOBAL_ALL','tropical','northern_midlatitude',
             'northern_highlatitude','southern_extratropical')
        and row['geology_pair']=='all' and row['subset'] in selected]
    print('ANONYMOUS_GENUS_DISJOINT_COLOUR_FORM='+json.dumps(
        highlights,sort_keys=True,allow_nan=False))
    print(json.dumps(dict(status=result['status'],
        n_original_islands=result['n_original_islands'],
        n_endpoints=len(result['endpoint_coverage']),
        n_summaries=len(result['pair_summaries']),
        pair_pseudoreplication_not_inferred_independent=True)))

if __name__=='__main__':
    main()
