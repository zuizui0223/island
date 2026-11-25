"""Island floral-composition beta diversity; NOT intraspecific adaptation.

Every original flora stays in coverage denominators. Each single/pair trait
uses its own observed-species mask. Distribution affinity is sum sqrt(p*q)
= 1 - normalized squared Hellinger distance. Profile-set beta partitions
are descriptive category replacement/nestedness, not historical change.
Only anonymous endpoint/region/geology summaries leave protected source CI.
"""
from __future__ import annotations
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree

PAIRS=(('flower_primary_color','floral_form'),('flower_primary_color','floral_symmetry'),
       ('floral_form','floral_symmetry'))
PRIMARY=('flower_primary_color','floral_form','flower_primary_color__floral_form')
MIN_SPECIES=5


def _seed(*parts):
    return int.from_bytes(sha256('|'.join(map(str,parts)).encode()).digest()[:8],'big')


def random_edges(names,k=16,seed=20261115):
    """Every eligible island is a focal; no category value is consulted."""
    n=len(names)
    if len(set(names))!=n or k<1:raise ValueError('Unique island keys and positive degree required')
    if n<2:return np.array([],dtype=int),np.array([],dtype=int)
    aa=[];bb=[]
    for i,name in enumerate(names):
        rng=np.random.default_rng(_seed(seed,name))
        j=rng.choice(n-1,size=min(k,n-1),replace=False)
        j=j+(j>=i)
        aa.extend([i]*len(j));bb.extend(j.tolist())
    return np.array(aa,dtype=int),np.array(bb,dtype=int)


def composition_metrics(p,a,b):
    p=np.asarray(p,dtype=float);a=np.asarray(a,dtype=int);b=np.asarray(b,dtype=int)
    if p.ndim!=2 or np.any(p<0) or not np.isfinite(p).all() or not np.allclose(p.sum(axis=1),1):
        raise ValueError('Each composition must be a finite normalized distribution')
    x,y=p[a],p[b]
    affinity=np.sqrt(x*y).sum(axis=1)
    xx,yy=x>0,y>0
    shared=(xx&yy).sum(axis=1)
    left=(xx&~yy).sum(axis=1);right=(yy&~xx).sum(axis=1)
    sor=(left+right)/(2*shared+left+right)
    small=np.minimum(left,right)
    turnover=np.divide(small,shared+small,out=np.zeros(len(a),dtype=float),where=(shared+small)>0)
    return dict(affinity=np.clip(affinity,0,1),distribution_overlap=np.minimum(x,y).sum(axis=1),
        profile_jaccard=shared/(shared+left+right),profile_sorensen=sor,
        profile_turnover=turnover,profile_nestedness=np.maximum(0,sor-turnover))


def exact_all_pair_affinity(p):
    n=len(p)
    if n<2:return None
    return float(np.clip((np.square(np.sqrt(p).sum(axis=0)).sum()-n)/(n*(n-1)),0,1))


def proportions(x,codes,n_categories):
    labels=csr_matrix((np.ones(len(codes)),(np.arange(len(codes)),codes)),
                      shape=(len(codes),n_categories))
    count=(x@labels).toarray()
    n=count.sum(axis=1)
    if np.any(n<MIN_SPECIES):raise ValueError('An unobserved/insufficient trait was treated as an outcome')
    return count/n[:,None]


def shuffle_codes(codes,groups,rng):
    if groups is None:return rng.permutation(codes)
    out=codes.copy()
    for ids in groups if isinstance(groups,list) else [np.flatnonzero(groups==g) for g in np.unique(groups)]:
        if len(ids)>1:out[ids]=rng.permutation(codes[ids])
    return out


def build_endpoint(states,flora,scope,geography):
    if scope not in ('all_recorded','original_native'):raise ValueError('Invalid original flora scope')
    names=sorted(states);pos={n:i for i,n in enumerate(names)}
    categories=sorted(set(states.values()));code={v:i for i,v in enumerate(categories)}
    codes=np.array([code[states[n]] for n in names],dtype=int)
    gen={};gids=[]
    for n in names:
        g=n.split()[0] if len(n.split())>=2 else '__fixed_'+n
        if g not in gen:gen[g]=len(gen)
        gids.append(gen[g])
    all_names=flora.groupby('island_id',sort=True).accepted_species.agg(set).to_dict()
    if scope=='original_native':
        obs=flora.loc[flora.origin_status.eq('native')].groupby('island_id').accepted_species.agg(set).to_dict()
    else:obs=all_names
    source_set=set(names);ids=[];rosters=[];trait_rosters=[];row=[];col=[]
    all_counts=[];coverage=[];qvals=[]
    for island in sorted(all_names):
        own=obs.get(island,set());known=own&source_set
        all_counts.append((len(own),len(known)))
        if len(known)<MIN_SPECIES:continue
        i=len(ids);ids.append(island);rosters.append(own);trait_rosters.append(known)
        selected=[pos[n] for n in sorted(known)]
        row.extend([i]*len(selected));col.extend(selected)
        count=np.bincount(codes[selected],minlength=len(categories));n=len(selected)
        qvals.append(float(.8*(1-np.sum(count*(count-1))/(n*(n-1)))))
        coverage.append(len(known)/len(own))
    x=csr_matrix((np.ones(len(row)),(row,col)),shape=(len(ids),len(names)))
    geo=geography.set_index('island_id').reindex(ids)
    for c in ('analysis_regime','geological_origin','spatial_block'):
        geo[c]=geo[c].fillna('unresolved')
    return dict(ids=ids,incidence=x,codes=codes,categories=categories,genera=np.array(gids),
        flora_sets=rosters,trait_sets=trait_rosters,geo=geo,counts=np.array(all_counts),
        coverage=np.array(coverage),raw_q=np.array(qvals))


def _taxonomic(sets,a,b):
    intersections=np.array([len(sets[i]&sets[j]) for i,j in zip(a,b)],dtype=int)
    sizes=np.array([len(s) for s in sets],dtype=int)
    return intersections/(sizes[a]+sizes[b]-intersections),intersections


def _focal_mean(v,a,n):
    v=np.asarray(v,dtype=float)
    sums=np.bincount(a,weights=v,minlength=n);count=np.bincount(a,minlength=n)
    return float(np.mean(sums[count>0]/count[count>0])) if count.any() else None


def summarize(p,a,b,full_sets,known_sets,mask=None):
    if mask is not None:a,b=a[mask],b[mask]
    if not len(a):return dict(n_directed_comparisons=0,n_unique_unordered_pairs=0,n_focal_islands=0)
    values=composition_metrics(p,a,b)
    jf,inter=_taxonomic(full_sets,a,b);jk,_=_taxonomic(known_sets,a,b)
    values['full_flora_taxonomic_jaccard']=jf;values['trait_mask_taxonomic_jaccard']=jk
    out=dict(n_directed_comparisons=len(a),n_unique_unordered_pairs=len(set(zip(np.minimum(a,b),np.maximum(a,b)))),
             n_focal_islands=len(set(a)),n_participating_islands=len(set(a)|set(b)),
             n_no_shared_full_flora_directed_comparisons=int((inter==0).sum()))
    out.update({'mean_'+key:_focal_mean(value,a,len(p)) for key,value in values.items()})
    return out


def nearest_edges(e,k=3):
    g=e['geo'];lat=pd.to_numeric(g.island_latitude,errors='coerce').to_numpy()
    lon=pd.to_numeric(g.island_longitude,errors='coerce').to_numpy()
    good=np.flatnonzero(np.isfinite(lat)&np.isfinite(lon)&(np.abs(lat)<=90)&(np.abs(lon)<=180))
    if len(good)<2:return np.array([],dtype=int),np.array([],dtype=int)
    la,lo=np.radians(lat[good]),np.radians(lon[good])
    xyz=np.column_stack((np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)))
    tree=cKDTree(xyz);d,idx=tree.query(xyz,k=min(k+1,len(good)))
    aa=[];bb=[]
    for local,i in enumerate(good):
        candidates=sorted([(float(dist),int(good[j])) for dist,j in zip(np.atleast_1d(d[local]),np.atleast_1d(idx[local])) if j!=local])[:k]
        for _,j in candidates:aa.append(int(i));bb.append(j)
    return np.array(aa,dtype=int),np.array(bb,dtype=int)


def _null(e,p,a,b,permutations,seed):
    _,inter=_taxonomic(e['flora_sets'],a,b);mask=inter==0
    aa,bb=a[mask],b[mask]
    if len(set(aa))<20:return [dict(reference=ref,testable=False,n_focal_islands=len(set(aa))) for ref in ('global','within_genus')]
    observed=_focal_mean(composition_metrics(p,aa,bb)['affinity'],aa,len(p))
    groups=[np.flatnonzero(e['genera']==g) for g in np.unique(e['genera'])]
    variable=sum(len(set(e['codes'][z]))>1 for z in groups)
    result=[]
    for ref in ('global','within_genus'):
        rng=np.random.default_rng(_seed(seed,ref));means=[]
        for it in range(permutations):
            shuffled=shuffle_codes(e['codes'],None if ref=='global' else groups,rng)
            pp=proportions(e['incidence'],shuffled,len(e['categories']))
            v=np.sqrt(pp[aa]*pp[bb]).sum(axis=1)
            means.append(_focal_mean(v,aa,len(pp)))
        v=np.array(means)
        result.append(dict(reference=ref,testable=True,n_permutations=permutations,
            n_focal_islands=len(set(aa)),n_directed_disjoint_flora_comparisons=len(aa),
            n_variable_global_genus_pools=variable,observed_affinity=observed,
            mean_null_affinity=float(v.mean()),observed_minus_null=float(observed-v.mean()),
            null95=[float(z) for z in np.quantile(v,[.025,.975])],
            p_excess_similarity=float((1+np.sum(v>=observed-1e-12))/(len(v)+1)),
            p_lower_similarity_diagnostic=float((1+np.sum(v<=observed+1e-12))/(len(v)+1)),
            degenerate_reference=bool(np.ptp(v)<1e-12),
            whole_species_profile_shared_across_islands=True))
    return result


def _bh(v):
    v=np.array(v,dtype=float);order=np.argsort(v);out=np.ones(len(v));last=1.
    for rank in range(len(v)-1,-1,-1):
        i=order[rank];last=min(last,v[i]*len(v)/(rank+1));out[i]=last
    return out


def analyse(states,flora,geography,*,pairs=PAIRS,permutations=199,seed=20261115,primary=PRIMARY,partners=16):
    required={'island_id','accepted_species','origin_status'}
    if not required<=set(flora) or flora.duplicated(['island_id','accepted_species']).any():
        raise ValueError('Unique original island-species membership required')
    if flora[list(required)].isna().any().any() or not flora.origin_status.isin(('native','introduced','unresolved')).all():
        raise ValueError('Invalid source native status/identity')
    needed={'island_id','analysis_regime','geological_origin','spatial_block','island_latitude','island_longitude'}
    if not needed<=set(geography) or geography.island_id.duplicated().any():
        raise ValueError('Unique original geographical metadata required')
    endpoints=list(states.items())
    for left,right in pairs:
        shared=set(states[left])&set(states[right])
        endpoints.append((left+'__'+right,{n:(states[left][n],states[right][n]) for n in shared}))
    result=[];tests=[];n_original=int(flora.island_id.nunique())
    for scope in ('all_recorded','original_native'):
        for tag,ss in endpoints:
            e=build_endpoint(ss,flora,scope,geography);n=len(e['ids'])
            out=dict(scope=scope,endpoint=tag,n_original_islands=n_original,n_eligible_islands=n,
                     n_with_any_trait=int((e['counts'][:,1]>0).sum()),
                     n_empty_original_scope=int((e['counts'][:,0]==0).sum()),
                     n_all_unordered_pairs=n*(n-1)//2,
                     mean_raw_finite5_rao=float(e['raw_q'].mean()) if n else None,
                     mean_recorded_fraction_among_eligible=float(e['coverage'].mean()) if n else None)
            if n<2:result.append(out);continue
            p=proportions(e['incidence'],e['codes'],len(e['categories']))
            a,b=random_edges(e['ids'],partners,_seed(seed,scope,tag,'edges'))
            out['n_focal_islands_in_uniform_graph']=len(set(a))
            out['exact_all_unordered_pair_mean_affinity']=exact_all_pair_affinity(p)
            out['uniform_pairs']=summarize(p,a,b,e['flora_sets'],e['trait_sets'])
            jf,inter=_taxonomic(e['flora_sets'],a,b)
            out['uniform_no_shared_flora']=summarize(p,a,b,e['flora_sets'],e['trait_sets'],inter==0)
            out['uniform_taxonomic_jaccard_at_most_0_1']=summarize(p,a,b,e['flora_sets'],e['trait_sets'],jf<=.1)
            near_a,near_b=nearest_edges(e)
            out['geographical_nearest_three']=summarize(p,near_a,near_b,e['flora_sets'],e['trait_sets'])
            # Descriptive geography: all pairs within each broad-region typed subset.
            # Symmetrize then focal-average; no falsely independent pair p-values.
            g=e['geo'];origin=g.geological_origin.to_numpy();region=g.analysis_regime.to_numpy()
            panels=[]
            for reg in sorted(set(region)-{'unresolved','geography_unresolved'}):
                loc=np.flatnonzero((region==reg)&np.isin(origin,['oceanic','continental']))
                if len(loc)<2:continue
                ia,ib=np.triu_indices(len(loc),1);left,right=loc[ia],loc[ib]
                for kind in ('oceanic_oceanic','continental_continental','oceanic_continental'):
                    if kind=='oceanic_continental':mask=origin[left]!=origin[right]
                    else:mask=(origin[left]==kind.split('_')[0])&(origin[right]==kind.split('_')[0])
                    x,y=left[mask],right[mask]
                    if not len(x):continue
                    aa,bb=np.r_[x,y],np.r_[y,x]
                    panels.append(dict(region=reg,geology_pair=kind,
                        **summarize(p,aa,bb,e['flora_sets'],e['trait_sets'])))
            out['within_region_geological_pair_descriptions']=panels
            if scope=='all_recorded' and tag in primary and permutations:
                null=_null(e,p,a,b,permutations,_seed(seed,tag,'null'))
                for z in null:tests.append(dict(endpoint=tag,scope=scope,**z))
            assert out['n_focal_islands_in_uniform_graph']==n
            result.append(out)
    if tests:
        q=_bh([t.get('p_excess_similarity',1) if not t.get('degenerate_reference',False) else 1 for t in tests])
        for t,qq in zip(tests,q):t['BH_excess_similarity_fixed_family']=float(qq) if t['testable'] else None
    return dict(status='GLOBAL_ISLAND_FLORAL_COMPOSITION_BETA_NOT_ADAPTATION',
        n_original_flora_islands=n_original,n_endpoints_per_scope=len(endpoints),endpoints=result,
        disjoint_flora_profile_permutation_tests=tests,
        claims=dict(no_intraspecific_change_estimated=True,no_effective_pollinator_assigned=True,
            all_original_flora_kept_in_denominators=True,only_endpoint_specific_two_trait_completion=True,
            equal_species_weights_not_plant_abundance=True,no_complete_three_trait_or_same_genus_filter=True,
            random_comparators_selected_without_flower_state_values=True,
            all_eligible_islands_are_uniform_graph_focals=True,
            exact_global_pair_affinity_but_other_full_population_distributions_sampled=True,
            no_pair_independence_pvalues_or_geology_causal_test=True,
            profile_permutations_preserve_species_incidence_and_global_joint_category_counts=True,
            reference_pool_observed_global_species_not_ancestral_source=True,
            genus_reference_not_phylogeny_or_missingness_correction=True,
            categorical_nestedness_not_historical_loss=True,not_KDE_hypervolume_overlap=True,
            only_anonymous_aggregate_no_source_taxon_or_island_identifiers=True))


def reconcile(result,baseline):
    if result['n_original_flora_islands']!=baseline['n_original_global_distinct_island_floras']:
        raise ValueError('Original 4505 flora denominator changed')
    errors=[]
    for e in result['endpoints']:
        old=[r for r in baseline['marginal_trait_and_pairwise_coverage']
             if r['source_plant_status_scope']==e['scope'] and r['trait_name']==e['endpoint']]
        n=sum(r['n_islands_with_at_least_five_direct_species'] for r in old)
        if n!=e['n_eligible_islands'] or sum(r['n_all_original_flora_islands'] for r in old)!=e['n_original_islands']:
            raise ValueError('Original endpoint mask changed')
        if n:
            mean=sum(r['n_islands_with_at_least_five_direct_species']*r.get('raw_observed_island_equal_finite5_mean',0) for r in old)/n
            error=abs(mean-e['mean_raw_finite5_rao']);errors.append(error)
            if error>1e-10:raise ValueError('Original Rao baseline changed')
    return dict(n_reconciled_scope_endpoints=len(result['endpoints']),max_raw_rao_difference=max(errors,default=0.))


def main():
    import argparse,hashlib
    from .global_partial_trait_flora_fd import strict_single_trait_original
    ap=argparse.ArgumentParser()
    for name in ('original-root','geography','baseline','output'):ap.add_argument('--'+name,type=Path,required=True)
    ap.add_argument('--permutations',type=int,default=199)
    args=ap.parse_args()
    paths={'traits':args.original_root/'fixed/source/current/evidence/direct_only_trait_ledger.csv.gz',
           'flora':args.original_root/'fixed/canonical/input/chapter1_status_flora.csv.gz'}
    sha={'traits':'9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745',
         'flora':'6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a'}
    for key,path in paths.items():
        h=hashlib.sha256()
        with path.open('rb') as f:
            for part in iter(lambda:f.read(1<<20),b''):h.update(part)
        if h.hexdigest()!=sha[key]:raise ValueError('Original immutable source changed: '+key)
    ledger=pd.read_csv(paths['traits'],dtype=str,usecols=['accepted_species','trait_name','normalized_value','evidence_scope','evidence_quality'])
    flora=pd.read_csv(paths['flora'],dtype=str,usecols=['island_id','accepted_species','origin_status'])
    geo=pd.read_csv(args.geography,dtype=str,usecols=['island_id','analysis_regime','geological_origin','spatial_block','island_latitude','island_longitude'])
    if flora.island_id.nunique()!=4505:raise ValueError('Original flora not 4505')
    states=strict_single_trait_original(ledger)
    out=analyse(states,flora,geo,permutations=args.permutations)
    out['baseline_reconciliation']=reconcile(out,json.loads(args.baseline.read_text()))
    out['source_sha256']=sha;out['seed']=20261115
    assert len(out['endpoints'])==26 and len(out['disjoint_flora_profile_permutation_tests'])==6
    text=json.dumps(out,sort_keys=True,indent=2,allow_nan=False)
    forbidden=('island_id','accepted_species','entity_ID','work_species','spatial_block','source_signature','latitude','longitude','species_name','genus_name')
    if any('"'+x+'"' in text for x in forbidden):raise ValueError('Protected identities in aggregate')
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(text+'\n')
    print(json.dumps(dict(status=out['status'],n_flora=4505,n_endpoints=26,baseline=out['baseline_reconciliation'])))

if __name__=='__main__':main()
