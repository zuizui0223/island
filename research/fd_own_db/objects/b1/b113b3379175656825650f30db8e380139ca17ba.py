"""Global per-trait FD with the ORIGINAL recording mask and taxonomic counts.

Exact expectation under globally permuting the FIRST trait among recorded
species within named genus/family strata; no trait imputation and no matching-
induced removal of islands. Family-unmapped species are fixed singleton
strata. A categorical taxon-conditioned residual is not adaptation, and these
ranks are NOT a dated phylogeny. Output contains anonymous aggregates only.
"""
from __future__ import annotations
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
import pandas as pd

N_MIN = 5
FACTOR = 0.8
PAIRS = (("flower_primary_color", "floral_form"),
         ("flower_primary_color", "floral_symmetry"),
         ("floral_form", "floral_symmetry"))


def finite_q(labels: list) -> float:
    """Expected empirical Rao of five DISTINCT species, categorical distance."""
    n = len(labels)
    if n < N_MIN:
        return float('nan')
    c = Counter(labels)
    return FACTOR * (1 - sum(x*(x-1) for x in c.values())/(n*(n-1)))


def missing_bounds(counts, total: int) -> tuple[float | None, float | None]:
    """Sharp PER-ISLAND bounds over a declared finite category codebook.

    All missing states in a modal category minimize mismatch. Integer water
    filling maximizes mismatch. These are NOT filled observations, not joint
    cross-island sharp bounds, and allow no unrecorded/new codebook category.
    """
    x = np.asarray(counts, dtype=float)
    if x.ndim != 1 or not len(x) or not np.isfinite(x).all():
        raise ValueError('Finite nonempty category counts required')
    if np.any(x < 0) or np.any(x != np.floor(x)) or int(total) != total or total < x.sum():
        raise ValueError('Invalid counts or full-flora denominator')
    x = x.astype(np.int64)
    total = int(total)
    if total < N_MIN:
        return None, None
    miss = total-int(x.sum())
    low = x.copy(); low[np.argmax(low)] += miss
    lo, hi = 0, total+1
    while hi-lo > 1:
        level = (lo+hi)//2
        if np.maximum(x, level).sum() <= total:
            lo = level
        else:
            hi = level
    high = np.maximum(x, lo)
    left = total-int(high.sum())
    indices = np.flatnonzero(high == lo)
    high[indices[:left]] += 1
    assert high.sum() == total and np.all(high >= x)
    def q(v):
        return float(FACTOR*(1-np.sum(v*(v-1))/(total*(total-1))))
    return q(low), q(high)


class TaxonNull:
    """First-trait global pool, with exact finite-population group moments."""
    def __init__(self, states: dict[str,str], groups: dict[str,str]):
        if not states or any(v is None for v in states.values()):
            raise ValueError('Nonempty actually observed first-trait pool required')
        self.states = states
        self.names = sorted(states)
        self.pos = {n:i for i,n in enumerate(self.names)}
        self.universe = set(self.names)
        categories = {c:i for i,c in enumerate(sorted(set(states.values())))}
        group_ids = {}
        memberships = []
        mapped = []
        for name in self.names:
            value = groups.get(name)
            valid = value is not None and str(value).strip() != ''
            key = ('mapped', str(value)) if valid else ('fixed_unmapped_species', name)
            if key not in group_ids:
                group_ids[key] = len(group_ids)
            memberships.append(group_ids[key]); mapped.append(valid)
        self.gid = np.array(memberships, dtype=int)
        self.mapped = np.array(mapped, dtype=bool)
        self.code = np.array([categories[states[n]] for n in self.names], dtype=int)
        table = np.zeros((len(group_ids), len(categories)), dtype=np.int64)
        np.add.at(table, (self.gid, self.code), 1)
        self.N = table.sum(axis=1)
        self.P = table/self.N[:,None]
        self.H = 1-np.sum(self.P*self.P, axis=1)
        self.variable = np.count_nonzero(table, axis=1) > 1

    def evaluate(self, selected: list[str], second: dict[str,str] | None = None) -> dict:
        if len(selected) != len(set(selected)) or not set(selected) <= self.universe:
            raise ValueError('Distinct selected species must belong to the observed pool')
        idx = np.array([self.pos[n] for n in selected], dtype=int)
        n = len(idx)
        if n < N_MIN:
            raise ValueError('Five observed species required for a finite5 score')
        raw_gid = self.gid[idx]
        groups, inv, m = np.unique(raw_gid, return_inverse=True, return_counts=True)
        N, P, H = self.N[groups], self.P[groups], self.H[groups]
        den = n*(n-1)
        def expected_same(hist):
            mu = hist.T @ P
            var_terms = np.zeros_like(hist, dtype=float)
            good = N > 1
            var_terms[good] = (hist[good]*(N[good,None]-hist[good]) /
                               (N[good,None]-1)) * H[good,None]
            return float(np.square(mu).sum()+var_terms.sum()-hist.sum())
        expected_first = expected_same(m[:,None]) / den
        first = [self.states[k] for k in selected]
        if second is None:
            labels = first
            expected_joint = expected_first
            movable = self.variable[groups] & (m < N)
            expected_pairing = None
        else:
            if not set(selected) <= set(second):
                raise ValueError('Paired endpoints require both observed traits')
            b = [second[k] for k in selected]
            blevels = {v:i for i,v in enumerate(sorted(set(b)))}
            mat = np.zeros((len(groups),len(blevels)),dtype=np.int64)
            np.add.at(mat,(inv,np.array([blevels[v] for v in b])),1)
            expected_joint = expected_same(mat) / den
            same_b = sum(c*(c-1) for c in Counter(b).values())/den
            expected_pairing = FACTOR*(expected_first*same_b-expected_joint)
            labels = list(zip(first,b))
            movable = self.variable[groups] & ((m < N) | (np.count_nonzero(mat,axis=1)>1))
        obs_q = finite_q(labels)
        equal_within = sum(c*(c-1) for c in Counter(zip(raw_gid,labels)).values())
        within = FACTOR*(np.sum(m*(m-1))-equal_within)/den
        out = dict(expected_q=float(FACTOR*(1-expected_joint)),
                   q_within=float(within),q_between=float(obs_q-within),
                   n_movable_groups=int(movable.sum()),
                   n_mapped_species=int(self.mapped[idx].sum()))
        if expected_pairing is not None:
            out['expected_pairing'] = float(expected_pairing)
        if not -1e-10 <= out['expected_q'] <= FACTOR+1e-10:
            raise ArithmeticError('Expected finite5 Q outside its mathematical range')
        return out


def _seed(*x):
    return int.from_bytes(sha256('|'.join(map(str,x)).encode()).digest()[:8],'big')


def endpoint_summary(q:pd.DataFrame) -> dict:
    known = q.loc[q.n_obs >= N_MIN]
    valid_bounds = q.loc[q.bound_lo.notna()]
    out = dict(n_original_islands=len(q),n_zero_native_or_empty_scope=int((q.n_total==0).sum()),
               n_with_any_observed_trait=int((q.n_obs>0).sum()),n_n5_observed=len(known),
               n_n5_genus_null=int(q.genus_expected_q.notna().sum()),
               n_n5_family_null=int(q.family_expected_q.notna().sum()),
               n_full_flora_N5_for_missingness_bounds=len(valid_bounds),
               mean_recording_fraction=(float(q.coverage.mean()) if q.coverage.notna().any() else None))
    if len(known):
        for c in ('raw_q','genus_expected_q','family_expected_q','genus_residual_q',
                  'family_residual_q','genus_q_within','genus_q_between',
                  'family_q_within','family_q_between','raw_pairing',
                  'genus_residual_pairing','family_residual_pairing','family_mapped_fraction'):
            if known[c].notna().any():out['mean_'+c]=float(known[c].mean())
        for rank in ('genus','family'):
            out['n_'+rank+'_zero_shufflable_groups'] = int((known[rank+'_movable_groups']==0).sum())
            denom=float(known.raw_q.sum())
            out[rank+'_between_pair_contribution_fraction']=(
                float(known[rank+'_q_between'].sum()/denom) if denom>1e-12 else None)
        assert out['n_n5_genus_null']==out['n_n5_family_null']==out['n_n5_observed']
    if len(valid_bounds):
        out['mean_missing_bound_lo']=float(valid_bounds.bound_lo.mean())
        out['mean_missing_bound_hi']=float(valid_bounds.bound_hi.mean())
        out['median_missing_bound_width']=float((valid_bounds.bound_hi-valid_bounds.bound_lo).median())
    return out


def _origin_contrast(q, metric, draws, seed):
    z=q.loc[q.origin.isin(('oceanic','continental')) & q[metric].notna() & q.block.notna()]
    counts=z.groupby('origin').size().to_dict()
    tab=z.groupby(['block','origin'])[metric].mean().unstack().reindex(columns=['oceanic','continental'])
    k=tab.notna().sum()
    out={'metric':metric,'n_oceanic':counts.get('oceanic',0),
         'n_continental':counts.get('continental',0),
         'n_oceanic_blocks':int(k.oceanic),'n_continental_blocks':int(k.continental),
         'descriptive_interval_available':False}
    if min(counts.get('oceanic',0),counts.get('continental',0))<10 or min(k)<5:
        return out
    v=tab.to_numpy(float); mask=np.isfinite(v)
    point=np.nanmean(v,axis=0)
    rng=np.random.default_rng(seed)
    w=rng.multinomial(len(v),np.full(len(v),1/len(v)),size=draws)
    denominator=w@mask.astype(float)
    ok=np.all(denominator>0,axis=1)
    bs=(w[ok]@np.nan_to_num(v))/denominator[ok]
    delta=bs[:,0]-bs[:,1]
    out.update(descriptive_interval_available=True,
               block_equal_oceanic_minus_continental=float(point[0]-point[1]),
               conditional_block_bootstrap95=[float(a) for a in np.quantile(delta,[.025,.975])],
               n_valid_bootstrap_draws=int(ok.sum()),
               shared_block_multiplicity_preserved=True)
    return out


def analyse(states:dict, flora:pd.DataFrame, geography:pd.DataFrame,
            family_map:dict, *,pairs=PAIRS,draws=999,seed=20261114)->dict:
    need={'island_id','accepted_species','origin_status'}
    if not need<=set(flora) or flora.duplicated(['island_id','accepted_species']).any():
        raise ValueError('Unique original island-species membership required')
    if flora[list(need)].isna().any().any() or not flora.origin_status.isin(('native','introduced','unresolved')).all():
        raise ValueError('Missing identity or original status')
    geo_cols={'island_id','analysis_regime','geological_origin','spatial_block'}
    if not geo_cols<=set(geography) or geography.island_id.duplicated().any():
        raise ValueError('Original unique geographic crosswalk required')
    endpoints=[]
    for trait,s in states.items():
        if not s:raise ValueError('Empty global trait codebook')
        endpoints.append((trait,s,None))
    for a,b in pairs:
        shared=set(states[a])&set(states[b])
        if not shared:raise ValueError('No original paired trait evidence')
        endpoints.append((a+'__'+b,{n:states[a][n] for n in shared},{n:states[b][n] for n in shared}))
    pools={}
    for tag,first,second in endpoints:
        genus={n:n.split()[0] for n in first if len(n.split())>=2}
        codebook=(sorted(set(first.values())) if second is None else
                  [(a,b) for a in sorted(set(first.values())) for b in sorted(set(second.values()))])
        pools[tag]=(TaxonNull(first,genus),TaxonNull(first,family_map),set(first),codebook)
    geo=geography.set_index('island_id')
    missing_geo=set(flora.island_id)-set(geo.index)
    private=[]
    for island,roster in flora.groupby('island_id',sort=True):
        if island in geo.index:
            row=geo.loc[island]
            region=row.analysis_regime if pd.notna(row.analysis_regime) else 'geography_unresolved'
            origin=row.geological_origin if pd.notna(row.geological_origin) else 'unresolved'
            block=row.spatial_block if pd.notna(row.spatial_block) else None
        else:
            region,origin,block='geography_unresolved','unresolved',None
        for scope,own in (('all_recorded',set(roster.accepted_species)),
                          ('original_native',set(roster.loc[roster.origin_status.eq('native'),'accepted_species']))):
            for tag,first,second in endpoints:
                pg,pf,known_names,codebook=pools[tag]
                selected=sorted(own&known_names); n=len(selected); N=len(own)
                labels=[first[x] for x in selected] if second is None else [(first[x],second[x]) for x in selected]
                counts=Counter(labels)
                lower,upper=missing_bounds([counts[k] for k in codebook],N)
                qrow=dict(island_id=island,region=region,origin=origin,block=block,
                          scope=scope,endpoint=tag,n_total=N,n_obs=n,
                          coverage=n/N if N else np.nan,bound_lo=lower,bound_hi=upper,
                          raw_q=np.nan,raw_pairing=np.nan,family_mapped_fraction=np.nan)
                for rank in ('genus','family'):
                    for key in ('expected_q','residual_q','q_within','q_between','movable_groups','residual_pairing'):
                        qrow[rank+'_'+key]=np.nan
                if n>=N_MIN:
                    qrow['raw_q']=finite_q(labels)
                    if second is not None:
                        pa=1-finite_q([first[x] for x in selected])/FACTOR
                        pb=1-finite_q([second[x] for x in selected])/FACTOR
                        pab=1-qrow['raw_q']/FACTOR
                        qrow['raw_pairing']=FACTOR*(pa*pb-pab)
                    for rank,pool in (('genus',pg),('family',pf)):
                        a=pool.evaluate(selected,second=second)
                        for key in ('expected_q','q_within','q_between'):
                            qrow[rank+'_'+key]=a[key]
                        qrow[rank+'_residual_q']=qrow['raw_q']-a['expected_q']
                        qrow[rank+'_movable_groups']=a['n_movable_groups']
                        if second is not None:
                            qrow[rank+'_residual_pairing']=qrow['raw_pairing']-a['expected_pairing']
                        if rank=='family':qrow['family_mapped_fraction']=a['n_mapped_species']/n
                private.append(qrow)
    data=pd.DataFrame(private)
    n_islands=flora.island_id.nunique()
    assert len(data)==n_islands*2*len(endpoints)
    global_support=[]
    for (scope,tag),q in data.groupby(['scope','endpoint'],sort=True):
        global_support.append(dict(scope=scope,endpoint=tag,**endpoint_summary(q)))
    regional=[];contrasts=[]
    metrics=('raw_q','genus_residual_q','family_residual_q','raw_pairing',
             'genus_residual_pairing','family_residual_pairing','coverage')
    for (region,scope,tag),q in data.groupby(['region','scope','endpoint'],sort=True):
        for origin,grp in [('all',q),*[(o,q.loc[q.origin.eq(o)]) for o in ('oceanic','continental')]]:
            if len(grp):regional.append(dict(region=region,origin=origin,scope=scope,endpoint=tag,**endpoint_summary(grp)))
        for metric in metrics:
            if q[metric].notna().any():
                contrasts.append(dict(region=region,scope=scope,endpoint=tag,
                    **_origin_contrast(q,metric,draws,_seed(seed,region,scope,tag,metric))))
    return dict(status='GLOBAL_TAXONOMIC_MASK_FD_EXPECTATION_AND_MISSING_BOUNDS',
        n_original_flora_islands=int(n_islands),n_private_island_scope_endpoint_rows=len(data),
        n_flora_islands_missing_geography=len(missing_geo),
        n_geologically_typed_actual_flora_islands=int(data[['island_id','origin']].drop_duplicates().origin.isin(('oceanic','continental')).sum()),
        original_direct_species_by_trait={k:len(v) for k,v in states.items()},
        family_candidate_map_global_species=len(family_map),
        global_endpoint_support=global_support,regional_taxonomic_decomposition=regional,
        geological_origin_same_mask_contrasts=contrasts,
        claims=dict(all_original_flora_islands_retained=True,
          no_three_trait_completion_required=True,no_missing_trait_filled_as_observed=True,
          raw_and_genus_and_family_use_identical_trait_species_mask=True,
          unknown_family_fixed_as_singleton_not_inferred_or_dropped=True,
          family_map_source_candidates_not_reviewed_taxonomy=True,
          null_reference_global_recorded_species_not_ancestral_mainland=True,
          paired_null_permutates_first_trait_and_holds_second_trait=True,
          bounds_assume_observed_global_category_codebook_closed=True,
          pair_bounds_relax_partially_observed_single_trait_information=True,
          per_island_bounds_not_joint_cross_island_sharp_bounds=True,
          zero_residual_in_invariant_groups_not_adaptive_stasis=True,
          rank_conditioning_not_full_phylogenetic_correction=True,
          block_intervals_descriptive_not_multiple_test_or_species_dependence_corrected=True,
          no_pollinator_shift_or_causal_adaptation_measured=True,
          source_taxon_island_or_family_labels_not_exported=True))


def reconcile_baseline(result:dict, baseline:dict) -> dict:
    """Independently saved original aggregate: SAME observed Q and sample."""
    if baseline['n_original_global_distinct_island_floras'] != result['n_original_flora_islands']:
        raise ValueError('Original flora denominator changed')
    errors=[]
    for x in result['global_endpoint_support']:
        old=[r for r in baseline['marginal_trait_and_pairwise_coverage']
             if r['source_plant_status_scope']==x['scope'] and r['trait_name']==x['endpoint']]
        n=sum(r['n_islands_with_at_least_five_direct_species'] for r in old)
        if n != x['n_n5_observed'] or sum(r['n_all_original_flora_islands'] for r in old)!=x['n_original_islands']:
            raise ValueError('Per-trait observed original island mask changed')
        if n:
            mean=sum(r['n_islands_with_at_least_five_direct_species']*
                r.get('raw_observed_island_equal_finite5_mean',0.) for r in old)/n
            error=abs(mean-x['mean_raw_q']);errors.append(error)
            if error>1e-10:raise ValueError('Raw observed Rao changed during taxonomic conditioning')
    return dict(n_reconciled_scope_endpoints=len(result['global_endpoint_support']),
                max_absolute_raw_mean_difference=max(errors,default=0.),
                historical_baseline_is_independently_saved_aggregate=True)


def main():
    import argparse,hashlib
    from .global_partial_trait_flora_fd import strict_single_trait_original
    from .colxr_family_opportunity import strict_candidate_family_map
    ap=argparse.ArgumentParser()
    for flag in ('original-root','geology','crosswalk','baseline','output'):
        ap.add_argument('--'+flag,type=Path,required=True)
    ap.add_argument('--seed',type=int,default=20261114)
    args=ap.parse_args()
    paths={
      'traits':args.original_root/'fixed/source/current/evidence/direct_only_trait_ledger.csv.gz',
      'flora':args.original_root/'fixed/canonical/input/chapter1_status_flora.csv.gz',
      'taxonomy':args.crosswalk}
    expected={
      'traits':'9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745',
      'flora':'6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a',
      'taxonomy':'079d9c27b0ca3aa51f0c7310f43ad4901b07a319d6063eab5a7d62b291bdd9e3'}
    observed={}
    for key,path in paths.items():
        h=hashlib.sha256()
        with path.open('rb') as f:
            for part in iter(lambda:f.read(1<<20),b''):h.update(part)
        observed[key]=h.hexdigest()
        if observed[key]!=expected[key]:raise ValueError('Frozen source hash mismatch: '+key)
    ledger=pd.read_csv(paths['traits'],dtype=str,usecols=['accepted_species','trait_name',
        'normalized_value','evidence_scope','evidence_quality'])
    flora=pd.read_csv(paths['flora'],dtype=str,usecols=['island_id','accepted_species','origin_status'])
    geography=pd.read_csv(args.geology,dtype=str,usecols=['island_id','analysis_regime',
        'geological_origin','spatial_block'])
    states=strict_single_trait_original(ledger)
    if flora.island_id.nunique()!=4505 or len(states)!=10:
        raise ValueError('Full 4505 flora / ten-trait source frame changed')
    if geography.geological_origin.isin(('oceanic','continental')).sum()!=434:
        raise ValueError('Independent geological source population changed')
    cross=pd.read_csv(paths['taxonomy'],dtype=str)
    union=set().union(*(set(s) for s in states.values()))
    families,audit=strict_candidate_family_map(cross,union)
    out=analyse(states,flora,geography,families,seed=args.seed,draws=999)
    baseline=json.loads(args.baseline.read_text())
    out['original_baseline_reconciliation']=reconcile_baseline(out,baseline)
    out['original_family_source_audit']=audit
    out['original_frozen_sha256']=observed
    out['seed']=args.seed
    assert out['n_private_island_scope_endpoint_rows']==117130
    assert out['n_geologically_typed_actual_flora_islands']==425
    text=json.dumps(out,sort_keys=True,indent=2,allow_nan=False)
    private=('accepted_species','island_id','entity_ID','work_species','source_signature',
             'spatial_block','candidate_family','source_family','genus','latitude','longitude')
    if any('"'+k+'"' in text for k in private):
        raise ValueError('Protected source identifiers in aggregate output')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(text+'\n')
    print(json.dumps({'status':out['status'],'n_original_flora':4505,
        'n_endpoints':len(out['global_endpoint_support']),
        'baseline_check':out['original_baseline_reconciliation'],
        'n_family_candidates':len(families)},sort_keys=True))


if __name__=='__main__':
    main()
