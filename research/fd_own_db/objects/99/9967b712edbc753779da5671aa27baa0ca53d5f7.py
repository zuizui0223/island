"""Same-mask marginal vs joint floral similarity and anonymous profile contributions.

The product-of-marginals is an analytical comparison, NEVER imputed observed
flowers. Its difference from joint affinity is signed, not a mutual information,
selection coefficient or proof of a particular pollinator syndrome.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .island_floral_beta import (PAIRS,build_endpoint,proportions,random_edges,
    _seed,_taxonomic,_focal_mean,composition_metrics)


def summarize_structure(p,categories,a,b,minimum_islands=10):
    if not len(a):return dict(n_directed_comparisons=0,n_focal_islands=0)
    if any(len(c)!=2 for c in categories):raise ValueError('Actual two-trait profiles required')
    first=sorted({c[0] for c in categories});second=sorted({c[1] for c in categories})
    map_a=np.array([[float(c[0]==x) for x in first] for c in categories])
    map_b=np.array([[float(c[1]==x) for x in second] for c in categories])
    pa,pb=p@map_a,p@map_b
    ba=composition_metrics(pa,a,b)['affinity'];bb=composition_metrics(pb,a,b)['affinity']
    contrib=np.sqrt(p[a]*p[b]);bj=contrib.sum(axis=1)
    if np.any(bj>np.minimum(ba,bb)+1e-10):raise ArithmeticError('Coarsening inequality violated')
    n=len(p);counts=np.bincount(a,minlength=n)
    weights=1/counts[a]/np.count_nonzero(counts)
    means=weights@contrib
    support=np.count_nonzero(p[np.unique(np.r_[a,b])]>0,axis=0)
    order=[int(i) for i in np.argsort(-means) if support[i]>=minimum_islands][:8]
    listed=[dict(first_state=categories[i][0],second_state=categories[i][1],
                 affinity_contribution=float(means[i]),n_observed_profile_islands=int(support[i])) for i in order]
    joint=_focal_mean(bj,a,n)
    return dict(n_directed_comparisons=len(a),n_focal_islands=len(set(a)),
        observed_joint_affinity=joint,same_mask_first_trait_affinity=_focal_mean(ba,a,n),
        same_mask_second_trait_affinity=_focal_mean(bb,a,n),
        marginal_product_reference_affinity=_focal_mean(ba*bb,a,n),
        joint_minus_product_reference=_focal_mean(bj-ba*bb,a,n),
        leading_observed_profile_contributions=listed,
        unlisted_contribution=float(joint-sum(x['affinity_contribution'] for x in listed)),
        no_marginal_cartesian_profile_promoted_to_observation=True)


def analyse_structure(states,flora,geo,baseline):
    out=[]
    for scope in ('all_recorded','original_native'):
        for first,second in PAIRS:
            tag=first+'__'+second
            names=set(states[first])&set(states[second])
            joint={n:(states[first][n],states[second][n]) for n in names}
            e=build_endpoint(joint,flora,scope,geo)
            old=next(x for x in baseline['endpoints'] if x['scope']==scope and x['endpoint']==tag)
            if len(e['ids'])!=old['n_eligible_islands']:raise ValueError('Original pair-specific island mask changed')
            p=proportions(e['incidence'],e['codes'],len(e['categories']))
            a,b=random_edges(e['ids'],16,_seed(20261115,scope,tag,'edges'))
            _,inter=_taxonomic(e['flora_sets'],a,b)
            allp=summarize_structure(p,e['categories'],a,b)
            dis=summarize_structure(p,e['categories'],a[inter==0],b[inter==0])
            if abs(allp['observed_joint_affinity']-old['uniform_pairs']['mean_affinity'])>1e-10:
                raise ValueError('Original comparator graph/raw joint affinity changed')
            if 'observed_joint_affinity' in dis:
                if abs(dis['observed_joint_affinity']-old['uniform_no_shared_flora']['mean_affinity'])>1e-10:
                    raise ValueError('Original disjoint flora mask changed')
            out.append(dict(scope=scope,endpoint=tag,n_original_flora_islands=int(flora.island_id.nunique()),
                n_eligible_islands=len(e['ids']),uniform_comparison=allp,disjoint_flora_comparison=dis))
    return dict(status='SAME_MASK_FLORAL_MARGINAL_JOINT_STRUCTURE_DIAGNOSTIC',
        n_original_flora_islands=int(flora.island_id.nunique()),same_mask_pair_endpoints=out,
        claims=dict(marginals_recomputed_on_identical_two_trait_species_and_pairs=True,
            structural_dimensionality_not_evidence_of_adaptation=True,
            product_reference_is_not_imputation_or_a_sampled_null_distribution=True,
            category_contributions_are_observed_traits_not_assigned_pollinator_identity=True,
            no_new_source_taxon_or_island_identifiers_exported=True))


def main():
    import argparse,hashlib
    from .global_partial_trait_flora_fd import strict_single_trait_original
    ap=argparse.ArgumentParser()
    for name in ('original-root','geography','beta-baseline','output'):ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args()
    tp=args.original_root/'fixed/source/current/evidence/direct_only_trait_ledger.csv.gz'
    fp=args.original_root/'fixed/canonical/input/chapter1_status_flora.csv.gz'
    for path,expected in [(tp,'9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745'),
                          (fp,'6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a')]:
        h=hashlib.sha256()
        with path.open('rb') as f:
            for part in iter(lambda:f.read(1<<20),b''):h.update(part)
        if h.hexdigest()!=expected:raise ValueError('Immutable original source hash changed')
    traits=pd.read_csv(tp,dtype=str,usecols=['accepted_species','trait_name','normalized_value','evidence_scope','evidence_quality'])
    flora=pd.read_csv(fp,dtype=str,usecols=['island_id','accepted_species','origin_status'])
    if flora.island_id.nunique()!=4505:raise ValueError('Original 4505 flora frame changed')
    geo=pd.read_csv(args.geography,dtype=str,usecols=['island_id','analysis_regime','geological_origin','spatial_block','island_latitude','island_longitude'])
    baseline=json.loads(args.beta_baseline.read_text())
    out=analyse_structure(strict_single_trait_original(traits),flora,geo,baseline)
    text=json.dumps(out,indent=2,sort_keys=True,allow_nan=False)
    for key in ('island_id','accepted_species','entity_ID','work_species','latitude','longitude','genus_name'):
        if '"'+key+'"' in text:raise ValueError('Restricted source identifier exported')
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(text+'\n')
    print('SAME original island pairs and observed joint affinity reconciled for all six endpoint/scope panels')

if __name__=='__main__':main()
