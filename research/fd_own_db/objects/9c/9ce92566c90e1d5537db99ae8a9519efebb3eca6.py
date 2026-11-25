import itertools
import numpy as np
import pandas as pd
import pytest
from fd.global_taxonomic_mask_fd import TaxonNull, finite_q, missing_bounds, analyse, reconcile_baseline


def compositions(n,k):
    if k==1:
        yield (n,)
        return
    for a in range(n+1):
        for b in compositions(n-a,k-1):yield (a,)+b


def test_exact_single_null_matches_every_label_permutation():
    states=dict(zip('abcdefg',['r','r','b','r','b','g','g']))
    groups=dict(zip('abcdefg',['A','A','A','B','B','B','C']))
    selected=list('abdeg')
    got=TaxonNull(states,groups).evaluate(selected)
    alternatives=[]
    for a in set(itertools.permutations(['r','r','b'])):
        for b in set(itertools.permutations(['r','b','g'])):
            z=dict(zip('abcdefg',a+b+('g',)))
            alternatives.append(finite_q([z[x] for x in selected]))
    assert got['expected_q']==pytest.approx(np.mean(alternatives))
    assert got['q_within']+got['q_between']==pytest.approx(finite_q([states[x] for x in selected]))
    assert got['n_movable_groups']==2


def test_pair_null_preserves_second_trait_and_two_trait_mask():
    a=dict(zip('abcdefg',['r','r','b','r','b','g','g']))
    b=dict(zip('abcdefg',['x','y','x','y','x','y','x']))
    groups=dict(zip('abcdefg',['A','A','A','B','B','B','C']))
    selected=list('abdeg')
    got=TaxonNull(a,groups).evaluate(selected,second=b)
    joint=[];coupling=[]
    for u in set(itertools.permutations(['r','r','b'])):
        for v in set(itertools.permutations(['r','b','g'])):
            first=dict(zip('abcdefg',u+v+('g',)))
            av=[first[x] for x in selected];bv=[b[x] for x in selected]
            qj=finite_q(list(zip(av,bv)))
            qa=finite_q(av);qb=finite_q(bv)
            d=.8*((1-qa/.8)*(1-qb/.8)-(1-qj/.8))
            joint.append(qj);coupling.append(d)
    assert got['expected_q']==pytest.approx(np.mean(joint))
    assert got['expected_pairing']==pytest.approx(np.mean(coupling))
    assert got['q_within']+got['q_between']==pytest.approx(finite_q([(a[x],b[x]) for x in selected]))


def test_complete_group_or_constant_group_is_uninformative():
    states=dict(zip('abcdef',['r','r','r','b','b','b']))
    fixed=dict(zip('abcdef',['A','A','A','B','B','B']))
    z=TaxonNull(states,fixed).evaluate(list('abcde'))
    assert z['n_movable_groups']==0
    assert z['expected_q']==pytest.approx(finite_q([states[x] for x in 'abcde']))
    varied=dict(zip('abcdef',['r','b','g','r','b','g']))
    z=TaxonNull(varied,fixed).evaluate(list(varied))
    assert z['n_movable_groups']==0
    assert z['expected_q']==pytest.approx(finite_q(list(varied.values())))


@pytest.mark.parametrize('counts,missing',[([2,1,0],3),([0,0,0],5),([4,0],2),([2,2],1),([3,2],0)])
def test_missing_value_bounds_are_sharp_by_exhaustive_enumeration(counts,missing):
    n=sum(counts)+missing
    lo,hi=missing_bounds(counts,n)
    qs=[]
    for addition in compositions(missing,len(counts)):
        labels=[j for j,(a,b) in enumerate(zip(counts,addition)) for _ in range(a+b)]
        qs.append(finite_q(labels))
    assert lo==pytest.approx(min(qs))
    assert hi==pytest.approx(max(qs))


def test_missing_native_scope_and_unknown_taxonomy_never_zero_filled():
    assert missing_bounds([0,0],0)==(None,None)
    assert missing_bounds([2,1],3)==(None,None)
    with pytest.raises(ValueError):missing_bounds([3,4],6)
    states=dict(zip('abcdef',['r','b','r','b','r','b']))
    z=TaxonNull(states,{'a':'A','b':'A'}).evaluate(list('abcde'))
    assert z['n_mapped_species']==2
    assert z['n_movable_groups']==0
    assert z['expected_q']==pytest.approx(finite_q([states[x] for x in 'abcde']))


def test_full_cohort_stays_when_only_one_trait_is_observed():
    states={'flower_primary_color':{f'G{i} sp':('w' if i%2 else 'r') for i in range(7)},
            'floral_form':{f'G{i} sp':('open' if i%2 else 'tube') for i in range(5)}}
    rows=[]
    for i in range(6):
        for j in range(7 if i<3 else 3):
            rows.append({'island_id':str(i),'accepted_species':f'G{j} sp',
                         'origin_status':'native' if i%2==0 else 'unresolved'})
    flora=pd.DataFrame(rows)
    geo=pd.DataFrame([{'island_id':str(i),'analysis_regime':'tropical',
        'geological_origin':'oceanic' if i<3 else 'continental','spatial_block':str(i)} for i in range(5)])
    out=analyse(states,flora,geo,{},pairs=(('flower_primary_color','floral_form'),),draws=99)
    assert out['n_original_flora_islands']==6
    assert out['n_private_island_scope_endpoint_rows']==6*2*3
    assert out['n_flora_islands_missing_geography']==1
    color=next(x for x in out['global_endpoint_support'] if x['scope']=='all_recorded' and x['endpoint']=='flower_primary_color')
    assert color['n_original_islands']==6
    assert color['n_n5_observed']==3
    assert color['n_with_any_observed_trait']==6
    assert color['n_n5_family_null']==color['n_n5_genus_null']==color['n_n5_observed']
    native=next(x for x in out['global_endpoint_support'] if x['scope']=='original_native' and x['endpoint']=='flower_primary_color')
    assert native['n_zero_native_or_empty_scope']==3
    assert out['claims']['no_missing_trait_filled_as_observed']


def test_baseline_reconciliation_catches_changed_mask():
    baseline=dict(n_original_global_distinct_island_floras=10,
        marginal_trait_and_pairwise_coverage=[dict(source_plant_status_scope='all_recorded',
            trait_name='color',n_islands_with_at_least_five_direct_species=7,
            n_all_original_flora_islands=10,raw_observed_island_equal_finite5_mean=.3)])
    result=dict(n_original_flora_islands=10,global_endpoint_support=[dict(scope='all_recorded',
        endpoint='color',n_n5_observed=7,n_original_islands=10,mean_raw_q=.3)])
    assert reconcile_baseline(result,baseline)['max_absolute_raw_mean_difference']==0.
    result['global_endpoint_support'][0]['n_n5_observed']=6
    with pytest.raises(ValueError,match='mask changed'):reconcile_baseline(result,baseline)
