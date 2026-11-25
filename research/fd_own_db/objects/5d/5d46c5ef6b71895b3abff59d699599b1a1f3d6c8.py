"""Original island pairwise composition source and mask integrity tests.

The math and source reproducibility exercises have no protected real plant IDs.
"""
import itertools
import json

import numpy as np
import pandas as pd
import pytest

from fd.island_functional_beta import (
    incidence_partition,prepare_endpoint,compare_block,
    analyse,strict_traits
)

def test_presence_turnover_and_nestedness_algebra():
    sor,turn,nest=incidence_partition(
        np.array([2,0,1]),np.array([2,2,3]),np.array([4,2,3]))
    assert sor.tolist()==pytest.approx([1/3,1,2/3])
    assert turn.tolist()==pytest.approx([0,1,2/3])
    assert nest.tolist()==pytest.approx([1/3,0,0])
    assert np.allclose(sor,turn+nest)
    with pytest.raises(ValueError):
        incidence_partition(np.array([3]),np.array([2]),np.array([4]))

def test_equal_alpha_can_have_no_shared_categories():
    values={'G a':'white','G b':'violet','G c':'yellow','G d':'red'}
    b=compare_block(prepare_endpoint([
        {'G a','G b'},{'G c','G d'}],values),0,2)
    assert b['overlap'][0,1]==pytest.approx(0)
    assert b['observed_taxonomic_sorensen'][0,1]==1
    assert b['category_turnover'][0,1]==1

def test_different_original_species_can_have_identical_flower_composition():
    values={'A a':'white','A b':'violet','B c':'white','B d':'violet'}
    b=compare_block(prepare_endpoint([
        {'A a','A b'},{'B c','B d'}],values),0,2)
    assert b['overlap'][0,1]==pytest.approx(1)
    assert b['observed_taxonomic_sorensen'][0,1]==1
    assert b['squared_frequency_distance'][0,1]==0
    assert b['global_null_squared_distance'][0,1]>0
    assert b['genus_null_squared_distance'][0,1]==pytest.approx(0)

def test_exact_shared_species_have_a_single_original_trait_state():
    values={'G a':'w','G b':'v','G c':'w','G d':'r'}
    b=compare_block(prepare_endpoint([
        {'G a','G b'},{'G a','G b'}],values),0,2)
    assert b['overlap'][0,1]==1
    assert b['global_null_squared_distance'][0,1]==pytest.approx(0)
    assert b['genus_null_squared_distance'][0,1]==pytest.approx(0)

def test_analytic_nulls_agree_with_exhaustive_original_state_permutations():
    names=['A a','A b','A c','B d','B e','C f']
    values=['r','r','b','b','g','g']
    states=dict(zip(names,values))
    rosters=[set(names[:3]),set(names[2:5]),set(names[3:])]
    out=compare_block(prepare_endpoint(rosters,states),0,3)
    def distance(mapping):
        cats=sorted(set(values))
        P=np.array([[sum(mapping[x]==c for x in island)/len(island)
           for c in cats] for island in rosters])
        return ((P[:,None,:]-P[None,:,:])**2).sum(2)
    null=np.mean([distance(dict(zip(names,v)))
        for v in set(itertools.permutations(values))],axis=0)
    assert np.allclose(out['global_null_squared_distance'],null,atol=1e-12)
    null2=[]
    for a in set(itertools.permutations(values[:3])):
        for b in set(itertools.permutations(values[3:5])):
            null2.append(distance(dict(zip(names,a+b+(values[5],)))))
    assert np.allclose(out['genus_null_squared_distance'],
                       np.mean(null2,axis=0),atol=1e-12)

def test_joint_profiles_only_use_real_observed_two_trait_combinations():
    values={'A a':('white','tube'),'A b':('white','open'),
            'B c':('white','tube'),'B d':('white','open')}
    prepared=prepare_endpoint([{'A a','A b'},{'B c','B d'}],values)
    assert len(prepared['categories'])==2
    assert compare_block(prepared,0,2)['overlap'][0,1]==1

def _fixture():
    names=[f'G{i} sp' for i in range(8)]
    traits={
      'flower_primary_color':{n:('white' if i%2 else 'violet')
                             for i,n in enumerate(names)},
      'floral_form':{n:('open' if i%2 else 'tube')
                    for i,n in enumerate(names[:6])}}
    flora=pd.DataFrame([dict(
        island_id=f'I{i}',accepted_species=name,
        origin_status='native' if i%2==0 else 'unresolved')
        for i in range(6) for name in
        (names[:6] if i<4 else names[4:])])
    geo=pd.DataFrame([dict(
        island_id=f'I{i}',analysis_regime='tropical',
        geological_origin='oceanic' if i<3 else 'continental',
        spatial_block=f'b{i//2}') for i in range(5)])
    return traits,flora,geo

def test_full_original_population_and_pairwise_only_recorded_traits():
    states,flora,geo=_fixture()
    out=analyse(states,flora,geo,
        pairs=(('flower_primary_color','floral_form'),),block_size=2)
    assert out['n_original_islands']==6
    assert out['n_without_geography']==1
    assert len(out['endpoint_coverage'])==6
    color=next(x for x in out['endpoint_coverage']
        if x['scope']=='all_recorded' and x['endpoint']=='flower_primary_color')
    assert color['n_original_islands']==6
    assert color['n_n5_islands']==4
    assert color['n_any_direct_islands']==6
    assert color['n_unordered_n5_pairs']==6
    row=next(x for x in out['pair_summaries']
        if x['scope']=='all_recorded' and x['endpoint']=='flower_primary_color'
        and x['region']=='GLOBAL_ALL' and
        x['geology_pair']=='all' and x['subset']=='all')
    assert row['n_pairs']==6 and row['n_unique_islands']==4
    assert row['mean_functional_overlap']==pytest.approx(1)
    assert row['island_equal_mean_functional_overlap']==pytest.approx(1)
    assert out['claims']['no_pair_count_treated_as_independent_replications']
    assert out['claims']['no_intraspecific_trait_change_measured']
    json.dumps(out,allow_nan=False)
    assert not any(k in json.dumps(out) for k in [
        '"I0"','"G0 sp"','"island_id"','"accepted_species"'])

def test_unknown_native_status_is_not_plant_absence():
    states,flora,geo=_fixture()
    flora['origin_status']='unresolved'
    out=analyse(states,flora,geo,pairs=(),block_size=3)
    native=[r for r in out['endpoint_coverage']
            if r['scope']=='original_native']
    assert all(r['n_n5_islands']==0 for r in native)
    assert not any(r['scope']=='original_native' for r in out['pair_summaries'])

def test_reordering_and_block_size_do_not_change_descriptive_result():
    states,flora,geo=_fixture()
    a=analyse(states,flora,geo,pairs=(),block_size=2)
    b=analyse(states,flora.sample(frac=1,random_state=3),
              geo,pairs=(),block_size=5)
    assert a['endpoint_coverage']==b['endpoint_coverage']
    for x,y in zip(a['pair_summaries'],b['pair_summaries']):
        assert x.keys()==y.keys()
        for key,value in x.items():
            if isinstance(value,float):
                assert value==pytest.approx(y[key],abs=1e-12)
            else:
                assert value==y[key]

def test_one_trait_record_does_not_require_other_trait_state():
    d=pd.DataFrame([
       dict(accepted_species='A a',trait_name='flower_primary_color',
            normalized_value='white',evidence_scope='species_direct',
            evidence_quality='high'),
       dict(accepted_species='B b',trait_name='floral_form',
            normalized_value='tube',evidence_scope='species_direct',
            evidence_quality='medium'),
       dict(accepted_species='C c',trait_name='floral_form',
            normalized_value='tube|open',evidence_scope='species_direct',
            evidence_quality='high')
    ])
    got=strict_traits(d)
    assert got['flower_primary_color']=={'A a':'white'}
    assert got['floral_form']=={'B b':'tube'}
    with pytest.raises(ValueError):
        strict_traits(pd.concat([d,d.iloc[:1]]))

def test_missing_trait_coverage_is_observed_information_not_zero_trait():
    states={f'G{i} sp':('white' if i%2 else 'blue') for i in range(6)}
    flora=pd.DataFrame([dict(
        island_id=i,accepted_species=f'G{k} sp',origin_status='native')
        for i in ('A','B') for k in range(8)])
    geo=pd.DataFrame([dict(
        island_id=i,analysis_regime='tropical',
        geological_origin='oceanic',spatial_block=i) for i in ('A','B')])
    out=analyse({'flower_primary_color':states},flora,geo,
                pairs=(),block_size=1)
    q=next(x for x in out['pair_summaries']
        if x['scope']=='all_recorded' and
        x['region']=='GLOBAL_ALL' and
        x['geology_pair']=='all' and x['subset']=='all')
    assert q['mean_pair_recorded_species_fraction']==pytest.approx(.75)
    assert q['mean_absolute_recording_fraction_difference']==0


def test_strict_full_flora_genus_disjoint_vs_species_disjoint():
    """Synthetic identities: two distinct species-sharing-free comparisons,
    one genus-sharing pair, and two zero-genus pairs."""
    from fd.island_functional_beta import accepted_botanical_genus, certified_genus_sets
    assert accepted_botanical_genus('A one')=='A'
    assert accepted_botanical_genus('Unresolved') is None
    assert accepted_botanical_genus('A sp.') is None
    names={
        'I0':[f'A a{k}' for k in range(5)],
        'I1':[f'A b{k}' for k in range(5)],
        'I2':[f'B c{k}' for k in range(5)]}
    st={'flower_primary_color':{sp:'white'
        for plants in names.values() for sp in plants}}
    flora=pd.DataFrame([dict(island_id=i,accepted_species=sp,
        origin_status='native') for i,plants in names.items() for sp in plants])
    geo=pd.DataFrame([dict(island_id=i,analysis_regime='tropical',
        geological_origin='unresolved',spatial_block=f'block{k}')
        for k,i in enumerate(names)])
    r=analyse(st,flora,geo,pairs=(),block_size=1)
    rows={x['subset']:x for x in r['pair_summaries'] if
        x['scope']=='all_recorded' and x['endpoint']=='flower_primary_color'
        and x['region']=='GLOBAL_ALL' and x['geology_pair']=='all'}
    assert rows['all']['n_pairs']==3
    assert rows['no_shared_full_flora_species']['n_pairs']==3
    assert rows['no_shared_full_flora_species_genus_shared']['n_pairs']==1
    assert rows['no_shared_full_flora_genera_strict']['n_pairs']==2
    assert rows['no_shared_direct_genera']['n_pairs']==2
    assert rows['genus_disjoint_different_blocks']['n_pairs']==2
    assert rows['no_shared_full_flora_genera_strict']['mean_functional_overlap']==1
    assert r['claims']['unresolved_accepted_genus_quarantined_not_assumed_distinct']
    gens,qual=certified_genus_sets([{'A one','Unknown'},{'B two'}])
    assert gens==[{'A'},{'B'}] and qual.tolist()==[False,True]

def test_missing_full_roster_genus_cannot_inflate_strict_disjoint_pairs():
    names={'I0':[f'A a{k}' for k in range(5)],
        'I1':[f'B b{k}' for k in range(5)]}
    st={'flower_primary_color':{sp:'white'
        for plants in names.values() for sp in plants}}
    flora=pd.DataFrame([dict(island_id=i,accepted_species=sp,
        origin_status='native') for i,plants in names.items() for sp in plants])
    flora=pd.concat([flora,pd.DataFrame([dict(
        island_id='I1',accepted_species='Unresolved',origin_status='native')])],
        ignore_index=True)
    geo=pd.DataFrame([dict(island_id=i,analysis_regime='tropical',
        geological_origin='unresolved',spatial_block=i)
        for i in names])
    r=analyse(st,flora,geo,pairs=(),block_size=1)
    rows={x['subset']:x for x in r['pair_summaries'] if
        x['scope']=='all_recorded' and x['endpoint']=='flower_primary_color'
        and x['region']=='GLOBAL_ALL' and x['geology_pair']=='all'}
    assert rows['no_shared_full_flora_species']['n_pairs']==1
    assert 'no_shared_full_flora_genera_strict' not in rows
    assert rows['no_shared_direct_genera']['n_pairs']==1
