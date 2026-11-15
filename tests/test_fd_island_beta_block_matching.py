"""Synthetic regressions for posthoc block-disjoint floral reference audit."""
import json
import numpy as np
import pandas as pd
import pytest
from fd_island_beta_block_matching import (
    match_blocks, graph_support, summarize, analyze,
)


def test_maximum_cardinality_never_reuses_any_spatial_block():
    a=[0,0,1,2,3,3,5]
    b=[1,2,3,3,4,5,6]
    x=match_blocks(a,b,seed=11)
    assert len(x)==3
    selected=np.r_[np.array(a)[x],np.array(b)[x]]
    assert len(set(selected))==2*len(x)
    assert np.array_equal(x,match_blocks(a,b,seed=11))
    assert graph_support(a,b)['n_candidate_blocks']==7
    with pytest.raises(ValueError):
        match_blocks([0],[0])


def test_reselecting_dyads_does_not_use_phenotype_outcome():
    a=np.array([0,0,0,2,2,3])
    b=np.array([1,1,3,3,4,4])
    x=match_blocks(a,b,seed=128)
    y=match_blocks(a,b,seed=128)
    assert np.array_equal(x,y)
    score=np.array([.1,-.2,.5,.7,1.,-1.])
    score_reversed=score[::-1]
    assert np.array_equal(x,y)
    assert np.mean(score[x])!=np.mean(score_reversed[x])


def test_matching_aggregates_are_not_biological_replications():
    edges={
        'i':np.array([0,1,2,3,4,5]),
        'j':np.array([6,7,8,9,10,11]),
        'block_i':np.array([0,0,1,2,3,3]),
        'block_j':np.array([1,2,2,3,4,5]),
        'obs':np.array([.4,.3,.6,.7,.8,.5]),
        'expected_g':np.array([.2]*6),
        'expected_r':np.array([.3]*6),
    }
    x=summarize(edges,draws=5,min_pairs=1,min_islands=2,min_blocks=2)
    assert x['n_candidate_pairs']==6
    assert x['matched_block_disjoint_pairs_min']>=1
    assert x['matched_block_disjoint_pairs_max']<=3
    assert x['n_outcome_blind_matching_variants']==5
    assert x['candidate_pair_mean_region_minus_global']==pytest.approx(-.1)
    assert x['selected_dyad_margins_not_statistical_independent_replications']
    assert summarize(edges,draws=2,min_pairs=100) is None


def _source():
    states={'flower_primary_color':{},'floral_form':{}}
    rows=[];geo=[]
    for i in range(12):
        trop=i<6
        region='tropical' if trop else 'northern_midlatitude'
        j=i%6
        for k in range(6):
            s=f'Gen{k%3} sp{i}_{k}'
            states['flower_primary_color'][s]='white' if (i+k)%3 else 'yellow'
            states['floral_form'][s]='tube' if (i+2*k)%4 else 'open'
            rows.append(dict(island_id=f'I{i}',accepted_species=s,
                             origin_status='native' if k<5 else 'introduced'))
        geo.append(dict(island_id=f'I{i}',analysis_regime=region,
                        spatial_block=f'{region}_block{j//2}',
                        island_latitude=0 if trop else 40,
                        island_longitude=j*9))
    return states,pd.DataFrame(rows),pd.DataFrame(geo)


def test_no_triple_gate_and_private_names_removed():
    states,flora,geo=_source()
    a=analyze(states,flora,geo,draws=5,min_pairs=1,min_islands=2,min_blocks=2)
    assert a['n_original_botanical_islands']==12
    assert a['eligible_n5_direct_islands_by_scope']=={
        'all_recorded':12,'original_native':12}
    assert any(x['scope']=='all_recorded' and x['region']=='tropical'
               for x in a['comparisons'])
    assert all(x['n_candidate_pairs']>=1 for x in a['comparisons'])
    assert a['claims']['comparison_requires_full_flora_species_disjointness']
    q=json.dumps(a,allow_nan=False)
    assert all(s not in q for s in (
        '"island_id"','"accepted_species"','"spatial_block"','"I0"'))


def test_full_original_species_overlap_masks_native_too():
    states,flora,geo=_source()
    rows=[dict(island_id='I1',accepted_species='Alien shared',
               origin_status='unresolved'),
          dict(island_id='I2',accepted_species='Alien shared',
               origin_status='unresolved')]
    change=pd.concat([flora,pd.DataFrame(rows)],ignore_index=True)
    a=analyze(states,flora,geo,draws=2,min_pairs=1,min_islands=2,min_blocks=2)
    b=analyze(states,change,geo,draws=2,min_pairs=1,min_islands=2,min_blocks=2)
    av=next(x for x in a['comparisons'] if x['scope']=='original_native'
            and x['region']=='tropical' and x['support']=='n5')
    bv=next(x for x in b['comparisons'] if x['scope']=='original_native'
            and x['region']=='tropical' and x['support']=='n5')
    assert bv['n_candidate_pairs']==av['n_candidate_pairs']-1


def test_invalid_inputs_fail_closed():
    states,flora,geo=_source()
    with pytest.raises(ValueError):
        analyze(states,pd.concat([flora,flora.iloc[:1]]),geo)
    with pytest.raises(ValueError):
        analyze(states,flora,geo.drop(columns=['spatial_block']))
    with pytest.raises(ValueError):
        summarize({'i':np.array([1])},draws=1)
