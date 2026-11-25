import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

PATH=Path(__file__).parents[1]/'src/fd/island_floral_beta.py'

def module():
    assert PATH.exists(), 'Island functional-beta analysis has not been implemented'
    s=importlib.util.spec_from_file_location('island_floral_beta',PATH)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def test_similarity_distinguishes_composition_from_equal_alpha():
    m=module()
    x=np.array([[.5,.5,0,0],[0,0,.5,.5],[.5,.5,0,0]])
    v=m.composition_metrics(x,np.array([0,0]),np.array([1,2]))
    assert v['affinity'].tolist()==pytest.approx([0,1])
    assert v['profile_sorensen'].tolist()==[1,0]
    assert m.exact_all_pair_affinity(x)==pytest.approx(1/3)

def test_turnover_nestedness_identity_and_subset():
    m=module()
    x=np.array([[1,0,0],[.5,.5,0],[0,0,1]])
    d=m.composition_metrics(x,np.array([0,0]),np.array([1,2]))
    assert d['profile_turnover'].tolist()==[0,1]
    assert d['profile_nestedness'].tolist()==pytest.approx([1/3,0])
    assert np.allclose(d['profile_sorensen'],d['profile_turnover']+d['profile_nestedness'])

def test_every_eligible_island_is_focal_not_selected_by_flower_states():
    m=module();a,b=m.random_edges([str(i) for i in range(21)],k=6,seed=3)
    a2,b2=m.random_edges([str(i) for i in range(21)],k=6,seed=3)
    assert np.array_equal(a,a2) and np.array_equal(b,b2)
    assert len(a)==126 and np.all(np.bincount(a)==6) and not (a==b).any()
    assert len(set(zip(a,b)))==126

def test_whole_profiles_permuted_without_fabricated_pairs():
    m=module();codes=np.array([0,1,2,0,1,2]);groups=np.array([0,0,1,1,2,2])
    rng=np.random.default_rng(7);p=m.shuffle_codes(codes,groups,rng)
    assert sorted(p)==sorted(codes)
    for g in np.unique(groups):assert sorted(p[groups==g])==sorted(codes[groups==g])

def inputs():
    rows=[];states={'color':{},'form':{}}
    for i in range(8):
        for j in range(6):
            n=f'G{j%3} species{i}_{j}'
            rows.append(dict(island_id=str(i),accepted_species=n,origin_status='native' if i%2==0 else 'unresolved'))
            if i<6:states['color'][n]='white' if j%2 else 'purple'
            if i<4:states['form'][n]='open' if j%2 else 'tube'
    geo=pd.DataFrame([dict(island_id=str(i),analysis_regime='tropical',
        geological_origin='oceanic' if i<4 else 'continental',spatial_block=str(i),
        island_latitude=float(i),island_longitude=float(i)) for i in range(7)])
    return states,pd.DataFrame(rows),geo

def test_partial_masks_all_islands_same_species_disjoint_similarity():
    m=module();s,f,g=inputs()
    out=m.analyse(s,f,g,pairs=(('color','form'),),permutations=9,seed=2,primary=('color',),partners=3)
    assert out['n_original_flora_islands']==8
    a=next(x for x in out['endpoints'] if x['scope']=='all_recorded' and x['endpoint']=='color')
    assert a['n_original_islands']==8 and a['n_eligible_islands']==6
    assert a['n_focal_islands_in_uniform_graph']==6 and a['n_all_unordered_pairs']==15
    assert a['uniform_pairs']['mean_affinity']==pytest.approx(1)
    assert a['uniform_no_shared_flora']['n_focal_islands']==6
    assert a['uniform_no_shared_flora']['mean_affinity']==pytest.approx(1)
    b=next(x for x in out['endpoints'] if x['scope']=='original_native' and x['endpoint']=='color')
    assert b['n_eligible_islands']==3
    assert b['n_empty_original_scope']==4
    joint=next(x for x in out['endpoints'] if x['scope']=='all_recorded' and x['endpoint']=='color__form')
    assert joint['n_eligible_islands']==4
    assert out['claims']['no_intraspecific_change_estimated']

def test_same_species_is_never_relabelled_independently_across_islands():
    m=module();s,f,g=inputs()
    e=m.build_endpoint(s['color'],f,'all_recorded',g)
    assert e['incidence'].shape[1]==len(s['color'])
    codes=m.shuffle_codes(e['codes'],e['genera'],np.random.default_rng(4))
    freq=m.proportions(e['incidence'],codes,len(e['categories']))
    assert np.allclose(freq.sum(axis=1),1)

def test_baseline_reconciliation_rejects_changed_eligible_mask():
    m=module();s,f,g=inputs();o=m.analyse(s,f,g,pairs=(('color','form'),),permutations=0,seed=2,primary=(),partners=2)
    old=dict(n_original_global_distinct_island_floras=8,marginal_trait_and_pairwise_coverage=[])
    for e in o['endpoints']:
        old['marginal_trait_and_pairwise_coverage'].append(dict(source_plant_status_scope=e['scope'],trait_name=e['endpoint'],
            n_all_original_flora_islands=8,n_islands_with_at_least_five_direct_species=e['n_eligible_islands'],
            raw_observed_island_equal_finite5_mean=e['mean_raw_finite5_rao']))
    a=m.reconcile(o,old);assert a['max_raw_rao_difference']<1e-12
    old['marginal_trait_and_pairwise_coverage'][0]['n_islands_with_at_least_five_direct_species']=0
    with pytest.raises(ValueError):m.reconcile(o,old)
