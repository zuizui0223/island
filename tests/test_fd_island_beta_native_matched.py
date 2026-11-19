"""Matched native/all-source comparison, synthetic only."""
import json
import pandas as pd
import pytest

from fd_island_beta_native_matched import analyze


def _data():
    states={'flower_primary_color':{},'floral_form':{}}
    flora=[];geo=[]
    for i in range(12):
        region='northern_midlatitude' if i<6 else 'tropical'
        block=f'{region}_b{i//3}'
        for k in range(8):
            plant=f'Gen{k%3} sp{i}_{k}'
            states['flower_primary_color'][plant]=(
                'white' if ((i+k)%4) else 'yellow')
            states['floral_form'][plant]=(
                'tube' if ((2*i+k)%3) else 'open')
            flora.append(dict(island_id=f'I{i}', accepted_species=plant,
                              origin_status='native' if k<6 else 'introduced'))
        geo.append(dict(island_id=f'I{i}', analysis_regime=region,
                        spatial_block=block, island_latitude=0 if i<6 else 20,
                        island_longitude=(i%6)*2))
    return states,pd.DataFrame(flora),pd.DataFrame(geo)


def test_reuses_exact_same_full_species_disjoint_pairs_in_both_scopes():
    s,f,g=_data()
    out=analyze(s,f,g,draws=2,min_pairs=1,min_islands=2,min_blocks=2)
    assert out['n_original_botanical_islands']==12
    assert out['n_original_native_direct_pair_eligible_islands']==12
    assert out['n_matched_eligible_same_region_other_block_full_species_disjoint_pairs']==18
    v=next(x for x in out['comparisons'] if (
        x['protocol']=='full_direct' and x['region']=='tropical'
        and x['distance_class']=='all_distances'))
    assert v['n_pairs']==9
    assert v['n_islands']==6
    assert v['n_physical_blocks']==2
    assert v['all_minus_native_paired_mean_residual_D2']==pytest.approx(
        v['all_recorded_mean_residual_D2']-v['native_mean_residual_D2'])
    assert out['claims']['identical_physical_island_pairs_within_each_stratum']
    assert not any(key in json.dumps(out) for key in [
        '"island_id"','"accepted_species"','"spatial_block"','"I0"'])


def test_same_original_full_species_records_are_excluded_even_if_geo_disjoint():
    s,f,g=_data()
    r=f.loc[f.island_id=='I0'].copy()
    r['island_id']='I3'
    f=pd.concat([f.loc[f.island_id!='I3'],r],ignore_index=True)
    out=analyze(s,f,g,draws=1,min_pairs=1,min_islands=2,min_blocks=2)
    assert out['n_matched_eligible_same_region_other_block_full_species_disjoint_pairs']==17


def test_equal_five_draws_are_reproducible_not_new_islands():
    s,f,g=_data()
    a=analyze(s,f,g,seed=17,draws=3,min_pairs=1,min_islands=2,min_blocks=2)
    b=analyze(s,f.sample(frac=1,random_state=2),g,seed=17,draws=3,
              min_pairs=1,min_islands=2,min_blocks=2)
    assert a==b
    q=[x for x in a['comparisons'] if x['protocol']=='equal_5_direct']
    assert q and all(x['n_montecarlo_thinnings']==3 for x in q)
    assert all(x['n_islands']<=12 and x['n_pairs']<=18 for x in q)


def test_full_same_plant_scope_gives_zero_matched_delta():
    s,f,g=_data()
    f['origin_status']='native'
    q=analyze(s,f,g,draws=1,min_pairs=1,min_islands=2,min_blocks=2)
    for row in q['comparisons']:
        if row['protocol']=='full_direct':
            assert row['all_minus_native_paired_mean_residual_D2']==pytest.approx(0)
            assert row['all_recorded_observed_D2']==pytest.approx(
                row['native_observed_D2'])


def test_invalid_source_rejected():
    s,f,g=_data()
    with pytest.raises(ValueError):
        analyze(s,pd.concat([f,f.iloc[:1]]),g)
    with pytest.raises(ValueError):
        analyze(s,f,g.drop(columns=['island_longitude']))
    with pytest.raises(ValueError):
        analyze(s,f,g,draws=0)
