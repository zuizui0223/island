"""Pure synthetic regressions; never require protected original botanical rows."""
import json
import numpy as np
import pandas as pd
import pytest

from fd_island_beta_geographic_null import Aggregate, analyze


def test_block_deletion_uses_each_dyad_exactly_once_for_each_endpoint():
    blocks = np.array([0, 1, 2, 3])
    a = Aggregate(4, blocks, 4)
    a.add(np.array([0, 0, 1]), np.array([1, 2, 3]),
          np.array([0.6, 0.7, 0.9]), np.array([0.5, 0.4, 0.4]))
    z = a.export(1, 2, 2)
    assert z['n_pairs'] == 3
    assert z['n_islands'] == 4
    assert z['pair_equal_mean_residual_D2'] == pytest.approx(.3)
    assert z['island_equal_mean_residual_D2'] == pytest.approx(
        ((.1+.3)/2 + (.1+.5)/2 + .3 + .5)/4)
    # Deleting block 0 leaves the 1--3 dyad, not double-subtracted.
    assert z['leave_one_physical_block_out_pair_mean_residual_D2_max'] == pytest.approx(.5)


def _fixture():
    states = {'flower_primary_color': {}, 'floral_form': {}}
    flora = []
    geo = []
    for i in range(12):
        region = 'northern_midlatitude' if i < 6 else 'tropical'
        b = f'{region}_b{i//3}'
        for k in range(6):
            name = f'Gen{i%4} sp{i}_{k}'
            states['flower_primary_color'][name] = 'white' if (i+k)%2 else 'blue'
            states['floral_form'][name] = 'open' if (i+2*k)%3 else 'tube'
            flora.append(dict(island_id=f'I{i}', accepted_species=name,
                              origin_status='native'))
        geo.append(dict(island_id=f'I{i}', analysis_regime=region,
                        spatial_block=b, geological_origin='unresolved',
                        island_latitude=0 if i < 6 else 20,
                        island_longitude=(i%6)*2))
    return states, pd.DataFrame(flora), pd.DataFrame(geo)


def test_regions_blocks_full_species_and_same_endpoint_source():
    states, flora, geo = _fixture()
    result = analyze(states, flora, geo, block_size=3,
                     min_pairs=1, min_islands=2, min_blocks=2)
    assert result['n_original_botanical_islands'] == 12
    assert result['eligible_islands_by_scope'] == {
        'all_recorded': 12, 'original_native': 12}
    assert all(row['region'] in ('northern_midlatitude', 'tropical')
               for row in result['comparisons'])
    source = [x for x in result['comparisons']
              if x['scope']=='all_recorded' and
                 x['region']=='tropical' and
                 x['support']=='n5' and
                 x['distance_class']=='all_distances']
    assert len(source) == 1
    assert source[0]['n_pairs'] == 9
    assert source[0]['n_physical_blocks'] == 2
    assert source[0]['n_islands'] == 6
    assert source[0]['pair_equal_mean_residual_D2'] == pytest.approx(
        source[0]['pair_equal_mean_observed_D2'] -
        source[0]['pair_equal_mean_genus_null_D2'])
    assert not any(x['support']=='both_n10' for x in result['comparisons'])
    assert json.dumps(result, allow_nan=False)
    for s in ['"island_id"', '"accepted_species"', '"spatial_block"']:
        assert s not in json.dumps(result)


def test_species_shared_and_same_block_pair_not_counted():
    states, flora, geo = _fixture()
    # Two island floras in different blocks but identical original taxa:
    # they must not enter any no-shared-accepted-species comparison.
    f2 = flora.loc[flora.island_id.eq('I0')].copy()
    f2['island_id'] = 'I3'
    flora = pd.concat([flora.loc[~flora.island_id.eq('I3')], f2],
                      ignore_index=True)
    result = analyze(states, flora, geo, min_pairs=1, min_islands=2,
                     min_blocks=2, block_size=1)
    main = next(x for x in result['comparisons']
                if x['scope']=='all_recorded' and
                   x['region']=='northern_midlatitude' and
                   x['support']=='n5' and
                   x['distance_class']=='all_distances')
    assert main['n_pairs'] == 8


def test_coverage_gate_is_outcome_blind():
    states, flora, geo = _fixture()
    original = analyze(states, flora, geo, min_pairs=1,
                       min_islands=2, min_blocks=2)
    # Add 60 non-trait-coded original plant names to every island.
    # The source species set changes; recorded outcome categories do not.
    extra = pd.DataFrame([
        dict(island_id=f'I{i}', accepted_species=f'Unknown{k}_is{i}',
             origin_status='unresolved')
        for i in range(12) for k in range(60)])
    more = pd.concat([flora, extra], ignore_index=True)
    changed = analyze(states, more, geo, min_pairs=1,
                      min_islands=2, min_blocks=2)
    assert not any(x['support']=='both_direct_fraction_ge_0p10'
                   for x in changed['comparisons'] if x['scope']=='all_recorded')
    assert any(x['support']=='both_direct_fraction_ge_0p10'
               for x in original['comparisons'] if x['scope']=='all_recorded')
    assert changed['eligible_islands_by_scope']==original['eligible_islands_by_scope']


def test_invalid_input_rejected():
    states, flora, geo = _fixture()
    with pytest.raises(ValueError):
        analyze(states, pd.concat([flora, flora.iloc[:1]]), geo)
    with pytest.raises(ValueError):
        analyze(states, flora, geo.drop(columns=['island_longitude']))
    with pytest.raises(ValueError):
        analyze(states, flora, geo, block_size=0)
