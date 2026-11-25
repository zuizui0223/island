import pandas as pd
import pytest
from fd.functional_matched import validate_matches, matched_retention, safe_matched_summary
from fd.functional_turnover import TRAITS


def _traits():
    out = []
    for island, taxon, genus, value in [
        ('A','A1','GA',1),('A','A2','GA',1),('A','A3','GA',1),
        ('N','N1','GN',1),('N','N2','GN',1),('N','N3','GN',1),
        ('C','C1','GC',0),('C','C2','GC',0),('C','C3','GC',0),
    ]:
        out.append({'island_id':island,'taxon_name':taxon,'genus':genus,
                    **{t:value for t in TRAITS}})
    return pd.DataFrame(out)


def _pairs():
    return pd.DataFrame([
      {'match_id':'m1','island_a':'A','island_b':'N','source_group':'S1','donor_group':'S1',
       'role':'local','distance_km':50},
      {'match_id':'m1','island_a':'A','island_b':'C','source_group':'S1','donor_group':'S2',
       'role':'distance_control','distance_km':70},
    ])


def test_local_advantage_and_independent_domain_computation():
    z=matched_retention(_traits(), _pairs(), min_carriers=2)
    for scope in ['all_taxa','species_exclusive','genus_exclusive']:
        row=z[(z.trait_scope==scope)&(z.trait=='DOMAIN_MEAN')]
        assert len(row)==2
        assert (row.local_advantage==1).all()
        assert (row.n_component_traits==3).all()


def test_distance_matching_is_mandatory():
    p=_pairs();p.loc[1,'distance_km']=250
    with pytest.raises(ValueError,match='distance'):
        validate_matches(p)


def test_repeated_focal_is_not_independent():
    p=_pairs()
    more=p.copy()
    more['match_id']='m2'
    more['island_b']=['N2','C2']
    with pytest.raises(ValueError,match='Each focal'):
        validate_matches(pd.concat([p,more],ignore_index=True))


def test_missing_domain_support_cannot_be_imputed():
    traits=_traits()
    traits.loc[traits.island_id=='N','autonomous_selfing']=None
    z=matched_retention(traits, _pairs(), min_carriers=2)
    r=z[(z.trait_scope=='genus_exclusive') & (z.domain=='reproductive_assurance') & (z.trait=='DOMAIN_MEAN')]
    assert len(r)==1
    assert r.iloc[0].n_component_traits==2
    traits.loc[traits.island_id=='N','selfing_mating_system']=None
    z=matched_retention(traits, _pairs(), min_carriers=2)
    assert z[(z.domain=='reproductive_assurance') & (z.trait=='DOMAIN_MEAN')].empty


def test_group_privacy_gate_and_weighting():
    rows=pd.DataFrame([{'match_id':f'm{i}', 'source_group':f'g{i//2}',
                       'trait_scope':'genus_exclusive', 'domain':'reproductive_assurance',
                       'trait':'DOMAIN_MEAN','local_advantage':float(i//2)} for i in range(10)])
    assert safe_matched_summary(rows,min_groups=6).empty
    result=safe_matched_summary(rows,min_groups=5)
    assert len(result)==1
    assert result.iloc[0].n_source_groups==5
    assert result.iloc[0].mean_equal_group_local_advantage==2.
    assert not {'match_id','source_group','island_a','island_b'} & set(result.columns)

def test_knowable_source_identity_is_enforced():
    p=_pairs();p.loc[1,'donor_group']='S1'
    with pytest.raises(ValueError,match='another source_group'):
        validate_matches(p)
    p=_pairs();p.loc[0,'donor_group']='UNKNOWN'
    with pytest.raises(ValueError,match='Unknown donor_group'):
        validate_matches(p)


def test_cli_writes_only_safe_empty_aggregate(tmp_path):
    from fd.functional_matched import main
    traits=tmp_path/'traits.csv';matches=tmp_path/'matches.csv';output=tmp_path/'out'
    _traits().to_csv(traits,index=False)
    _pairs().to_csv(matches,index=False)
    main(['--private-traits',str(traits),'--private-matches',str(matches),
          '--output',str(output),'--min-carriers','2'])
    import json
    assert sorted(p.name for p in output.iterdir())==['functional_matched_grouped.csv','status.json']
    assert json.loads((output/'status.json').read_text())['released_comparison_rows']==0
    assert (output/'functional_matched_grouped.csv').read_text().startswith('trait_scope,domain,trait,')
