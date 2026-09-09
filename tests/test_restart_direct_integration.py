import importlib.util
import json
import copy
from pathlib import Path
import pandas as pd
import pytest

root = Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location('integration',root/'scripts/integrate_reviewed_restart.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def sample():
    records=json.loads((root/'data/v2/staging/traits/restart_reviewed_20260907.json').read_text())
    records=[r for r in records if r['accepted_species']=='Neolitsea sericea']
    assert len(records)==1
    c=pd.DataFrame([dict(accepted_species='Neolitsea sericea',axis='reproductive_assurance',quality='',trait_names='',trait_composition='',source_groups='',source_lineages='')])
    d=pd.DataFrame(columns=['accepted_species','trait_name','quality'])
    o={'traits':{'mating_system':{'domain':'reproductive_assurance','allowed_values':['predominantly_outcrossing']}}}
    return c,d,records,o


def test_new_direct_cell():
    c,d,r,o=sample()
    after,ledger,change=module.integrate(c,d,r,o)
    assert after.quality.tolist()==['medium']
    assert after.trait_composition.iloc[0]=='mating_system=["predominantly_outcrossing"]'
    assert c.quality.iloc[0]==''
    assert len(ledger)==len(change)==1


def test_all_base_records_have_distinct_species_trait_cells():
    import yaml
    records=json.loads((root/'data/v2/staging/traits/restart_reviewed_20260907.json').read_text(encoding='utf-8'))
    c=pd.DataFrame([dict(accepted_species=r['accepted_species'],axis=r['axis'],quality='',trait_names='',trait_composition='',source_groups='',source_lineages='') for r in records])
    d=pd.DataFrame(columns=['accepted_species','trait_name','quality'])
    o=yaml.safe_load((root/'config/trait_ontology.yml').read_text(encoding='utf-8'))
    after,ledger,change=module.integrate(c,d,records,o)
    assert len(ledger)==len(change)==len(records)
    assert not ledger.duplicated(['accepted_species','trait_name']).any()
    assert after.quality.ne('').all()


def test_existing_low_is_not_overwritten():
    c,d,r,o=sample()
    c.loc[0,'quality']='low'
    with pytest.raises(ValueError,match='empty baseline'):
        module.integrate(c,d,r,o)


def test_duplicate_or_ontology_change_is_rejected():
    c,d,r,o=sample()
    with pytest.raises(ValueError,match='duplicate trait'):
        module.integrate(c,d,r+r,o)
    r[0]['axis']='floral_structural_complexity'
    with pytest.raises(ValueError,match='ontology'):
        module.integrate(c,d,r,o)


def test_two_traits_share_one_cell_without_losing_multistate():
    c,d,r,o=sample()
    r2=copy.deepcopy(r[0])
    r2.update(trait_name='self_incompatibility',normalized_value='SC|SI',state_set=['SC','SI'],quality='high',record_id='second')
    o['traits']['self_incompatibility']={'domain':'reproductive_assurance','allowed_values':['SC','SI']}
    after,ledger,changes=module.integrate(c,d,r+[r2],o)
    assert len(after)==1 and len(ledger)==2
    assert after.quality.iloc[0]=='high'
    assert 'self_incompatibility=["SC","SI"]' in after.trait_composition.iloc[0]
    assert len(changes[['accepted_species','axis']].drop_duplicates())==1


def test_reward_never_projects_to_structure():
    c,d,r,o=sample()
    r[0].update(trait_name='reward_type',axis='floral_structural_complexity',normalized_value='nectar')
    o['traits']['reward_type']={'domain':'floral_architecture','allowed_values':['nectar']}
    with pytest.raises(ValueError,match='ontology'):
        module.integrate(c,d,r,o)


def correction_sample():
    receipt=json.loads((root/'data/v2/staging/traits/restart_harrisia_correction_20260907.json').read_text())
    c=pd.DataFrame([dict(accepted_species='Harrisia portoricensis',axis='reproductive_assurance',quality='high',trait_composition='self_incompatibility=["SC"]')])
    d=pd.DataFrame([dict(accepted_species='Harrisia portoricensis',trait_name='self_incompatibility',quality='high',normalized_value='SC',source_lineages='url:https://europepmc.org/article/MED/21622342')])
    d=pd.concat([d,pd.DataFrame([dict(accepted_species='Harrisia portoricensis',trait_name='autonomous_selfing_capacity',quality='medium',normalized_value='absent',source_lineages='doi:10.1038/ncomms13313')])],ignore_index=True)
    o={'traits':{'self_incompatibility':{'allowed_values':['SC','mixed_or_variable']},'autonomous_selfing_capacity':{'allowed_values':['absent','autonomous']}}}
    return c,d,receipt,o


def test_partial_compatibility_correction_does_not_claim_new_axis():
    c,d,r,o=correction_sample()
    after,ledger=module.correct_harrisia(c,d,r,o)
    assert after.quality.equals(c.quality)
    assert ledger.set_index('trait_name').normalized_value.to_dict()=={'self_incompatibility':'mixed_or_variable','autonomous_selfing_capacity':'absent'}
    assert c.trait_composition.iloc[0]=='self_incompatibility=["SC"]'
    with pytest.raises(ValueError,match='precondition'):
        module.correct_harrisia(after,ledger,r,o)


def test_correction_rejects_changed_source_and_genus_training():
    c,d,r,o=correction_sample()
    d.loc[0,'source_lineages']='unrelated'
    with pytest.raises(ValueError,match='precondition'):
        module.correct_harrisia(c,d,r,o)
    r['genus_rule_training_allowed']=True
    with pytest.raises(ValueError,match='Unapproved'):
        module.correct_harrisia(c,d,r,o)
