import importlib.util
import json
from pathlib import Path
import pandas as pd
import pytest

root = Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location('integration',root/'scripts/integrate_reviewed_restart.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def sample():
    records=json.loads((root/'data/v2/staging/traits/restart_reviewed_20260907.json').read_text())
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
