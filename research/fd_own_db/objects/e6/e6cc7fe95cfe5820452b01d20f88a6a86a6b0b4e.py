import importlib.util
import pathlib
import sys
import pandas as pd
import pytest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"))
from range_map_partner_incidence import validate_config, traits_for_mode

def config():
    return dict(provider="synthetic",release="test",codebook_url="https://example.invalid/codebook",authorized_use_note="synthetic only",
        fields=dict(scientific="sci_name",presence="PRESENCE",origin="ORIGIN",seasonal="SEASONAL"),
        keep=dict(presence=[1,2],origin=[1,2],seasonal=[1,2]),diet_mode="fruit",
        hashes={k:"a"*64 for k in ("ranges","islands","traits")})

def test_provider_schema():
    assert validate_config(config())["diet_mode"]=="fruit"
    x=config();x["fields"]["seasonal"]="bad OR 1=1"
    with pytest.raises(ValueError):validate_config(x)
    x=config();x["keep"].pop("seasonal")
    with pytest.raises(ValueError):validate_config(x)

def test_trait_modes(tmp_path):
    f=tmp_path/"traits.csv"
    pd.DataFrame(dict(scientific=["Abc def","Abc ghi","Jkl mno"],family=["F"]*3,clade=["A","A","B"],diet_fraction=[0,10,90])).to_csv(f,index=False)
    assert traits_for_mode(f,"fruit").diet_nect.tolist()==[0,100,100]
    assert traits_for_mode(f,"nectar").diet_nect.tolist()==[0,10,90]
