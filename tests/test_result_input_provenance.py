import json
from pathlib import Path
import importlib.util
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("result_input_provenance",ROOT/"scripts/check_result_input_provenance.py")
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
REG=json.loads((ROOT/"research/result_input_provenance_20261125.json").read_text())
S=REG["sources"]

def test_current_registry():
    assert mod.validate(REG)==len(REG["results"])
def test_source_not_count():
    assert mod.classify(S,["own_gbif","own_traits"])==("island","own")
    assert mod.classify(S,["gift_island","own_traits"])==("island","mixed")
    assert mod.classify(S,["gift_island","gift_traits"])==("fd","external_only")
    assert mod.classify(S,["euppollnet","own_plant_aggregate"])==("island","mixed")
    assert mod.classify(S,["glopl","geography"])==("fd","external_only")
    assert mod.classify(S,["gift_island","unknown"])==("hold","unverified")
def test_fails_closed():
    with pytest.raises(ValueError): mod.classify(S,["unregistered"])
    with pytest.raises(ValueError): mod.classify({"a":{"type":"derived","parents":["a"]}},["a"])
