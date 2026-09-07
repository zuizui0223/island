import importlib.util
import json
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

spec = importlib.util.spec_from_file_location('recovery', Path(__file__).parents[1] / 'scripts/recover_verified_trait_checkpoint.py')
recovery = importlib.util.module_from_spec(spec)
spec.loader.exec_module(recovery)


def inputs(tmp_path, names, states, ontology):
    coverage = pd.DataFrame([dict(accepted_species='Example species', axis='reproductive_assurance', quality='low', trait_composition='', trait_names='', source_groups='', source_lineages='')])
    sidecar = tmp_path / 'low.csv.gz'
    pd.DataFrame([dict(accepted_species='Example species', axis='reproductive_assurance', trait_names=json.dumps(names), predicted_state_sets=json.dumps([json.dumps(s) for s in states]), support_source_lineages='["doi:source1","doi:source2"]', family_inference='false', global_fallback='false')]).to_csv(sidecar,index=False)
    allowed=tmp_path/'ontology.json'
    allowed.write_text(json.dumps({'traits': {k: {'allowed_values': v} for k,v in ontology.items()}}))
    return coverage,sidecar,allowed


def test_independently_sorted_lists_do_not_swap_reproductive_traits(tmp_path):
    args=inputs(tmp_path,['autonomous_selfing_capacity','self_incompatibility'],[['SC'],['autonomous']],{'self_incompatibility':['SC','SI'],'autonomous_selfing_capacity':['autonomous','absent']})
    with patch.object(recovery,'LOW_SIDECAR_SHA256',recovery.sha256(args[1])):
        fixed,audit=recovery.restore_low_values(*args)
    assert fixed.iloc[0].trait_composition == 'autonomous_selfing_capacity=["autonomous"]|self_incompatibility=["SC"]'
    assert fixed.iloc[0].quality == 'low'
    assert len(audit)==1


def test_ambiguous_state_association_fails_closed(tmp_path):
    args=inputs(tmp_path,['trait_a','trait_b'],[['x'],['y']],{'trait_a':['x','y'],'trait_b':['x','y']})
    with patch.object(recovery,'LOW_SIDECAR_SHA256',recovery.sha256(args[1])):
        with pytest.raises(ValueError,match='ambiguous'):
            recovery.restore_low_values(*args)


def test_existing_checkpoint_is_not_overwritten(tmp_path):
    with pytest.raises(ValueError,match='never overwritten'):
        recovery.recover(tmp_path/'missing-source',tmp_path)


def test_deduplicated_state_can_belong_to_two_traits(tmp_path):
    args=inputs(tmp_path,['autonomous_selfing_capacity','cleistogamy','self_incompatibility'],[['SC','SI'],['absent']],{'self_incompatibility':['SC','SI'],'autonomous_selfing_capacity':['autonomous','absent'],'cleistogamy':['present','absent']})
    with patch.object(recovery,'LOW_SIDECAR_SHA256',recovery.sha256(args[1])):
        fixed,_=recovery.restore_low_values(*args)
    assert fixed.iloc[0].trait_composition == 'autonomous_selfing_capacity=["absent"]|cleistogamy=["absent"]|self_incompatibility=["SC","SI"]'
