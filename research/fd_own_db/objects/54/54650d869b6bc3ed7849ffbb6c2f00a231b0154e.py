import hashlib
import pandas as pd
import pytest

from fd.functional_ledger_harmonization import (
    strict_states, ledger_overlap, paired_source_ledger_sensitivity,
)


def _ledger(names, values):
    return pd.DataFrame({
        "accepted_species": names,
        "trait_name": ["floral_symmetry"]*len(names),
        "normalized_value": values,
    })


def test_state_counts_and_ambiguous_record_exclusion():
    a=_ledger(["s1","s2","s3","s4"], ["actinomorphic","zygomorphic",
                                     "actinomorphic|zygomorphic","actinomorphic"])
    b=_ledger(["s1","s2","s4","s5"],["actinomorphic","actinomorphic",
                                    "actinomorphic","zygomorphic"])
    z=ledger_overlap(a,b).iloc[0]
    assert z.snapshot_strict == 3
    assert z.source_strict == 4
    assert z.both_strict == 3
    assert z.conflicting_states == 1
    assert z.only_snapshot == 0
    assert z.only_source == 1
    assert strict_states(a).shape[0] == 3


def _predictions(n=45):
    out=[]
    for i in range(n):
        for j in range(4):
            name=f"Species_{i}_{j}"
            fold=hashlib.sha256(name.encode()).digest()[0]%5
            for ledger, p_source in (("snapshot",.65),("source",.80)):
                out.append({
                    "scope": "All", "trait": "floral_symmetry", "region": "tropical",
                    "island_id": f"private_island_{i}", "accepted_species": name,
                    "spatial_block": f"private_block_{i//3}", "ledger": ledger,
                    "fold": fold, "state": 1, "p_global": .50, "p_source": p_source,
                })
    return pd.DataFrame(out)


def test_same_target_spatial_block_bootstrap_and_privacy():
    x=paired_source_ledger_sensitivity(_predictions(), bootstrap_draws=399)
    assert len(x)==1
    assert x.iloc[0].n_identical_target_islands==45
    assert x.iloc[0].n_spatial_blocks==15
    assert x.iloc[0].source_minus_snapshot_gain>0
    assert x.iloc[0].block_bootstrap_lower>0
    assert not {"island_id","accepted_species","spatial_block","genus"}&set(x.columns)


def test_changed_ledger_state_is_rejected():
    x=_predictions()
    x.loc[1,"state"]=0
    with pytest.raises(ValueError,match="species/state"):
        paired_source_ledger_sensitivity(x,bootstrap_draws=399)


def test_missing_name_in_one_ledger_is_rejected():
    x=_predictions()
    with pytest.raises(ValueError,match="species/state"):
        paired_source_ledger_sensitivity(x.drop(index=0),bootstrap_draws=399)


def test_duplicate_and_fold_leakage_are_rejected():
    x=_predictions()
    with pytest.raises(ValueError,match="Duplicate"):
        paired_source_ledger_sensitivity(
            pd.concat([x,x.iloc[[0]]],ignore_index=True),bootstrap_draws=399)
    x=_predictions()
    x.loc[0,"fold"]=(x.loc[0,"fold"]+1)%5
    with pytest.raises(ValueError,match="fold mismatch"):
        paired_source_ledger_sensitivity(x,bootstrap_draws=399)


def test_small_or_nonindependent_groups_are_suppressed():
    x=_predictions()
    x["spatial_block"]="same_block"
    assert paired_source_ledger_sensitivity(x,bootstrap_draws=399).empty
