"""Synthetic strict SC vs autonomous selfing evidence and provenance tests."""
import numpy as np
import pandas as pd
import pytest
from fd_reproductive_parallel_isolation import (
    direct_binary_states,expected_from_genus,
)


def test_autonomous_is_not_equated_with_sc_and_variable_is_missing():
    states={
       "self_incompatibility":{
           "A first":"SC","A second":"SI","B third":"mixed_or_variable"},
       "autonomous_selfing_capacity":{
           "A first":"absent","A second":"autonomous","B third":"delayed",
           "C fourth":"mixed_or_variable"}}
    out,audit=direct_binary_states(states)
    assert out["self_incompatibility"]=={"A first":1,"A second":0}
    assert out["autonomous_selfing_capacity"]=={
        "A first":0,"A second":1}
    assert len(audit)==2
    assert audit[1]["n_other_direct_status_excluded"]==2
    assert audit[0]["n_other_direct_status_excluded"]==1
    assert out["self_incompatibility"]["A first"]==1
    assert out["autonomous_selfing_capacity"]["A first"]==0


def test_in_region_genus_expected_share_not_evolution_or_exact_ancestors():
    status={"GenA first":1,"GenA second":0,
            "GenB first":1,"GenB second":1}
    original=[{"GenA first","GenB first"},{"GenA second","GenB second"}]
    focal=[{"GenA first","GenB first"},{"GenA second","GenB second"}]
    ex=expected_from_genus(focal,original,status)
    assert ex==pytest.approx([.75,.75])
    assert [1.-ex[0],.5-ex[1]]==pytest.approx([.25,-.25])


def test_invalid_reproductive_ontology_not_silently_imputed():
    with pytest.raises(ValueError):
        direct_binary_states({"autonomous_selfing_capacity":{}})
    with pytest.raises(ValueError):
        direct_binary_states({
          "autonomous_selfing_capacity":{"A first":"absent"},
          "self_incompatibility":{"A first":"SC"}})
