"""Synthetic tests for marginal-preserving floral colour/form coupling scores."""
import numpy as np
import pandas as pd
import pytest

from fd.marginal_fixed_coupling import (
    margins, match_marginals, marginal_fixed_source_gain,
    marginal_fixed_placebo_gain,
    profile_axis_codes,
)
from fd.trait_kernels import phenotype_kernels


def profiles():
    return {
        0: (frozenset({"red"}), frozenset({"tube"})),
        1: (frozenset({"red"}), frozenset({"open"})),
        2: (frozenset({"blue"}), frozenset({"tube"})),
        3: (frozenset({"blue"}), frozenset({"open"})),
    }


def test_raking_exactly_fixes_both_marginals_without_losing_coupling():
    c,f=profile_axis_codes(profiles())
    g=np.array([.25,.25,.25,.25])
    local=np.array([.49,.01,.01,.49])
    matched, info=match_marginals(g,local,c,f)
    assert np.allclose(matched,local)
    assert np.allclose(margins(g,c,f)[0],margins(matched,c,f)[0],atol=1e-8)
    assert np.allclose(margins(g,c,f)[1],margins(matched,c,f)[1],atol=1e-8)
    assert info["max_abs_marginal_error"] < 1e-8


def test_pure_marginal_shift_is_removed_even_when_probabilities_change():
    c,f=profile_axis_codes(profiles())
    g=np.array([.2,.3,.1,.4])
    # Multiply by row and column effects only: odds ratios are invariant,
    # so IPF must return the original joint table.
    local=g*np.array([4.,4.,1.,1.])*np.array([1.,3.,1.,3.])
    local/=local.sum()
    assert not np.allclose(g,local)
    matched, info=match_marginals(g,local,c,f)
    assert np.allclose(matched,g,atol=1e-8)
    assert info["max_abs_marginal_error"] < 1e-8


class Toy:
    n_profiles=4
    def probabilities(self,source,genus):
        assert genus=="Gen"
        prior=np.ones(4)/4
        if source=="real":
            local=np.array([.49,.01,.01,.49])
        elif source=="fake1":
            local=np.array([.01,.49,.49,.01])
        elif source=="fake2":
            local=np.array([.2,.3,.3,.2])
        else:
            raise ValueError("Unknown source")
        return np.stack([prior,prior,prior,local])


def test_coupling_gain_with_identical_margins_outcome_blind_sources():
    data=pd.DataFrame({
        "accepted_species":["Gen a","Gen b"],
        "profile_id":[0,3], "genus":["Gen","Gen"],
    })
    actual=marginal_fixed_source_gain(data,Toy(),"real",profiles())
    wrong=marginal_fixed_source_gain(data,Toy(),"fake1",profiles())
    assert actual["gain_marginal_fixed_joint"] > .2
    assert wrong["gain_marginal_fixed_joint"] < 0
    compare=marginal_fixed_placebo_gain(
        data,Toy(),"real",["fake1","fake2"],profiles())
    assert compare["actual_minus_placebo_marginal_fixed_joint"]>0
    assert compare["maximum_ipf_error"]<1e-8
    assert set(compare)=={
        "actual_minus_placebo_marginal_fixed_joint",
        "actual_minus_placebo_marginal_fixed_brier",
        "maximum_ipf_error","mean_ipf_iterations"}
    with pytest.raises(ValueError):
        marginal_fixed_placebo_gain(data,Toy(),"real",["fake1","real"],profiles())


def test_no_false_coupling_gain_when_local_changes_only_marginals():
    c,f=profile_axis_codes(profiles())
    g=np.array([.2,.3,.1,.4])
    shifted=g*np.array([2,2,1,1])
    shifted/=shifted.sum()
    got,_=match_marginals(g,shifted,c,f)
    assert np.allclose(got,g,atol=1e-8)


def test_invalid_probability_and_no_renormalization_pseudocounts():
    c,f=profile_axis_codes(profiles())
    g=np.ones(4)/4
    with pytest.raises(ValueError):
        match_marginals(g,np.array([.4,.6,0.,0.]),c,f)
    with pytest.raises(ValueError):
        match_marginals(g,np.array([.1,.1,.1,.1]),c,f)
    with pytest.raises(ValueError):
        match_marginals(g,g,c,f,max_iter=0)
