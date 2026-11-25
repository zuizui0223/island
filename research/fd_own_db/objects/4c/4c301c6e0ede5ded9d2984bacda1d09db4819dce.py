"""Pure synthetic tests. Original licensed taxa and island IDs never enter CI."""
import itertools

import numpy as np
import pandas as pd
import pytest

from fd.pollinator_functional_diversity import (
    exact_colour_architecture_excess, sample_morphology, direct_triple,
)


def test_exact_shuffle_correction_averages_zero_over_all_colour_permutations():
    # Exact finite-sample null is independently confirmed by exhaustive labels.
    form = np.array(["open_radial", "open_radial", "closed", "closed"],dtype=str)
    sym = np.array(["actinomorphic", "actinomorphic", "zygomorphic", "zygomorphic"],dtype=str)
    perms = sorted(set(itertools.permutations(("red","red","blue","blue"))))
    scores = [exact_colour_architecture_excess(np.column_stack((p,form,sym)))
              for p in perms]
    assert len(scores)==6
    assert abs(np.mean(scores))<1e-12
    assert max(scores)>0
    assert min(scores)<0


def test_same_marginal_rao_can_hide_different_joint_functional_diversity():
    n = 30
    forms = np.array(["open_radial"]*15+["closed"]*15,dtype=str)
    symmetry = np.array(["actinomorphic"]*30,dtype=str)
    aligned = np.column_stack((["red"]*15+["blue"]*15,forms,symmetry))
    mixed = np.column_stack((
        ["red"]*8+["blue"]*7+["red"]*7+["blue"]*8,
        forms,symmetry))
    rng1=np.random.default_rng(1)
    rng2=np.random.default_rng(1)
    a=sample_morphology(aligned,repeats=5,rng=rng1)
    b=sample_morphology(mixed,repeats=5,rng=rng2)
    assert a["marginal_rao"] == pytest.approx(b["marginal_rao"])
    assert a["joint_rao"] != pytest.approx(b["joint_rao"])
    assert a["joint_hill2"] != pytest.approx(b["joint_hill2"])
    assert a["colour_architecture_excess_mi"] > b["colour_architecture_excess_mi"]


def test_direct_guard_rejects_genus_imputation_and_conflicts():
    inp=pd.DataFrame([
        ("A","flower_primary_color","red","species_direct","high"),
        ("A","floral_form","open_radial","species_direct","high"),
        ("A","floral_symmetry","actinomorphic","species_direct","high"),
        ("B","flower_primary_color","red","genus_consensus","medium"),
        ("B","floral_form","open_radial","genus_consensus","medium"),
        ("B","floral_symmetry","actinomorphic","genus_consensus","medium"),
        ("C","flower_primary_color","red|pink","species_direct","high"),
        ("C","floral_form","open_radial","species_direct","high"),
        ("C","floral_symmetry","actinomorphic","species_direct","high"),
    ],columns=("accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"))
    x=direct_triple(inp)
    assert list(x.index)==["A"]
    with pytest.raises(ValueError):
        direct_triple(pd.concat([inp,inp.iloc[[0]]],ignore_index=True))


def test_exact_mi_requires_real_profiles():
    with pytest.raises(ValueError):
        exact_colour_architecture_excess(np.array([["red","open_radial","actinomorphic"]]))
