"""Synthetic full same-species floral colour×form coupling invariance tests."""
import numpy as np
import pytest

from fd.paired_nominal_trait_diversity import nominal_two_axis_diversity


STATES=[
    ("colour0","form0"),
    ("colour0","form1"),
    ("colour1","form0"),
    ("colour1","form1"),
]


def test_additive_gower_cannot_detect_true_same_species_pairing():
    independent=nominal_two_axis_diversity([.25]*4,STATES)
    correlated=nominal_two_axis_diversity([.5,0,0,.5],STATES)
    anticorrelated=nominal_two_axis_diversity([0,.5,.5,0],STATES)
    for r in (independent,correlated,anticorrelated):
        assert r["nominal_first_axis_Rao_Q"] == pytest.approx(.5)
        assert r["nominal_second_axis_Rao_Q"] == pytest.approx(.5)
        assert r["nominal_Gower_additive_Rao_Q"] == pytest.approx(.5)
        assert r["joint_Q_if_source_axes_independent"] == pytest.approx(.75)
    assert independent["true_joint_category_Rao_Q"] == pytest.approx(.75)
    assert correlated["true_joint_category_Rao_Q"] == pytest.approx(.5)
    assert anticorrelated["true_joint_category_Rao_Q"] == pytest.approx(.5)
    assert independent["true_joint_Q_minus_independent_marginal_Q"] == pytest.approx(0)
    assert correlated["true_joint_Q_minus_independent_marginal_Q"] == pytest.approx(-.25)
    assert anticorrelated["true_joint_Q_minus_independent_marginal_Q"] == pytest.approx(-.25)
    assert correlated["additive_Gower_cannot_detect_cross_axis_association"]


def test_joint_diversity_excess_can_be_positive_at_fixed_marginals():
    source=nominal_two_axis_diversity([.6,.2,.2,0],STATES)
    expected_independent=nominal_two_axis_diversity([.64,.16,.16,.04],STATES)
    assert source["nominal_first_axis_Rao_Q"] == pytest.approx(.32)
    assert source["nominal_second_axis_Rao_Q"] == pytest.approx(.32)
    assert source["nominal_Gower_additive_Rao_Q"] == pytest.approx(.32)
    assert expected_independent["nominal_Gower_additive_Rao_Q"] == pytest.approx(.32)
    assert source["true_joint_Q_minus_independent_marginal_Q"] == pytest.approx(.0224)
    assert expected_independent["true_joint_Q_minus_independent_marginal_Q"] == pytest.approx(0)


def test_additive_category_gower_raoq_equals_weighted_marginal_simpson_exactly():
    rng=np.random.default_rng(20261120)
    for n in (4,9):
        pairs=[(f"a{i}",f"b{j}") for i in range(n) for j in range(n)]
        for _ in range(5):
            p=rng.dirichlet(np.ones(len(pairs)))
            observed=nominal_two_axis_diversity(p,pairs,weight_first=.2,weight_second=.8)
            D=np.array([
                [.2*(int(a!=aa))+.8*(int(b!=bb))
                 for aa,bb in pairs]
                for a,b in pairs],float)
            assert observed["nominal_Gower_additive_Rao_Q"] == pytest.approx(
                p@D@p,abs=1e-12)
            assert observed["nominal_Gower_additive_Rao_Q"] == pytest.approx(
                .2*observed["nominal_first_axis_Rao_Q"]+
                .8*observed["nominal_second_axis_Rao_Q"],abs=1e-12)


@pytest.mark.parametrize("p,categories,weight",[
    ([.4,.4],STATES,.5),
    ([.5,.5],[("c0","f0"),("c0","f0")],.5),
    ([.5,.5],["c0","f1"],.5),
    ([-.1,1.1],STATES[:2],.5),
    ([.3,.3],STATES[:2],.5),
    ([float("nan"),.5],STATES[:2],.5),
    ([.5,.5],STATES[:2],-1),
    ([.5,.5],STATES[:2],.8),
])
def test_pair_category_ontology_simplex_and_distance_weight_guards(p,categories,weight):
    with pytest.raises(ValueError):
        nominal_two_axis_diversity(p,categories,weight_first=weight)
