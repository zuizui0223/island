"""Tiny real-branch-length rooted tree, no false topology or imputed tips."""
import numpy as np
import pytest

from fd.true_backbone_phylo_distance import true_tip_patristic_distance


def tree():
    edge=np.array([[7,5],[7,6],[5,1],[5,2],[6,3],[6,4]])
    lengths=np.array([2.,3.,1.,1.5,2.,1.])
    return edge,lengths


def test_exact_backbone_lca_patristic_distances():
    edge,weights=tree()
    tips=np.array([1,2,3,4])
    D=true_tip_patristic_distance(edge,weights,tips,n_original_tree_tips=4)
    assert D.shape==(4,4)
    assert D[0,1]==pytest.approx(2.5)
    assert D[2,3]==pytest.approx(3.)
    assert D[0,2]==pytest.approx(8.)
    assert D[1,3]==pytest.approx(7.5)
    assert np.array_equal(D,D.T)
    assert np.all(np.diag(D)==0)


def test_missing_taxa_and_unsupported_tree_edges_fail_closed():
    edge,w=tree()
    with pytest.raises(ValueError,match="not a measured"):
        true_tip_patristic_distance(edge,w,np.array([1,5]),n_original_tree_tips=4)
    with pytest.raises(ValueError,match="branch length"):
        true_tip_patristic_distance(edge,np.array([-1.,2.,1.,1.,1.,1.]),
                                    np.array([1,2]),n_original_tree_tips=4)
    with pytest.raises(ValueError):
        true_tip_patristic_distance(edge,w,np.array([1,1]),n_original_tree_tips=4)
