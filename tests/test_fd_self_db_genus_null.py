import numpy as np
from scripts.fd_self_db_genus_null import random_assignments, all_metrics, parent_draws, perm_p
from scripts.fd_self_db_isolation import morphology, seed_for, METRICS

def test_whole_profiles_and_genus_margins():
    codes=np.array([0,1,1,2,3,3])
    g=["A","A","A","B","B","C"]
    a,v=random_assignments(codes,g,30,2)
    assert np.array_equal(a[0],codes)
    assert np.array_equal(a[:,5],np.full(31,3))
    for row in a:
        assert sorted(row[:3])==[0,1,1] and sorted(row[3:5])==[2,3]
    assert v.tolist()==[True,True,True,True,True,False]

def test_vectorized_metrics_reproduce_parent():
    profiles=np.array([["a","x","r"],["b","x","r"],["a","y","s"],["b","y","s"]])
    rng=np.random.default_rng(52)
    p=rng.integers(0,4,size=(3,50))
    draws=parent_draws(50,"synthetic","all_recorded")
    m=all_metrics(p,profiles,draws)
    for i in range(3):
        expected=morphology(profiles[p[i]],30,seed_for("synthetic","triple_cfs_n30","all_recorded"))
        assert np.allclose(m[i],[expected[k] for k in METRICS],atol=1e-12)

def test_exact_30_and_constant_null():
    profiles=np.array([["a","x","r"],["b","y","s"]])
    p=np.tile(np.arange(30)%2,(3,1))
    a=all_metrics(p,profiles,parent_draws(30,"i","all_recorded"))
    assert np.allclose(a,a[0])
    assert perm_p(0,np.zeros(399))==1

def test_permutation_tail_plus_one():
    assert perm_p(10,np.zeros(399))==.005
