from itertools import combinations, permutations
import numpy as np
import pandas as pd
from scripts.fd_self_db_isolation import expected_q, excess_mi, morphology, fit_cluster, adjust_fixed, DIST, COV

def test_exact_q_matches_exhaustive_draws():
    a = np.array(["A","A","A","B","C"])
    values = []
    for subset in combinations(a,3):
        _, counts = np.unique(subset,return_counts=True)
        values.append(1-np.square(counts/3).sum())
    assert np.isclose(expected_q([3,1,1],3),np.mean(values))

def test_mi_null_exact():
    a=np.array([["a","x","r"],["a","x","r"],["b","y","z"],["b","y","z"]])
    v=[]
    for labels in permutations(a[:,0]):
        b=a.copy()
        b[:,0]=labels
        v.append(excess_mi(b))
    assert abs(np.mean(v))<1e-12
    assert excess_mi(a)>0

def test_joint_not_marginal():
    a=np.array([["a","x"],["a","x"],["b","y"],["b","y"]])
    b=np.array([["a","x"],["a","y"],["b","x"],["b","y"]])
    x=morphology(a,4,1); y=morphology(b,4,1)
    assert np.isclose(x["marginal_q"],y["marginal_q"])
    assert x["joint_q"]<y["joint_q"]

def test_cluster_and_support():
    rng=np.random.default_rng(4)
    cols=[DIST,*COV]
    p=pd.DataFrame(rng.normal(size=(160,len(cols))),columns=cols)
    p["spatial_block"]=np.repeat(np.arange(16),10)
    p["y"]=2*p[DIST]+rng.normal(0,.1,160)
    x=fit_cluster(p,"y",cols)
    assert x["status"]=="fit" and x["lo"]>0
    assert fit_cluster(p.iloc[:40],"y",cols)["status"]=="not_testable"

def test_fixed_slots():
    rows=[dict(scope="a",variant="b",family="c",status="fit",p=.01),dict(scope="a",variant="b",family="c",status="not_testable")]
    adjust_fixed(rows)
    assert np.isclose(rows[0]["q_bh"],.02)
    assert np.isclose(rows[0]["q_by"],.03)
    assert rows[1]["family_slots"]==2 and "q_bh" not in rows[1]
