import numpy as np
from scripts.fd_self_db_capacity_space import calc_metrics, coverage_bounds, p_no, normalize_roster
import pandas as pd

def test_all_or_none_capacity():
    near=np.array([[1]*5+[0]*5,[0]*5+[1]*5],bool)
    for v in (0,1):
        d=calc_metrics(near,np.full(10,v),["a"]*10)
        assert np.isclose(d["share_n10"],v)
        assert abs(d["excess_random_n10"])<1e-10
        assert abs(d["excess_genus_full"])<1e-10

def test_rare_profiles_coverage_and_genus():
    near=np.array([[1]*5+[0]*5,[0]*5+[1]*5],bool)
    labels=np.array([1]*5+[0]*5)
    x=calc_metrics(near,labels,["a"]*5+["b"]*5)
    assert np.isclose(x["share_n10"],.5)
    assert x["excess_random_n10"]<0
    assert abs(x["excess_genus_full"])<1e-10

def test_unknown_bounds():
    near=np.array([[1,0,0],[0,1,0],[0,0,1]],bool)
    lo,hi=coverage_bounds(near,[1,0,np.nan])
    assert np.allclose([lo,hi],[1/3,2/3])

def test_no_hit():
    assert p_no(10,5,10)==0 and p_no(10,0,10)==1
    assert np.isclose(p_no(10,1,5),.5)

def test_roster_ignores_status():
    r=pd.DataFrame({"island_id":["a","a","b"],"species":["X y","X y","Y z"],"origin_status":["native","introduced","unknown"]})
    x=normalize_roster(r)
    assert len(x)==2 and list(x.columns)==["island_id","accepted_species"]
