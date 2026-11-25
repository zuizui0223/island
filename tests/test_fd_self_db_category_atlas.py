import numpy as np
from scripts.fd_self_db_category_atlas import no_hit,partition,exact_genus_count

def test_exact_volume_partition():
    near=np.array([[1]*5+[0]*5,[0]*5+[1]*5,[1]*10],bool)
    d=partition(10,near,np.array([1]*5+[0]*5,bool),np.array([0]*5+[1]*5,bool),
                np.ones(3)/3,np.array([0,2,1]))
    assert np.allclose([d["positive_only_mass"],d["negative_only_mass"],d["shared_mass"]],[1/3]*3)
    assert np.isclose(d["total_mass"],sum(d[k] for k in ("positive_only_mass","negative_only_mass","shared_mass")))
    assert np.isclose(d["positive_share"],2/3)

def test_no_hit_finite():
    assert np.isclose(no_hit(10,1,5),.5)
    assert np.allclose(no_hit(10,np.array([0,5,10]),10),[1,0,0])

def test_exact_genus_distribution():
    y=np.array([1,1,0,0])
    mask=np.array([1,1,0,0],bool)
    e,p,m=exact_genus_count(y,mask,np.array(["a"]*4))
    assert np.isclose(e,1) and np.isclose(p,1/3) and m==2
    e,p,m=exact_genus_count(y,mask,np.array(["a","a","b","b"]))
    assert e==2 and p==1 and m==0
