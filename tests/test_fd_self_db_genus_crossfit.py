import numpy as np
from scripts.fd_self_db_genus_crossfit import probabilities,fold,simulate,proper_scores

def test_fold_outcome_changes_do_not_change_own_fold_predictions():
    names=[f"Genus species{i}" for i in range(35)]
    codes=np.arange(len(names))%3
    genus=["Genus"]*len(names)
    x,n,f=probabilities(names,codes,genus,3,3)
    change=f==f[0]
    changed=codes.copy();changed[change]=(changed[change]+1)%3
    y,m,g=probabilities(names,changed,genus,3,3)
    assert np.allclose(x[:,change],y[:,change])
    assert np.array_equal(f,g)

def test_no_congener_is_exact_blind_fallback():
    names=[f"Species {i}" for i in range(25)]
    p,n,f=probabilities(names,np.arange(25)%3,names,3,4)
    assert np.all(n==0)
    assert np.allclose(p[0],p[1]) and np.allclose(p[0],p[2])

def test_coherent_sampling_and_proper_score():
    p=np.array([[1.,0.],[0.,1.]])
    u=np.random.default_rng(2).random((100,2))
    a=simulate(p,u)
    assert np.array_equal(a,np.tile([0,1],(100,1)))
    assert np.allclose(proper_scores(p,[0,1]),0)
