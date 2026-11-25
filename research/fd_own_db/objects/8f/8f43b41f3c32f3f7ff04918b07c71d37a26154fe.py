import numpy as np
import pandas as pd
from scipy import sparse
import pytest
from fd.partner_flora_e2 import run_partner_test, run_partner_clades, holm_fixed

def synthetic():
    rng=np.random.default_rng(19)
    rows=[];plants=[];inc=[]
    # Same regional genus pool: repeated species are not counted twice in expectation.
    for i in range(100):
        rows.append(dict(island=i,group=i//5,area_km2=2+rng.uniform(0,20),dist_km=rng.uniform(1,100),lat=rng.uniform(25,60)))
        v=rng.integers(0,2,14);v[0]=1
        inc.append(v)
        for g in range(8):
            k=int(rng.random()<(.65 if v[13] else .35))
            plants.append(dict(island=i,species=f"P{g} s{k}",genus=f"P{g}",grade="high",trait=k))
    return pd.DataFrame(rows),pd.DataFrame(plants),sparse.csr_matrix(inc),np.array([0]*13+[2])

def test_gate_and_bombus():
    d,p,m,l=synthetic()
    r=run_partner_test(d,p,m,l,mode="bombus",n_boot=19,n_perm=19)
    assert r["status"]=="insufficient_support" and r["p_primary"] is None and "E2" not in r
    r=run_partner_test(d,p,m,l,mode="bombus",min_groups=4,min_islands=20,n_boot=39,n_perm=39)
    # n is constant by construction, so frozen design must fail rank rather than silently delete control.
    assert r["status"]=="insufficient_rank"
    p=p[~((p.genus=="P7")&(p.island%3==0))].copy()
    r=run_partner_test(d,p,m,l,mode="bombus",min_groups=4,min_islands=20,n_boot=49,n_perm=49)
    assert r["status"]=="ok" and r["E2"]>0 and "animal_genus" not in r["permutation_p"]
    assert abs(r["E1"]-r["E2"]-r["E_expected"])<1e-10

def test_diet_and_fixed_holm():
    d,p,m,l=synthetic()
    p=p[~((p.genus=="P7")&(p.island%3==0))]
    l[2:7]=1;l[7:12]=2
    g=np.repeat(np.arange(7),2)
    r=run_partner_test(d,p,m,l,g,mode="diet",min_groups=4,min_islands=20,n_boot=39,n_perm=39)
    assert r["status"] in {"ok","insufficient_null_support"}
    if r["status"]=="ok":assert set(r["permutation_p"])=={"plant","animal_genus","joint"}
    r=run_partner_test(d,p,m,l,np.arange(14),mode="diet",min_groups=4,min_islands=20,n_boot=19,n_perm=19)
    assert r["status"]=="insufficient_animal_null_support" and r["p_primary"] is None
    assert holm_fixed({"BOMB_NORTH":.02})=={"NECTAR_BIRD":1.,"BOMB_NORTH":.06,"NECTAR_BAT":1.}
    with pytest.raises(ValueError):holm_fixed({"posthoc":.001})

def test_clade_gate_and_fruit_binary():
    d,p,m,l=synthetic()
    p=p[~((p.genus=="P7")&(p.island%3==0))]
    l[2:7]=1;l[7:12]=2
    g=np.repeat(np.arange(7),2)
    cl=np.repeat(["A","A","A","A","B","B","B"],2)
    z=run_partner_clades(d,p,m,l,g,cl,["A","B","unobserved"],n_boot=19,n_perm=19)
    assert len(z)==3 and all(v["status"]=="insufficient_support" and v["holm_p_clade_family"]==1 for v in z.values())
    with pytest.raises(ValueError):run_partner_clades(d,p,m,l,g,np.arange(14),["A"])
    # When only the arbitrary 1/2 fruit bins vary within genus, the binary fruit-use
    # permutation must not claim a supported movable animal null.
    binary_blocks=np.where(l==0,0,1)
    r=run_partner_test(d,p,m,l,binary_blocks,mode="fruit",min_groups=4,min_islands=20,n_boot=19,n_perm=19)
    assert r["status"] in {"insufficient_animal_null_support","insufficient_null_support"}
