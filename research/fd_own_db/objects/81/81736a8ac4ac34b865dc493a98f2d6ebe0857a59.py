import numpy as np
import pandas as pd
import pytest
from fd.globi_family_guild_audit import audit_guilds

def _fixture():
    rng=np.random.default_rng(234)
    n=80
    idx=[f"G{i}" for i in range(n)]
    family=pd.Series([f"F{i//4}" for i in range(n)],index=idx)
    z=rng.uniform(0,1,n)
    X=pd.DataFrame({"ref_effort":rng.uniform(1,3,n),
                    "ref_plant":rng.uniform(1,2,n),
                    "form_open":z},index=idx)
    y=np.column_stack([.2+.18*z,.3-.12*z,.15+.02*z,.15-.04*z,.2-.04*z])
    y=y/y.sum(axis=1)[:,None]
    return X,family,y

def test_anonymous_family_equal_and_genus_equal_outputs():
    x,fam,y=_fixture()
    out=audit_guilds(x,fam,y,source_fields=("ref_effort","ref_plant"),
                     bootstrap=299,shuffles=99)
    assert len(out)==6
    assert set(out.n_families)=={20}
    assert out.loc[out.channel.eq("five_channel_mean"),"genus_equal_MSE_gain"].iloc[0]>0
    assert out.n_morphology_permutable_genera.eq(80).all()
    assert not {"genus","family","accepted_species","island_id","provider_id"} & set(out)
    assert not out.local_pollination_measured.any()
    assert not out.provider_independent.any()

def test_unresolved_taxonomy_and_missing_morphology_rejected():
    x,fam,y=_fixture()
    fam.iloc[0]="AMBIG:FamilyA&FamilyB"
    with pytest.raises(ValueError,match="Unresolved"):
        audit_guilds(x,fam,y,source_fields=("ref_effort","ref_plant"),
                     bootstrap=299,shuffles=99)
    x,fam,y=_fixture()
    x.iloc[0,-1]=np.nan
    with pytest.raises(ValueError,match="Unknown"):
        audit_guilds(x,fam,y,source_fields=("ref_effort","ref_plant"),
                     bootstrap=299,shuffles=99)

def test_compositional_invalid_and_genus_duplicate_fail_closed():
    x,fam,y=_fixture()
    y[0,0]=.95
    with pytest.raises(ValueError,match="Invalid"):
        audit_guilds(x,fam,y,source_fields=("ref_effort","ref_plant"),
                     bootstrap=299,shuffles=99)
    x,fam,y=_fixture()
    x.index=["G0"]*len(x)
    with pytest.raises(ValueError,match="unique genus"):
        audit_guilds(x,fam,y,source_fields=("ref_effort","ref_plant"),
                     bootstrap=299,shuffles=99)
