import pandas as pd
import pytest
from fd.source_network_universe import audit_full_cohort


def data(n=45):
    islands=pd.DataFrame({
        "island_id":[f"I{i}" for i in range(n)],
        "analysis_regime":["north"]*n,
    })
    src=[]
    for i in range(n):
        donor=["shared", f"D{i%9}a", f"D{i%9}b", f"D{i%9}c", f"D{i%9}d"]
        for rank,e in enumerate(donor,1):
            src.append({"island_id":f"I{i}","source_mode":"geo_k5",
                        "source_rank":rank,"entity_ID":e})
    return islands,pd.DataFrame(src)


def test_all_islands_audited_before_trait_support_filters():
    g,s=data()
    out=audit_full_cohort(g,s)
    assert out.shape[0]==1
    x=out.iloc[0]
    assert x.n_total_islands==45
    assert x.n_source_signatures==9
    assert x.components_any_shared_source==1
    assert x.largest_component_islands_any_shared==45
    assert x.components_three_shared_sources==9
    assert not {"island_id","entity_ID","source_id","source_signature"} & set(out.columns)


def test_missing_label_and_duplicate_focal_rejected():
    g,s=data()
    with pytest.raises(ValueError,match="Unlabeled"):
        audit_full_cohort(g,s[s.island_id!="I0"])
    with pytest.raises(ValueError,match="exactly one regime"):
        audit_full_cohort(pd.concat([g,g.iloc[[0]]]),s)


def test_duplicate_source_entity_rejected():
    g,s=data()
    s.loc[0,"entity_ID"]=s.loc[1,"entity_ID"]
    with pytest.raises(ValueError,match="Duplicate"):
        audit_full_cohort(g,s)


def test_small_group_suppressed():
    g,s=data()
    assert audit_full_cohort(g.iloc[:20],s[s.island_id.isin(g.iloc[:20].island_id)]).empty


def test_privacy_threshold_enforced():
    g,s=data()
    with pytest.raises(ValueError,match="Privacy"):
        audit_full_cohort(g,s,min_islands=2)
