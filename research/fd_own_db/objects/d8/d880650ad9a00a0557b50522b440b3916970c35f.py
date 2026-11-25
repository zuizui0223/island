import pandas as pd
import pytest
from fd.source_signature_dependency import audit_source_overlap, _components


def synthetic(n=48):
    rows, src=[], []
    for i in range(n):
        sources=["common",f"v{i%8}_1",f"v{i%8}_2",f"v{i%8}_3",f"v{i%8}_4"]
        for rank,entity in enumerate(sources,1):
            src.append(dict(island_id=f"Island_{i}",source_mode="geo_k5",
                            source_rank=rank,entity_ID=entity))
        rows.append(dict(region="north",scope="All",trait="SC",
                         ledger="snapshot",island_id=f"Island_{i}",
                         spatial_block=f"block_{i//3}",seed=20261110,
                         brier_advantage=.001))
    return pd.DataFrame(rows),pd.DataFrame(src)


def test_one_shared_donor_links_distinct_source_signatures():
    a,s=synthetic()
    z=audit_source_overlap(a,s,bootstrap_draws=299)
    assert len(z)==1
    q=z.iloc[0]
    assert q.n_islands==48
    assert q.n_spatial_blocks==16
    assert q.n_source_signatures==8
    assert q.components_any_donor_overlap==1
    assert q.components_three_of_five_overlap==8
    assert q.top_two_source_signature_islands==12
    assert q.largest_donor_component_islands==48
    assert q.block_equal_gain_drop_top_two==pytest.approx(.001)
    assert q.component_equal_gain==pytest.approx(.001)
    assert not {"island_id","source_signature","entity_ID","spatial_block"} & set(z.columns)


def test_repeated_seed_is_not_second_island_replication():
    a,s=synthetic()
    b=a.copy();b.seed=20261111
    z=audit_source_overlap(pd.concat([a,b],ignore_index=True),s,bootstrap_draws=299)
    assert z.iloc[0].n_islands==48


def test_duplicate_and_missing_donors_raise():
    a,s=synthetic()
    with pytest.raises(ValueError,match="Duplicate"):
        audit_source_overlap(pd.concat([a,a.iloc[[0]]],ignore_index=True),s,bootstrap_draws=299)
    with pytest.raises(ValueError,match="five donors"):
        audit_source_overlap(a,s.iloc[1:],bootstrap_draws=299)
    t=s.copy();t.loc[0,"entity_ID"]=t.loc[1,"entity_ID"]
    with pytest.raises(ValueError,match="Duplicate source donor"):
        audit_source_overlap(a,t,bootstrap_draws=299)


def test_suppress_when_only_one_spatial_block():
    a,s=synthetic()
    a.spatial_block="same_block"
    z=audit_source_overlap(a,s,bootstrap_draws=299)
    assert z.empty
    assert "components_any_donor_overlap" in z.columns


def test_independent_components_differ_from_signatures():
    assert _components([
        ("a","b","c","d","e"),
        ("a","f","g","h","i"),
        ("j","k","l","m","n"),
    ],1)==(2,2)
