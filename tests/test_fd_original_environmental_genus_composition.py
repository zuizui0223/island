"""Source-private GIFT island genus-source environment predictions, synthetic only."""
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from fd_original_environmental_genus_composition import (
    _floats,_xyz,_jaccard,genus_source_matrix,
    _nearest_probabilities,_proper_genus_score,
    _group_summary,OUTPUT_COLS,
)


def test_predeclared_climate_area_isolation_features_and_longitude_wrap():
    df=pd.DataFrame([
        dict(climate_pc1=0,climate_pc2=1,climate_pc3=2,climate_pc4=3,
             area_km2=1,distance_to_continent_km=0),
        dict(climate_pc1=1,climate_pc2=0,climate_pc3=1,climate_pc4=2,
             area_km2=100,distance_to_continent_km=100),
        dict(climate_pc1=2,climate_pc2=3,climate_pc3=0,climate_pc4=1,
             area_km2=5000,distance_to_continent_km=1000),
    ])
    z=_floats(df)
    assert z.shape==(3,6) and np.isfinite(z).all()
    assert z[0,4]<z[2,4]
    pts=_xyz(np.array([0,0]),np.array([179.9,-179.9]))
    assert np.linalg.norm(pts[0]-pts[1])<.01
    df.loc[0,"area_km2"]=-1
    with pytest.raises(ValueError,match="negative"):
        _floats(df)


def test_missing_and_source_entity_jaccard_not_genus_roster():
    assert _jaccard({"a","b"},{"b","c"})==pytest.approx(1/3)
    assert _jaccard({"a"},{"z"})==0
    assert _jaccard(set(),set())==0


def test_full_4505_GIFT_genus_source_uses_only_positive_accepted_flora():
    islands=[f"Island_{k:04d}" for k in range(4505)]
    records=[]
    for i,name in enumerate(islands):
        for j in range(10):
            records.append(dict(
                island_id=name,
                accepted_species=f"{'Alpha' if j<7 else 'Beta'} species{i}_{j}",
                origin_status="native" if j<8 else "introduced"))
    flora=pd.DataFrame(records)
    matrix,nsp,ng,kn=genus_source_matrix(flora,islands,"all_recorded")
    assert matrix.shape==(4505,2)
    assert kn==2 and nsp.min()==nsp.max()==10
    assert ng.min()==ng.max()==2
    assert np.allclose(np.asarray(matrix.sum(axis=1)).ravel(),1)
    assert matrix[0,0]==pytest.approx(.7)
    q,n,genus_total,names=genus_source_matrix(flora,islands,"original_native")
    assert names==2 and n.min()==8
    assert q[0,0]==pytest.approx(7/8)
    bad=flora.iloc[:-1]
    with pytest.raises(ValueError,match="parent"):
        genus_source_matrix(bad,islands,"all_recorded")


def test_reference_and_environment_use_same_source_genus_support_and_unknown():
    G=csr_matrix([
        [.8,.2,0],
        [.9,.1,0],
        [.2,.8,0],
        [.1,.9,0],
        [0.,0.,1.],
    ])
    donor=np.array([0,1,2,3])
    uniform=np.asarray(G[donor].mean(axis=0)).ravel()
    scores=np.array([.01,.04,9.,12.])
    model=_nearest_probabilities(G,donor,scores,uniform,k=2)
    assert model[0]>uniform[0]
    assert np.isclose(model.sum(),1)
    truth=G.getrow(1)
    known=uniform>0
    base=_proper_genus_score(uniform,truth,known)
    near=_proper_genus_score(model,truth,known)
    assert near[0]<base[0] and near[1]<base[1]
    unknown=G.getrow(4)
    u=_proper_genus_score(model,unknown,known)
    assert np.isfinite(u).all()
    with pytest.raises(ValueError,match="probability"):
        _proper_genus_score(np.array([.3,.3,.3]),truth,known)


def test_negative_environmental_predictive_gain_retained_in_block_bootstrap():
    rows=[]
    for i in range(75):
        d=dict(spatial_block=f"b{i//3}",first_source=f"s{i//3}")
        for j,col in enumerate(OUTPUT_COLS):
            d[col]=-.025 if col.startswith("environment_source") else .015
        rows.append(d)
    frame=pd.DataFrame(rows)
    a=_group_summary(frame,"spatial_block",seed=44,nboot=99)
    assert a["status"]=="ORIGINAL_GIFT_TAXONOMIC_GENUS_SOURCE_COMPOSITION_GROUP_PREDICTION"
    assert a["n_source_groups"]==25
    assert a["environment_source_nearest_minus_region_equal_Brier_group_equal_mean"]==pytest.approx(-.025)
    assert a["environment_source_nearest_minus_region_equal_Brier_group_resample_95"]==pytest.approx([-.025,-.025])
    assert a["climate_nearest_minus_region_equal_logloss_group_equal_mean"]==pytest.approx(.015)
    b=_group_summary(frame.iloc[:12],"spatial_block",seed=44,nboot=99)
    assert b["status"]=="NOT_IDENTIFIABLE_SOURCE_GEOGRAPHIC_GROUP_SUPPORT"
