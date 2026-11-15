"""Original directional sector isolation regression tests, synthetic only."""
import numpy as np
import pandas as pd
import pytest

from fd_directional_flower_isolation_gradients import (
    sector,holm,weighted_cluster_fit,expected_same_region_genus_sector,
    build_island_rows, MAIN, EXTRA
)


def make_regression(n=120):
    rows=[]
    for i in range(n):
        x=np.log2(1.+i*9.0)
        block="block"+str(i//3)
        rows.append(dict(
            response=.18-.023*x+.011*np.log1p(10+i)+.07*(.08+.0002*i)
                     +.003*np.log(5+i%12)+.001*(30+i%25)+.002*np.sin(i*1.3),
            log2_isolation_plus1=x,log1p_area=np.log1p(10+i),
            direct_coverage=.08+.0002*i, log_direct_n=np.log(5+i%12),
            absolute_latitude=30+i%25,
            spatial_block=block,
            source_cluster="source"+str(i//6)))
    return pd.DataFrame(rows)


def test_exact_cluster_regression_direction_and_spatial_independence():
    x=make_regression()
    r=weighted_cluster_fit(x,cluster="spatial_block",n_boot=99,seed=7)
    assert r["status"]=="ORIGINAL_ISLAND_BLOCK_EQUAL_DESCRIPTIVE_REGRESSION"
    assert r["n_clusters"]==40 and r["n_islands"]==120
    assert r["estimate_per_doubling"]==pytest.approx(-.023,abs=.002)
    assert r["n_valid_block_resamples"]>=80
    assert r["delete_one_block_min"]==pytest.approx(-.023,abs=.002)
    assert r["block_resample_quantile_95"][1]<0
    y=weighted_cluster_fit(x,cluster="source_cluster",n_boot=99,seed=7)
    assert y["n_clusters"]==20
    assert y["estimate_per_doubling"]==pytest.approx(-.023,abs=.002)


def test_holm_adjusts_all_six_named_outcomes_not_only_positive():
    raw=[.001,.02,.03,.05,.1,.8]
    adj=holm(raw)
    assert adj[0]==pytest.approx(.006)
    assert adj[1]==pytest.approx(.1)
    assert adj[2]==pytest.approx(.12)
    assert adj[5]==pytest.approx(.8)
    assert all(a>=b for a,b in zip(adj,raw))
    with pytest.raises(ValueError):holm([.1,float("nan")])


def test_predeclared_botanical_state_sectors_exact_not_pollinator_guild():
    assert sector("floral_symmetry","zygomorphic","bilateral_zygomorphic")
    assert not sector("floral_symmetry","actinomorphic","bilateral_zygomorphic")
    assert sector("floral_form","bell_campanulate","tube_or_bell")
    assert sector("floral_form","tubular","tube_or_bell")
    assert not sector("floral_form","funnel_trumpet","tube_or_bell")
    assert sector("tube_depth_class","deep","deepest")
    assert not sector("tube_depth_class","intermediate","deepest")
    with pytest.raises(ValueError):
        sector("floral_symmetry","actinomorphic","bee_pollination")
    assert len(MAIN)==6


def test_genus_expected_residual_controls_documented_taxa_mix():
    states={"A first":"tubular","A second":"open_radial",
            "B first":"tubular","B second":"tubular"}
    region=[{"A first","B first"},{"A second","B second"}]
    target=[{"A first","B first"},{"A second","B second"}]
    z=expected_same_region_genus_sector(
        target,states,region,trait="floral_form",sector_name="tube_or_bell")
    assert z==pytest.approx([.75,.75])
    assert (np.array([1.,.5])-z)==pytest.approx([.25,-.25])


def _synthetic_original_labels():
    s={"floral_form":{},"floral_symmetry":{},"tube_depth_class":{}}
    for k,v in (("floral_form",["bell_campanulate","bilabiate","brush_puff",
            "composite_head","funnel_trumpet","open_radial","other_described",
            "papilionaceous","salverform","spurred","tubular","urn_urceolate"]),
           ("floral_symmetry",["actinomorphic","asymmetric","zygomorphic"]),
           ("tube_depth_class",["absent_or_open","shallow","intermediate","deep"])):
        for j,state in enumerate(v):
            s[k][f"Canon{j} plant{j}"]=state
    f=[];geo=[]
    for i in range(90):
        island=f"I{i}"
        geo.append(dict(
            island_id=island,
            analysis_regime="northern_midlatitude" if i<45 else "tropical",
            spatial_block=f"Block{i//3}",
            distance_to_continent_km=3+i*10,area_km2=5+i,
            island_latitude=35+i%20))
        for k in range(12):
            sp=f"Gen{k%4} name{i}_{k}"
            f.append(dict(island_id=island,accepted_species=sp,
                          origin_status="native" if k<7 else "introduced"))
            s["floral_form"][sp]=("open_radial" if k%3 else "tubular")
            s["floral_symmetry"][sp]=("actinomorphic" if k%3 else "zygomorphic")
            s["tube_depth_class"][sp]=("shallow" if k%3 else "deep")
    assigned=pd.DataFrame([
        dict(island_id=f"I{i}",source_mode="geo_k5",source_rank=1,
             entity_ID=f"S{i//5}") for i in range(90)])
    return s,pd.DataFrame(f),pd.DataFrame(geo),assigned


def test_region_only_congener_source_no_triple_completion_required():
    s,f,g,a=_synthetic_original_labels()
    q=build_island_rows(s,f,g,a)
    assert len(q)>len(MAIN)
    y=q[("all_recorded","tropical","floral_form","tube_or_bell")]
    assert len(y)==45
    assert y.spatial_block.nunique()==15
    assert y.n_direct.min()==12
    assert y.genus_relative_response.notna().all()
    assert all("pollination" not in x for x in y.columns)


def test_rank_deficiency_and_missing_source_checks():
    z=make_regression(12)
    q=weighted_cluster_fit(z,cluster="spatial_block",n_boot=99)
    assert q["status"]=="INSUFFICIENT_ORIGINAL_GEOGRAPHIC_SUPPORT"
    s,f,g,a=_synthetic_original_labels()
    with pytest.raises(ValueError):
        build_island_rows(s,f,g.drop(columns=["distance_to_continent_km"]),a)
    with pytest.raises(ValueError):
        build_island_rows(s,pd.concat([f,f.iloc[:1]]),g,a)
