"""Synthetic checks of the fixed-dyad, source-pool dependent genus null."""
import json
import numpy as np
import pandas as pd
import pytest
from fd_island_beta_regional_genus import analyze, _subset_prepared, _regional_donor_audit
from fd_island_functional_beta import prepare_endpoint, compare_block


def _fixture():
    states={"flower_primary_color":{},"floral_form":{}}
    rows=[];g=[]
    for i in range(12):
        region="northern_midlatitude" if i<6 else "tropical"
        for k in range(6):
            species=f"Gen{k%3} sp{i}_{k}"
            states["flower_primary_color"][species]=(
                "white" if (i+k)%3 else "violet")
            states["floral_form"][species]=(
                "open_radial" if (i*2+k)%4 else "bell")
            rows.append(dict(island_id=f"I{i}", accepted_species=species,
                             origin_status="native" if k<5 else "introduced"))
        g.append(dict(island_id=f"I{i}",
                      analysis_regime=region,
                      spatial_block=f"{region}_b{i//2}",
                      island_latitude=10 if i<6 else -10,
                      island_longitude=(i%6)*3))
    return states,pd.DataFrame(rows),pd.DataFrame(g)


def _find(out,scope,region,gate="n5",distance="all_distances"):
    return next(x for x in out["comparisons"]
                if x["scope"]==scope and x["region"]==region and
                   x["support"]==gate and x["distance_class"]==distance)


def test_regional_and_global_references_share_exact_same_pairs():
    states,flora,geo=_fixture()
    out=analyze(states,flora,geo,batch=2,min_pairs=1,min_islands=2,min_blocks=2)
    assert out["n_original_botanical_flora"]==12
    assert out["eligible_by_scope"]["all_recorded"]["n_eligible_islands"]==12
    assert out["eligible_by_scope"]["original_native"]["n_eligible_islands"]==12
    q=_find(out,"all_recorded","northern_midlatitude")
    assert q["n_pairs"]==12 # six islands, three pairs within blocks excluded
    assert q["n_active_islands"]==6 and q["n_physical_blocks"]==3
    assert q["region_minus_global_residual_D2"]==pytest.approx(
        q["region_genus_residual_D2"]-q["global_genus_residual_D2"])
    assert q["global_genus_residual_D2"]==pytest.approx(
        q["mean_observed_D2"]-q["global_genus_expected_D2"])
    assert q["region_genus_residual_D2"]==pytest.approx(
        q["mean_observed_D2"]-q["region_genus_expected_D2"])
    assert q["n_region_reference_species"]==36
    assert out["claims"]["exclude_shared_species_on_full_original_botanical_rosters"]
    assert not any(s in json.dumps(out) for s in (
        '"island_id"','"accepted_species"','"spatial_block"','"I0"'))


def test_unused_donor_changes_global_but_not_region_local_reference():
    states,flora,geo=_fixture()
    baseline=analyze(states,flora,geo,batch=3,min_pairs=1,min_islands=2,min_blocks=2)
    # A globally measured Direct species not recorded on ANY focal island
    # modifies the global null, not the regional in-sample exchange support.
    states2={k:dict(v) for k,v in states.items()}
    states2["flower_primary_color"]["Gen0 extra_donor"]="red"
    states2["floral_form"]["Gen0 extra_donor"]="tubular"
    added=analyze(states2,flora,geo,batch=3,min_pairs=1,min_islands=2,min_blocks=2)
    a=_find(baseline,"all_recorded","northern_midlatitude")
    b=_find(added,"all_recorded","northern_midlatitude")
    assert a["n_pairs"]==b["n_pairs"]
    assert a["mean_observed_D2"]==pytest.approx(b["mean_observed_D2"])
    assert a["region_genus_expected_D2"]==pytest.approx(
        b["region_genus_expected_D2"])
    assert a["global_genus_expected_D2"]!=pytest.approx(
        b["global_genus_expected_D2"],abs=1e-10)


def test_full_flora_overlap_excluded_even_when_native_endpoints_disjoint():
    states,flora,geo=_fixture()
    one=flora.loc[flora.island_id=="I0"].iloc[:1].copy()
    # A non-trait-coded accepted plant shared across target blocks must remove pair.
    # Add to both I0 and I2 while preserving original native trait carriers.
    for island in ("I0","I2"):
        flora=pd.concat([flora,pd.DataFrame([dict(
            island_id=island,accepted_species="X common_non_trait",
            origin_status="unresolved")])],ignore_index=True)
    out=analyze(states,flora,geo,batch=2,min_pairs=1,min_islands=2,min_blocks=2)
    assert _find(out,"all_recorded","northern_midlatitude")["n_pairs"]==11
    assert _find(out,"original_native","northern_midlatitude")["n_pairs"]==11


def test_local_reference_degree_zero_cannot_be_independent_validation():
    # With no regional within-genus trait variation, the in-sample null
    # is deterministic and reproduces observed category frequency distances.
    states={"flower_primary_color":{},"floral_form":{}}
    rosters=[]
    for i in range(4):
        z=set()
        for k in range(5):
            sp=f"UniqueGen{i}_{k} plant"
            z.add(sp)
            states["flower_primary_color"][sp]="white" if i<2 else "blue"
            states["floral_form"][sp]="open" if i<2 else "bell"
        rosters.append(z)
    joint={sp:(states["flower_primary_color"][sp],states["floral_form"][sp])
           for z in rosters for sp in z}
    p=prepare_endpoint(rosters,joint)
    m=compare_block(p,0,4)
    assert np.max(np.abs(
        m["squared_frequency_distance"]-
        m["genus_null_squared_distance"]))<1e-12
    audit=_regional_donor_audit(joint,set(joint))
    assert audit["n_region_genera_with_multiple_direct_profiles"]==0
    assert audit["fraction_distinct_region_species_in_variable_genus"]==0


def test_stable_against_geographic_batch_and_original_flora_reorder():
    states,flora,geo=_fixture()
    one=analyze(states,flora,geo,batch=2,min_pairs=1,min_islands=2,min_blocks=2)
    two=analyze(states,flora.sample(frac=1,random_state=19),
                geo.sample(frac=1,random_state=5),batch=5,
                min_pairs=1,min_islands=2,min_blocks=2)
    assert one["eligible_by_scope"]==two["eligible_by_scope"]
    a=sorted(one["comparisons"],key=lambda x:(
        x["scope"],x["region"],x["support"],x["distance_class"]))
    b=sorted(two["comparisons"],key=lambda x:(
        x["scope"],x["region"],x["support"],x["distance_class"]))
    for x,y in zip(a,b):
        for key in x:
            assert x[key]==pytest.approx(y[key]) if isinstance(
                x[key],float) else x[key]==y[key]


def test_invalid_sources_fail_closed():
    s,f,g=_fixture()
    with pytest.raises(ValueError):
        analyze(s,pd.concat([f,f.iloc[:1]]),g)
    with pytest.raises(ValueError):
        analyze(s,f,g.drop(columns=["island_longitude"]))
    with pytest.raises(ValueError):
        analyze(s,f,g,batch=0)
