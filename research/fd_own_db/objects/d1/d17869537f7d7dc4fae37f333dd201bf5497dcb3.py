"""Original-source network geometry and positive-guild support, synthetic only."""
import json
import pandas as pd
import pytest

from fd.fullframe_wol_network_bridge import (
    full_botanical_exact_geometry_gate,
    anonymous_network_extension_counts,
    source_scoped_region_network_audit,
)


def _fixture():
    ids=[f"gshhg_2.3.7_h_{i:020x}" for i in range(4505)]
    historical=set(ids[:928])
    flora=pd.DataFrame([{
        "island_id":x,"accepted_species":f"Genus species{k}",
        "origin_status":"native","endemic_status":"not_endemic",
    } for k,x in enumerate(ids)])
    # The old 928 had one unreproducible polygon; all other 3577
    # full-GIFT floras now have potentially verifiable polygon IDs.
    polygons=pd.DataFrame([{
        "island_id":x,"geometry_sha256":x.rsplit("_",1)[1],
        "geometry":"synthetic polygon placeholder",
    } for x in ids if x!=ids[927]])
    receipt={
        "GSHHG_archive_SHA256":
           "8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41",
        "GSHHG_version":"2.3.7","GSHHG_resolution":"h",
        "global_full_universe_identity_verified":False,
        "n_original_focal_islands":21,
        "n_focal_islands_with_exact_historical_polygon":21,
        "n_IDs_missing_from_geometry":9,
    }
    return ids,historical,flora,polygons,receipt


def test_verified_full_flora_geometry_preserves_928_history_and_missingness():
    ids,historical,flora,polygons,receipt=_fixture()
    verified,audit=full_botanical_exact_geometry_gate(
        polygons,flora,historical,receipt)
    assert len(verified)==4504
    assert audit["n_original_accepted_GIFT_botanical_island_floras"]==4505
    assert audit["n_historical_928_direct_floras_with_exact_polygon"]==927
    assert audit["n_original_island_floras_without_verified_polygon"]==1
    assert not audit["all_8265_original_GSHHG_hashes_rebuilt"]
    assert audit["n_source_original_universe_ids_missing_from_rebuild"]==9
    corrupted=polygons.copy()
    corrupted.loc[0,"geometry_sha256"]="abcdef0123456789abcd"
    with pytest.raises(ValueError,match="incompatible"):
        full_botanical_exact_geometry_gate(corrupted,flora,historical,receipt)
    broken=dict(receipt,global_full_universe_identity_verified=True)
    with pytest.raises(ValueError,match="global GSHHG"):
        full_botanical_exact_geometry_gate(polygons,flora,historical,broken)


def test_outside_historical_928_positive_sites_are_new_opportunity_not_absence():
    ids,historical,_,_,_=_fixture()
    old=[dict(network_id=f"old{k}",island_id=ids[k%16]) for k in range(75)]
    new=[dict(network_id=f"new{k}",island_id=ids[1000+k]) for k in range(3)]
    q=anonymous_network_extension_counts(pd.DataFrame(old+new),historical)
    assert q["n_source_point_networks_inside_verified_full_botanical_GIFT_polygons"]==78
    assert q["n_original_927_polygon_network_sites_reproduced"]==75
    assert q["n_additional_verified_botanical_island_network_sites"]==3
    assert q["n_additional_verified_botanical_islands_with_network_sites"]==3
    duplicate=pd.DataFrame(old+[old[0]])
    with pytest.raises(ValueError,match="ambiguous"):
        anonymous_network_extension_counts(duplicate,historical)


def test_regional_genuine_flower_pair_with_exact_accepted_plant_only():
    ids,historical,flora,_,_=_fixture()
    # Assign one exact accepted plant on one original tropical island.
    plant=flora.loc[1000,"accepted_species"]
    legitimate=ids[1000]
    edges=pd.DataFrame([
        dict(network_id="net1",island_id=legitimate,plant=plant,
             guild="bird_family_supported"),
        dict(network_id="net1",island_id=legitimate,
             plant="Another species",guild="bees_family_supported"),
    ])
    floral={"floral_symmetry":{plant:"actinomorphic"},
            "tube_depth_class":{plant:"deep"},
            "floral_form":{plant:"tubular"}}
    regs={i:"tropical" for i in ids}
    out=source_scoped_region_network_audit(
        edges,flora,floral,regs,source_audit={})
    assert out["n_global_parent_botanical_islands"]==4505
    assert len(out["regional_original_source_opportunity"])==1
    x=out["regional_original_source_opportunity"][0]
    assert x["n_original_published_networks_with_same_accepted_plant_and_guild"]==1
    assert x["n_original_islands_with_source_positive_plant_guild_documentation"]==1
    cells=x["true_direct_paired_flower_guild_documentation_cells"]
    yes=[q for q in cells if
         q["botanical_scope"]=="native_label" and
         q["verified_real_joint_morphology"]=="actinomorphic_deep" and
         q["source_network_visitor_taxonomy"]=="bird_family_supported"]
    assert len(yes)==1
    assert yes[0]["n_plant_guild_network_documented_units"]==1
    assert not any('"'+k+'"' in json.dumps(out) for k in (
        "island_id","accepted_species","visitor_species","network_id"))
