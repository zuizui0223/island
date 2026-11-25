"""Synthetic-only plant species × exact island visitor evidence gates."""
import pandas as pd
import pytest

from fd.globi_local_real_syndrome_gate import (
    GUILDS,PAIRS,real_same_species_island_visitor_audit
)


def fixtures():
    names=["Alpha plant1","Alpha plant2","Beta plant1"]
    rows=[dict(island_id=f"s{i}",accepted_species=names[i%3],
               origin_status=("native" if i%3<2 else "unresolved"))
          for i in range(928)]
    f=pd.DataFrame(rows)
    s={
        "floral_symmetry":{
            "Alpha plant1":"actinomorphic",
            "Alpha plant2":"zygomorphic",
            "Beta plant1":"actinomorphic",
        },
        "tube_depth_class":{
            "Alpha plant1":"deep",
            "Alpha plant2":"deep",
            "Beta plant1":"shallow",
        },
        "floral_form":{
            "Alpha plant1":"tubular",
            "Alpha plant2":"tubular",
            "Beta plant1":"open_radial",
        },
    }
    links=[
        dict(island_id="s0",accepted_species="Alpha plant1",
             origin_status="native",guild="flower_visiting_birds",
             reference_id="r1",provider_id="archive:dominant",
             evidence_class="flower_visit",latitude=1.0,longitude=2.0),
        dict(island_id="s0",accepted_species="Alpha plant1",
             origin_status="native",guild="lepidoptera",
             reference_id="r2",provider_id="archive:minor",
             evidence_class="pollination",latitude=1.0,longitude=2.0),
        dict(island_id="s1",accepted_species="Alpha plant2",
             origin_status="native",guild="lepidoptera",
             reference_id="r3",provider_id="archive:minor",
             evidence_class="flower_visit",latitude=1.2,longitude=2.2),
    ]
    return f,pd.DataFrame(links),s


def row(out,trait,scope,state,source,guild):
    return next(x for x in out["support"]
                if x["real_same_accepted_species_morphological_pair"]==trait
                and x["original_plant_origin_scope"]==scope
                and x["actual_direct_pair_state"]==state
                and x["GloBI_source_variant"]==source
                and x["documented_original_visitor_guild"]==guild)


def test_real_direct_pair_original_island_guild_documentation():
    f,v,states=fixtures()
    x=real_same_species_island_visitor_audit(
        f,v,states,"archive:dominant")
    assert x["n_original_GIFT_botanical_island_population"]==4505
    assert x["n_original_selected_colour_form_symmetry_complete_islands"]==928
    assert x["n_source_guild_trait_plant_origin_support_cells"]==3*3*2*2*5
    q=row(x,"actinomorphic_deep","strict_native","joint_true",
          "original_citations","flower_visiting_birds")
    assert q["n_original_island_plant_pairs_with_both_direct_traits"]>0
    assert q["n_cited_island_plant_pairs"]==1
    assert q["n_islands_with_source_cited_exact_plant_and_guild"]==1
    assert q["n_source_archive_labels"]==1
    assert q["n_island_plant_pairs_with_2plus_archive_labels"]==0
    d=row(x,"actinomorphic_deep","strict_native","joint_true",
          "largest_archive_omitted","flower_visiting_birds")
    assert d["n_cited_island_plant_pairs"]==0
    m=row(x,"actinomorphic_deep","strict_native","joint_true",
          "largest_archive_omitted","lepidoptera")
    assert m["n_cited_island_plant_pairs"]==1
    assert m["n_island_plant_pairs_with_pollination_predicate"]==1


def test_morphology_negative_is_real_both_axes_and_not_animal_absence():
    f,v,states=fixtures()
    x=real_same_species_island_visitor_audit(f,v,states,"archive:dominant")
    q=row(x,"actinomorphic_deep","strict_native",
          "both_axes_known_not_joint","original_citations","lepidoptera")
    assert q["n_cited_island_plant_pairs"]==1
    assert q["unrecorded_guild_is_not_absent"]


def test_unobserved_second_axis_not_completed_by_marginal_plant_trait():
    f,v,states=fixtures()
    states["tube_depth_class"].pop("Alpha plant1")
    x=real_same_species_island_visitor_audit(f,v,states,"archive:dominant")
    q=row(x,"actinomorphic_deep","strict_native","joint_true",
          "original_citations","flower_visiting_birds")
    assert q["n_cited_island_plant_pairs"]==0


def test_false_locality_or_unknown_animal_guild_not_silently_inferred():
    f,v,states=fixtures()
    v.loc[0,"island_id"]="non_focal"
    with pytest.raises(ValueError,match="NOT an exact original island"):
        real_same_species_island_visitor_audit(
            f,v,states,"archive:dominant")
    f,v,states=fixtures()
    v.loc[0,"guild"]="unknown_insect"
    with pytest.raises(ValueError,match="guilds invalid"):
        real_same_species_island_visitor_audit(
            f,v,states,"archive:dominant")


def test_global_universe_and_original_928_drift_fails_closed():
    f,v,states=fixtures()
    with pytest.raises(ValueError,match="original 4505"):
        real_same_species_island_visitor_audit(
            f,v,states,"archive:dominant",original_global_flora_n=100)
    with pytest.raises(ValueError,match="927 polygon"):
        real_same_species_island_visitor_audit(
            f,v,states,"archive:dominant",verified_polygons_n=400)
    with pytest.raises(ValueError,match="cohort drift"):
        real_same_species_island_visitor_audit(
            f.iloc[:-1],v,states,"archive:dominant")


def test_no_documented_guild_does_not_create_negative_biological_visitor_state():
    f,v,s=fixtures()
    x=real_same_species_island_visitor_audit(
        f,v.iloc[0:0],s,"archive:dominant")
    assert all(p["n_cited_island_plant_pairs"]==0 for p in x["support"])
    assert all(p["unrecorded_guild_is_not_absent"] for p in x["support"])
