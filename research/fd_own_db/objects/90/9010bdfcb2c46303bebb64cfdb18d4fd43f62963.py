"""All original GIFT islands remain present in guild co-occurrence evidence power audit."""
import pandas as pd

from fd.global_guild_cooccurrence_trait_support import guild_cooccurrence_global
from fd.global_partial_pollinator_guild_bridge import CHANNELS


def _fixture():
    n=4001
    names=["Genus one","Genus two","Genus three",
           "Other a","Third a","Fourth a","Fifth a"]
    ledger=pd.DataFrame([
       dict(accepted_species=name,trait_name=trait,
            normalized_value=("white" if i%2 else "purple")
                if trait=="flower_primary_color"
                else ("radial" if i%2 else "tube") if trait=="floral_form"
                else "actinomorphic" if trait=="floral_symmetry" else "self_capable",
            evidence_scope="species_direct",evidence_quality="high")
       for i,name in enumerate(names)
       for trait in ("flower_primary_color","floral_form",
                     "floral_symmetry","autonomous_selfing_capacity")
    ])
    flora=[]
    geography=[]
    channels=[]
    for i in range(n):
        island=f"island{i:04d}"
        roster=names if i<50 else names[:2]
        for name in roster:
            flora.append(dict(island_id=island,accepted_species=name,
                              origin_status="native"))
        geography.append(dict(
            island_id=island,analysis_regime="tropical",
            geological_origin="oceanic" if i<250 else
                              "continental" if i<500 else "unresolved",
            spatial_block=f"block{i%30}"))
        if i<500:
            for guild in CHANNELS:
                if guild=="bombus":
                    state="disrupted" if i%7==0 else "retained"
                else:
                    state="retained" if i%2==0 else "unresolved"
                channels.append(dict(island_id=island,channel_id=guild,
                         channel_state=state,source_state="available"))
    return pd.DataFrame(flora),ledger,pd.DataFrame(geography),pd.DataFrame(channels)


def test_global_source_guild_assembly_retains_original_4001_island_denominator():
    flora,ledger,geo,channels=_fixture()
    x=guild_cooccurrence_global(flora,ledger,geo,channels)
    assert x["n_original_global_plant_flora_islands"]==4001
    assert x["n_original_guilds"]==5
    panel=next(p for p in x["original_global_guild_co_observation_by_region_and_geology"]
               if p["region"]=="GLOBAL_ALL" and p["geological_origin"]=="all")
    assert panel["n_ALL_islands_in_original_botanical_flora"]==4001
    assert panel["n_at_least_one_recorded_guild"]>0
    assert panel["n_any_guild_unresolved"]>3501
    target=next(t for t in panel["source_guild_nondetection_plus_other_recorded_pairs"]
        if t["originally_qualified_NONDETECTION_guild"]=="bombus"
        and t["other_guild_GBIF_RECORDED"]=="flower_visiting_birds")
    assert target["n_both_source_qualified_guild_opportunity"]>0
    assert target["n_candidate_with_5_direct__all_recorded__flower_primary_color"]>0
    assert target["n_candidate_with_5_direct__all_recorded__flower_primary_color__floral_form"]>0
    assert x["source_scope"]["independent_pollinator_GBIF_nondetection_requires_source_and_effort"]
    assert "island_id" not in x["original_global_guild_co_observation_by_region_and_geology"][0]


def test_no_original_three_trait_complete_case_gate():
    flora,ledger,geo,channels=_fixture()
    # Most original source islands have only 2 Direct plants, and ALL
    # 4001 must be retained; individual five-plant outcomes just have
    # their own smaller evidence denominator.
    result=guild_cooccurrence_global(flora,ledger,geo,channels)
    assert result["original_GBIF_evaluable_all_five_guilds"]<4001
    assert result["n_original_global_plant_flora_islands"]==4001
    assert result["source_scope"]["observed_plant_species_with_fewer_than_five_trait_states_still_in_data"]
