"""GBIF source-background effort missingness and floral recording sensitivity."""
import pandas as pd
import numpy as np
import pytest

from fd.global_guild_effort_trait_coverage import original_global_flora_guild_effort_audit


def _source():
    accepted=[f"Genus{z} plant" for z in range(12)]
    ledger=[]
    for t in ("flower_primary_color","floral_form","floral_symmetry",
              "autonomous_selfing_capacity","inflorescence_display"):
        for idx,sp in enumerate(accepted[:(8 if t!="floral_form" else 5)]):
            ledger.append(dict(accepted_species=sp,trait_name=t,
                normalized_value="white" if idx%2 else "purple",
                evidence_scope="species_direct",evidence_quality="high"))
    flora=[]
    geography=[]
    channel=[]
    for i in range(4001):
        island=f"i{i:05}"
        roster=accepted[:i%12+1]
        for sp in roster:
            flora.append(dict(island_id=island,accepted_species=sp,
                              origin_status="native" if i%2==0 else "unresolved"))
        geography.append(dict(island_id=island,
                analysis_regime="tropical" if i<900 else "northern_midlatitude",
                geological_origin="oceanic" if i<300 else "continental" if i<600 else "unresolved",
                spatial_block=f"block{i%20}"))
        if i<500:
            for c in ("bombus","non_bombus_bees","lepidoptera",
                      "diptera","flower_visiting_birds"):
                state="retained" if i%2 else "disrupted"
                channel.append(dict(
                    island_id=island,channel_id=c,
                    channel_state=state,source_state="available",
                    background_record_count=10+i%50,
                    background_spatial_units=1+i%15,
                    background_temporal_units=1+i%5,
                    distinct_dataset_count=1+i%8))
    return pd.DataFrame(flora),pd.DataFrame(ledger),pd.DataFrame(geography),pd.DataFrame(channel)


def test_effort_fields_original_and_unresolved_global_plant_sources():
    a=original_global_flora_guild_effort_audit(*_source())
    assert a["n_original_botanical_islands"]==4001
    assert a["n_original_GBIF_island_guild_cells"]==4001*5
    g=next(z for z in a["original_plant_flora_and_GBIF_background_effort_by_group"]
        if z["region"]=="GLOBAL_ALL" and z["geological_origin"]=="all"
        and z["guild"]=="bombus")
    assert g["n_GBIF_guild_recorded"]==250
    assert g["n_guild_source_qualified_NONdetection"]==250
    assert g["n_guild_unresolved"]==3501
    assert g["n_guild_state_source_missing_not_a_nondetection"]==3501
    assert g["recorded_background_record_count"]["n_nonmissing"]==250
    assert g["qualified_nondetection_background_spatial_units"]["n_nonmissing"]==250
    assert g["n_background_record_effort_bins_with_both_qualified_states"]>=1
    assert a["scientific_limits"]["missing_background_value_not_inferred_zero"]


def test_background_negative_effort_fails_closed():
    args=list(_source())
    args[3].loc[0,"background_record_count"]=-1
    with pytest.raises(ValueError,match="negative"):
        original_global_flora_guild_effort_audit(*args)
