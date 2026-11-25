import pandas as pd
import pytest

from fd.pollinator_shift_opportunity import (
    GUILDS, audit_shift_opportunity,
)


def _synthetic():
    islands = pd.DataFrame({
        "island_id": [f"secret{i}" for i in range(42)],
        "analysis_regime": ["tropical"]*21+["northern_midlatitude"]*21,
        "spatial_block": [f"block{i//3}" for i in range(42)],
        "n_direct_complete_species": [30+i%2 for i in range(42)],
    })
    rows=[]
    for i in range(42):
        for guild in GUILDS:
            lost=(i<11 and guild=="bombus")
            rows.append({
                "island_id":f"secret{i}",
                "channel_id":guild,
                "channel_state":"disrupted" if lost else "retained",
                "source_state":"available",
            })
    return islands,pd.DataFrame(rows)


def test_observed_bird_with_bombus_nondetection_is_not_certified_shift():
    x,y=_synthetic()
    data,receipt=audit_shift_opportunity(x,y)
    assert len(data)==20
    a=data.query("adequately_undetected_guild == 'bombus' and co_recorded_candidate_guild == 'flower_visiting_birds'").iloc[0]
    assert a.n_islands_candidate_co_record==11
    assert a.n_with_all_other_guilds_retained==11
    assert a.n_spatial_blocks==4
    assert not a.actual_species_level_pollinator_shift_certified
    assert receipt["n_certified_evolutionary_pollinator_shifts"]==0
    assert not {"island_id","spatial_block","genus","provider_id"} & set(data.columns)


def test_unresolved_is_not_absence_and_must_not_be_promoted():
    x,y=_synthetic()
    idx=y.index[(y.channel_id=="bombus") & (y.island_id=="secret0")][0]
    y.loc[idx,"channel_state"]="unresolved"
    z,_=audit_shift_opportunity(x,y)
    row=z.query("adequately_undetected_guild == 'bombus' and co_recorded_candidate_guild == 'flower_visiting_birds'").iloc[0]
    assert row.n_islands_candidate_co_record==10
    y.loc[idx,"source_state"]="unresolved"
    y.loc[idx,"channel_state"]="disrupted"
    with pytest.raises(ValueError,match="source opportunity"):
        audit_shift_opportunity(x,y)


def test_missing_and_duplicate_channel_records_fail_closed():
    x,y=_synthetic()
    with pytest.raises(ValueError,match="does not match"):
        audit_shift_opportunity(x,y.iloc[1:])
    with pytest.raises(ValueError,match="duplicate"):
        audit_shift_opportunity(x,pd.concat([y,y.iloc[[0]]],ignore_index=True))


def test_below_threshold_flora_excluded_not_padded():
    x,y=_synthetic()
    x.loc[0,"n_direct_complete_species"]=29
    with pytest.raises(ValueError,match="species count"):
        audit_shift_opportunity(x,y)
