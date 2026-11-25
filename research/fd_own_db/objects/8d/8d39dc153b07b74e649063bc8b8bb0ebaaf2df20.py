import pandas as pd
import pytest
from fd.pollinator_shift_source_pool import source_pool_audit,GUILDS


def fixture(n=12):
    flora,traits,sources,gift,guilds,geo,globi=[],[],[],[],[],[],[]
    # Six genera, 40 species; 30 Direct species on every island.
    for k in range(40):
        sp=f"Plant{k%6} species{k}"
        for name,value in (
            ("flower_primary_color","white"),
            ("floral_form","open_radial" if k%2 else "tubular"),
            ("floral_symmetry","actinomorphic" if k%2 else "zygomorphic"),
        ):
            traits.append(dict(accepted_species=sp,trait_name=name,normalized_value=value,
                               evidence_scope="species_direct",evidence_quality="high"))
        for j in range(5):
            gift.append(dict(entity_ID=100+j,work_species=sp))
    for k in range(n):
        island=f"Island{k}"
        geo.append(dict(island_id=island,analysis_regime="tropical",
                        spatial_block=f"block{k//2}"))
        for i in range(5):
            sources.append(dict(island_id=island,source_mode="geo_k5",
                                source_rank=i+1,entity_ID=100+i))
        for j in range(30):
            flora.append(dict(island_id=island,
                              accepted_species=f"Plant{j%6} species{j}",
                              origin_status="native"))
        for group in GUILDS:
            guilds.append(dict(
                island_id=island,channel_id=group,source_state="available",
                channel_state="disrupted" if group=="bombus" else "retained",
            ))
    for j in range(6):
        globi.append(dict(genus=f"Plant{j}",n_independent_references=8,
            **dict(zip(("refs__"+name for name in GUILDS),(2,3,1,4,2)))))
    return tuple(pd.DataFrame(rows) for rows in
                 (flora,traits,sources,gift,guilds,geo,globi))


def test_cohort_and_GIFT_species_overlap_but_no_shift_inferred():
    out,receipt=source_pool_audit(*fixture())
    assert len(out)==1 and out.n_islands.iloc[0]==12
    assert receipt["n_islands_atleast_30_direct_complete_species"]==12
    assert receipt["n_bombus_nondetected_bird_recorded"]==12
    assert receipt["n_candidate_island_species_records"]==360
    assert receipt["n_native_same_named_source_species_records"]==360
    assert receipt["confirmed_same_species_evolutionary_pollinator_shift_in_data"] is False
    assert receipt["direct_species_traits_have_independent_island_populations"] is False
    assert not {"island_id","genus","accepted_species","entity_ID","spatial_block"} & set(out.columns)


def test_qualified_Bombus_nondetection_cannot_be_inferred_without_source():
    v=list(fixture())
    ix=v[4].index[v[4].channel_id.eq("bombus")][0]
    v[4].loc[ix,"source_state"]="unresolved"
    with pytest.raises(ValueError,match="source opportunity"):
        source_pool_audit(*v)


def test_incomplete_five_source_assignment_fails():
    v=list(fixture())
    v[2]=v[2].iloc[1:]
    with pytest.raises(ValueError,match="Five distinct"):
        source_pool_audit(*v)


def test_no_auto_promote_genus_consensus_to_Direct():
    v=list(fixture())
    v[1].loc[0,"evidence_scope"]="genus_consensus"
    out,receipt=source_pool_audit(*v)
    assert out.empty
    assert receipt["n_islands_atleast_30_direct_complete_species"]==0


def test_suppress_small_region_cells_without_exporting_original_island_rows():
    out,receipt=source_pool_audit(*fixture(n=8))
    assert out.empty
    assert receipt["n_bombus_nondetected_bird_recorded"]==8


def test_duplicate_GloBI_genus_and_trait_claim_fail_closed():
    v=list(fixture())
    v[-1]=pd.concat([v[-1],v[-1].iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError,match="Duplicate GloBI"):
        source_pool_audit(*v)
    v=list(fixture())
    v[1]=pd.concat([v[1],v[1].iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError,match="Duplicate original"):
        source_pool_audit(*v)
