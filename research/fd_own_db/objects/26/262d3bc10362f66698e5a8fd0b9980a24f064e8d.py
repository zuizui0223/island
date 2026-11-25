"""Photo-measured source overlap is an opportunity audit, not local adaptation."""
import numpy as np
import pandas as pd
import pytest

from fd.fcp_geological_endemic_opportunity import (
    audit_fcp_endemic_geology_opportunity
)


def setup():
    all_names=["Aex islandensis","Bex native","Aex mainland"]+[
        f"Fillergenus species{i}" for i in range(4472)]
    n_class=np.full(4475,50,dtype=int)
    n_class[:3357]+=1
    fcp=pd.DataFrame({
        "species":all_names,
        "n_source_photo":np.full(4475,100),
        "n_classified_photo":n_class,
        "D_photo_pairwise_unbiased":np.full(4475,0.4),
        "W_axis_sq_alignment":np.full(4475,0.5),
        "n_local50km_colour_photo_pairs":np.full(4475,15),
    })
    flora=[];geo=[];assign=[]
    for i in range(4):
        island=f"island{i}"
        geo.append(dict(island_id=island,
            geological_origin="oceanic" if i%2 else "continental",
            analysis_regime="tropical",spatial_block=f"block{i}"))
        for sp,origin,endemic,floristic in (
            ("Aex islandensis","native","endemic","endemic"),
            ("Bex native","native","nonendemic","native_nonendemic"),
        ):
            flora.append(dict(island_id=island,accepted_species=sp,
                origin_status=origin,endemic_status=endemic,
                floristic_status=floristic))
        for rank in range(1,6):
            assign.append(dict(island_id=island,entity_ID=f"s{rank}",
                               source_mode="geo_k5",source_rank=str(rank)))
    donor=pd.DataFrame([dict(entity_ID=f"s{k}",work_species="Aex mainland")
                        for k in range(1,6)])
    return fcp,pd.DataFrame(flora),pd.DataFrame(geo),pd.DataFrame(assign),donor


def test_exact_fcp_endemic_source_congener_and_repeat_photo_guard():
    r=audit_fcp_endemic_geology_opportunity(*setup())
    assert r["n_original_geologically_typed_islands"]==4
    assert r["n_FCP_species_original"]==4475
    assert r["n_FCP_original_photo_measurements"]==447500
    assert r["n_FCP_colour_classifiable_original"]==227107
    a=next(x for x in r["original_regional_photo_opportunity"] if
           x["region"]=="tropical" and
           x["geological_origin"]=="oceanic" and
           x["native_endemic_status"]=="endemic")
    assert a["n_original_exact_FCP_measured_species_rows"]==2
    assert a["n_original_exact_FCP_distinct_species"]==1
    assert a["n_species_with_GIFT_candidate_mainland_FCP_same_genus"]==1
    assert a["n_independent_geographic_blocks_with_FCP_species"]==2
    assert r["scope"]["FCP_species_wide_photo_colour_is_not_island_local_photo_trait"]


def test_endemic_status_contradiction_fail_closed():
    args=list(setup())
    args[1].loc[0,"origin_status"]="introduced"
    with pytest.raises(ValueError,match="contradictory"):
        audit_fcp_endemic_geology_opportunity(*args)
