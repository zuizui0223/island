import pandas as pd
import pytest
from fd.globi_provider_provenance_gate import (
    ORIGINAL_GLOBI_SHA256, source_coverage_gate, drop_largest_provider,
)


def evidence():
    rows=[]
    for i in range(60):
        genus=f"Genus{i}"
        for provider in range(5):
            for j in range(4):
                rows.append({
                    "plant_genus":genus,
                    "visitor_channel":("bombus","non_bombus_bees","lepidoptera","diptera","flower_visiting_birds")[provider],
                    "provider_id":f"Dataset{provider}",
                    "reference_id":f"Reference_{i}_{provider}_{j}",
                    "argument_type":"SUPPORTS", "evidence_class":"flower_visit",
                })
    return pd.DataFrame(rows)


def test_provider_feasibility_anonymous_and_supported_omission():
    q=source_coverage_gate(evidence(),source_sha256=ORIGINAL_GLOBI_SHA256)
    assert q["n_independent_source_providers"]==5
    assert q["n_observed_plant_genera"]==60
    assert not q["provider_holdout_completed"]
    assert not {"provider_id","plant_genus","reference_id"} & set(q)
    raw,meta=drop_largest_provider(evidence(),source_sha256=ORIGINAL_GLOBI_SHA256)
    assert raw.shape==(60,5)
    assert meta["n_provider_removed"]==1
    assert not meta["heldout_family_prediction_recomputed"]


def test_existing_genera_aggregate_cannot_pass_source_independence():
    summary=pd.DataFrame({"genus":["Genus1"],"effective_channel_number":[1.]})
    with pytest.raises(ValueError,match="provider/ref keys absent"):
        source_coverage_gate(summary,source_sha256=ORIGINAL_GLOBI_SHA256)


def test_sha_and_claim_semantics_are_strict():
    with pytest.raises(ValueError,match="SHA256"):
        source_coverage_gate(evidence(),source_sha256="unverified")
    q=evidence()
    q.loc[0,"argument_type"]="REFUTES"
    with pytest.raises(ValueError,match="REFUTES"):
        source_coverage_gate(q,source_sha256=ORIGINAL_GLOBI_SHA256)


def test_shared_reference_id_dropped_from_both_dataset_groups():
    q=evidence()
    q.loc[0,"reference_id"]=q.loc[4,"reference_id"]
    r=source_coverage_gate(q,source_sha256=ORIGINAL_GLOBI_SHA256)
    assert r["n_cross_provider_reference_collisions"]==1
    _,meta=drop_largest_provider(q,source_sha256=ORIGINAL_GLOBI_SHA256)
    assert meta["n_cross_provider_references_discarded"]==1
