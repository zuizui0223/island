"""Synthetic study-method gate: metadata description != negative guild survey."""
import pandas as pd
import pytest

from fd.euppollnet_study_protocol_audit import study_year_original_protocol_gate


def _fixture():
    pos=pd.DataFrame([
        dict(Study_id="S1",Network_id="N1",Date="2020-06-18",
             island_id="I1",plant="Alpha alba",
             Pollinator_order="Hymenoptera",
             Pollinator_family="Apidae",
             Pollinator_accepted_name="Bombus terrestris"),
        dict(Study_id="S1",Network_id="N1",Date="2020-06-18",
             island_id="I1",plant="Beta rubra",
             Pollinator_order="Lepidoptera",
             Pollinator_family="Papilionidae",
             Pollinator_accepted_name="Papilio machaon"),
        dict(Study_id="S2",Network_id="N2",Date="2021-05-09",
             island_id="I2",plant="Gamma alba",
             Pollinator_order="Hymenoptera",
             Pollinator_family="Apidae",
             Pollinator_accepted_name="Apis mellifera"),
    ])
    meta=pd.DataFrame([
        dict(Study_id="S1",Year="2020",Sampling_method="Transect",
             Sampling_days="5",Sampling_period="30",
             Taxa_recorded="Hymenoptera, Lepidoptera",
             Hymenoptera_percent="0.7",Lepidoptera_percent="0.3"),
        dict(Study_id="S3",Year="2020",Sampling_method="Transect",
             Sampling_days="5",Sampling_period="30",
             Taxa_recorded="Bees",
             Hymenoptera_percent="1.0",Lepidoptera_percent=None),
        dict(Study_id="S3",Year="2020",Sampling_method="Plot",
             Sampling_days="3",Sampling_period="30",
             Taxa_recorded="Bees",
             Hymenoptera_percent="1.0",Lepidoptera_percent=None),
    ])
    return pos,meta


def test_source_documented_methods_do_not_infer_unrecorded_guild_absence():
    visits,meta=_fixture()
    q=study_year_original_protocol_gate(visits,meta)
    assert q["n_original_parent_island_floras"]==4505
    assert q["n_original_positive_same_plant_source_visit_rows"]==3
    assert q["n_original_distinct_source_study_ids"]==2
    assert q["n_original_distinct_source_study_years"]==2
    assert q["n_study_years_with_exact_unique_public_metadata"]==1
    assert q["n_study_years_without_unique_original_metadata"]==1
    assert q["n_original_source_metadata_study_year_keys_ambiguous"]==1
    assert q["n_matched_study_years_with_both_Hymenoptera_and_Lepidoptera_reported_percentages"]==1
    assert q["n_source_study_years_with_positive_links_from_two_animal_guild_classes"]==1
    assert q["original_method_labels_frequency"]=={"Transect":1}
    assert q["original_study_method_is_not_complete_guild_nondetection_protocol"]
    assert q["no_inferred_missing_visits_or_pollen_delivery"]


def test_missing_original_source_method_fields_and_false_parent_rejected():
    visits,meta=_fixture()
    with pytest.raises(ValueError,match="4505"):
        study_year_original_protocol_gate(visits,meta,
                                         n_original_global_GIFT_floras=928)
    with pytest.raises(ValueError,match="metadata columns missing"):
        study_year_original_protocol_gate(visits,
                                          meta.drop(columns=["Taxa_recorded"]))
