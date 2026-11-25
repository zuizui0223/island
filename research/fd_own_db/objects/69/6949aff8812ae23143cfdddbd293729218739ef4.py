"""Source-pinned original European flower count availability: synthetic control."""
import gzip
import pandas as pd
import pytest
from fd.euppollnet_flower_availability import (
    audit_public_flower_source,source_original_flower_receipt,
)


def _fixture():
    f=pd.DataFrame([
        dict(Flower_data_merger="X",Flower_count="10",
             Plant_accepted_name="A alba"),
        dict(Flower_data_merger="Y",Flower_count="5",
             Plant_accepted_name="B alba"),
        dict(Flower_data_merger="Y",Flower_count="7",
             Plant_accepted_name="C alba"),
        dict(Flower_data_merger="Z",Flower_count=None,
             Plant_accepted_name="D alba"),
    ]).to_csv(index=False).encode()
    i=gzip.compress(pd.DataFrame([
        dict(Flower_data_merger="X",Study_id="One"),
        dict(Flower_data_merger="X",Study_id="One"),
        dict(Flower_data_merger="Y",Study_id="Two"),
        dict(Flower_data_merger="Q",Study_id="Three"),
    ]).to_csv(index=False).encode())
    return f,i


def test_ambiguous_public_flower_count_keys_not_counted_as_verified_availability():
    f,i=_fixture()
    q=audit_public_flower_source(f,i,allow_synthetic_source=True)
    assert q["n_original_flowers_raw_rows"]==4
    assert q["n_original_interaction_raw_rows"]==4
    assert q["n_original_distinct_join_keys_shared"]==2
    assert q["n_distinct_flower_join_keys_with_multiple_source_rows"]==1
    assert q["n_original_join_keys_shared_and_unique_flower_count_rows"]==1
    assert q["n_original_interaction_rows_on_any_shared_join_key"]==3
    assert q["n_original_interaction_rows_with_unique_flower_key"]==2
    assert q["original_quantitative_flower_abundance_column_candidates"][
        "Flower_count"]["n_nonnegative_numeric"]==3
    assert q["conclusions"]["flower_source_duplicate_join_keys_not_auto_aggregated"]


def test_missing_published_flower_merger_key_fails_closed():
    f,i=_fixture()
    wrong=f.replace(b"Flower_data_merger",b"unknown_field")
    with pytest.raises(ValueError,match="Flower_data_merger"):
        audit_public_flower_source(wrong,i,allow_synthetic_source=True)


def test_unpinned_flower_count_bytes_not_accepted():
    f,_=_fixture()
    with pytest.raises(ValueError,match="source changed"):
        source_original_flower_receipt(f)


def _synthetic_gift_positive_and_flower_availability():
    source=pd.DataFrame([
      dict(Study_id="S1",Network_id="N1",island_id="I1",
           plant="Alpha alba",origin_status="native",endemic_status="endemic",
           Date="2020-06-18",Flower_data_merger="A",
           Pollinator_order="Hymenoptera",Pollinator_family="Apidae",
           Pollinator_accepted_name="Bombus terrestris"),
      dict(Study_id="S1",Network_id="N1",island_id="I1",
           plant="Alpha alba",origin_status="native",endemic_status="endemic",
           Date="2020-06-18",Flower_data_merger="A",
           Pollinator_order="Lepidoptera",Pollinator_family="Papilionidae",
           Pollinator_accepted_name="Papilio machaon"),
      dict(Study_id="S2",Network_id="N2",island_id="I2",
           plant="Beta rubra",origin_status="native",endemic_status="not_endemic",
           Date="2021-07-03",Flower_data_merger="B",
           Pollinator_order="Lepidoptera",Pollinator_family="Sphingidae",
           Pollinator_accepted_name="Macroglossum stellatarum"),
      dict(Study_id="S2",Network_id="N2",island_id="I2",
           plant="Beta rubra",origin_status="native",endemic_status="not_endemic",
           Date="2021-07-03",Flower_data_merger="C",
           Pollinator_order="Lepidoptera",Pollinator_family="Sphingidae",
           Pollinator_accepted_name="Macroglossum stellatarum"),
      dict(Study_id="S2",Network_id="N2",island_id="I2",
           plant="Beta rubra",origin_status="native",endemic_status="not_endemic",
           Date="2021-07-03",Flower_data_merger="D",
           Pollinator_order="Lepidoptera",Pollinator_family="Sphingidae",
           Pollinator_accepted_name="Macroglossum stellatarum"),
    ])
    flowers=pd.DataFrame([
      dict(Study_id="S1",Network_id="N1",Day="18",Month="6",Year="2020",
           Flower_data_merger="A",Flower_count="20",Units="Flowers"),
      dict(Study_id="S2",Network_id="N2",Day="3",Month="7",Year="2021",
           Flower_data_merger="B",Flower_count="9",Units="Flower heads"),
      dict(Study_id="S2",Network_id="N2",Day="3",Month="7",Year="2021",
           Flower_data_merger="B",Flower_count="8",Units="Flower heads"),
      dict(Study_id="S2",Network_id="N2",Day="4",Month="7",Year="2021",
           Flower_data_merger="C",Flower_count="8",Units="Flower heads"),
      dict(Study_id="S2",Network_id="N2",Day="3",Month="7",Year="2021",
           Flower_data_merger="D",Flower_count="0",Units="Flower heads"),
    ])
    traits=dict(floral_symmetry={"Alpha alba":"actinomorphic",
                  "Beta rubra":"zygomorphic"},
        tube_depth_class={"Alpha alba":"deep","Beta rubra":"deep"},
        floral_form={"Alpha alba":"tubular","Beta rubra":"bell_campanulate"},
        flower_primary_color={"Alpha alba":"white","Beta rubra":"red_pink"})
    return source,flowers,traits


def test_exact_site_day_study_gift_flower_supply_rejects_duplicate_zero_wrong_day():
    from fd.euppollnet_flower_availability import (
      qualified_gift_flower_supply_on_positive_visits)
    visits,counts,traits=_synthetic_gift_positive_and_flower_availability()
    q=qualified_gift_flower_supply_on_positive_visits(visits,counts,traits)
    assert q["parent_original_GIFT_botanical_floras"]==4505
    assert q["n_original_source_positive_same_island_botanical_visits"]==5
    assert q["n_original_source_flower_keys_ambiguous_not_automatically_summed"]==1
    assert q["n_original_visitor_rows_with_unique_original_flower_key"]==4
    assert q["n_original_visitors_with_same_published_study_network_calendar_day"]==3
    assert q["n_original_positive_visitors_with_nonzero_numeric_flower_count_and_unit"]==2
    assert q["n_original_GIFT_islands_with_qualified_flower_exposure"]==1
    assert q["n_original_source_studies_with_qualified_flower_exposure"]==1
    assert q["original_qualified_source_flower_unit_labels_counts"]=={"Flowers":2}
    matched=[r for r in q["true_source_same_species_flower_x_animal_positive_exposure_cells"] if
       r["botanical_scope"]=="native_endemic" and
       r["genuine_same_species_flower_pair"]=="actinomorphic_deep" and
       r["original_animal_taxonomic_guild"]=="butterfly_family"]
    assert len(matched)==1 and matched[0]["n_original_qualified_flower_count_linked_visitor_rows"]==1
    assert q["limitations"]["no_calibrated_pollen_deposition_or_seed_set_estimated"]


def test_gift_original_flower_denominator_requires_all_source_components():
    from fd.euppollnet_flower_availability import (
      qualified_gift_flower_supply_on_positive_visits)
    visits,counts,traits=_synthetic_gift_positive_and_flower_availability()
    with pytest.raises(ValueError,match="4505"):
        qualified_gift_flower_supply_on_positive_visits(
            visits,counts,traits,original_GIFT_botanical_floras=927)
    with pytest.raises(ValueError,match="missing"):
        qualified_gift_flower_supply_on_positive_visits(
            visits,counts.drop(columns=["Units"]),traits)
