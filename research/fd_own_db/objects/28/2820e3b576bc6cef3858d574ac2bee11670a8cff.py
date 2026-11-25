import pytest
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_dopi_v2_positive_visitors import (
    audit_positive_visitor_rows,identify_visitor_group,source_metadata,DOI
)


def test_morphology_never_infers_effective_bird_butterfly_pollination():
    assert identify_visitor_group("Hymenoptera","Apidae","Bombus terrestris")=="Bombus"
    assert identify_visitor_group("Hymenoptera","Apidae","Apis mellifera")=="other_bee"
    assert identify_visitor_group("Lepidoptera","Sphingidae","Macroglossum stellatarum")=="Sphingidae_hawkmoth"
    assert identify_visitor_group("Lepidoptera","Nymphalidae","Vanessa atalanta")=="butterfly_family"
    assert identify_visitor_group("Lepidoptera",None,"unknown butterfly")=="unresolved_Lepidoptera_family"
    assert identify_visitor_group("Diptera","Syrphidae","Episyrphus balteatus")=="Syrphidae_hoverfly"
    assert identify_visitor_group("Passeriformes",None,"bird")=="bird_order_not_effective_pollinator"


def metadata():
    return {"metadata":{"title":"DoPI The Database of Pollinator Interactions",
                        "doi":DOI,"version":"2.0"},
            "files":[{"key":"DoPI_dataset_2.0.csv","size":1000,
                      "checksum":"md5:xxxx",
                      "links":{"self":"https://zenodo.org/api/records/22876365/files/DoPI_dataset.csv/content"}}]}


def test_metadata_is_version_qualified_and_exact_download():
    q=source_metadata(metadata())
    assert q["doi"]==DOI
    assert q["dataset_file"]=="DoPI_dataset_2.0.csv"
    assert q["version"]=="2.0"


def test_wrong_dopi_release_or_missing_dataset_not_silently_substituted():
    d=metadata()
    d["metadata"]["doi"]="10.5281/zenodo.999"
    with pytest.raises(ValueError,match="not immutable"):
        source_metadata(d)
    d=metadata()
    d["files"]=[]
    with pytest.raises(ValueError,match="single"):
        source_metadata(d)


def test_positive_links_do_not_turn_missing_study_guild_into_absence():
    rows=[
     dict(plant_species="Plant a",pollinator_species="Bombus terrestris",
          pollinator_family="Apidae",pollinator_order="Hymenoptera",
          doi="study1",year="2022",pollinator_survey="All",plant_survey="All"),
     dict(plant_species="Plant a",pollinator_species="Bombus terrestris",
          pollinator_family="Apidae",pollinator_order="Hymenoptera",
          doi="study1",year="2022",pollinator_survey="All",plant_survey="All"),
     dict(plant_species="Plant a",pollinator_species="Vanessa atalanta",
          pollinator_family="Nymphalidae",pollinator_order="Lepidoptera",
          doi="study1",year="unknown",pollinator_survey="All",plant_survey="All"),
     dict(plant_species="Plant b",pollinator_species="Episyrphus balteatus",
          pollinator_family="Syrphidae",pollinator_order="Diptera",
          doi="study2",year="2019",pollinator_survey="Group",plant_survey="Group"),
     dict(plant_species="Plant b",pollinator_species="missing",
          pollinator_family="",pollinator_order="Lepidoptera",
          doi="",publication="",year="2019"),
    ]
    z=audit_positive_visitor_rows(rows,fields=set().union(*(x.keys() for x in rows)))
    assert z["n_original_raw_records"]==5
    assert z["n_positive_taxon_study_rows"]==4
    assert z["n_study_plant_positive_units"]["Bombus"]==1
    assert z["n_study_plant_positive_units"]["butterfly_family"]==1
    assert z["n_studies_by_positive_guild"]["Syrphidae_hoverfly"]==1
    assert z["year_min"]==2019 and z["year_max"]==2022
    assert z["flags"]["no_zero_absence_derived"] is True


def test_missing_plant_or_pub_fields_fail_closed():
    with pytest.raises(ValueError,match="taxonomic"):
        audit_positive_visitor_rows([],fields={"plant_species","publication"})
    with pytest.raises(ValueError,match="study"):
        audit_positive_visitor_rows([],fields={
            "plant_species","pollinator_species","pollinator_family","pollinator_order"})

def test_pollination_evidence_grades_and_source_study_deduplication():
    from audit_dopi_v2_positive_visitors import pollination_evidence_grade
    assert pollination_evidence_grade("1")=="confirmed_pollen_transfer_or_seed_set"
    assert pollination_evidence_grade("2")=="pollen_carriage_inferred"
    assert pollination_evidence_grade("3")=="flower_visit_only_unconfirmed_pollen"
    assert pollination_evidence_grade("4")=="robbery_herbivory_predation_or_nonpollinating"
    assert pollination_evidence_grade("unknown")=="missing_or_unrecognized_evidence_grade"
    rows=[]
    for study,plant,animal,fam,order,grade in [
        ("ref1","Plant a","Bombus terrestris","Apidae","Hymenoptera","1"),
        ("ref1","Plant a","Bombus terrestris","Apidae","Hymenoptera","3"),
        ("ref2","Plant b","Bombus lucorum","Apidae","Hymenoptera","2"),
        ("ref2","Plant c","Apis mellifera","Apidae","Hymenoptera","4"),
        ("ref3","Plant d","Vanessa atalanta","Nymphalidae","Lepidoptera","3"),
        ("ref3","Plant e","Macroglossum stellatarum","Sphingidae","Lepidoptera",""),
    ]:
        rows.append(dict(doi=study,publication=study,plant_species=plant,
                         pollinator_species=animal,pollinator_family=fam,
                         pollinator_order=order,pollination=grade))
    from audit_dopi_v2_positive_visitors import audit_positive_visitor_rows
    z=audit_positive_visitor_rows(rows,fields=set(rows[0]))
    gr=z["n_study_plant_guild_units_by_best_documented_evidence"]
    assert gr["confirmed_pollen_transfer_or_seed_set"]==1
    assert gr["pollen_carriage_inferred"]==1
    assert gr["flower_visit_only_unconfirmed_pollen"]==1
    assert gr["robbery_herbivory_predation_or_nonpollinating"]==1
    assert gr["missing_or_unrecognized_evidence_grade"]==1
    assert z["n_original_pollination_raw_code_counts"]["1"]==1
    assert z["n_original_pollination_raw_code_counts"]["3"]==2
    assert z["flags"]["category_4_not_a_pollination_event"]
