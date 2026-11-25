"""Bounded iNaturalist photo discovery: deterministic source cohorts, fake API tests."""
import pandas as pd

from fd.oceanic_endemic_inat_photo_pilot import (
    select_original_endemic_colour_pairs,query_exact_source_taxon,
    original_photo_pair_pilot,photo_observation_qualifies
)


def source():
    genera=("Gena","Genb","Genc","Gend")
    rows=[];geo=[];assign=[];donors=[];traits=[]
    for genus in genera:
        traits.append(dict(accepted_species=f"{genus} mainland",
            trait_name="flower_primary_color",normalized_value="red",
            evidence_scope="species_direct",evidence_quality="high"))
    for origin in ("oceanic","continental"):
        for k in range(5):
            island=f"{origin}_{k}"
            geo.append(dict(
                island_id=island,analysis_regime="tropical",
                geological_origin=origin,spatial_block=f"b{origin}_{k}",
                weigelt_archipelago=f"arch{origin}_{k}"))
            for genus in genera:
                taxon=f"{genus} {origin}_{k}"
                rows.append(dict(
                    island_id=island,accepted_species=taxon,
                    origin_status="native",endemic_status="endemic",
                    floristic_status="endemic"))
                traits.append(dict(accepted_species=taxon,
                    trait_name="flower_primary_color",normalized_value="white",
                    evidence_scope="species_direct",evidence_quality="high"))
            for rank in range(1,6):
                src=f"s{rank}"
                assign.append(dict(island_id=island,source_mode="geo_k5",
                    source_rank=rank,entity_ID=src))
    for rank in range(1,6):
        for genus in genera:
            donors.append(dict(entity_ID=f"s{rank}",work_species=f"{genus} mainland"))
    return [pd.DataFrame(traits),pd.DataFrame(rows),pd.DataFrame(geo),
            pd.DataFrame(assign),pd.DataFrame(donors)]


def fake(endpoint,params):
    if endpoint=="taxa":
        name=params["q"]
        return {"results":[dict(id=100+len(name),name=name,rank="species",
                                is_active=True)]}
    if endpoint=="observations":
        taxon_id=params["taxon_id"]
        obs=dict(
            taxon=dict(id=taxon_id,name=f"Fake{taxon_id} species",rank="species"),
            quality_grade="research",obscured=False,geoprivacy=None,
            geojson={"coordinates":[140.,35.]},positional_accuracy=100,
            observed_on="2026-05-01",photos=[dict(
                id=333,license_code="cc-by",url="https://images.example.test/photo.jpg")])
        return {"total_results":1,"results":[obs]}
    raise ValueError("unexpected endpoint")


def fake_exact(end,params):
    if end=="taxa":
        return fake(end,params)
    raise RuntimeError("observations disabled for unit test")


def test_original_source_pairing_uses_both_geologies_without_colour_value():
    source_data=source()
    x,a=select_original_endemic_colour_pairs(*source_data,max_pairs=3)
    assert len(x)==6 and a["n_original_distinct_source_pair_eligible_genera"]==4
    assert a["n_selected_original_accepted_endemic_plant_names"]==6
    assert len({r["genus"] for r in x})==3
    assert {r["origin"] for r in x}=={"oceanic","continental"}
    copied=source_data[0].copy()
    copied["normalized_value"]="violet"
    y,b=select_original_endemic_colour_pairs(copied,*source_data[1:],max_pairs=3)
    assert x==y and a==b


def test_provider_taxon_record_exact_and_observed_flower_not_assumed():
    rec=query_exact_source_taxon("Gena oceanic_0",fake_exact)
    assert rec["provider_status"]=="provider_error"
    def no_api(end,params):
        return {"results":[]} if end=="taxa" else {}
    rec=query_exact_source_taxon("Gena oceanic_0",no_api)
    assert rec["provider_status"]=="no_unambiguous_exact_public_taxon"


def test_genuine_public_photo_requires_exact_species_and_allowed_license():
    obs=dict(taxon=dict(id=44,name="Gena alba",rank="species"),
        quality_grade="research",obscured=False,geoprivacy=None,
        geojson={"coordinates":[140.,35.]},positional_accuracy=100,
        observed_on="2026-05-01",photos=[dict(
            id=333,license_code="cc-by",url="https://images.example.test/photo.jpg")])
    assert photo_observation_qualifies(obs,44,"Gena alba")
    obs["photos"][0]["license_code"]="all-rights-reserved"
    assert not photo_observation_qualifies(obs,44,"Gena alba")
    obs["photos"][0]["license_code"]="cc-by"
    obs["obscured"]=True
    assert not photo_observation_qualifies(obs,44,"Gena alba")


def test_bounded_taxon_source_retry_not_called_from_unit_test():
    x,a=select_original_endemic_colour_pairs(*source(),max_pairs=2)
    def empty_exact(end,params):
        return {"results":[]} if end=="taxa" else {}
    out=original_photo_pair_pilot(x,a,client=empty_exact)
    assert out["provider_query_counts_by_geology"]["oceanic"][
        "source_selected_endemic_accepted_taxa"]==2
    assert out["n_shared_genera_with_eligible_public_exact_photo_on_both_geologies"]==0
    assert out["source_query_protocol"]["no_original_photo_pixels_downloaded"]
