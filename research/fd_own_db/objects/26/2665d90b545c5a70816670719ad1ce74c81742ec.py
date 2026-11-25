"""Exact public location vs island-centroid radius is NOT island polygon presence."""
from fd.oceanic_endemic_near_island_photo import (
    source_island_centroid_radius,
    public_near_observation_source_gate,
    original_near_island_photo_provider_pilot,
)


def test_source_centroid_radius_clipped_from_area():
    assert source_island_centroid_radius(1)==5.
    assert source_island_centroid_radius(1000000)==75.


def test_public_open_and_taxon_privacy_required_even_if_near():
    sample=dict(
        taxon={"id":44,"name":"Species alba","rank":"species"},
        quality_grade="research",obscured=False,geoprivacy=None,
        geojson={"coordinates":[140.,20.]},positional_accuracy=100,
        observed_on="2026-04-01",taxon_geoprivacy=None,
        photos=[dict(id=333,license_code="cc-by",url="https://i.example/photo.jpg")])
    assert public_near_observation_source_gate(sample,44,"Species alba",
        latitude=20.,longitude=140.,radius_km=5.)
    sample["taxon_geoprivacy"]="obscured"
    assert not public_near_observation_source_gate(sample,44,"Species alba",
        latitude=20.,longitude=140.,radius_km=5.)
    sample["taxon_geoprivacy"]=None
    assert not public_near_observation_source_gate(sample,44,"Species alba",
        latitude=25.,longitude=140.,radius_km=5.)


def test_private_original_thirteen_genus_pilot_mock_provider():
    names={}
    selected=[]
    for k in range(13):
        for typ in ("oceanic","continental"):
            name=f"Genus{k} {typ}"
            selected.append(dict(
                plant=name,genus=f"Genus{k}",origin=typ,
                spatial_block=f"block{k}_{typ}",
                island_id=f"island{k}_{typ}",
                island_area_km2=50.,
                island_latitude=20.,
                island_longitude=140.,
            ))
    def mock(endpoint,params):
        if endpoint=="taxa":
            n=len(names)+100
            names[n]=params["q"]
            return {"results":[dict(id=n,name=params["q"],
                rank="species",is_active=True)]}
        tid=params["taxon_id"]
        observation=dict(
            taxon=dict(id=tid,name=names[tid],rank="species"),
            quality_grade="research",obscured=False,geoprivacy=None,
            taxon_geoprivacy=None,
            geojson={"coordinates":[140.,20.]},
            positional_accuracy=50,observed_on="2026-04-01",
            photos=[dict(id=222,license_code="cc-by",
                url="https://i.example/photo.jpg")],
        )
        return {"results":[observation],"total_results":1}
    result=original_near_island_photo_provider_pilot(
        selected,dict(n_selected_original_distinct_genera=13),client=mock)
    assert result["n_matched_genera_with_centroid_near_public_photo_in_both_geologies"]==13
    assert result["original_public_geo_photo_provider_by_island_type"][
        "oceanic"]["n_public_near_island_centroid_qualified_photo_taxa"]==13
    assert result["provider_query_safeguards"]["centroid_circle_not_evidence_of_island_presence"]
