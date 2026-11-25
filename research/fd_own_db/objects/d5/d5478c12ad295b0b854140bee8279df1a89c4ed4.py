"""Verify locality/source photo AND gating and taxon-level privacy."""
import pandas as pd

from fd.oceanic_endemic_vs_GIFT_photo_locality_join import (
    original_endemic_centroid_vs_GIFT_mainland_global_photo,
)


def input_sample():
    pairs=[]
    for i in range(13):
        for typ in ("oceanic","continental"):
            pairs.append(dict(
                genus=f"G{i}",origin=typ,
                plant=f"G{i} {typ}",botanical_source_congener=f"G{i} mainland",
                island_id=f"island_{i}_{typ}",
                spatial_block=f"geo_{i}_{typ}",archipelago=f"arch_{i}_{typ}",
                island_latitude=20.,island_longitude=140.,island_area_km2=75.,
                gift_rank=1,
            ))
    return pairs,dict(n_original_same_genus_pairings=26)


def test_both_local_island_source_and_congener_global_photo():
    pairs,audit=input_sample()
    name_ids={}
    reverse={}
    def fake(endpoint,params):
        if endpoint=="taxa":
            name=params["q"]
            if name not in name_ids:
                nid=len(name_ids)+100
                name_ids[name]=nid; reverse[nid]=name
            return {"results":[dict(name=name,id=name_ids[name],
                                    rank="species",is_active=True)]}
        if endpoint=="observations":
            taxon_id=int(params["taxon_id"])
            name=reverse[taxon_id]
            genus=int(name.split()[0][1:])
            if "lat" in params and name.endswith("oceanic"):
                permitted=genus<3
            elif "lat" in params and name.endswith("continental"):
                permitted=genus<2
            else:
                permitted=genus<12
            photos=([dict(id=123456,license_code="cc-by",
                          url="https://example.org/accepted.jpg")]
                    if permitted else [])
            return {"results":[dict(
                taxon=dict(id=taxon_id,name=name,rank="species"),
                photos=photos,geojson={"coordinates":[140.,20.]},
                positional_accuracy=100,observed_on="2026-06-01",
                quality_grade="research",obscured=False,
                geoprivacy=None,taxon_geoprivacy=None,
            )],"total_results":1}
        raise AssertionError("Unrecognized original provider endpoint")
    out=original_endemic_centroid_vs_GIFT_mainland_global_photo(
        pairs,audit,client=fake)
    assert out["n_original_different_accepted_GIFT_source_congener_species_queried"]==13
    a=out["support_by_geological_island_class"]
    assert a["oceanic"][
        "n_pairs_with_centroid_near_endemic_photo_AND_global_GIFT_congener_photo"]==3
    assert a["continental"][
        "n_pairs_with_centroid_near_endemic_photo_AND_global_GIFT_congener_photo"]==2
    assert out[
        "n_shared_genera_with_centroid_near_photo_and_source_congener_global_photo_in_BOTH_classes"]==2
    assert out["qualification"]["source_congener_photo_GLOBALLY_not_proven_mainland_local"]
    assert out["qualification"]["no_flower_pixels_or_image_colour_ROI_measured"]


def test_exact_original_geo_support_cannot_be_dropped():
    pairs,audit=input_sample()
    del pairs[0]["island_area_km2"]
    try:
        original_endemic_centroid_vs_GIFT_mainland_global_photo(
            pairs,audit,client=lambda *args,**kwargs:{"results":[]})
    except ValueError as e:
        assert "geographic data missing" in str(e)
    else:
        raise AssertionError("Missing island GSHHG area should fail closed")
