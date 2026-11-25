"""Original island-endemic/GIFT-mainland candidate photo support: source-only tests."""
import pandas as pd

from fd.oceanic_endemic_mainland_congener_photo import (
    source_independent_mainland_controls,
    original_mainland_congener_public_photo_gate,
)


def _source():
    selected=[]
    flora=[]
    sources=[]
    mainland=[]
    ledger=[]
    for g in range(13):
        genus=f"Genus{g}"
        accepted=f"{genus} mainland"
        ledger.append(dict(accepted_species=accepted,trait_name="flower_primary_color",
            normalized_value="yellow",evidence_scope="species_direct",
            evidence_quality="high"))
        for rank in range(1,6):
            mainland.append(dict(entity_ID=f"S{rank}",work_species=accepted))
        for typ in ("oceanic","continental"):
            plant=f"{genus} {typ}"
            island=f"{typ}-{g}"
            selected.append(dict(genus=genus,origin=typ,plant=plant,
                                 island_id=island))
            flora.append(dict(island_id=island,accepted_species=plant))
            ledger.append(dict(accepted_species=plant,trait_name="flower_primary_color",
                normalized_value="red" if typ=="oceanic" else "white",
                evidence_scope="species_direct",evidence_quality="high"))
            for rank in range(1,6):
                sources.append(dict(island_id=island,source_mode="geo_k5",
                                    entity_ID=f"S{rank}",source_rank=rank))
    return selected,pd.DataFrame(ledger),pd.DataFrame(flora),pd.DataFrame(sources),pd.DataFrame(mainland)


def test_26_focal_endemics_preselected_13_source_genera_and_rank():
    s,l,f,a,m=_source()
    pairs,audit=source_independent_mainland_controls(s,l,f,a,m)
    assert len(pairs)==26
    assert audit["n_distinct_mainland_candidate_species"]==13
    assert all(x["gift_rank"]==1 and x["botanical_source_congener"].endswith(
        "mainland") for x in pairs)
    assert audit["original_selection_source_only_pre_image_pre_colour"]
    scrambled=l.copy()
    scrambled["normalized_value"]="purple"  # no trait outcome used in source choice
    p2,a2=source_independent_mainland_controls(s,scrambled,f,a,m)
    assert p2==pairs and a2==audit


def test_photo_provider_results_not_used_to_reselect_source():
    selected,l,flora,assign,mainland=_source()
    paired,audit=source_independent_mainland_controls(
        selected,l,flora,assign,mainland)
    key={}
    def fake(endpoint,params):
        if endpoint=="taxa":
            n=str(params["q"])
            if n not in key:key[n]=len(key)+101
            return {"results":[{"name":n,"rank":"species",
                                "id":key[n],"is_active":True}]}
        if endpoint=="observations":
            tid=int(params["taxon_id"])
            name=next(k for k,v in key.items() if v==tid)
            qualifies=("mainland" in name or name.endswith("oceanic"))
            photo=[{"id":12345,"license_code":"cc-by",
                    "url":"https://example.test/flower.jpg"}] if qualifies else []
            return {"total_results":1,"results":[{
                "taxon":{"id":tid,"name":name,"rank":"species"},
                "quality_grade":"research","obscured":False,
                "geoprivacy":None,"observed_on":"2026-05-01",
                "positional_accuracy":250,
                "geojson":{"coordinates":[138.0,35.0]},
                "photos":photo
            }]}
        raise AssertionError(endpoint)
    z=original_mainland_congener_public_photo_gate(
        paired,audit,client=fake)
    assert z["n_unique_source_mainland_accepted_congener_names_queried"]==13
    assert z["n_unique_original_endemic_accepted_names_queried"]==26
    x=z["public_photo_opportunity_by_geological_island_type"]
    assert x["oceanic"]["n_source_genus_pairs_with_endemic_and_mainland_global_photo"]==13
    assert x["continental"]["n_source_genus_pairs_with_endemic_and_mainland_global_photo"]==0
    assert x["continental"]["n_mainland_GIFT_source_congeners_with_global_public_photo"]==13
    assert z["source_scope"]["source_photo_is_GLOBAL_not_confirmed_mainland_locality"]
