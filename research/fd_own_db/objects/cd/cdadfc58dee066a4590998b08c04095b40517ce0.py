"""Synthetic geological-origin observational source versus endemic proxy safeguards."""
import numpy as np
import pandas as pd

from fd.geological_establishment_adaptation import (
    candidate_endemic_divergence, geology_origin_flora_pilot
)


def test_endemic_distance_is_not_imputed_as_adaptation():
    names=["Ab genus","Ab endem","Cd native","Cd endem"]
    source={"Ab source","Cd source"}
    traits={
        "Ab genus":("red","tubular","actinomorphic"),
        "Ab endem":("white","tubular","zygomorphic"),
        "Cd native":("white","open_radial","zygomorphic"),
        "Cd endem":("red","tubular","actinomorphic"),
        "Ab source":("red","tubular","actinomorphic"),
        "Cd source":("white","open_radial","zygomorphic")
    }
    flora=pd.DataFrame({
        "accepted_species":names,
        "origin_status":["native"]*4,
        "endemic_status":["nonendemic","endemic","nonendemic","endemic"],
        "floristic_status":["native_nonendemic","endemic","native_nonendemic","endemic"],
    })
    z=candidate_endemic_divergence(flora,set(names),source,traits)
    assert z["n_endemic_original_direct"]==2
    assert z["n_endemic_source_genus_comparable"]==2
    assert np.isclose(z["endemic_source_congener_category_mismatch_mean"],2.5)
    assert np.isnan(z["endemic_minus_native_source_congener_mismatch"])


def test_928_geology_joins_never_promote_unknown_origin_or_local_native():
    names=[f"Genus{i%40} species{i}" for i in range(2234)]
    direct=pd.DataFrame([
        ("white","open_radial","actinomorphic")
        if i%2 else ("red","tubular","zygomorphic")
        for i in range(len(names))
    ],index=names,columns=[
        "flower_primary_color","floral_form","floral_symmetry"])
    island_ids=[f"i{k:04d}" for k in range(928)]
    flora=[]
    for i in range(928):
        for j,name in enumerate(names[:35]):
            endemic=(j<5)
            flora.append({
                "island_id":island_ids[i],"accepted_species":name,
                "origin_status":"native",
                "endemic_status":"endemic" if endemic else "nonendemic",
                "floristic_status":"endemic" if endemic else "native_nonendemic",
                "strict_source_native":True,
                "regional_native_compatible":True
            })
    classified=pd.DataFrame(flora)
    geo=pd.DataFrame([{"island_id":t,"analysis_regime":"tropical",
        "spatial_block":f"b{i%40}"} for i,t in enumerate(island_ids)])
    original_geo=[]
    for i,t in enumerate(island_ids):
        geological_origin=(
            "oceanic" if i<220 else
            "continental" if i<440 else "unresolved"
        )
        original_geo.append({
            "island_id":t,"geological_origin":geological_origin,
            "weigelt_archipelago":f"arch{i%80}",
            "geological_class_source":"Roeble_2024_510_island_geology",
            "LGM_connectivity":"lgm_unconnected" if i%2 else "lgm_connected",
            "area_km2":str(100+((i%7)*10)),
            "distance_to_continent_km":str(60+((i%7)*20)),
            **{f"climate_pc{k}":str((i%9)*0.1+k*0.02) for k in range(1,5)},
        })
    geological=pd.DataFrame(original_geo)
    first={f"s{k}":set(names[35:]) for k in range(35)}
    near={t:tuple(f"s{(i+j)%35}" for j in range(5))
          for i,t in enumerate(island_ids)}
    result=geology_origin_flora_pilot(
        direct,classified,geo,geological,near,first,seed=16,boots=1999)
    assert result["status"].startswith("GEOLOGICAL_OCEANIC_CONTINENTAL")
    assert result["n_original_all_928_Direct_island_floras"]==928
    assert result["n_geology_labelled_source_comparable_all_flora"]==440
    assert result["independent_archived_geology_label_counts"]["oceanic"]==220
    assert result["independent_archived_geology_label_counts"]["continental"]==220
    assert result["analysis"]["establishment_survival_not_observable_from_species_checklists"]
    assert result["analysis"]["endemic_source_congener_comparison_is_not_adaptation_or_sister_phylogeny"]
    assert all("island_id" not in q for q in result["groups"])
    assert all("source_signature" not in q for q in result["geological_origin_associations"])
    assert sum(x["origin_contrast_supported_sampling"]
               for x in result["geological_origin_associations"])>=1
