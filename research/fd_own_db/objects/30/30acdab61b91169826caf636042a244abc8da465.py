"""GIFT source congeners are descriptive controls, not evolutionary ancestors."""
import pandas as pd
import numpy as np
import pytest
from fd.geological_endemic_trait_divergence import (
    plant_source_discordance, geological_endemic_state_discordance,
    _status_contrast_within_identical_genus
)

def test_source_congener_frequency_not_favorable_match():
    out=plant_source_discordance({"Genus foo"},
        {"Genus bar":"red","Genus baz":"blue","Genus qux":"blue"},
        {"Genus foo":"red"})
    assert len(out)==1
    assert out[0]["source_genus_congener_count"]==3
    assert np.isclose(out[0]["discordance"],2/3)


def synthetic():
    geo=[];flora=[];nearest=[];donor=[]
    base=["Agen","Bgen","Cgen"]
    names=[]
    for origin in ("oceanic","continental"):
        for i in range(2):
            island=f"{origin}_{i}"
            geo.append(dict(island_id=island,geological_origin=origin,
                analysis_regime="tropical",spatial_block=f"{origin}_block_{i}"))
            for genus in base:
                for status in ("endemic","native_nonendemic"):
                    plant=f"{genus} {status}"
                    if plant not in names:names.append(plant)
                    flora.append(dict(island_id=island,accepted_species=plant,
                         origin_status="native",
                         endemic_status="endemic" if status=="endemic" else "nonendemic",
                         floristic_status=status))
            for rank in range(1,6):
                nearest.append(dict(island_id=island,source_mode="geo_k5",
                   source_rank=str(rank),entity_ID=f"s{rank}"))
    for rank in range(1,6):
        for genus in base:
            donor.append(dict(entity_ID=f"s{rank}",work_species=f"{genus} mainland"))
    x=[]
    for nm in names+[f"{g} mainland" for g in base]:
        x.append(dict(accepted_species=nm,trait_name="flower_primary_color",
            normalized_value="red" if nm.endswith(" endemic") else "blue",
            evidence_scope="species_direct",evidence_quality="high"))
    return pd.DataFrame(x),pd.DataFrame(flora),pd.DataFrame(geo),pd.DataFrame(nearest),pd.DataFrame(donor)


def test_provenance_and_paired_same_island_support_and_small_n_suppression():
    r=geological_endemic_state_discordance(*synthetic(),boots=999)
    assert r["n_geology_typed_islands"]==4
    assert r["n_original_islands_with_5_source_region_trait_eligibility"]==4
    assert r["source_provenance"]["original_endemic_status_separate_from_native"]
    for origin in ("continental","oceanic"):
        z=next(x for x in r["endemic_vs_native_same_island_congener_support"]
            if x["geological_origin"]==origin and x["trait"]=="flower_primary_color")
        assert z["n_islands_with_both_ge3_congener_comparisons"]==2
        assert z["n_islands_with_at_least_two_same_genus_endemic_nonendemic_pairs"]==2
        assert z["n_common_endemic_nonendemic_genus_units_in_this_group"]==6
        assert not z["same_genus_endemic_native_comparison_reportable"]
        assert not z["paired_status_comparison_reportable"]
    x=next(x for x in r["endemic_vs_native_same_island_congener_support"]
        if x["geological_origin"]=="oceanic" and x["trait"]=="flower_primary_color")
    assert x["n_status_variable_endemic_nonendemic_genus_island_units"]==6
    assert x["n_genus_island_units_with_nonzero_source_discordance"]==6
    assert x["n_islands_with_any_within_genus_status_flower_variation"]==2
    assert all(not x["origin_comparison_identifiable_in_supported_cohort"]
               for x in r["geological_origin_difference_of_endemic_status_effect"])


def test_endemic_native_origin_conflict_refuses_source():
    args=list(synthetic())
    args[1].loc[0,"origin_status"]="introduced"
    with pytest.raises(ValueError,match="Endemic source"):
        geological_endemic_state_discordance(*args,boots=999)

def test_categorical_zero_needs_focal_variation_support_audit():
    e=[dict(accepted_species="A endemic",discordance=0.2,_private_original_state="radial"),
       dict(accepted_species="B endemic",discordance=0.7,_private_original_state="bilateral")]
    n=[dict(accepted_species="A native",discordance=0.2,_private_original_state="radial"),
       dict(accepted_species="B native",discordance=0.7,_private_original_state="bilateral")]
    number,varies,nonzero,delta=_status_contrast_within_identical_genus(e,n)
    assert (number,varies,nonzero,delta)==(2,0,0,0.0)
    n[1]["_private_original_state"]="radial"
    # If donor reference has equal frequencies, status state variation need
    # not produce a nonzero discordance: separate those two audit concepts.
    number,varies,nonzero,delta=_status_contrast_within_identical_genus(e,n)
    assert number==2 and varies==1 and nonzero==0 and delta==0.0
