"""Exact native + global single-trait rather than triple-complete floral diversity.

Synthetic tests contain no confidential original island botanical evidence.
The original full 8k restored source runner asserts the 928 historic
complete-triple subset only as a backward-compatibility audit.
"""
import numpy as np
import pandas as pd
import pytest

from fd.global_partial_trait_flora_fd import (
    TRAITS,PAIR_TRAITS,strict_single_trait_original,
    single_trait_fd,pair_trait_fd,original_global_partial_traits
)


def test_independent_flower_colour_only_no_form_symmetry_required():
    names=["A alba","B alba","C alba","D alba","E alba",
           "F alba","G alba"]
    colour={p:("white" if i<3 else "purple") for i,p in enumerate(names[:5])}
    z=single_trait_fd(colour,set(names))
    assert z["n_direct_species"]==5
    # P(two different distinct plant names have different colours)=3*2*2/(5*4)=0.6
    assert np.isclose(z["finite_5_species_rao"],0.8*0.6)
    assert z["n_distinct_recorded_genera"]==5
    assert np.isclose(z["finite_5_genus_equal_rao"],z["finite_5_species_rao"])
    assert np.isnan(single_trait_fd({"A alba":"white"},set(names))[
        "finite_5_species_rao"])


def test_unbalanced_species_versus_genus_weights_are_separate():
    categories={"A one":"white","A two":"white","A three":"white",
                "A four":"white","B one":"black","C one":"black"}
    r=single_trait_fd(categories,set(categories))
    assert r["n_direct_species"]==6 and r["n_distinct_recorded_genera"]==3
    assert np.isnan(r["finite_5_genus_equal_rao"])
    assert np.isfinite(r["finite_5_species_rao"])


def test_two_trait_coupling_uses_pair_not_third_trait():
    names=[f"G{k} sp" for k in range(5)]
    color={n:("white" if i<=2 else "blue") for i,n in enumerate(names)}
    form={n:("tubular" if i<=2 else "open_radial")
          for i,n in enumerate(names)}
    x=pair_trait_fd(color,form,set(names))
    assert x["n_two_trait_original_species"]==5
    assert x["finite_5_joint_rao"]>0
    assert x["finite_5_pairing_departure"]<0
    assert np.isnan(pair_trait_fd(color,{"G0 sp":"tubular"},set(names))[
        "finite_5_pairing_departure"])


def test_reject_ambiguous_direct_and_preserve_other_trait():
    d=pd.DataFrame([
       dict(accepted_species="Alpha a",trait_name="flower_primary_color",
            normalized_value="white",evidence_scope="species_direct",
            evidence_quality="high"),
       dict(accepted_species="Beta a",trait_name="floral_form",
            normalized_value="open_radial",evidence_scope="species_direct",
            evidence_quality="medium"),
       dict(accepted_species="Gamma a",trait_name="floral_symmetry",
            normalized_value="zygomorphic|actinomorphic",evidence_scope="species_direct",
            evidence_quality="high"),
    ])
    r=strict_single_trait_original(d)
    assert r["flower_primary_color"]=={"Alpha a":"white"}
    assert r["floral_form"]=={"Beta a":"open_radial"}
    assert r["floral_symmetry"]=={}
    with pytest.raises(ValueError):
        strict_single_trait_original(pd.concat([d,d.iloc[:1]],ignore_index=True))


def _global_source():
    # 4000 islands, 928 with >=30 colour×form×symmetry Direct recorded
    # original real-world gate; source 2,3,4 traits can have additional
    # incomplete observations without gaining complete-case eligibility.
    names=[f"G{k%37} plant{k}" for k in range(68)]
    fl=[]
    for k in range(4000):
        island=f"i{k:04}"
        roster=names[:35] if k<928 else names[35:42]
        for n in roster:
            fl.append(dict(island_id=island,accepted_species=n,
                origin_status="native",endemic_status="nonendemic",
                floristic_status="native_nonendemic"))
    flora=pd.DataFrame(fl)
    # First 35 species have triples, next 7 colour only, so individual
    # trait coverage is wider than three-trait completion.
    direct=[]
    for k,n in enumerate(names[:42]):
        for trait in (TRAITS[:3] if k<35 else (TRAITS[0],)):
            direct.append(dict(accepted_species=n,trait_name=trait,
                normalized_value=("white" if k%2 else "violet") if trait==TRAITS[0]
                    else ("open_radial" if k%2 else "tubular") if trait==TRAITS[1]
                    else "actinomorphic",
                evidence_scope="species_direct",evidence_quality="high"))
    geo=[]
    assignments=[]
    for k in range(4000):
        island=f"i{k:04}"
        if k<440:
            origin="oceanic" if k%2 else "continental"
        else:
            origin="unresolved"
        geo.append(dict(island_id=island,geological_origin=origin,
            analysis_regime="tropical",
            spatial_block=f"block{k%33}",LGM_connectivity="lgm_connected"))
        if k<440:
            for rank in range(1,6):
                assignments.append(dict(island_id=island,
                    source_mode="geo_k5",source_rank=str(rank),entity_ID=f"s{rank}"))
    mainland=pd.DataFrame([
        dict(entity_ID=f"s{k}",work_species=n)
        for k in range(1,6) for n in names[40:42]
    ])
    return [pd.DataFrame(direct),flora,pd.DataFrame(geo),
            pd.DataFrame(assignments),mainland]


def test_full_4000_islands_main_trait_gate_not_full_three_evidence():
    data=_global_source()
    result=original_global_partial_traits(*data,seed=20261114,bootstrap_draws=999)
    assert result["n_original_global_distinct_island_floras"]==4000
    assert result["n_original_three_trait_complete30_islands_retained_for_SENSITIVITY"]==928
    assert result["n_original_distinct_islands_with_geological_origin_typed"]==440
    assert result["n_original_partial_island_scope_endpoint_rows"]==4000*2*13
    records=[r for r in result["marginal_trait_and_pairwise_coverage"]
        if r["region"]=="tropical" and r["geological_origin"]=="unresolved"
        and r["source_plant_status_scope"]=="original_native"]
    one=next(x for x in records if x["trait_name"]=="flower_primary_color")
    pair=next(x for x in records if x["trait_name"]=="flower_primary_color__floral_form")
    assert one["n_islands_with_at_least_five_direct_species"]==3560
    # Of the 3560 unresolved-geology islands, 488 have all three traits,
    # but 3072 have colour ONLY; pairwise cannot invent their missing form.
    assert pair["n_islands_with_at_least_five_direct_species"]==488
    missing=[x for x in result["oceanic_minus_continental_conditional_traits"]
             if x["metric"]=="observed_direct_coverage"]
    assert missing
    assert result["source_design"]["original_direct_trait_observation_fraction_by_island_and_scope"]
    assert result["source_design"]["native_status_scope_without_original_native_species_is_NA_not_zero"]
    assert result["source_design"]["original_island_complete_case_gate_removed_for_main_analysis"]
