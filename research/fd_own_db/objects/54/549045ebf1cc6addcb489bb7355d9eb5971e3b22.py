"""Synthetic taxonomy/metadata-safe cross-geology endemic source comparison."""
import pandas as pd
import numpy as np
import pytest

from fd.geological_endemic_shared_genus import (
    calculate_endemic_shared_genus_geology, _exact_signed_genus_sensitivity,
)


def example():
    names=[]
    geo=[]
    flora=[]
    source=[]
    assign=[]
    genera=["Aa","Bb","Cc","Dd","Ee","Ff"]
    for kind in ("oceanic","continental"):
        for island_nr in range(3):
            island=f"{kind}_{island_nr}"
            geo.append(dict(island_id=island,geological_origin=kind,
                            analysis_regime="tropical",spatial_block=island))
            for genus in genera:
                plant=f"{genus} {kind}_{island_nr}"
                names.append(plant)
                flora.append(dict(
                    island_id=island,accepted_species=plant,origin_status="native",
                    endemic_status="endemic",floristic_status="endemic"))
            for rank in range(1,6):
                assign.append(dict(island_id=island,entity_ID=f"s{rank}",
                                   source_mode="geo_k5",source_rank=str(rank)))
    for rank in range(1,6):
        for genus in genera:
            source.append(dict(entity_ID=f"s{rank}",work_species=f"{genus} mainland"))
    allnames=names+[f"{genus} mainland" for genus in genera]
    ledger=pd.DataFrame([dict(
        accepted_species=n,trait_name="flower_primary_color",
        normalized_value=("red" if "oceanic" in n else "yellow" if "continental" in n else "red"),
        evidence_scope="species_direct",evidence_quality="high")
        for n in allnames])
    return [ledger,pd.DataFrame(flora),pd.DataFrame(geo),
            pd.DataFrame(assign),pd.DataFrame(source)]


def test_shared_genus_opportunity_and_state_variability_are_separate():
    d=calculate_endemic_shared_genus_geology(*example(),draws=999)
    assert d["original_geological_islands_with_valid_labels"]==6
    row=next(x for x in d["original_shared_genus_geological_comparisons"]
             if x["region"]=="tropical" and x["trait"]=="flower_primary_color")
    assert row["n_same_genus_shared_across_geological_types"]==6
    assert row["n_shared_genera_with_different_original_endemic_trait_categories"]==6
    # Merely 3 blocks per geological type: fail original independence gate.
    assert row["n_oceanic_blocks_covered_by_shared_genera"]==3
    assert not row["qualified_to_compare_source_endemic_discordance_within_shared_genera"]
    assert "paired_genus_oceanic_minus_continental_discordance" not in row
    assert d["method"]["accepted_plant_species_not_island_records_are_taxon_units"]


def test_source_taxon_status_provenance_conflict_fails_closed():
    inp=example()
    inp[1].loc[0,"origin_status"]="introduced"
    with pytest.raises(ValueError,match="endemic plant native origin"):
        calculate_endemic_shared_genus_geology(*inp,draws=999)


def test_no_outcome_driven_decision_when_nonvariable():
    inp=example()
    inp[0]["normalized_value"]="red"
    r=calculate_endemic_shared_genus_geology(*inp,draws=999)
    x=next(z for z in r["original_shared_genus_geological_comparisons"]
           if z["region"]=="tropical" and z["trait"]=="flower_primary_color")
    assert x["n_shared_genera_with_different_original_endemic_trait_categories"]==0
    assert not x["qualified_to_compare_source_endemic_discordance_within_shared_genera"]


def test_exact_2_pow_n_signed_genus_and_loo_are_finite():
    z=_exact_signed_genus_sensitivity(np.array([1., 1., 1., 1., 1.]))
    assert z["n_exact_sign_reassignments"]==32
    assert z["exact_two_sided_genus_label_sign_flip_p"]==2/32
    assert z["leave_one_genus_out_effect_range"]==[1.0,1.0]
    assert z["leave_one_genus_out_same_direction_fraction"]==1.0
    zero=_exact_signed_genus_sensitivity(np.zeros(5))
    assert zero["exact_two_sided_genus_label_sign_flip_p"]==1.0
    assert zero["n_genus_pairs_with_nonzero_difference"]==0
