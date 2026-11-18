import pandas as pd
import numpy as np
import pytest

from island_v2.fd_globi_provider_family_transfer import (
    direct_morphology, family_crosswalk, family_transfer,
)

def inputs(n=72):
    flora=[]
    families=[]
    fixed=[]
    before=[]
    after=[]
    for i in range(n):
        genus=f"G{i:03d}"
        fam=f"Family{i//3}"
        open_rate=(i%3)/2
        families.append({"genus":genus,"family":fam})
        fixed.append({"genus":genus,"n_matched_plant_species":2+i%4})
        for j in range(4):
            for trait, state in [
                ("flower_primary_color","white"),
                ("floral_form","open_radial" if j < (i%3)+1 else "tubular"),
                ("floral_symmetry","actinomorphic"),
            ]:
                flora.append({
                    "accepted_species":f"{genus} species{j}",
                    "trait_name":trait,"normalized_value":state,
                    "evidence_scope":"species_direct","evidence_quality":"high",
                })
        before.append([12+int(10*open_rate),12,9,5+int(10*open_rate),8])
        after.append([4+int(5*open_rate),4,4,3+int(5*open_rate),3])
    bb=pd.DataFrame(before,columns=("bombus","non_bombus_bees","lepidoptera","flower_visiting_birds","diptera"))
    aa=pd.DataFrame(after,columns=bb.columns)
    idx=[r["genus"] for r in families]
    bb.index=idx;aa.index=idx
    return bb,aa,pd.DataFrame(flora),pd.DataFrame(fixed),pd.DataFrame(families)


def test_frozen_source_matched_genus_and_family_training():
    bb,aa,flora,fixed,family=inputs()
    out=family_transfer(bb,aa,flora,fixed,family,bootstrap=299)
    assert out["status"]=="COMPLETED_EXPLORATORY_PROVIDER_OMITTED_FAMILY_TRANSFER"
    assert out["n_same_genera_before_after_with_direct_traits_and_unambiguous_family"]==72
    assert out["n_same_families"]==24
    assert out["after"]["family_equal_open_form_gain"]>0
    assert out["pre"]["genus_equal_baseline_MSE"]>0
    assert not out["local_pollinator_presence_or_effectiveness_measured"]
    assert not out["original_taxon_or_provider_IDs_released"]
    assert "genus" not in out and "family" not in out
    assert out["before_after_comparison_is_same_genus_same_family"]
    assert out["bootstrap_is_paired_at_family_level"]
    assert out["family_equal_delta_ci025"] <= out["family_equal_after_minus_before_gain"]
    assert out["family_equal_delta_ci975"] >= out["family_equal_after_minus_before_gain"]
    assert set(out["per_guild_after_minus_before"]) == {
        "bombus", "non_bombus_bees", "lepidoptera",
        "flower_visiting_birds", "diptera",
    }
    assert all(q["family_equal_delta_ci025"] <=
               q["family_equal_after_minus_before_gain"] <=
               q["family_equal_delta_ci975"]
               for q in out["per_guild_after_minus_before"].values())
    assert set(out["pre"]["individual_guilds"]) == set(
        out["after"]["individual_guilds"])
    assert abs(out["genus_equal_after_minus_before_gain"] -
               (out["after"]["genus_equal_open_form_gain"] -
                out["pre"]["genus_equal_open_form_gain"])) < 1e-12


def test_loss_of_coverage_gates_are_not_negative_ecological_effect():
    before,after,flora,fixed,tax=inputs()
    after.iloc[12:]=0
    result=family_transfer(before,after,flora,fixed,tax,bootstrap=299)
    assert result["status"]=="NOT_IDENTIFIABLE_AFTER_PROVIDER_OMISSION"
    assert not result["heldout_original_provider_prediction_completed"]
    assert "after" not in result


def test_competing_original_family_labels_excluded():
    before,after,flora,fixed,tax=inputs()
    tax=pd.concat([tax,pd.DataFrame([{"genus":"G000","family":"OtherFamily"}])])
    result=family_transfer(before,after,flora,fixed,tax,bootstrap=299)
    assert result["n_taxonomy_conflicting_genera"]==1
    assert result["n_same_genera_before_after_with_direct_traits_and_unambiguous_family"]==71


def test_duplicate_species_direct_measurements_are_rejected():
    before,after,flora,fixed,tax=inputs()
    flora=pd.concat([flora,flora.iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError,match="Duplicated"):
        family_transfer(before,after,flora,fixed,tax,bootstrap=299)
