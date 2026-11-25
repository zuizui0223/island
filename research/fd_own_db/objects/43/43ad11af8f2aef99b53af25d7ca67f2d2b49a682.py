import pandas as pd
import pytest
from fd.independent_reassembly_cv import evaluate

def cohort(n=40):
    return pd.DataFrame({
        "island_id":[f"I{i}" for i in range(n)],
        "archipelago":[f"A{i//2}" for i in range(n)],
        "spatial_block":[f"B{i//2}" for i in range(n)],
        "first_source":[f"S{i//2}" for i in range(n)],
        "trait_provider":[f"P{i%2}" for i in range(n)],
        "n_target":[10]*n,"n_positive":[8]*n,
        "p_baseline":[.2]*n,"p_genus":[.7]*n,
        "p_guild":[.8]*n,"guild_supported":[True]*n,
        "target_species_heldout":[True]*n,"archipelago_heldout":[True]*n,
        "source_provider_heldout":[True]*n})

def test_exact_brier_and_grouped_cv():
    z=evaluate(cohort(),draws=299)
    assert z["source_genus"]["gain_genus"]==pytest.approx(.35)
    assert z["guild_subset"]["gain_guild"]==pytest.approx(.01)

def test_guild_nondetection_not_zero_effect():
    d=cohort();d["guild_supported"]=False;d["p_guild"]=float("nan")
    z=evaluate(d,draws=299)
    assert z["guild_subset"]["status"]=="NOT_IDENTIFIABLE_INDEPENDENT_SUPPORT"
    assert "gain_guild" not in z["guild_subset"]

def test_group_pseudoreplication():
    d=cohort();d["archipelago"]="one"
    assert evaluate(d,draws=299)["source_genus"]["status"]=="NOT_IDENTIFIABLE_INDEPENDENT_SUPPORT"

@pytest.mark.parametrize("col,val",[
    ("p_genus",1.2),("guild_supported","nondetected"),
    ("n_positive",11),("archipelago",None)])
def test_invalid_inputs(col,val):
    d=cohort();d[col]=val
    with pytest.raises(ValueError):
        evaluate(d,draws=299)

def test_duplicate_islands():
    d=cohort();d.loc[1,"island_id"]=d.loc[0,"island_id"]
    with pytest.raises(ValueError):
        evaluate(d,draws=299)

@pytest.mark.parametrize("flag",["target_species_heldout","archipelago_heldout","source_provider_heldout"])
def test_unproven_out_of_fold_is_rejected(flag):
    d=cohort();d.loc[0,flag]=False
    with pytest.raises(ValueError):
        evaluate(d,draws=299)

def test_unsupported_guild_predictions_fail_closed():
    d=cohort();d["guild_supported"]=False
    with pytest.raises(ValueError):
        evaluate(d,draws=299)

def test_partly_supported_guild_no_fabricated_absences():
    d=cohort()
    d.loc[0:19,"guild_supported"]=False
    d.loc[0:19,"p_guild"]=float("nan")
    z=evaluate(d,draws=299)
    assert z["source_genus"]["status"]=="DESCRIPTIVE_HELDOUT_PREDICTIVE_COMPARISON"
    assert z["guild_subset"]["status"]=="NOT_IDENTIFIABLE_INDEPENDENT_SUPPORT"
