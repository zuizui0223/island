import pandas as pd
import pytest
from fd.paired_source_genus_cv import paired_source_genus_gain

def frame(n=40):
    rows=[]
    for i in range(n):
        for sector,pp,pg in [("actinomorphic_deep",.2,.7),("tubular_bell_deep",.7,.2)]:
            rows.append(dict(island_id=f"I{i}",archipelago=f"A{i//2}",
                spatial_block=f"B{i//2}",first_source=f"S{i//2}",
                trait_provider=f"P{i%2}",sector=sector,n_target=10,
                n_positive=8,p_pooled=pp,p_genus=pg,
                target_species_heldout=True,archipelago_heldout=True,
                source_provider_heldout=True))
    return pd.DataFrame(rows)

def test_common_island_paired_gain():
    z=paired_source_genus_gain(frame(),draws=299)
    assert z["n_common_islands"]==40
    assert z["delta_gain_radial_minus_tubular"]==pytest.approx(.70)

def test_no_unpaired_substitution():
    d=frame()
    d=d.loc[~((d.island_id=="I0")&(d.sector=="tubular_bell_deep"))]
    z=paired_source_genus_gain(d,draws=299)
    assert z["n_common_islands"]==39

def test_insufficient_groups_not_effect_zero():
    d=frame();d["archipelago"]="A"
    z=paired_source_genus_gain(d,draws=299)
    assert z["status"]=="NOT_IDENTIFIABLE_PAIRED_INDEPENDENT_SUPPORT"
    assert "delta_gain_radial_minus_tubular" not in z

def test_unverified_holdout_rejected():
    d=frame();d.loc[0,"archipelago_heldout"]=False
    with pytest.raises(ValueError):
        paired_source_genus_gain(d,draws=299)

def test_endpoint_source_mismatch_rejected():
    d=frame();d.loc[0,"first_source"]="OTHER"
    with pytest.raises(ValueError):
        paired_source_genus_gain(d,draws=299)

def test_duplicate_sector_rejected():
    d=frame();d=pd.concat([d,d.iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError):
        paired_source_genus_gain(d,draws=299)
