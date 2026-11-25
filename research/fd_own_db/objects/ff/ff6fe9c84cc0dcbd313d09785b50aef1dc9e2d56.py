import pandas as pd
import pytest
from fd.functional_support_null import validate_counts, score_pairs, safe_aggregate


def row(scope="All",region="tropical",a="I1",b="I2",trait="generalized_form"):
    x={"scope":scope,"region":region,"island_a":a,"island_b":b,"trait":trait}
    for stage,n in (("all",100),("species",80),("genus",60)):
        for side,k in (("a",40),("b",30 if stage=="all" else 25)):
            x[f"{stage}_n_{side}"]=n
            x[f"{stage}_k_{side}"]=k
    x["genus_k_b"]=5
    return x


def test_missing_and_non_nested_input_fails():
    x=row();x["genus_n_a"]=101
    with pytest.raises(ValueError,match="Non-nested support"):
        validate_counts(pd.DataFrame([x]))
    x=row();x["genus_k_a"]=61
    with pytest.raises(ValueError,match="exceeds"):
        validate_counts(pd.DataFrame([x]))
    x=row();x["genus_k_a"]=41
    with pytest.raises(ValueError,match="Non-nested positive"):
        validate_counts(pd.DataFrame([x]))
    x=row();x["trait"]="made_up"
    with pytest.raises(ValueError,match="Unknown"):
        validate_counts(pd.DataFrame([x]))


def test_sequential_genus_exclusion_exceeds_random_support_loss():
    s=score_pairs(pd.DataFrame([row()]),draws=399,seed=11)
    assert s.extra_genus_gap_beyond_support.iloc[0]>.20
    pd.testing.assert_frame_equal(s,score_pairs(pd.DataFrame([row()]),draws=399,seed=11))


def test_no_extra_genus_removal_is_exact_null():
    x=row()
    for side in "ab":
        x[f"genus_n_{side}"]=x[f"species_n_{side}"]
        x[f"genus_k_{side}"]=x[f"species_k_{side}"]
    s=score_pairs(pd.DataFrame([x]),draws=399,seed=9)
    assert abs(s.extra_genus_gap_beyond_support.iloc[0])<1e-12


def test_aggregate_never_emits_private_identifiers():
    x=[]
    for i in range(20):
        for trait in ("self_compatibility","autonomous_selfing"):
            x.append(row(a=f"A{i}",b=f"B{i}",trait=trait))
    scored=score_pairs(pd.DataFrame(x),draws=399,seed=2)
    out=safe_aggregate(scored)
    assert len(out)==4  # two atomic + regional domain + pooled domain
    assert not set(("island_a","island_b","taxon_name","genus"))&set(out.columns)
    assert (out.n_pairs==20).all()
    assert "reproductive_assurance" in set(out.measure)
    with pytest.raises(ValueError,match="Privacy"):
        safe_aggregate(scored,min_pairs=5)
