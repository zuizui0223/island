"""Synthetic original GIFT genus-overlap count-matching and leakage tests."""
import numpy as np
import pytest

from fd_original_stage2_name_disjoint_flower_forecast import (
    paired_genus_matched_source_priors,
    GENUS_MATCHED_MAX_DIRECT_PER_GENUS,
    GENUS_MATCHED_MIN_SHARED_GENERA,
    GENUS_MATCHED_MIN_TOTAL_SPECIES,
    GENUS_MATCHED_COLUMNS,
    BALANCED_SOURCE_REPLICATE_SEEDS,
)
from fd.two_stage_floral_forecast import proper_losses


def synthetic_donors():
    counts_l=(30,8,6,5,5)
    counts_e=(6,6,6,6,6)
    local=set()
    external=set()
    states={}
    for i, g in enumerate(("Alpha", "Beta", "Gamma", "Delta", "Epsilon")):
        for j in range(counts_l[i]):
            name=f"{g} localsp{j:03d}"
            local.add(name);states[name]=g
        for j in range(counts_e[i]):
            name=f"{g} extsp{j:03d}"
            external.add(name);states[name]=g
    forbidden={"Alpha localsp000","Beta extsp000"}
    indices={g:i for i,g in enumerate(("Alpha","Beta","Gamma","Delta","Epsilon"))}
    return local,external,forbidden,states,indices


def test_genus_matching_removes_source_genus_composition_differences():
    local,external,forbidden,states,indices=synthetic_donors()
    assert len(local)!=len(external)
    matched=paired_genus_matched_source_priors(local,external,forbidden,states,indices)
    assert matched is not None
    priors,meta=matched
    assert set(priors)==set(BALANCED_SOURCE_REPLICATE_SEEDS)=={0,1,2}
    assert meta["n_common_direct_source_genera"]==GENUS_MATCHED_MIN_SHARED_GENERA==5
    assert meta["n_genus_matched_direct_source_names_per_side"]>=GENUS_MATCHED_MIN_TOTAL_SPECIES==20
    for seed,(pl,pe) in priors.items():
        np.testing.assert_array_equal(pl,pe)
        assert np.isclose(pl.sum(),1) and np.isclose(pe.sum(),1)
        assert np.isfinite(pl).all() and (pl>0).all()
    # Neither source has an inherent preferred category after equal
    # genus coverage and equal within-genus N; no fictitious "effect".
    results,n2=paired_genus_matched_source_priors(
        local,external,forbidden,states,indices)
    assert meta==n2
    for seed in priors:
        np.testing.assert_array_equal(priors[seed][0],results[seed][0])


def test_real_within_genus_difference_survives_matched_source_taxon_weights():
    local,external,forbidden,states,indices=synthetic_donors()
    # Only Alpha is different on the EXTERNAL donor side; the genus list
    # and donor numbers in each matched genus stay fixed.
    for name in external:
        if name.startswith("Alpha "):
            states[name]="Beta"
    priors,meta=paired_genus_matched_source_priors(
        local,external,forbidden,states,indices)
    for seed,(pl,pe) in priors.items():
        assert pl[indices["Alpha"]]>pe[indices["Alpha"]]
        assert pe[indices["Beta"]]>pl[indices["Beta"]]
        bl,ll=proper_losses(pl,[indices["Alpha"]]*5)
        be,le=proper_losses(pe,[indices["Alpha"]]*5)
        assert bl<be and ll<le


def test_full_accepted_name_blacklist_applied_before_matching_both_sources():
    local,external,forbidden,states,indices=synthetic_donors()
    m,_=paired_genus_matched_source_priors(local,external,forbidden,states,indices)
    changed=dict(states)
    for s in forbidden:
        changed[s]="Gamma" if states[s]!="Gamma" else "Epsilon"
    new,_=paired_genus_matched_source_priors(local,external,forbidden,changed,indices)
    for seed in m:
        np.testing.assert_array_equal(m[seed][0],new[seed][0])
        np.testing.assert_array_equal(m[seed][1],new[seed][1])
    with pytest.raises(ValueError,match="ambiguous"):
        paired_genus_matched_source_priors(
            local|{"Beta extsp000"},external,forbidden,states,indices)


def test_unsupported_source_overlap_not_inferred_as_biological_negative():
    local,external,forbidden,states,indices=synthetic_donors()
    only_three_genus={x for x in external if x.split()[0] in {"Alpha","Beta","Gamma"}}
    assert paired_genus_matched_source_priors(
        local,only_three_genus,forbidden,states,indices) is None
    a,b={"Alpha local_a"}, {"Alpha ext_a"}
    one={"Alpha local_a":"Alpha","Alpha ext_a":"Beta"}
    assert paired_genus_matched_source_priors(
        a,b,{"Alpha focal"},one,{"Alpha":0,"Beta":1}) is None


def test_full_source_three_seed_two_loss_columns_immutable():
    for seed in BALANCED_SOURCE_REPLICATE_SEEDS:
        for direction in ("genusmatched_external","genusmatched_local"):
            prefix=f"{direction}_seed{seed}_gain_vs_"
            assert any(k.startswith(prefix) and k.endswith("_Brier")
                       for k in GENUS_MATCHED_COLUMNS)
            assert any(k.startswith(prefix) and k.endswith("_logloss")
                       for k in GENUS_MATCHED_COLUMNS)
        assert (f"genusmatched_minus_equalN_external_seed{seed}_Brier"
                in GENUS_MATCHED_COLUMNS)
