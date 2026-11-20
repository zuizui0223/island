"""Synthetic only: test original GIFT source balanced categorical reference control."""
import numpy as np
import pytest

from fd_original_stage2_name_disjoint_flower_forecast import (
    BALANCED_SOURCE_REPLICATE_SEEDS,
    BALANCED_SOURCE_GAIN_COLUMNS,
    BALANCED_SOURCE_MAX_DIRECT_NAMES,
    balanced_source_prior,
    paired_balanced_source_priors,
)
from fd.two_stage_floral_forecast import proper_losses


def sample_original():
    # Fully synthetic botanical accepted names and true categorical Direct states.
    local = {f"Locgen species{i:04d}" for i in range(42)}
    external = {f"Extgen species{i:04d}" for i in range(31)}
    states = {name:("open" if i%3 else "deep")
              for i,name in enumerate(sorted(local|external))}
    forbidden = {"Locgen species0001", "Extgen species0001"}
    return local,external,states,forbidden


def test_same_actual_accepted_name_count_and_frozen_donor_sample_repeats():
    a,b,states,forbidden=sample_original()
    idx={"open":0,"deep":1}
    models,n=paired_balanced_source_priors(a,b,forbidden,states,idx)
    assert n==30
    assert set(models)==set(BALANCED_SOURCE_REPLICATE_SEEDS)=={0,1,2}
    for seed,(local,ext) in models.items():
        assert local.shape==(2,) and ext.shape==(2,)
        assert (local>0).all() and (ext>0).all()
        np.testing.assert_allclose(local.sum(),1)
        np.testing.assert_allclose(ext.sum(),1)
        local_2=balanced_source_prior(
            a-forbidden,states,idx,n_names=n,seed=seed)
        np.testing.assert_array_equal(local,local_2)
    again,n2=paired_balanced_source_priors(a,b,forbidden,states,idx)
    assert n2==n
    for seed in models:
        np.testing.assert_array_equal(models[seed][0],again[seed][0])
        np.testing.assert_array_equal(models[seed][1],again[seed][1])


def test_blocked_species_recurrence_and_region_overlap_not_leaking():
    a,b,states,forbidden=sample_original()
    idx={"open":0,"deep":1}
    models,n=paired_balanced_source_priors(a,b,forbidden,states,idx)
    assert n==30
    # Matched donor trait state of forbidden accepted species cannot matter.
    altered=dict(states)
    for name in forbidden:
        altered[name]="deep" if states[name]=="open" else "open"
    models2,n2=paired_balanced_source_priors(a,b,forbidden,altered,idx)
    assert n2==n
    for seed in models:
        np.testing.assert_array_equal(models[seed][0],models2[seed][0])
        np.testing.assert_array_equal(models[seed][1],models2[seed][1])
    with pytest.raises(ValueError,match="region exclusive"):
        paired_balanced_source_priors(a|{"Extgen species0000"},b,
                                     forbidden,states,idx)


def test_missing_source_and_same_name_balance_boundaries():
    a,b,states,forbidden=sample_original()
    with pytest.raises(ValueError,match=">=20"):
        balanced_source_prior(a-forbidden,states,{"open":0,"deep":1},
                              n_names=19,seed=0)
    with pytest.raises(ValueError,match=">=20"):
        balanced_source_prior(a-forbidden,states,{"open":0,"deep":1},
                              n_names=100,seed=0)
    with pytest.raises(ValueError,match=">=20"):
        balanced_source_prior(a-forbidden,states,{"open":0,"deep":1},
                              n_names=BALANCED_SOURCE_MAX_DIRECT_NAMES+1,seed=0)
    with pytest.raises(ValueError,match=">=20"):
        paired_balanced_source_priors(a,set(sorted(b)[:10]),forbidden,
                                     states,{"open":0,"deep":1})


def test_balanced_reference_delta_sign_is_exactly_paired_on_same_real_codes():
    a,b,states,forbidden=sample_original()
    refs,n=paired_balanced_source_priors(
        a,b,forbidden,states,{"open":0,"deep":1})
    codes=[0,0,0,1,0,1]
    for seed,(local,external) in refs.items():
        lb,ll=proper_losses(local,codes)
        eb,el=proper_losses(external,codes)
        assert np.isfinite([lb,ll,eb,el]).all()
        assert (lb-eb)==pytest.approx(-(eb-lb))
        assert (ll-el)==pytest.approx(-(el-ll))
        for sign in ("balanced_external_seed"+str(seed),
                     "balanced_local_seed"+str(seed)):
            assert any(c.startswith(sign) and c.endswith("_Brier")
                       for c in BALANCED_SOURCE_GAIN_COLUMNS)
