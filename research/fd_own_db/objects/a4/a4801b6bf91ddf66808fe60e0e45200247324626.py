"""Finite categorical volume must be sensitive to geometry, not just k profiles."""
from itertools import product

import numpy as np
import pandas as pd
import pytest

from fd.discrete_functional_volume import (
    TRAITS,FiniteProfileReference,MEASURES,_resample_group,
)


def _triple_species(n_repeats=10):
    combos=list(product(("white","purple","yellow"),
                        ("open_radial","tubular","composite_head"),
                        ("actinomorphic","zygomorphic")))
    assert len(combos)==18
    names=[f"synthetic_{j}" for j in range(n_repeats*len(combos))]
    frames=[combos[j%len(combos)] for j in range(len(names))]
    return pd.DataFrame(frames,index=names,columns=TRAITS)


def test_categorical_one_mismatch_volume_is_not_a_recorded_profile_count():
    x=_triple_species()
    space=FiniteProfileReference.from_direct_species(x,null_draws=199,seed=7)
    assert len(space.profiles)==18
    compact=np.asarray([
        ("white","open_radial","actinomorphic"),
        ("white","tubular","actinomorphic"),
        ("white","composite_head","actinomorphic")
    ])
    spread=np.asarray([
        ("white","open_radial","actinomorphic"),
        ("purple","tubular","zygomorphic"),
        ("yellow","composite_head","actinomorphic")
    ])
    a=space.coverage(compact)
    b=space.coverage(spread)
    assert a[4]==b[4]==3
    assert a[0]==b[0]==pytest.approx(3/18)
    assert b[1]!=pytest.approx(a[1])
    assert a[1]>a[0] and b[1]>b[0]
    assert a[2]>=a[1] and b[2]>=b[1]
    assert 0<=a[3]<=1


def test_profile_k_null_separate_from_colour_architecture_shuffle():
    x=_triple_species()
    a=FiniteProfileReference.from_direct_species(x,null_draws=199,seed=9)
    assert sorted(a.expected_volume_r1_given_k)==list(range(1,19))
    r=a.sample_scores(x.to_numpy(dtype=str),n=30,repeats=15,
                      null_shuffles=19,rng=np.random.default_rng(12))
    assert set(r)==set(MEASURES)
    assert 0<r["occupancy_r0"]<=r["metric_volume_r1"]<=r["metric_volume_r2"]<=1
    assert -1<=r["metric_volume_r1_richness_residual"]<=1
    assert -1<=r["colour_architecture_pairing_excess_r1"]<=1
    b=a.sample_scores(x.to_numpy(dtype=str),n=30,repeats=15,
                      null_shuffles=19,rng=np.random.default_rng(12))
    assert r==b


def test_null_uses_counterfactual_centers_not_fake_reference_occupancy():
    x=_triple_species()
    # Eliminate observed purple/tubular/zygomorphic while keeping enough
    # real Direct profiles, then allow a null colour permutation to form it.
    x=x.loc[~(
       x.flower_primary_color.eq("purple")&
       x.floral_form.eq("tubular")&
       x.floral_symmetry.eq("zygomorphic")
    )]
    a=FiniteProfileReference.from_direct_species(x,null_draws=199)
    assert len(a.profiles)==17
    missing=("purple","tubular","zygomorphic")
    assert missing not in a.profiles
    assert missing in a.masks
    x2=np.asarray([missing],dtype=str)
    r=a.coverage(x2)
    assert r[0]==0 # exact reference occupancy; hypothetical null centers carry no mass
    assert r[1]>0


def test_original_multistate_profile_is_rejected():
    x=_triple_species()
    x.iloc[0,0]="white|yellow"
    with pytest.raises(ValueError,match="Ambiguous"):
        FiniteProfileReference.from_direct_species(x,null_draws=199)


def test_block_and_mainland_source_aggregation_independent_guard():
    x=pd.DataFrame({
        "spatial_block":[f"block{i%8}" for i in range(40)],
        "source_signature":[f"source{i%5}" for i in range(40)],
        "metric":np.linspace(-.1,.1,40)})
    for key in ("spatial_block","source_signature"):
        z=_resample_group(x,"metric",key,seed=22,draws=999)
        assert len(z)==3
        assert all(np.isfinite(i) for i in z)
    x.source_signature="same"
    assert np.isnan(_resample_group(x,"metric","source_signature",seed=2)[0])
