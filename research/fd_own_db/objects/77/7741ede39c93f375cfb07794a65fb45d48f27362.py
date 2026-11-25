"""Source native labels and exact expected finite-n floral profile coverage tests."""
from itertools import product
import numpy as np
import pandas as pd
import pytest
from fd.regional_native_provenance_volume import (
    TRAITS,regional_native_candidates,ExactFiniteProfileSupport,
    cluster_equal_bootstrap
)


def source_direct():
    combos=list(product(
        ("white","purple","yellow","red","blue","pink","green","orange","brown","cream"),
        ("open_radial","tubular","brush","composite","zygomorphic"),
        ("actinomorphic","zygomorphic")))
    # Synthetic original-pattern reference: exactly 29 rare (1–2),
    # 46 middle (3–25) and 25 common (>=26), totaling 2234 species.
    counts=np.asarray([2]*29+[10]*46+[68]*25,dtype=int)
    counts[75:91]+=1
    assert counts.sum()==2234
    rows=[c for c,n in zip(combos,counts) for _ in range(int(n))]
    assert len(rows)==2234
    return pd.DataFrame(
        rows,index=[f"Genus{i%30} plant_{i}" for i in range(2234)],
        columns=TRAITS)


def test_wcvp_regional_native_does_not_override_introduced_and_unknown():
    original=pd.DataFrame([
      ("i1","Alpha a","native"),
      ("i1","Alpha b","introduced"),
      ("i1","Alpha c","unresolved"),
      ("i2","Alpha c","unresolved"),
      ("i1","Alpha d","unresolved"),
    ],columns=["island_id","accepted_species","origin_status"])
    mapping=pd.DataFrame([
      ("i1","SEA","accepted"),("i2","SEA","ambiguous")
    ],columns=["island_id","tdwg_l3_code","tdwg_match_status"])
    ranges=pd.DataFrame([
      ("Alpha b","SEA"),("Alpha c","SEA|IND"),("Alpha d","IND")
    ],columns=["accepted_species","native_l3_codes"])
    result=regional_native_candidates(original,mapping,ranges)
    sel=result.set_index(["island_id","accepted_species"])
    assert bool(sel.at[("i1","Alpha a"),"regional_native_compatible"])
    assert not bool(sel.at[("i1","Alpha b"),"regional_native_compatible"])
    assert bool(sel.at[("i1","Alpha c"),"regional_WCVP_upgrade_only"])
    assert not bool(sel.at[("i2","Alpha c"),"regional_WCVP_upgrade_only"])
    assert not bool(sel.at[("i1","Alpha d"),"regional_native_compatible"])
    with pytest.raises(ValueError,match="contradict"):
        regional_native_candidates(original,mapping,pd.concat([
            ranges,ranges.iloc[[2]].assign(native_l3_codes="SEA")]))


def test_analytic_30_species_coverage_exact_at_n_equals_pool_size():
    t=source_direct()
    model=ExactFiniteProfileSupport(t,reference_draws=1999,reference_seed=1)
    assert len(model.rare)==29 and len(model.common)==25
    subset=list(t.index[:30])
    r=model.source_scores(subset,seed=7,k_distribution_draws=2000)
    selected={tuple(t.loc[name].values) for name in subset}
    expected=np.array([
      any(sum(a!=b for a,b in zip(q,p))<=1 for p in selected)
      for q in model.profiles
    ])
    assert r["rare_coverage"]==pytest.approx(expected[model.rare].mean(),abs=1e-10)
    assert r["common_coverage"]==pytest.approx(expected[model.common].mean(),abs=1e-10)
    assert r["estimated_mean_distinct_k"]==pytest.approx(len(selected))
    assert r==model.source_scores(subset,seed=7,k_distribution_draws=2000)


def test_true_independent_geographic_bootstrap_does_not_use_flora_replicates():
    rows=pd.DataFrame({
        "block":[f"block{i%10}" for i in range(40)],
        "source":[f"s{i%15}" for i in range(40)],
        "gap":np.linspace(-.05,.08,40)
    })
    for x in ("block","source"):
        a=cluster_equal_bootstrap(rows,"gap",x,draws=999,seed=7)
        assert np.isfinite(a["mean"])
        assert len(a["ci95"])==2
        assert a["n_independent_clusters"]>=10
        assert a==cluster_equal_bootstrap(rows,"gap",x,draws=999,seed=7)
