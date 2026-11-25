"""Synthetic checks: private name-level data are never used in CI."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from fd.core import ProfileSpace
from fd.genus_source_null import (
    _adjust, _distance_statistics, _permutable_groups, _source_input,
    infer_for_scope, load_private_maps, selected_native_pairs,
    validate_island_counts,
)

PROFILES = {
    0: (frozenset({"red"}), frozenset({"bell"})),
    1: (frozenset({"blue"}), frozenset({"bell"})),
    2: (frozenset({"red"}), frozenset({"tubular"})),
    3: (frozenset({"blue"}), frozenset({"tubular"})),
}


def synthetic():
    rows = []
    # Three groups, each containing distinct profiles. Some taxa occur in
    # both regions: all islands must reuse that taxon's same assigned profile.
    specs = [("Abc one", "Abc", 0),
             ("Abc two", "Abc", 1),
             ("Abc three", "Abc", 2),
             ("Abc four", "Abc", 3),
             ("Def one", "Def", 0),
             ("Def two", "Def", 1)]
    membership = {
        "C1": ["Abc one", "Abc two", "Def one"],
        "O1": ["Abc two", "Abc three", "Def two"],
        "C2": ["Abc one", "Abc four", "Def two"],
        "O2": ["Abc two", "Abc four", "Def one"],
    }
    for island, names in membership.items():
        region = "northern_midlatitude" if island.endswith("1") else "tropical"
        for name in names:
            genus, profile = next((g,p) for n,g,p in specs if n==name)
            rows.append({"island_id": island, "accepted_species": name,
                         "genus": genus, "profile_id": profile,
                         "analysis_regime": region})
    pairs = [
        {"continental_id": "C1", "oceanic_id": "O1",
         "region": "northern_midlatitude"},
        {"continental_id": "C2", "oceanic_id": "O2",
         "region": "tropical"},
    ]
    return pd.DataFrame(rows), pairs


def test_genus_signatures_preserve_profile_multisets():
    data,pairs=synthetic()
    table=_source_input(data,pairs,4)
    group_sig,n_groups,active=_permutable_groups(
        table["names"],"genus_region_signature")
    assert n_groups >= 2
    assert all(len(g)>=2 for g in group_sig)
    rng=np.random.default_rng(31)
    original=table["original"].copy()
    perm=original.copy()
    for group in group_sig:
        perm[group]=original[group][rng.permutation(len(group))]
        assert np.array_equal(np.sort(perm[group]),np.sort(original[group]))
    # Taxon identity and its occurrences on all islands remain untouched.
    assert len(table["islands"]) == 4
    assert all(len(index)==3 for index in table["islands"])
    assert active >= 2


def test_finite_permutations_replay_for_source_signature():
    data,pairs=synthetic()
    space=ProfileSpace(PROFILES)
    a=infer_for_scope("all",data,pairs,space,policy="genus_region_signature",
                      draws=99,seed=782)
    b=infer_for_scope("all",data,pairs,space,policy="genus_region_signature",
                      draws=99,seed=782)
    assert a==b
    assert len(a)==6
    assert all(0 < row["p_upper_approx"] <= 1 for row in a)
    assert len([row for row in a if row["region"]=="global"]) == 2
    assert all(row["n_islands"]==4 for row in a)
    with pytest.raises(ValueError):
        infer_for_scope("all",data,pairs,space,policy="not_a_null",draws=99,seed=782)


def test_published_counts_fail_closed():
    data,_=synthetic()
    expected=data.groupby("island_id").agg(
        n_names=("accepted_species","nunique"),
        n_profiles=("profile_id","nunique"),
    ).reset_index()
    # If one island carries two identical profiles the expected number differs:
    # the checker must fail, not silently drop the mismatch.
    assert validate_island_counts(data,expected,"all")["mismatches"]==0
    broken=expected.copy()
    broken.loc[0,"n_profiles"]=4
    with pytest.raises(ValueError,match="published island count mismatches"):
        validate_island_counts(data,broken,"all")


def test_private_profile_map_identity_validation(tmp_path):
    path=tmp_path/"private.csv"
    pd.DataFrame({
        "accepted_species":["A alpha","A beta"],
        "genus":["A","A"],"profile_id":[0,1],
    }).to_csv(path,index=False)
    profile=load_private_maps(path,"all",set(PROFILES),strict=False)
    assert len(profile)==2
    with pytest.raises(ValueError,match="population changed"):
        load_private_maps(path,"all",set(PROFILES),strict=True)
    # Exact source name/genus population must be validated.
    pd.DataFrame({
        "accepted_species":["A alpha","A beta"],
        "genus":["X","A"],"profile_id":[0,1],
    }).to_csv(path,index=False)
    with pytest.raises(ValueError,match="genus/name"):
        load_private_maps(path,"all",set(PROFILES),strict=False)


def test_native_eligibility_is_based_on_taxon_support_only():
    matched=[
        {"continental_id":"C1","oceanic_id":"O1","region":"tropical"},
        {"continental_id":"C2","oceanic_id":"O2","region":"tropical"},
    ]
    first=pd.DataFrame({"island_id":["C1"]*21+["O1"]*23+["C2"]*2+["O2"]*25})
    second=pd.DataFrame({"island_id":["C1"]*22+["O1"]*22+["C2"]*21+["O2"]*22})
    selected=selected_native_pairs(matched,first,second,min_names=20)
    assert len(selected)==1
    assert selected[0]["continental_id"]=="C1"
