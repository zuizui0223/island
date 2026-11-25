"""Synthetic source/taxon identity matching guards (no private original data)."""
import numpy as np
import pytest

from fd.native_taxon_matched_pairing import (
    RANKS, sample_taxon_matched_pairing, species_genus, support_slots
)
from fd.regional_native_flora_fd import METRICS


def example():
    own = set()
    source = set()
    reference = {}
    family = {}
    for i in range(9):
        genus=f"G{i}"
        for j, side in enumerate(("local", "donor")):
            name=f"{genus} species{j}"
            (own if side=="local" else source).add(name)
            family[name]=f"Family{i}"
            reference[name]=(
                ["pale", "blue", "red"][i % 3],
                "open_radial" if i%2 else "tubular",
                "actinomorphic" if i%2 else "zygomorphic"
            )
    return own, source, reference, family


def test_exact_genus_species_disjoint_and_reproducible():
    own,donor,direct,family=example()
    a,b,shared=support_slots(own,donor,rank="genus",families=family)
    assert len(shared)==9
    assert set(a)==set(b)==set(shared)
    r1=sample_taxon_matched_pairing(own,donor,direct,rank="genus",
        families=family,n=5,repeats=40,seed=20261114)
    r2=sample_taxon_matched_pairing(own,donor,direct,rank="genus",
        families=family,n=5,repeats=40,seed=20261114)
    assert r1==r2
    assert all(np.isfinite(r1["matched"][m]) for m in METRICS)
    assert all(abs(r1["matched"][m])<1e-12 for m in METRICS)
    assert r1["n_identity_matchable_groups"]==9


def test_family_same_and_different_genus_are_distinct():
    own,donor,direct,family=example()
    assert len(support_slots(own,donor,rank="family",families=family)[2])==9
    assert not support_slots(own,donor,rank="family_different_genus",families=family)[2]
    for i in range(9):
        name=f"Other{i} species1"
        donor.add(name)
        family[name]=f"Family{i}"
        direct[name]=direct[f"G{i} species0"]
    a,b,support=support_slots(own,donor,rank="family_different_genus",families=family)
    assert len(support)==9
    for k, candidates in support.items():
        assert all(species_genus(lhs)!=species_genus(rhs)
                   for lhs,options in candidates.items() for rhs in options)
    r=sample_taxon_matched_pairing(own,donor,direct,
        rank="family_different_genus",families=family,n=8,repeats=40,seed=3)
    assert all(np.isfinite(r["matched"][m]) for m in METRICS)


def test_unsupported_rank_and_overlap_fail_closed():
    own,donor,direct,family=example()
    with pytest.raises(ValueError):
        support_slots(own,donor|{next(iter(own))},rank="genus",families=family)
    with pytest.raises(ValueError):
        support_slots(own,donor,rank="phylo_guess",families=family)
    with pytest.raises(ValueError):
        sample_taxon_matched_pairing(own,donor,direct,rank="genus",
            families=family,n=10,repeats=40,seed=1)
    with pytest.raises(ValueError):
        sample_taxon_matched_pairing(own,donor,direct,rank="genus",
            families=family,n=5,repeats=3,seed=1)
