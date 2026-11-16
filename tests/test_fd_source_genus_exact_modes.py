"""Synthetic equivalence to independently validated exact genus decomposition."""
import numpy as np
import pytest

from fd_direct_source_genus_accounting import exact_genus_bin_decomposition
from run_fd_source_genus_exact_modes import (
    source_genus_counts, exact_genus_from_source_counts,
)


def check(local,source,states):
    original=exact_genus_bin_decomposition(local,source,states)
    cached=exact_genus_from_source_counts(
        local,source_genus_counts(source,states),states)
    for field in (
        "original_direct_binary_sector_difference",
        "shared_genus_composition","shared_genus_within_species_states",
        "native_exclusive_genus","source_exclusive_genus_negative",
        "n_both_status_shared_genera",
        "n_native_source_carriers_from_shared_genera"
    ):
        assert np.isclose(cached[field],original[field],atol=1e-12)


def test_shared_and_unique_original_botanical_genera():
    local={"Alpha islandone","Alpha islandtwo","Beta islandone",
           "Gamma islandone"}
    source={"Alpha outsideone","Alpha outsidetwo","Beta outsideone",
            "Delta outsideone","Delta outsidetwo"}
    states={
        "Alpha islandone":1,"Alpha islandtwo":0,"Beta islandone":1,
        "Gamma islandone":1,"Alpha outsideone":0,"Alpha outsidetwo":1,
        "Beta outsideone":0,"Delta outsideone":0,"Delta outsidetwo":1
    }
    check(local,source,states)


def test_no_shared_genus_is_still_exactly_reconciled():
    local={"Alpha islandone","Alpha islandtwo","Beta islandone"}
    source={"Gamma otherone","Gamma othertwo","Delta otherone"}
    states={x: int("one" in x) for x in local|source}
    check(local,source,states)


def test_complete_shared_genera_but_different_category_frequencies():
    local={"Alpha isleone","Alpha isletwo","Beta isleone",
           "Beta isletwo"}
    source={"Alpha mainlandone","Alpha mainlandtwo","Beta mainlandone",
            "Beta mainlandtwo","Beta mainlandthree"}
    states={x:int(x.startswith("Alpha")) for x in local|source}
    states["Alpha mainlandtwo"]=0
    states["Beta isletwo"]=1
    check(local,source,states)


def test_unrecorded_donor_trait_cannot_be_silently_filled():
    with pytest.raises(ValueError,match="unrecorded"):
        source_genus_counts({"Alpha unseen"}, {"Beta seen":1})


def test_null_reference_from_empty_original_names_fails():
    with pytest.raises(ValueError,match="empty"):
        exact_genus_from_source_counts(set(),{"Alpha":[1,1]}, {})
