import importlib

import pandas as pd

mod = importlib.import_module("island_v2.native_status_trait_analysis")


def test_binary_mapping_excludes_cross_side_multistate_records():
    evidence = pd.DataFrame(
        {
            "accepted_species": ["A a", "B b", "C c", "D d"],
            "trait_name": ["flower_primary_color"] * 4,
            "normalized_value": [
                "white",
                "red_pink|yellow_orange",
                "white|red_pink",
                "multicolored_variable",
            ],
            "resolution_status": ["resolved"] * 4,
        }
    )
    out = mod.direct_binary_species(
        evidence,
        trait_name="flower_primary_color",
        positive={"white", "green_brown_inconspicuous"},
        negative={"yellow_orange", "red_pink", "blue_purple"},
    ).set_index("accepted_species")
    assert out.loc["A a", "state"] == 1
    assert out.loc["B b", "state"] == 0
    assert "C c" not in out.index
    assert "D d" not in out.index


def test_confirmed_native_excludes_introduced_and_unresolved():
    status = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i1"],
            "accepted_species": ["A a", "B b", "C c"],
            "origin_status": ["native", "introduced", "unresolved"],
        }
    )
    out = mod.confirmed_native_flora(status)
    assert out["accepted_species"].tolist() == ["A a"]


def test_native_counts_keep_species_trials():
    status = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i2"],
            "accepted_species": ["A a", "B b", "A a"],
            "origin_status": ["native", "native", "native"],
        }
    )
    evidence = pd.DataFrame(
        {
            "accepted_species": ["A a", "B b"],
            "trait_name": ["self_incompatibility", "self_incompatibility"],
            "normalized_value": ["SC", "SI"],
            "resolution_status": ["resolved", "resolved"],
        }
    )
    native = mod.confirmed_native_flora(status)
    counts = mod.build_native_counts(native, evidence)
    sc = counts.loc[counts["outcome"].eq("self_compatibility")].set_index("island_id")
    assert sc.loc["i1", "trials"] == 2
    assert sc.loc["i1", "successes"] == 1
    assert sc.loc["i2", "trials"] == 1
    assert sc.loc["i2", "successes"] == 1
