import importlib

import pandas as pd

mod = importlib.import_module("island_v2.floristic_status_resolve")


def test_endemism_conflict_does_not_erase_native_origin():
    raw = pd.DataFrame(
        {
            "island_id": ["i1", "i1"],
            "accepted_species": ["A a", "A a"],
            "origin_status": ["native", "native"],
            "endemic_status": ["endemic", "nonendemic"],
            "status_reference": ["r1", "r2"],
        }
    )
    out = mod.collapse_gift_origin(raw).iloc[0]
    assert out["origin_status"] == "native"
    assert bool(out["gift_endemic_claim"])
    assert bool(out["gift_nonendemic_claim"])


def test_wcvp_multiple_native_regions_is_nonendemic_even_if_gift_claims_endemic():
    origin = pd.DataFrame(
        {
            "island_id": ["i1"],
            "accepted_species": ["A a"],
            "origin_status": ["native"],
            "origin_conflict": [False],
            "gift_endemic_claim": [True],
            "gift_nonendemic_claim": [False],
            "gift_status_references": ["r1"],
        }
    )
    wcvp = pd.DataFrame(
        {"accepted_species": ["A a"], "n_native_l3": [2], "native_l3_codes": ["AAA|BBB"]}
    )
    mapping = pd.DataFrame(
        {"island_id": ["i1"], "tdwg_l3_code": ["AAA"], "tdwg_l3_name": ["A"], "tdwg_match_status": ["accepted"]}
    )
    out = mod.resolve_endemism(origin, wcvp, mapping).iloc[0]
    assert out["endemic_status"] == "nonendemic"
    assert out["endemism_scope"] == "tdwg_l3"


def test_endemic_requires_gift_claim_and_single_matching_wcvp_region():
    origin = pd.DataFrame(
        {
            "island_id": ["i1", "i2"],
            "accepted_species": ["A a", "B b"],
            "origin_status": ["native", "native"],
            "origin_conflict": [False, False],
            "gift_endemic_claim": [True, False],
            "gift_nonendemic_claim": [False, False],
            "gift_status_references": ["r1", "r2"],
        }
    )
    wcvp = pd.DataFrame(
        {
            "accepted_species": ["A a", "B b"],
            "n_native_l3": [1, 1],
            "native_l3_codes": ["AAA", "BBB"],
        }
    )
    mapping = pd.DataFrame(
        {
            "island_id": ["i1", "i2"],
            "tdwg_l3_code": ["AAA", "BBB"],
            "tdwg_l3_name": ["A", "B"],
            "tdwg_match_status": ["accepted", "accepted"],
        }
    )
    out = mod.resolve_endemism(origin, wcvp, mapping).set_index("accepted_species")
    assert out.loc["A a", "endemic_status"] == "endemic"
    assert out.loc["B b", "endemic_status"] == "unresolved"
