import pandas as pd

from island_v2.chapter1_wcvp_partition_diagnostic import classify_wcvp_partitions


def test_classify_wcvp_partitions_preserves_source_status_and_splits_unresolved():
    flora = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i1", "i1", "i2"],
            "accepted_species": ["native", "intro", "compat", "incompat", "unknown"],
            "origin_status": ["native", "introduced", "unresolved", "unresolved", "unresolved"],
            "floristic_status": [
                "native_nonendemic",
                "introduced",
                "unresolved",
                "unresolved",
                "unresolved",
            ],
        }
    )
    ranges = pd.DataFrame(
        {
            "accepted_species": ["native", "intro", "compat", "incompat"],
            "native_l3_codes": ["A", "A", "A|B", "B"],
        }
    )
    mapping = pd.DataFrame(
        {
            "island_id": ["i1", "i2"],
            "tdwg_l3_code": ["A", ""],
            "tdwg_match_status": ["accepted", "no_match"],
        }
    )
    out = classify_wcvp_partitions(flora, ranges, mapping).set_index("accepted_species")
    assert out.loc["native", "wcvp_partition"] == "source_native"
    assert out.loc["intro", "wcvp_partition"] == "source_introduced"
    assert out.loc["compat", "wcvp_partition"] == "wcvp_compatible_unresolved"
    assert out.loc["incompat", "wcvp_partition"] == "wcvp_incompatible_unresolved"
    assert out.loc["unknown", "wcvp_partition"] == "wcvp_unclassifiable_unresolved"
