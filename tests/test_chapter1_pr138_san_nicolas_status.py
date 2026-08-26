import pandas as pd

from island_v2.chapter1_pr138_san_nicolas_status import (
    build_status_ledger,
    parse_nps_san_nicolas_text,
)


def test_parse_and_build_status_ledger_fail_closed():
    text = """
    Achillea millefolium yarrow N
    Brassica nigra black mustard E
    Encelia californica bush sunflower N1
    Example conflict common name N
    Example conflict common name E
    Family Asteraceae
    """
    parsed = parse_nps_san_nicolas_text(text)
    by_key = parsed.set_index("species_key")
    assert by_key.loc["achillea millefolium", "origin_status"] == "native"
    assert by_key.loc["brassica nigra", "origin_status"] == "introduced"
    assert by_key.loc["encelia californica", "origin_status"] == "native"
    assert by_key.loc["example conflict", "origin_status"] == "unresolved"
    assert bool(by_key.loc["example conflict", "status_conflict"])

    flora = pd.DataFrame(
        {
            "island_id": ["sni", "sni", "sni", "other"],
            "accepted_species": [
                "Achillea millefolium",
                "Brassica nigra",
                "Unknown species",
                "Achillea millefolium",
            ],
        }
    )
    ledger = build_status_ledger(
        flora,
        parsed,
        island_id="sni",
        status_reference="https://example.org/sni.pdf",
    )
    got = ledger.set_index("accepted_species")
    assert got.loc["Achillea millefolium", "origin_status"] == "native"
    assert got.loc["Brassica nigra", "origin_status"] == "introduced"
    assert got.loc["Unknown species", "origin_status"] == "unresolved"
    assert set(ledger["endemic_status"]) == {"unresolved"}
