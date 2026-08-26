import pandas as pd

from island_v2.chapter1_pr138_san_nicolas_status import (
    build_status_ledger,
    parse_cedros_oberbauer_text,
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
    assert by_key.loc["encelia californica", "origin_status"] == "introduced"
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
        status_source="test source",
        status_reference="https://example.org/sni.pdf",
        evidence_prefix="sni-test",
    )
    got = ledger.set_index("accepted_species")
    assert got.loc["Achillea millefolium", "origin_status"] == "native"
    assert got.loc["Brassica nigra", "origin_status"] == "introduced"
    assert got.loc["Unknown species", "origin_status"] == "unresolved"
    assert set(ledger["endemic_status"]) == {"unresolved"}


def test_cedros_parser_uses_actual_last_appendix_and_asterisk_status():
    text = """
    The island supports about 224 species (Appendix 2.).
    Avena barbata appears in vegetation prose before the appendix.
    Appendix 2.
    Preliminary Annotated List of Vascular Plants of Isla de Cedros, Baja California, Mexico
    Abronia maritima Nutt. ex Wats. Grows on western dunes.
    *Avena fatua L. Found on Cerro de Cedros.
    *Bromus rubens L. Found on Cerro de Cedros.
    Endemic to north facing slopes in northern canyons.
    Typha latifolia L. Listed by Moran as occurring.
    """
    parsed = parse_cedros_oberbauer_text(text).set_index("species_key")
    assert parsed.loc["abronia maritima", "origin_status"] == "native"
    assert parsed.loc["avena fatua", "origin_status"] == "introduced"
    assert parsed.loc["bromus rubens", "origin_status"] == "introduced"
    assert parsed.loc["typha latifolia", "origin_status"] == "native"
    assert "endemic to" not in parsed.index
    assert "avena barbata" not in parsed.index
