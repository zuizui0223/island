import pandas as pd

from island_v2.chapter1_pr138_san_nicolas_status import (
    build_status_ledger,
    parse_cch2_san_nicolas_html,
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


def test_cedros_parser_uses_last_appendix_and_preserves_qualified_status():
    text = """
    The island supports about 224 species (Appendix 2.).
    Avena barbata appears in vegetation prose before the appendix.
    Appendix 2.
    Preliminary Annotated List of Vascular Plants of Isla de Cedros, Baja California, Mexico
    Abronia maritima Nutt. ex Wats. Grows on western dunes.  *Avena fatua L. Found on Cerro de Cedros.
    *Bromus rubens L. Found on Cerro de Cedros. Typha latifolia L. Listed by Moran as occurring.
    Endemic to north facing slopes in northern canyons.
    """
    parsed = parse_cedros_oberbauer_text(text).set_index("species_key")
    assert parsed.loc["abronia maritima", "origin_status"] == "unresolved"
    assert parsed.loc["avena fatua", "origin_status"] == "unresolved"
    assert parsed.loc["avena fatua", "status_basis"] == "presumably_introduced_asterisk"
    assert parsed.loc["bromus rubens", "origin_status"] == "unresolved"
    assert parsed.loc["typha latifolia", "origin_status"] == "unresolved"
    assert "endemic to" not in parsed.index
    assert "avena barbata" not in parsed.index


def test_cedros_ocr_match_is_unique_and_fail_closed():
    text = """
    Appendix 2.
    *Ceanothus verrucoslis Nutt. Found on upper slopes.
    Salvia cedroensis Greene. Known from the island.
    """
    parsed = parse_cedros_oberbauer_text(text)
    flora = pd.DataFrame(
        {
            "island_id": ["cedros", "cedros", "cedros"],
            "accepted_species": [
                "Ceanothus verrucosus",
                "Salvia cedrosensis",
                "Salvia cedronensis",
            ],
        }
    )
    ledger = build_status_ledger(
        flora,
        parsed,
        island_id="cedros",
        status_source="Oberbauer 1993",
        status_reference="https://example.org/cedros.pdf",
        evidence_prefix="cedros-test",
        allow_fuzzy_name_match=True,
    ).set_index("accepted_species")

    assert ledger.loc["Ceanothus verrucosus", "origin_status"] == "unresolved"
    assert ledger.loc["Ceanothus verrucosus", "name_match_method"] == "ocr_fuzzy_unique"
    assert (
        ledger.loc["Ceanothus verrucosus", "status_basis"]
        == "presumably_introduced_asterisk"
    )
    assert ledger.loc["Salvia cedrosensis", "origin_status"] == "unresolved"
    assert ledger.loc["Salvia cedronensis", "origin_status"] == "unresolved"


def test_cch2_parser_uses_terminal_taxa_and_island_status_notes():
    html = """
    <div class="taxon-container" id="parent">
      <div class="taxon-div"><span class="taxon-span">Example alpha</span></div>
    </div>
    <div class="taxon-container" id="native-child">
      <div class="taxon-div"><span class="taxon-span">Example alpha var. alpha</span></div>
      <div class="note-div">Channel Islands endemic; Voucher 1</div>
    </div>
    <div class="taxon-container" id="introduced-child">
      <div class="taxon-div"><span class="taxon-span">Example alpha var. beta</span></div>
      <div class="note-div">Native in CA, but presumably introduced to island</div>
    </div>
    <div class="taxon-container" id="introduced">
      <div class="taxon-div"><span class="taxon-span">Plantus secundus</span></div>
      <div class="note-div">Non-native, naturalized; <a>Voucher 2</a></div>
    </div>
    <div class="taxon-container" id="unclear">
      <div class="taxon-div"><span class="taxon-span">Plantus tertius</span></div>
      <div class="note-div">Native to CA, status on island unclear</div>
    </div>
    <div class="taxon-container" id="qualified-introduction">
      <div class="taxon-div"><span class="taxon-span">Plantus quartus</span></div>
      <div class="note-div">Native in CA, but possibly introduced to island</div>
    </div>
    """
    parsed = parse_cch2_san_nicolas_html(html).set_index("species_key")

    assert parsed.loc["example alpha", "origin_status"] == "unresolved"
    assert bool(parsed.loc["example alpha", "status_conflict"])
    assert parsed.loc["plantus secundus", "origin_status"] == "introduced"
    assert parsed.loc["plantus tertius", "origin_status"] == "unresolved"
    assert parsed.loc["plantus quartus", "origin_status"] == "unresolved"
    assert parsed.loc["plantus quartus", "status_basis"] == "island_introduction_qualified"
