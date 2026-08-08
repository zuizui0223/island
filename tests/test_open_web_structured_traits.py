from __future__ import annotations

from island_v2.open_web_evidence import Page
from island_v2.open_web_structured_traits import extract_structured_characteristics


def _page(text: str) -> Page:
    return Page(
        requested_url="https://example.org/species/alpha/beta",
        final_url="https://example.org/species/alpha/beta",
        status_code=200,
        content_type="text/html",
        title="Alpha beta",
        text=text,
        language="en",
        retrieved_at_utc="2026-07-26T00:00:00Z",
        content_sha256="abc",
    )


def test_extracts_complete_structured_flower_fields() -> None:
    page = _page(
        "Flower petal color\nyellow\n"
        "Flower symmetry\nthere are two or more ways to evenly divide the flower "
        "(the flower is radially symmetrical)\n"
        "Flower diameter\n20–30 mm\n"
        "Cleistogamous flowers\nthere are no cleistogamous flowers on the plant\n"
        "Corolla morphology\nNA"
    )
    assert extract_structured_characteristics(
        page, trait_name="flower_primary_color"
    ) == [("yellow", "yellow_orange", "Flower petal color: yellow")]
    assert extract_structured_characteristics(
        page, trait_name="floral_symmetry"
    )[0][1] == "actinomorphic"
    assert extract_structured_characteristics(
        page, trait_name="flower_size_class"
    )[0][1] == "medium"
    assert extract_structured_characteristics(page, trait_name="cleistogamy")[0][1] == "absent"


def test_multiple_structured_options_remain_multistate() -> None:
    page = _page(
        "Cleistogamous flowers\n"
        "the plant has some cleistogamous flowers\n"
        "there are no cleistogamous flowers on the plant\n"
        "Corolla morphology\nNA"
    )
    rows = extract_structured_characteristics(page, trait_name="cleistogamy")
    assert len(rows) == 1
    assert rows[0][1] == "multistate_variable"
    assert "has some cleistogamous" in rows[0][2]
    assert "no cleistogamous" in rows[0][2]


def test_multiple_colours_are_not_arbitrarily_collapsed() -> None:
    page = _page(
        "Flower petal color\nblue to purple\nwhite\nLeaf type\nsimple"
    )
    rows = extract_structured_characteristics(
        page, trait_name="flower_primary_color"
    )
    assert len(rows) == 1
    assert rows[0][1] == "multicolored_variable"


def test_does_not_infer_from_unlabelled_or_unmapped_fields() -> None:
    page = _page("Leaf color\nyellow\nFruit diameter\n20–30 mm")
    assert extract_structured_characteristics(page, trait_name="flower_primary_color") == []
    assert extract_structured_characteristics(page, trait_name="flower_size_class") == []


def test_extracts_colon_delimited_university_profile_fields() -> None:
    page = _page(
        "Flowers:\n"
        "Flower Color:\nPink\nPurple/Lavender\nWhite\n"
        "Flower Inflorescence:\nRaceme\nSpike\n"
        "Flower Shape:\nTubular\n"
        "Flower Size:\n1-3 inches\n"
        "Flower Description:\nTwo-inch tubular flowers.\n"
    )
    assert extract_structured_characteristics(
        page, trait_name="flower_primary_color"
    )[0][1] == "multicolored_variable"
    assert extract_structured_characteristics(
        page, trait_name="inflorescence_display"
    )[0][1] == "raceme_spike_panicle"
    assert extract_structured_characteristics(page, trait_name="floral_form")[0][1] == "tubular"
    assert extract_structured_characteristics(
        page, trait_name="flower_size_class"
    )[0][1] == "large"
