from __future__ import annotations

from island_v2.open_web_evidence import Page
from island_v2.open_web_structured_traits import (
    extract_botanical_description,
    extract_structured_characteristics,
)


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


def test_rejects_ambiguous_upper_bound_flower_size() -> None:
    page = _page("Flower Size:\n< 1 inch\nFlower Description:\nSmall flowers")
    assert extract_structured_characteristics(page, trait_name="flower_size_class") == []


def test_extracts_organ_scoped_plantnet_narrative_without_calyx_colour() -> None:
    page = _page(
        "Description:\n"
        "Inflorescence an erect spike arising from a basal leaf rosette. "
        "Calyx c. 4 cm long, purple. "
        "Corolla c. 4–4.5 cm long, white; tube thickened distally.\n"
        "Distribution and occurrence:\nCultivated and naturalised."
    )
    assert extract_botanical_description(
        page,
        trait_name="flower_primary_color",
        treatment_end_pattern="Text by",
    ) == [
        (
            "Corolla c. 4–4.5 cm long, white;",
            "white",
            "Corolla c. 4–4.5 cm long, white;",
        )
    ]
    assert extract_botanical_description(
        page, trait_name="flower_size_class"
    )[0][:2] == ("4–4.5 cm", "large")
    assert extract_botanical_description(
        page, trait_name="inflorescence_display"
    )[0][1] == "raceme_spike_panicle"


def test_narrative_tube_depth_requires_explicit_tube_measurement() -> None:
    overall = _page("Description:\nCorolla c. 4 cm long, white; tube thickened distally.")
    explicit = _page("Description:\nCorolla tube 18–22 mm long, white.")
    assert extract_botanical_description(overall, trait_name="tube_depth_class") == []
    assert (
        extract_botanical_description(explicit, trait_name="tube_depth_class")[0][1]
        == "intermediate"
    )


def test_narrative_size_rejects_calyx_and_inflorescence_measurements() -> None:
    page = _page(
        "Description:\n"
        "Calyx lobes form the upper lip of corolla, c. 4 cm long. "
        "Flowers borne in racemes to 10 cm long; "
        "corolla tubular-campanulate, 3–4 cm long."
    )
    assert extract_botanical_description(page, trait_name="flower_size_class") == [
        ("3–4 cm", "large", "corolla tubular-campanulate, 3–4 cm long.")
    ]


def test_narrative_colour_rejects_calyx_colour_near_corolla_reference() -> None:
    page = _page(
        "Description:\n"
        "Calyx lobes form the upper lip of corolla, purple. Corolla white."
    )
    assert extract_botanical_description(
        page, trait_name="flower_primary_color"
    ) == [("Corolla white.", "white", "Corolla white.")]


def test_narrative_rejects_component_lengths_as_whole_flower_size() -> None:
    page = _page(
        "Description:\nPetals 3–6 mm long. Corolla tube 4–8 mm long. "
        "Flowers 10–20 mm diam. Flowers solitary on stalks 5–40 cm long."
    )
    assert extract_botanical_description(page, trait_name="flower_size_class") == [
        ("10–20 mm", "small", "Flowers 10–20 mm diam")
    ]


def test_narrative_size_rejects_spathe_and_branch_measurements() -> None:
    page = _page(
        "Description:\nFlowers subtended by a spathe 20–30 mm long. "
        "Flowers on branches 4–6 cm long. Corolla 12–16 mm long."
    )
    assert extract_botanical_description(page, trait_name="flower_size_class") == [
        ("12–16 mm", "small", "Corolla 12–16 mm long.")
    ]


def test_narrative_tube_depth_rejects_bare_anther_and_width_measurements() -> None:
    page = _page(
        "Description:\nTube 6–12 mm long. Anther tube 1 mm long. "
        "Perianth tube 6–12 mm wide. Corolla tube slender, 7–14 mm long."
    )
    assert extract_botanical_description(page, trait_name="tube_depth_class") == [
        (
            "Corolla tube slender, 7–14 mm long.",
            "intermediate",
            "Corolla tube slender, 7–14 mm long.",
        )
    ]


def test_narrative_solitary_spikelets_are_not_solitary_flower_display() -> None:
    page = _page(
        "Description:\nInflorescence compound with many solitary spikelets. "
        "Flowers solitary or in pairs."
    )
    assert extract_botanical_description(
        page, trait_name="inflorescence_display"
    ) == [("Flowers solitary or in pairs.", "solitary", "Flowers solitary or in pairs.")]


def test_narrative_flower_head_shape_is_not_individual_floral_form() -> None:
    page = _page("Description:\nFlower heads cup-shaped. Corolla campanulate.")
    assert extract_botanical_description(page, trait_name="floral_form") == [
        ("Corolla campanulate.", "bell_campanulate", "Corolla campanulate.")
    ]


def test_narrative_colour_preserves_suffix_variants_as_multicolour() -> None:
    page = _page(
        "Description:\nCorolla white, greenish or yellowish, with maroon throat."
    )
    assert extract_botanical_description(
        page, trait_name="flower_primary_color"
    )[0][1] == "multicolored_variable"
