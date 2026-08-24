from __future__ import annotations

import hashlib
from pathlib import Path

import pytest
import yaml
from bs4 import BeautifulSoup

from island_v2.open_web_evidence import (
    DomainRegistry,
    Page,
    StrictNameResolver,
    canonical_url,
    dedupe_candidates,
    extract_rules,
    identity_gate,
    process_pages,
)

ROOT = Path(__file__).resolve().parents[1]


def _config() -> dict:
    return yaml.safe_load(
        (ROOT / "config" / "open_web_trait_acquisition.yml").read_text(encoding="utf-8")
    )


def _mikawa_page() -> Page:
    raw = (ROOT / "tests" / "fixtures" / "open_web" / "mikawa_species.html").read_text(
        encoding="utf-8"
    )
    soup = BeautifulSoup(raw, "html.parser")
    return Page(
        requested_url="https://mikawanoyasou.org/data/mikawainunohige.htm",
        final_url="https://mikawanoyasou.org/data/mikawainunohige.htm",
        status_code=200,
        content_type="text/html",
        title=soup.title.get_text(" ", strip=True),
        text="\n".join(soup.stripped_strings),
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256=hashlib.sha256(raw.encode()).hexdigest(),
    )


def test_mikawa_registry_example_is_tier_b_and_page_scoped() -> None:
    registry = DomainRegistry.load(ROOT / "config" / "open_web_domain_registry.yml")
    record = registry.lookup("https://mikawanoyasou.org/data/mikawainunohige.htm")
    assert record is not None
    assert record.recommended_evidence_tier == "B"
    assert record.source_type == "specialist_regional_flora"
    assert registry.lookup("https://mikawanoyasou.org/index.shtml") is None


def test_identity_gate_accepts_exact_species_and_family() -> None:
    resolver = StrictNameResolver(
        [
            {
                "accepted_species": "Eriocaulon mikawanum",
                "family": "Eriocaulaceae",
            }
        ]
    )
    resolution, page_name, reason = identity_gate(
        _mikawa_page(),
        searched_name="Eriocaulon mikawanum",
        expected_species="Eriocaulon mikawanum",
        resolver=resolver,
    )
    assert not reason
    assert page_name == "Eriocaulon mikawanum"
    assert resolution.method == "exact_accepted_name"


def test_inventory_identity_uses_page_subject_not_a_name_buried_in_body() -> None:
    resolver = StrictNameResolver(
        [
            {"accepted_species": "Plantus focalis", "family": "Plantaceae"},
            {"accepted_species": "Plantus referenced", "family": "Plantaceae"},
        ]
    )
    page = Page(
        requested_url="https://mikawanoyasou.org/data/focal.htm",
        final_url="https://mikawanoyasou.org/data/focal.htm",
        status_code=200,
        content_type="text/html",
        title="地域植物 Plantus focalis Plantaceae",
        text="Plantus focalis の花は白色。参考: Plantus referenced の花は赤色。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    resolution, page_name, reason = identity_gate(
        page,
        searched_name="",
        expected_species="",
        resolver=resolver,
    )
    assert not reason
    assert page_name == "Plantus focalis"
    assert resolution.accepted_species == "Plantus focalis"


def test_identity_gate_rejects_an_accepted_hybrid_even_without_x_in_title() -> None:
    resolver = StrictNameResolver(
        [
            {
                "accepted_species": "Erica arborea x Erica scoparia",
                "family": "Ericaceae",
            }
        ],
        [
            {
                "synonym": "Erica colorans",
                "accepted_species": "Erica arborea x Erica scoparia",
                "family": "Ericaceae",
                "backbone": "GBIF",
            },
            {
                "synonym": "Erica colorans",
                "accepted_species": "Erica arborea x Erica scoparia",
                "family": "Ericaceae",
                "backbone": "WFO",
            },
        ],
    )
    page = Page(
        requested_url="https://mikawanoyasou.org/data/erica_colorans.htm",
        final_url="https://mikawanoyasou.org/data/erica_colorans.htm",
        status_code=200,
        content_type="text/html",
        title="エリカ・コロランス Erica colorans Ericaceae",
        text="Erica colorans の花冠は鐘形。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    _, _, reason = identity_gate(
        page,
        searched_name="",
        expected_species="",
        resolver=resolver,
    )
    assert reason == "accepted_species_is_hybrid"


def test_identity_gate_rejects_infraspecific_page_subject() -> None:
    resolver = StrictNameResolver([{"accepted_species": "Plantus alba", "family": "Plantaceae"}])
    page = Page(
        requested_url="https://mikawanoyasou.org/data/infraspecies.htm",
        final_url="https://mikawanoyasou.org/data/infraspecies.htm",
        status_code=200,
        content_type="text/html",
        title="Plantus alba var. minor Plantaceae",
        text="Plantus alba var. minor の花は白色。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    _, _, reason = identity_gate(
        page,
        searched_name="",
        expected_species="",
        resolver=resolver,
    )
    assert reason == "page_subject_rank_is_infraspecific"


def test_identity_gate_rejects_form_dot_rank() -> None:
    resolver = StrictNameResolver([{"accepted_species": "Plantus alba", "family": "Plantaceae"}])
    page = Page(
        requested_url="https://mikawanoyasou.org/data/form.htm",
        final_url="https://mikawanoyasou.org/data/form.htm",
        status_code=200,
        content_type="text/html",
        title="Plantus alba form. normalis Plantaceae",
        text="Plantus alba form. normalis の花は白色。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    _, _, reason = identity_gate(
        page,
        searched_name="",
        expected_species="",
        resolver=resolver,
    )
    assert reason == "page_subject_rank_is_infraspecific"


def test_identity_gate_rejects_cultivar_scientific_name_subject() -> None:
    resolver = StrictNameResolver([{"accepted_species": "Plantus alba", "family": "Plantaceae"}])
    page = Page(
        requested_url="https://mikawanoyasou.org/data/cultivar.htm",
        final_url="https://mikawanoyasou.org/data/cultivar.htm",
        status_code=200,
        content_type="text/html",
        title="園芸植物 Plantus alba Plantaceae",
        text="学名 Plantus alba 'Snow White'.\n野生種の花は白色。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    _, _, reason = identity_gate(
        page,
        searched_name="",
        expected_species="",
        resolver=resolver,
    )
    assert reason == "labeled_scientific_name_is_cultivar"


def test_identity_gate_rejects_title_taxon_when_labeled_scientific_name_differs() -> None:
    resolver = StrictNameResolver(
        [
            {"accepted_species": "Plantus alba", "family": "Plantaceae"},
            {"accepted_species": "Alienus roseus", "family": "Plantaceae"},
        ]
    )
    page = Page(
        requested_url="https://example.org/wrong-title.htm",
        final_url="https://example.org/wrong-title.htm",
        status_code=200,
        content_type="text/html",
        title="Plantus alba Plantaceae",
        text="学名 Alienus roseus\n花は赤色。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    _, _, reason = identity_gate(
        page,
        searched_name="",
        expected_species="",
        resolver=resolver,
    )
    assert reason == "labeled_scientific_name_mismatch"


def test_synonym_requires_two_backbones_with_same_accepted_name_and_family() -> None:
    master = [{"accepted_species": "Plantus alba", "family": "Plantaceae"}]
    one = [
        {
            "synonym": "Planta candida",
            "accepted_species": "Plantus alba",
            "family": "Plantaceae",
            "backbone": "GBIF",
        }
    ]
    assert StrictNameResolver(master, one).resolve("Planta candida").status.startswith("rejected")
    two = one + [{**one[0], "backbone": "WFO"}]
    resolution = StrictNameResolver(master, two).resolve("Planta candida")
    assert resolution.status == "accepted"
    assert resolution.backbone_count == 2


def test_rules_extract_flower_traits_but_not_fruit_colour() -> None:
    page = _mikawa_page()
    form = extract_rules(page, trait_name="floral_form", config=_config())
    size = extract_rules(page, trait_name="flower_size_class", config=_config())
    colour = extract_rules(page, trait_name="flower_primary_color", config=_config())
    assert any(value == "tubular" for _, value, _ in form)
    assert any(value == "very_small" for _, value, _ in size)
    assert not any("褐色" in excerpt for _, _, excerpt in colour)


def test_head_inflorescence_is_not_mislabelled_as_individual_floral_form() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/juncus.htm",
        final_url="https://mikawanoyasou.org/data/juncus.htm",
        status_code=200,
        content_type="text/html",
        title="イグサ Juncus example Juncaceae",
        text="頭花は直径5mmの半球形で、1頭花に3～6個の花がつく。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert not extract_rules(page, trait_name="floral_form", config=_config())
    assert not extract_rules(page, trait_name="flower_size_class", config=_config())
    assert any(
        value == "composite_display"
        for _, value, _ in extract_rules(page, trait_name="inflorescence_display", config=_config())
    )


def test_flower_measurement_ignores_leaf_size_earlier_in_same_paragraph() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/measure.htm",
        final_url="https://mikawanoyasou.org/data/measure.htm",
        status_code=200,
        content_type="text/html",
        title="試験種 Plantus mensura Plantaceae",
        text="葉は長さ30cm。花冠は筒状、長さ12mm、白色。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert extract_rules(page, trait_name="flower_size_class", config=_config()) == [
        ("12mm", "small", "花冠は筒状、長さ12mm、白色。")
    ]


def test_flower_measurement_does_not_cross_sentence_into_fruit_size() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/fruit.htm",
        final_url="https://mikawanoyasou.org/data/fruit.htm",
        status_code=200,
        content_type="text/html",
        title="試験種 Plantus fructus Plantaceae",
        text="花は雌雄別株。果実は直径8mm。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert not extract_rules(page, trait_name="flower_size_class", config=_config())


def test_plant_height_after_flower_count_is_not_flower_size() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/height.htm",
        final_url="https://mikawanoyasou.org/data/height.htm",
        status_code=200,
        content_type="text/html",
        title="試験種 Plantus altus Plantaceae",
        text="花を1個つけ、高さ10～30cm。花冠は長さ12mm。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert extract_rules(page, trait_name="flower_size_class", config=_config()) == [
        ("12mm", "small", "花冠は長さ12mm。")
    ]


@pytest.mark.parametrize(
    "text",
    [
        "花は小花柄が長さ1.5～2mm。",
        "花期の葉は長さ2～7cm。",
        "花時に長さ約1.5～3.5cm、果時には後屈する。",
        "花が100～500個つき、花序は長さ4～14cm。",
        "花がつき、長さ6～9mm。",
        "花を数個つけ、長さ1cm以下。",
        "花を密生し、長さ1～2cm。",
        "花が少なく、長さ1cm以下、根際に雌小穂がつく。",
        "花托筒は円筒形、長さ6～8mm。",
        "花のつくシュートでは長さ3～5cm。",
        "別種は花が長さ1.3cmで対象種に似ているが、対象種は長さ2cm。",
        "地下部は花に比べて非常に大きく、半径50cm。",
        "開花受粉の葯は長さ0.4～0.7mm。",
        "花弁状の苞は長さ2～6cm。",
        "花の護頴は長さ8～15mmの芒をつける。",
        "花の下部で直径1.5～2mm。",
        "花被に包まれた果実は直径0.8～1mm。",
        "花被の筒部は長さ2～3mm。",
        "花冠は車形で、基部は長さ0.5～0.7mmが融着する。",
        "花は無く、または1～3個つき、小舌は長さ2～3.5mm。",
        "Flora of Chinaでは長さ3～8mm、刺針状花被片がある。",
    ],
)
def test_nested_or_nonflower_dimensions_are_not_flower_size(text: str) -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/nonflower.htm",
        final_url="https://mikawanoyasou.org/data/nonflower.htm",
        status_code=200,
        content_type="text/html",
        title="試験種 Plantus nonflorus Plantaceae",
        text=text,
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert not extract_rules(page, trait_name="flower_size_class", config=_config())


def test_corolla_tube_length_is_tube_depth_not_whole_flower_size() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/tube.htm",
        final_url="https://mikawanoyasou.org/data/tube.htm",
        status_code=200,
        content_type="text/html",
        title="試験種 Plantus tubus Plantaceae",
        text="花冠の筒部は長さ12mm。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert not extract_rules(page, trait_name="flower_size_class", config=_config())
    assert extract_rules(page, trait_name="tube_depth_class", config=_config()) == [
        ("12mm", "intermediate", "花冠の筒部は長さ12mm。")
    ]


def test_spathe_tube_is_not_flower_tube_depth() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/spathe.htm",
        final_url="https://mikawanoyasou.org/data/spathe.htm",
        status_code=200,
        content_type="text/html",
        title="試験種 Plantus spatha Plantaceae",
        text="仏炎苞は帯白色。筒部は長楕円形、長さ約3.5㎝。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert not extract_rules(page, trait_name="tube_depth_class", config=_config())


@pytest.mark.parametrize(
    "text",
    [
        "筒部が短い花冠は直径8～12mm。",
        "筒部に毛があり、花被片は長さ5.3mm。",
        "花冠筒部の内面に長さ1mmの鱗片がある。",
        "筒部と拡大部に毛があり、拡大部は長さ4mm。",
        "花冠筒部は萼より短く、花冠裂片は長さ5mm。",
        "花冠筒部は萼の1.5倍、花柱は長さ1mm。",
    ],
)
def test_nearby_non_tube_dimension_is_not_tube_depth(text: str) -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/not-tube-depth.htm",
        final_url="https://mikawanoyasou.org/data/not-tube-depth.htm",
        status_code=200,
        content_type="text/html",
        title="試験種 Plantus tubus Plantaceae",
        text=text,
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert not extract_rules(page, trait_name="tube_depth_class", config=_config())


def test_appended_genus_catalogue_cannot_supply_a_species_trait() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/focal.htm",
        final_url="https://mikawanoyasou.org/data/focal.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus focalis Plantaceae",
        text=(
            "Plantus focalis の主記載には花形の記述がない。\n"
            "対象属 family Plantaceae－ genus Plantus\n"
            "Plantus other の花冠は筒状。"
        ),
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert not extract_rules(page, trait_name="floral_form", config=_config())


def test_appended_colon_genus_catalogue_cannot_create_multistate_symmetry() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/focal-colon.htm",
        final_url="https://mikawanoyasou.org/data/focal-colon.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus focalis Plantaceae",
        text=(
            "花は左右相称。\n対象属 family: Plantaceae －genus Plantus\n属では花はまれに放射相称。"
        ),
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert extract_rules(page, trait_name="floral_symmetry", config=_config()) == [
        ("左右相称", "zygomorphic", "花は左右相称。")
    ]


def test_mikawa_family_marker_stops_japanese_genus_catalogue() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/focal-mikawa.htm",
        final_url="https://mikawanoyasou.org/data/focal-mikawa.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus focalis Plantaceae",
        text=(
            "対象種の花冠は白色。\n"
            "対象属 科family:シソ科 Plantaceae- 亜科 Subfamily- 属 genus: Plantus\n"
            "別種の花冠は筒状で赤色。"
        ),
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    boundary = r"fami(?:ly|liy)\s*:?"
    assert not extract_rules(
        page,
        trait_name="floral_form",
        config=_config(),
        treatment_end_pattern=boundary,
    )
    assert extract_rules(
        page,
        trait_name="flower_primary_color",
        config=_config(),
        treatment_end_pattern=boundary,
    )[0][1:] == (
        "white",
        "対象種の花冠は白色。",
    )


def test_mikawa_misspelled_familiy_marker_stops_genus_catalogue() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/focal-familiy.htm",
        final_url="https://mikawanoyasou.org/data/focal-familiy.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus focalis Plantaceae",
        text="対象種には形の記述なし。familiy Plantaceae－genus Plantus。別種の花冠は筒状。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert not extract_rules(
        page,
        trait_name="floral_form",
        config=_config(),
        treatment_end_pattern=r"fami(?:ly|liy)\s*:?",
    )


@pytest.mark.parametrize(
    "text",
    [
        "痩果は4稜形（舌状花では3稜形）。",
        "別の植物は黄色の蝶形花で対象種に似ている。",
        "花は両性、完全花で、多くの種は放射相称。",
    ],
)
def test_form_terms_in_non_species_context_are_rejected(text: str) -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/form-context.htm",
        final_url="https://mikawanoyasou.org/data/form-context.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus focalis Plantaceae",
        text=text,
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert not extract_rules(page, trait_name="floral_form", config=_config())


def test_green_ovary_is_not_primary_flower_colour() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/ovary.htm",
        final_url="https://mikawanoyasou.org/data/ovary.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus ovary Plantaceae",
        text="雌花は緑色の球形の子房がある。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert not extract_rules(page, trait_name="flower_primary_color", config=_config())


def test_japanese_multicolour_flower_is_retained_as_state_set() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/multicolour.htm",
        final_url="https://mikawanoyasou.org/data/multicolour.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus varius Plantaceae",
        text="花冠は白色（～淡紫色）。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    result = extract_rules(page, trait_name="flower_primary_color", config=_config())
    assert result == [("blue_purple;white", "multicolored_variable", "花冠は白色（～淡紫色）。")]


def test_colour_attached_to_flower_survives_inflorescence_context() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/flower-colour.htm",
        final_url="https://mikawanoyasou.org/data/flower-colour.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus purpureus Plantaceae",
        text="短い円錐花序に淡紅紫色の花を固まってつける。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    result = extract_rules(page, trait_name="flower_primary_color", config=_config())
    assert result == [("紫色", "blue_purple", "短い円錐花序に淡紅紫色の花を固まってつける。")]


def test_yellowish_white_flower_is_multicolour_not_white_only() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/yellowish-white.htm",
        final_url="https://mikawanoyasou.org/data/yellowish-white.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus pallidus Plantaceae",
        text="花冠は淡黄白色。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    result = extract_rules(page, trait_name="flower_primary_color", config=_config())
    assert result == [("white;yellow_orange", "multicolored_variable", "花冠は淡黄白色。")]


def test_corolla_lobe_spot_is_retained_in_multicolour_state_set() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/corolla-spots.htm",
        final_url="https://mikawanoyasou.org/data/corolla-spots.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus maculatus Plantaceae",
        text="花冠は青紫色。花冠裂片に白色の斑点がある。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert extract_rules(page, trait_name="flower_primary_color", config=_config()) == [
        (
            "blue_purple;white",
            "multicolored_variable",
            "花冠は青紫色。花冠裂片に白色の斑点がある。",
        )
    ]


def test_generic_japanese_colour_term_does_not_capture_nearby_fruit() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/fruit-colour.htm",
        final_url="https://mikawanoyasou.org/data/fruit-colour.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus fructus Plantaceae",
        text="花は小さい。果実は赤色。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert not extract_rules(page, trait_name="flower_primary_color", config=_config())


@pytest.mark.parametrize(
    "text",
    [
        "雌花に褐色の毛が密生する。",
        "雌花の花柱に白色の腺点がある。",
        "花穂は緑色。",
        "小穂は白緑色で、小花を5個もつ。",
        "雌花の鱗片の芒は褐色を帯びる。",
    ],
)
def test_generic_colour_term_does_not_capture_non_corolla_structure(text: str) -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/non-corolla-colour.htm",
        final_url="https://mikawanoyasou.org/data/non-corolla-colour.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus structura Plantaceae",
        text=text,
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert not extract_rules(page, trait_name="flower_primary_color", config=_config())


def test_comparison_species_cannot_expand_target_colour_state_set() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/target-colour.htm",
        final_url="https://mikawanoyasou.org/data/target-colour.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus target Plantaceae",
        text="花冠は緑色～紫色。Plantus other の花冠は白色。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert extract_rules(page, trait_name="flower_primary_color", config=_config()) == [
        (
            "blue_purple;green_brown_inconspicuous",
            "multicolored_variable",
            "花冠は緑色～紫色。",
        )
    ]


def test_dried_specimen_colour_does_not_expand_fresh_flower_state_set() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/dried-colour.htm",
        final_url="https://mikawanoyasou.org/data/dried-colour.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus caeruleus Plantaceae",
        text="花冠は青色、乾くと白色。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert extract_rules(page, trait_name="flower_primary_color", config=_config()) == [
        ("青色", "blue_purple", "花冠は青色、乾くと白色。")
    ]


def test_stamen_that_only_looks_like_flower_colour_is_not_primary_colour() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/stamen-colour.htm",
        final_url="https://mikawanoyasou.org/data/stamen-colour.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus stamineus Plantaceae",
        text="花は小さく、黄色の花に見えるのは5個の雄しべである。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert not extract_rules(page, trait_name="flower_primary_color", config=_config())


def test_colour_state_negated_after_term_is_rejected() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/absent-spot.htm",
        final_url="https://mikawanoyasou.org/data/absent-spot.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus immaculatus Plantaceae",
        text="花冠に黄色の斑紋がない。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert not extract_rules(page, trait_name="flower_primary_color", config=_config())


def test_colour_states_in_separate_paragraphs_are_not_globally_collapsed() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/paragraph-colour.htm",
        final_url="https://mikawanoyasou.org/data/paragraph-colour.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus paragraphus Plantaceae",
        text="花冠は白色。\n花冠は紫色。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert extract_rules(page, trait_name="flower_primary_color", config=_config()) == [
        ("白色", "white", "花冠は白色。"),
        ("紫色", "blue_purple", "花冠は紫色。"),
    ]


@pytest.mark.parametrize(
    "text",
    [
        "頭花は",
        "花冠が白色のものを",
        "は全体に明るい緑色で、頭花の白い舌状花がはっきり見える",
    ],
)
def test_incomplete_trait_excerpt_is_rejected(text: str) -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/incomplete.htm",
        final_url="https://mikawanoyasou.org/data/incomplete.htm",
        status_code=200,
        content_type="text/html",
        title="対象種 Plantus capitatus Plantaceae",
        text=text,
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert not extract_rules(page, trait_name="floral_form", config=_config())
    assert not extract_rules(page, trait_name="flower_primary_color", config=_config())
    assert not extract_rules(page, trait_name="inflorescence_display", config=_config())


def test_inflorescence_or_calyx_symmetry_is_not_flower_symmetry() -> None:
    page = Page(
        requested_url="https://mikawanoyasou.org/data/grass.htm",
        final_url="https://mikawanoyasou.org/data/grass.htm",
        status_code=200,
        content_type="text/html",
        title="試験草 Gramina example Poaceae",
        text="総状花序は左右相称。萼は放射相称。",
        language="ja",
        retrieved_at_utc="2026-07-25T00:00:00Z",
        content_sha256="fixture",
    )
    assert not extract_rules(page, trait_name="floral_symmetry", config=_config())


def test_unknown_domain_and_tier_trait_restrictions_fail_closed() -> None:
    registry = DomainRegistry.load(ROOT / "config" / "open_web_domain_registry.yml")
    resolver = StrictNameResolver(
        [{"accepted_species": "Eriocaulon mikawanum", "family": "Eriocaulaceae"}]
    )

    class Fetcher:
        def fetch(self, _url: str) -> Page:
            return _mikawa_page()

    unknown = process_pages(
        [{"result_url": "https://unknown.example/species", "trait_name": "floral_form"}],
        registry=registry,
        resolver=resolver,
        config=_config(),
        fetcher=Fetcher(),  # type: ignore[arg-type]
    )
    assert not unknown["candidates"]
    assert unknown["identity_audit"][0]["reason"].startswith("domain_not_approved")

    tier_b_reproduction = process_pages(
        [
            {
                "result_url": "https://mikawanoyasou.org/data/mikawainunohige.htm",
                "accepted_species": "Eriocaulon mikawanum",
                "searched_name": "Eriocaulon mikawanum",
                "trait_name": "self_incompatibility",
            }
        ],
        registry=registry,
        resolver=resolver,
        config=_config(),
        fetcher=Fetcher(),  # type: ignore[arg-type]
    )
    assert not tier_b_reproduction["candidates"]
    assert any(
        row["reason"] == "source_tier_trait_not_allowed"
        for row in tier_b_reproduction["identity_audit"]
    )


def test_lineage_deduplicator_does_not_count_republished_excerpt_twice() -> None:
    base = {
        "accepted_species": "Plantus alba",
        "trait_name": "floral_form",
        "normalized_value": "tubular",
        "supporting_excerpt": "Flowers are tubular.",
        "source_lineage": "url:https://one.example/page",
    }
    selected, duplicates = dedupe_candidates(
        [base, {**base, "source_lineage": "url:https://mirror.example/page"}]
    )
    assert len(selected) == 1
    assert len(duplicates) == 1
    assert duplicates[0]["duplicate_reason"] == "same_excerpt_republished_or_mirrored"


def test_lineage_deduplicator_uses_citation_and_content_fingerprint() -> None:
    base = {
        "accepted_species": "Plantus alba",
        "trait_name": "floral_form",
        "normalized_value": "tubular",
        "supporting_excerpt": "Flowers are tubular.",
        "source_citation": "Flora Example (2020), treatment 4",
        "content_fingerprint": "abc123",
        "source_lineage": "url:https://one.example/page",
    }
    citation_mirror = {
        **base,
        "supporting_excerpt": "Flowers  are tubular!",
        "source_lineage": "url:https://mirror.example/page",
        "content_fingerprint": "different",
    }
    selected, duplicates = dedupe_candidates([base, citation_mirror])
    assert len(selected) == 1
    assert duplicates[0]["duplicate_reason"] == ("same_citation_and_normalized_excerpt")

    content_mirror = {
        **base,
        "source_citation": "Different citation",
        "supporting_excerpt": "Tubular flowers occur.",
        "source_lineage": "url:https://other.example/page",
    }
    selected, duplicates = dedupe_candidates([base, content_mirror])
    assert len(selected) == 1
    assert duplicates[0]["duplicate_reason"] == "same_content_fingerprint"


def test_canonical_url_removes_tracking_without_changing_source_path() -> None:
    assert canonical_url("https://www.Example.org/a/?utm_source=x&id=4#frag") == (
        "https://example.org/a?id=4"
    )
