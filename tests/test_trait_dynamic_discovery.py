from island_v2.trait_dynamic_discovery import (
    _mentions_name,
    _parse_bing_rss,
    _parse_duck_html,
    _source_type,
    _unwrap_duck_url,
)


def test_parse_bing_rss() -> None:
    xml = """<rss><channel><item><title>Flora page</title><link>https://example.org/a</link><description>Flowers yellow.</description></item></channel></rss>"""
    rows = _parse_bing_rss(xml)
    assert rows == [{"title": "Flora page", "url": "https://example.org/a", "snippet": "Flowers yellow."}]


def test_parse_duck_html_and_unwrap() -> None:
    html = '<a class="result__a" href="//duckduckgo.com/l/?uddg=https%3A%2F%2Fexample.org%2Fflora">Result</a><a class="result__snippet">Flowers white.</a>'
    rows = _parse_duck_html(html)
    assert rows[0]["url"] == "https://example.org/flora"
    assert rows[0]["snippet"] == "Flowers white."
    assert _unwrap_duck_url(rows[0]["url"]) == rows[0]["url"]


def test_species_identity_check_accepts_synonym_name() -> None:
    assert _mentions_name("Draba brachycarpa has white petals.", ["Abdra brachycarpa", "Draba brachycarpa"])
    assert not _mentions_name("Another species has white petals.", ["Abdra brachycarpa"])


def test_source_type_prioritizes_flora_and_pdf() -> None:
    assert _source_type("https://www.efloras.org/florataxon.aspx?id=1") == "flora_or_monograph"
    assert _source_type("https://example.org/paper.pdf") == "paper_or_flora_pdf"
