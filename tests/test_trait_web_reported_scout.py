from pathlib import Path

import pandas as pd

from island_v2.trait_web_reported_scout import (
    PROHIBITED_OUTPUT_COLUMNS,
    extract_reported,
    load_keyword_layers,
    scout_web_reported,
)

CONFIG, TRAIT_LAYER = load_keyword_layers(
    {
        "M0": Path("config/m0_descriptive_keywords.yml"),
        "M1": Path("config/m1_reported_keywords.yml"),
    }
)


def test_layers_merge_tags_m0_and_m1():
    assert TRAIT_LAYER["flower_primary_color"] == "M0"
    assert TRAIT_LAYER["self_incompatibility_reported"] == "M1"


def test_extract_reported_keeps_full_provenance_schema():
    text = "The species is self-incompatible and bird-pollinated, with red flowers."
    rows = extract_reported(
        "Genus species", text, "wikipedia", "wikipedia",
        "https://en.wikipedia.org/wiki/Genus_species", "species_direct", CONFIG, TRAIT_LAYER,
    )
    by = {(r["trait_name"], r["provisional_candidate_value"]) for r in rows}
    assert ("self_incompatibility_reported", "self_incompatible") in by
    assert ("pollination_vector_reported", "bird") in by
    assert ("flower_primary_color", "red") in by
    r = rows[0]
    for col in ("candidate_class", "source_type", "source_url", "raw_description",
                "evidence_scope", "wild_or_cultivated", "review_status"):
        assert col in r
    assert all(r["candidate_class"] == "reported" for r in rows)
    assert all(r["review_status"] == "unreviewed" for r in rows)


def test_no_explicit_statement_yields_no_reported_candidate():
    rows = extract_reported(
        "Genus species", "A shrub found in montane forest.", "wikipedia", "wikipedia",
        "u", "species_direct", CONFIG, TRAIT_LAYER,
    )
    assert rows == []


def _fake_getter(qid="Q1", enwiki="Genus species", wd_desc="", extract="", missing=False):
    def getter(url, params):
        if "wikidata.org" in url and params.get("action") == "wbsearchentities":
            return {"search": [{"id": qid, "description": wd_desc}]} if qid else {"search": []}
        if "wikidata.org" in url and params.get("action") == "wbgetentities":
            return {"entities": {qid: {"sitelinks": {"enwiki": {"title": enwiki}} if enwiki else {}}}}
        if "wikipedia.org" in url:
            if missing:
                return {"query": {"pages": {"-1": {"missing": ""}}}}
            return {"query": {"pages": {"1": {"title": enwiki or "Genus species", "extract": extract}}}}
        return {}

    return getter


def test_scout_builds_reported_from_wikipedia_and_reports_coverage():
    getter = _fake_getter(
        enwiki="Genus species",
        extract="Flowers white; the plant is self-compatible and wind-pollinated.",
    )
    df = pd.DataFrame({"accepted_species": ["Genus species", "Other taxon"]})
    frame, report = scout_web_reported(df, getter, CONFIG, TRAIT_LAYER, max_taxa=50)
    assert PROHIBITED_OUTPUT_COLUMNS.isdisjoint(frame.columns)
    vals = set(frame["provisional_candidate_value"])
    assert {"white", "self_compatible", "wind"} <= vals
    assert (frame["candidate_class"] == "reported").all()
    assert report["n_species_with_reported"] >= 1
    assert set(report["n_by_source"]).issubset({"wikipedia", "wikidata"})


def test_missing_wikipedia_page_yields_zero_and_does_not_abort():
    getter = _fake_getter(enwiki="", missing=True)
    df = pd.DataFrame({"accepted_species": ["Nonexistent taxon"]})
    frame, report = scout_web_reported(df, getter, CONFIG, TRAIT_LAYER, max_taxa=50)
    assert len(frame) == 0
    assert report["n_species_zero_reported"] == 1


def test_evidence_scope_species_direct_when_title_matches():
    getter = _fake_getter(enwiki="Genus species", extract="Flowers red.")
    df = pd.DataFrame({"accepted_species": ["Genus species"]})
    frame, _ = scout_web_reported(df, getter, CONFIG, TRAIT_LAYER, max_taxa=50)
    assert (frame["evidence_scope"] == "species_direct").all()
