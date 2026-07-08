from pathlib import Path

import pandas as pd

from island_v2.trait_descriptive_scout import (
    PROHIBITED_OUTPUT_COLUMNS,
    extract_candidates,
    load_keywords,
    scout_descriptive,
)

CONFIG = load_keywords(Path("config/m0_descriptive_keywords.yml"))


def test_extract_color_and_symmetry_with_excerpt():
    text = "Flowers with yellow petals, corolla zygomorphic and two-lipped."
    rows = extract_candidates("Genus species", text, {}, CONFIG)
    values = {(r["trait_name"], r["provisional_candidate_value"]) for r in rows}
    assert ("flower_primary_color", "yellow") in values
    assert ("floral_symmetry", "zygomorphic") in values
    assert all(r["review_status"] == "unreviewed" for r in rows)
    assert all(r["candidate_kind"] == "descriptive_keyword_match" for r in rows)
    yellow = next(r for r in rows if r["provisional_candidate_value"] == "yellow")
    assert "yellow" in yellow["evidence_excerpt"].lower()


def test_synonym_maps_to_canonical_value():
    rows = extract_candidates("Genus species", "Petals violet to mauve.", {}, CONFIG)
    assert any(r["provisional_candidate_value"] == "purple" for r in rows)
    # canonical value, not the surface term
    assert all(r["provisional_candidate_value"] != "violet" for r in rows)


def test_multiple_colors_are_not_collapsed():
    rows = extract_candidates("Genus species", "Corolla white, sometimes pink.", {}, CONFIG)
    colors = {r["provisional_candidate_value"] for r in rows if r["trait_name"] == "flower_primary_color"}
    assert {"white", "pink"} <= colors


def test_whole_word_only_no_substring_false_match():
    # "redirect" must not trigger a 'red' colour candidate
    rows = extract_candidates("Genus species", "See redirect note; greenhouse grown.", {}, CONFIG)
    assert not any(r["provisional_candidate_value"] == "red" for r in rows)


def test_no_controlled_term_yields_no_candidate():
    rows = extract_candidates("Genus species", "A tall shrub of montane forest.", {}, CONFIG)
    assert rows == []


def _fake_getter(match_key, descriptions):
    def getter(url, params):
        if "/species/match" in url:
            return {"usageKey": match_key} if match_key else {}
        if url.endswith("/descriptions"):
            return {"results": descriptions}
        return {}

    return getter


def test_scout_descriptive_builds_candidates_and_report():
    getter = _fake_getter(
        123,
        [
            {"type": "morphology", "description": "Flowers red and actinomorphic.", "source": "Flora X"},
            {"type": "distribution", "description": "Native to island Y."},  # wrong type -> ignored
        ],
    )
    df = pd.DataFrame({"accepted_species": ["Genus species"]})
    frame, report = scout_descriptive(df, getter, CONFIG, max_taxa=50, max_descriptions=20)
    assert PROHIBITED_OUTPUT_COLUMNS.isdisjoint(frame.columns)
    vals = set(frame["provisional_candidate_value"])
    assert "red" in vals and "actinomorphic" in vals
    # distribution text was skipped -> no candidates from it
    assert report["n_species_with_candidates"] == 1
    assert report["primary_source"] == "gbif_species_descriptions"
    assert frame.iloc[0]["source_url"] == "https://www.gbif.org/species/123"


def test_species_without_key_yields_zero_candidates():
    getter = _fake_getter(None, [])
    df = pd.DataFrame({"accepted_species": ["Unknown taxon"]})
    frame, report = scout_descriptive(df, getter, CONFIG, max_taxa=50, max_descriptions=20)
    assert len(frame) == 0
    assert report["n_species_zero_candidates"] == 1


def test_source_failure_does_not_abort_run():
    def flaky(url, params):
        if "/species/match" in url:
            raise RuntimeError("gbif down")
        return {}

    df = pd.DataFrame({"accepted_species": ["Genus species"]})
    frame, report = scout_descriptive(df, flaky, CONFIG, max_taxa=50, max_descriptions=20)
    assert len(frame) == 0
    assert report["n_lookup_errors"] == 1
