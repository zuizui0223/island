import json

import pandas as pd
import pytest
import typer

from island_v2.trait_source_discovery import (
    GRADE_A,
    GRADE_B,
    GRADE_C,
    GRADE_D,
    OUTPUT_COLUMNS,
    PROHIBITED_OUTPUT_COLUMNS,
    _require_nominated_island,
    discover_sources,
    grade_source_hint,
    parse_crossref_works,
    parse_openalex_works,
)


def _fake_getter(responses):
    """Build a getter that returns canned JSON keyed by URL substring."""

    def getter(url, params):
        for key, payload in responses.items():
            if key in url:
                if callable(payload):
                    return payload(url, params)
                return payload
        return {}

    return getter


OPENALEX_PAYLOAD = {
    "results": [
        {
            "display_name": "Pollination biology of Bombus-visited island plants",
            "doi": "https://doi.org/10.1234/abcd",
            "publication_year": 2019,
            "open_access": {"is_oa": True, "oa_status": "green", "oa_url": "https://oa.example/1"},
            "primary_location": {
                "source": {"display_name": "Journal of Island Botany"},
                "pdf_url": "https://oa.example/1.pdf",
            },
            "authorships": [{"author": {"display_name": "A. Botanist"}}],
        }
    ]
}

CROSSREF_PAYLOAD = {
    "message": {
        "items": [
            {
                "title": ["A floral trait revision"],
                "DOI": "10.5678/efgh",
                "published": {"date-parts": [[2021, 3]]},
                "container-title": ["Phytotaxa"],
                "author": [{"given": "B.", "family": "Taxonomist"}],
            }
        ]
    }
}

UNPAYWALL_PAYLOAD = {
    "is_oa": True,
    "oa_status": "gold",
    "best_oa_location": {"url": "https://unpaywall.example/x", "url_for_pdf": "https://unpaywall.example/x.pdf"},
}


def test_parse_openalex_extracts_seed_without_trait_values():
    seeds = parse_openalex_works(OPENALEX_PAYLOAD, "Antidesma vogelianum", limit=5)
    assert len(seeds) == 1
    seed = seeds[0]
    assert seed.doi == "10.1234/abcd"
    assert seed.source_api == "openalex"
    assert seed.is_open_access == "true"
    assert seed.candidate_page_locator.startswith("HUMAN_REVIEW_REQUIRED")
    assert seed.lead_status == "unreviewed_literature_seed"


def test_parse_crossref_extracts_year_and_venue():
    seeds = parse_crossref_works(CROSSREF_PAYLOAD, "Adenia lobata", limit=5)
    assert seeds[0].doi == "10.5678/efgh"
    assert seeds[0].publication_year == "2021"
    assert seeds[0].venue == "Phytotaxa"


def test_discover_sources_produces_only_lead_columns():
    getter = _fake_getter(
        {
            "openalex": OPENALEX_PAYLOAD,
            "crossref": CROSSREF_PAYLOAD,
            "unpaywall": UNPAYWALL_PAYLOAD,
        }
    )
    taxa = pd.DataFrame(
        {
            "accepted_species": ["Antidesma vogelianum", "Adenia lobata"],
            "island_id": ["isl_jp", "isl_jp"],
            "genus": ["Antidesma", "Adenia"],
            "family": ["Phyllanthaceae", "Passifloraceae"],
        }
    )
    frame, report = discover_sources(
        taxa, getter, contact_email="test@example.org", max_taxa=5, max_seeds_per_taxon=5
    )
    # No finalized-value column may ever appear.
    assert PROHIBITED_OUTPUT_COLUMNS.isdisjoint(frame.columns)
    for col in OUTPUT_COLUMNS:
        assert col in frame.columns
    # Provenance passthrough is preserved.
    assert "island_id" in frame.columns
    assert set(frame["query_taxon"]) == {"Antidesma vogelianum", "Adenia lobata"}
    assert report["n_taxa_queried"] == 2
    assert report["n_literature_seeds"] == len(frame)
    assert report["n_open_access_seeds"] >= 1
    # Unpaywall receipt was attached to the DOI-bearing OpenAlex seed.
    oa_seed = frame.loc[frame["source_api"] == "openalex"].iloc[0]
    assert oa_seed["candidate_pdf_url"] == "https://unpaywall.example/x.pdf"


def test_discover_sources_is_bounded_by_max_taxa():
    getter = _fake_getter({"openalex": OPENALEX_PAYLOAD, "crossref": {}, "unpaywall": {}})
    taxa = pd.DataFrame({"accepted_species": [f"Species number{i}" for i in range(20)]})
    _frame, report = discover_sources(
        taxa, getter, contact_email="t@example.org", max_taxa=5, max_seeds_per_taxon=3
    )
    assert report["n_taxa_queried"] == 5


def test_discover_sources_rejects_prohibited_input_column():
    getter = _fake_getter({})
    taxa = pd.DataFrame({"accepted_species": ["X y"], "trait_value": ["white"]})
    with pytest.raises(typer.BadParameter):
        discover_sources(taxa, getter, contact_email="t@example.org", max_taxa=5, max_seeds_per_taxon=5)


def test_source_failure_does_not_abort_run():
    def flaky(url, params):
        if "openalex" in url:
            raise RuntimeError("network down")
        if "crossref" in url:
            return CROSSREF_PAYLOAD
        return UNPAYWALL_PAYLOAD

    taxa = pd.DataFrame({"accepted_species": ["Adenia lobata"]})
    frame, report = discover_sources(
        taxa, flaky, contact_email="t@example.org", max_taxa=5, max_seeds_per_taxon=5
    )
    assert report["n_lookup_errors"] >= 1
    assert report["n_literature_seeds"] >= 1  # crossref still produced a lead


def test_source_grade_hint_maps_to_study_vocabulary():
    # Scholarly work with venue + DOI -> A.
    assert grade_source_hint("openalex", "Phytotaxa", "10.1/x", "https://oa.example/x")[0] == GRADE_A
    # Curated/institutional host -> B (even if scholarly signal is weak).
    assert grade_source_hint("crossref", "", "10.2/y", "https://powo.science.kew.org/taxon/1")[0] == GRADE_B
    # Unvetted web (Wikipedia) -> D, regardless of a DOI.
    assert grade_source_hint("openalex", "Journal", "10.3/z", "https://es.wikipedia.org/wiki/X")[0] == GRADE_D
    # Scholarly DOI but no venue -> C.
    assert grade_source_hint("crossref", "", "10.4/w", "")[0] == GRADE_C
    # Nothing usable -> below C, sorted last.
    assert grade_source_hint("", "", "", "")[1] == 5


def test_discovery_ranks_abc_above_unvetted_web():
    # Two seeds for one taxon: an A-grade journal work and a D-grade wiki page.
    payload = {
        "results": [
            {
                "display_name": "A flora treatment",
                "doi": "https://doi.org/10.1/aaa",
                "publication_year": 2018,
                "open_access": {"is_oa": True, "oa_status": "green", "oa_url": "https://oa.example/a"},
                "primary_location": {"source": {"display_name": "Kew Bulletin"}, "pdf_url": ""},
                "authorships": [],
            },
            {
                "display_name": "Wikipedia summary",
                "doi": "",
                "publication_year": 2020,
                "open_access": {"is_oa": True, "oa_status": "bronze", "oa_url": "https://en.wikipedia.org/wiki/X"},
                "primary_location": {"source": {"display_name": ""}, "pdf_url": ""},
                "authorships": [],
            },
        ]
    }
    getter = _fake_getter({"openalex": payload, "crossref": {}, "unpaywall": {}})
    taxa = pd.DataFrame({"accepted_species": ["Genus species"]})
    frame, report = discover_sources(taxa, getter, "t@example.org", max_taxa=5, max_seeds_per_taxon=5)

    # A-grade lead is ranked first; D-grade wiki lead is last.
    assert frame.iloc[0]["provisional_source_reliability_hint"] == GRADE_A
    assert frame.iloc[-1]["provisional_source_reliability_hint"] == GRADE_D
    assert list(frame["priority_rank"]) == sorted(frame["priority_rank"])
    assert report["n_hint_track_eligible_A_B_C"] >= 1
    assert report["n_hint_D_unvetted_web"] == 1
    # Grade hint is about the source, never a finalized biological value.
    assert PROHIBITED_OUTPUT_COLUMNS.isdisjoint(frame.columns)


def test_nomination_gate_blocks_unnominated_island(tmp_path):
    report_path = tmp_path / "report.json"
    report_path.write_text(json.dumps({"eligible_island_ids": ["isl_jp"]}), encoding="utf-8")
    # Eligible island passes.
    _require_nominated_island(report_path, "isl_jp")
    # Non-eligible island is refused.
    with pytest.raises(typer.BadParameter):
        _require_nominated_island(report_path, "isl_cape")
