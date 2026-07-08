import json

import pandas as pd
import pytest
import typer

from island_v2.trait_source_discovery import (
    GRADE_A,
    GRADE_D,
    OUTPUT_COLUMNS,
    PROHIBITED_OUTPUT_COLUMNS,
    TRAIT_GROUP_QUERIES,
    _require_nominated_island,
    discover_sources,
    grade_source_hint,
    openalex_abstract,
    parse_crossref_works,
    parse_openalex_works,
)

GROUP_A = "A_floral_morphology_colour"
GROUP_C = "C_reproductive_assurance"
CFG_A = TRAIT_GROUP_QUERIES[GROUP_A]
CFG_C = TRAIT_GROUP_QUERIES[GROUP_C]


def _fake_getter(responses):
    """Return canned JSON keyed by URL substring (network-free)."""

    def getter(url, params):
        for key, payload in responses.items():
            if key in url:
                return payload(url, params) if callable(payload) else payload
        return {}

    return getter


def _openalex_work(title, doi="10.1/x", venue="Journal", oa=True, oa_url="https://oa/x", abstract_inv=None):
    return {
        "display_name": title,
        "doi": f"https://doi.org/{doi}" if doi else None,
        "publication_year": 2019,
        "open_access": {"is_oa": oa, "oa_status": "green", "oa_url": oa_url},
        "primary_location": {"source": {"display_name": venue}, "pdf_url": ""},
        "authorships": [{"author": {"display_name": "A. Botanist"}}],
        "abstract_inverted_index": abstract_inv,
    }


def test_parse_openalex_scores_trait_relevance_and_never_emits_trait_values():
    # A floral-morphology title is relevant to group A; the group is tagged.
    payload = {"results": [_openalex_work("Floral morphology and corolla of the genus")]}
    seeds = parse_openalex_works(payload, "Genus species", limit=5, trait_group=GROUP_A, group_cfg=CFG_A)
    assert len(seeds) == 1
    seed = seeds[0]
    assert seed.query_trait_group == GROUP_A
    assert seed.query_template.startswith('"Genus species"')
    assert seed.title_relevance_score >= 2  # "floral", "corolla", "morpholog"
    assert "corolla" in seed.trait_keywords_matched
    assert seed.likely_evidence_type == "floral_morphology_or_colour"
    assert seed.lead_status == "unreviewed_literature_seed"
    assert seed.candidate_page_locator.startswith("HUMAN_REVIEW_REQUIRED")


def test_offtopic_paper_scores_zero_relevance():
    # High-quality source, but irrelevant to floral traits -> title_relevance 0.
    payload = {"results": [_openalex_work("Micropropagation and tissue culture protocol")]}
    seeds = parse_openalex_works(payload, "Genus species", limit=5, trait_group=GROUP_A, group_cfg=CFG_A)
    assert seeds[0].title_relevance_score == 0
    assert seeds[0].trait_keywords_matched == ""


def test_openalex_abstract_reconstructs_from_inverted_index():
    inv = {"autonomous": [0], "selfing": [1], "observed": [2]}
    payload = {"results": [_openalex_work("A study", abstract_inv=inv)]}
    seeds = parse_openalex_works(payload, "Genus species", limit=5, trait_group=GROUP_C, group_cfg=CFG_C)
    # abstract keywords for group C matched even though the title has none.
    assert seeds[0].title_relevance_score == 0
    assert seeds[0].abstract_relevance_score >= 1
    assert openalex_abstract({"abstract_inverted_index": inv}) == "autonomous selfing observed"


def test_parse_crossref_extracts_year_and_venue_and_strips_abstract_markup():
    payload = {
        "message": {
            "items": [
                {
                    "title": ["Herkogamy and mating system"],
                    "DOI": "10.5/y",
                    "published": {"date-parts": [[2021, 3]]},
                    "container-title": ["Phytotaxa"],
                    "abstract": "<jats:p>self-incompatibility was tested</jats:p>",
                }
            ]
        }
    }
    seeds = parse_crossref_works(payload, "Adenia lobata", limit=5, trait_group=GROUP_C, group_cfg=CFG_C)
    assert seeds[0].publication_year == "2021"
    assert seeds[0].venue == "Phytotaxa"
    assert seeds[0].title_relevance_score >= 1  # "herkogam", "mating system"
    assert seeds[0].abstract_relevance_score >= 1  # "self-incompatib"


def test_source_grade_hint_maps_to_study_vocabulary():
    assert grade_source_hint("openalex", "Phytotaxa", "10.1/x", "https://oa/x")[0] == GRADE_A
    assert grade_source_hint("crossref", "", "10.2/y", "https://powo.science.kew.org/taxon/1")[0] == "B_curated_database_or_institution"
    assert grade_source_hint("openalex", "Journal", "10.3/z", "https://es.wikipedia.org/wiki/X")[0] == GRADE_D
    assert grade_source_hint("", "", "", "")[1] == 5


def test_discover_sources_runs_three_trait_groups_with_only_lead_columns():
    getter = _fake_getter(
        {
            "openalex": {"results": [_openalex_work("Floral morphology of the species")]},
            "crossref": {},
            "unpaywall": {"is_oa": True, "oa_status": "gold", "best_oa_location": {"url_for_pdf": "https://x.pdf"}},
        }
    )
    taxa = pd.DataFrame({"accepted_species": ["Genus species"], "island_id": ["isl_jp"]})
    frame, report = discover_sources(taxa, getter, "t@example.org", max_taxa=5, max_seeds_per_taxon=5)

    assert PROHIBITED_OUTPUT_COLUMNS.isdisjoint(frame.columns)
    for col in OUTPUT_COLUMNS:
        assert col in frame.columns
    # One openalex seed per trait group -> all three groups represented.
    assert set(frame["query_trait_group"]) == set(TRAIT_GROUP_QUERIES)
    assert report["n_seeds_group_A_floral_morphology_colour"] == 1
    assert report["n_seeds_group_C_reproductive_assurance"] == 1
    assert "island_id" in frame.columns


def test_zero_relevance_leads_are_not_ranked_first():
    # Same taxon, group A: a relevant floral paper and an off-topic genomics paper.
    payload = {
        "results": [
            _openalex_work("Genome assembly and SNP markers", doi="10.9/z", venue="Genomics J"),
            _openalex_work("Floral morphology, corolla and inflorescence", doi="10.1/a", venue="Kew Bulletin"),
        ]
    }
    getter = _fake_getter({"openalex": payload, "crossref": {}, "unpaywall": {}})
    taxa = pd.DataFrame({"accepted_species": ["Genus species"]})
    frame, report = discover_sources(taxa, getter, "t@example.org", max_taxa=5, max_seeds_per_taxon=5)

    group_a = frame.loc[frame["query_trait_group"] == "A_floral_morphology_colour"].reset_index(drop=True)
    # The floral paper (relevance > 0) heads the group; the genomics paper (0) is last.
    assert group_a.iloc[0]["title_relevance_score"] > 0
    assert "Floral morphology" in group_a.iloc[0]["title"]
    assert int(group_a.iloc[-1]["title_relevance_score"]) == 0
    assert report["n_zero_title_relevance_seeds"] >= 1


def test_discover_sources_rejects_prohibited_input_column():
    getter = _fake_getter({})
    taxa = pd.DataFrame({"accepted_species": ["X y"], "trait_value": ["white"]})
    with pytest.raises(typer.BadParameter):
        discover_sources(taxa, getter, "t@example.org", max_taxa=5, max_seeds_per_taxon=5)


def test_source_failure_does_not_abort_run():
    def flaky(url, params):
        if "openalex" in url:
            raise RuntimeError("network down")
        if "crossref" in url:
            return {"message": {"items": [{"title": ["Pollination biology"], "DOI": "10.5/y", "container-title": ["J"]}]}}
        return {}

    taxa = pd.DataFrame({"accepted_species": ["Adenia lobata"]})
    frame, report = discover_sources(taxa, flaky, "t@example.org", max_taxa=5, max_seeds_per_taxon=5)
    assert report["n_lookup_errors"] >= 1
    assert report["n_literature_seeds"] >= 1  # crossref still produced leads


def test_nomination_gate_blocks_unnominated_island(tmp_path):
    report_path = tmp_path / "report.json"
    report_path.write_text(json.dumps({"eligible_island_ids": ["isl_jp"]}), encoding="utf-8")
    _require_nominated_island(report_path, "isl_jp")
    with pytest.raises(typer.BadParameter):
        _require_nominated_island(report_path, "isl_cape")


def test_parse_scores_taxon_relevance_from_title():
    # Title naming the queried binomial -> strong taxon relevance.
    payload = {"results": [_openalex_work("Floral morphology of Genus species in situ")]}
    seeds = parse_openalex_works(payload, "Genus species", limit=5, trait_group=GROUP_A, group_cfg=CFG_A)
    assert seeds[0].taxon_relevance_score == 4
    assert seeds[0].taxon_match_kind == "binomial_title"


def test_trait_relevant_but_off_species_lead_scores_zero_taxon():
    # Matches floral trait keywords but is about a DIFFERENT taxon -> taxon gate = 0.
    payload = {"results": [_openalex_work("Floral morphology and corolla of Dendrobium nobile")]}
    seeds = parse_openalex_works(payload, "Genus species", limit=5, trait_group=GROUP_A, group_cfg=CFG_A)
    assert seeds[0].title_relevance_score > 0          # trait-relevant
    assert seeds[0].taxon_relevance_score == 0         # but not about the target species
    assert seeds[0].taxon_match_kind == "none"


def test_taxon_relevant_lead_outranks_trait_relevant_off_species():
    payload = {
        "results": [
            # trait-relevant but off-species (wrong taxon)
            _openalex_work("Floral morphology and corolla of Dendrobium nobile", doi="10.9/off"),
            # names the target species
            _openalex_work("Pollination and corolla of Genus species", doi="10.1/on"),
        ]
    }
    getter = _fake_getter({"openalex": payload, "crossref": {}, "unpaywall": {}})
    taxa = pd.DataFrame({"accepted_species": ["Genus species"]})
    frame, report = discover_sources(taxa, getter, "t@example.org", max_taxa=5, max_seeds_per_taxon=5)

    group_a = frame.loc[frame["query_trait_group"] == "A_floral_morphology_colour"].reset_index(drop=True)
    # The on-species lead heads the group even though both are trait-relevant.
    assert "Genus species" in group_a.iloc[0]["title"]
    assert int(group_a.iloc[0]["taxon_relevance_score"]) > 0
    assert int(group_a.iloc[-1]["taxon_relevance_score"]) == 0
    assert report["n_taxon_relevant_seeds"] >= 1
    assert report["n_zero_taxon_relevance_seeds"] >= 1
