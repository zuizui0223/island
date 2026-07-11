from __future__ import annotations

import pandas as pd

from island_v2.search_result_trait_candidates import extract_search_result_candidates


def test_extracts_explicit_reported_traits_and_likely_pollination() -> None:
    search = pd.DataFrame(
        [
            {
                "accepted_species": "Plantus alba",
                "searched_name": "Plantus alba",
                "search_query": '"Plantus alba" flowers',
                "search_engine": "test",
                "search_rank": "1",
                "result_title": "Plantus alba growing guide",
                "result_url": "https://example.org/alba",
                "result_snippet": (
                    "Plantus alba produces tubular white flowers and attracts pollinators."
                ),
            }
        ]
    )
    species = pd.DataFrame(
        [{"accepted_species": "Plantus alba", "family": "Lamiaceae"}]
    )

    frame, report = extract_search_result_candidates(search, species)

    values = {
        (row.trait_name, row.provisional_candidate_value, row.candidate_class)
        for row in frame.itertuples(index=False)
    }
    assert ("flower_primary_color", "white", "reported") in values
    assert ("floral_form", "tubular", "reported") in values
    assert ("pollen_vector_mode", "biotic", "likely") in values
    assert report["n_species_with_any_candidate"] == 1


def test_rejects_irrelevant_search_result_even_when_trait_words_occur() -> None:
    search = pd.DataFrame(
        [
            {
                "accepted_species": "Plantus alba",
                "searched_name": "Plantus alba",
                "search_query": '"Plantus alba" flower color',
                "search_engine": "test",
                "search_rank": "1",
                "result_title": "Orange dress ideas",
                "result_url": "https://example.org/dress",
                "result_snippet": "Beautiful orange flowers and purple patterns for summer outfits.",
            }
        ]
    )
    species = pd.DataFrame(
        [{"accepted_species": "Plantus alba", "family": "Unknownaceae"}]
    )

    frame, report = extract_search_result_candidates(search, species)

    assert frame.empty
    assert report["n_species_with_species_scoped_search_result"] == 0


def test_family_prior_is_likely_proxy_only() -> None:
    search = pd.DataFrame(
        columns=[
            "accepted_species",
            "searched_name",
            "search_query",
            "search_engine",
            "search_rank",
            "result_title",
            "result_url",
            "result_snippet",
        ]
    )
    species = pd.DataFrame(
        [
            {"accepted_species": "Gramina ventosa", "family": "Poaceae"},
            {"accepted_species": "Orchis ficta", "family": "Orchidaceae"},
        ]
    )

    frame, _ = extract_search_result_candidates(search, species)

    assert set(frame["candidate_class"]) == {"proxy"}
    assert set(frame["inference_status"]) == {"likely"}
    got = dict(zip(frame["accepted_species"], frame["provisional_candidate_value"]))
    assert got == {"Gramina ventosa": "abiotic_wind", "Orchis ficta": "biotic"}
