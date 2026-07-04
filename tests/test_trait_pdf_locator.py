import importlib

import pandas as pd

locator = importlib.import_module("island_v2.trait_pdf_locator")


def _leads() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "query_taxon": "Test plant",
                "island_id": "island_core",
                "doi": "10.1000/test",
                "source_api": "openalex",
                "is_open_access": "true",
                "candidate_pdf_url": "https://example.org/open.pdf",
            },
            {
                "query_taxon": "Closed plant",
                "island_id": "island_core",
                "doi": "10.1000/closed",
                "source_api": "crossref",
                "is_open_access": "false",
                "candidate_pdf_url": "https://example.org/closed.pdf",
            },
            {
                "query_taxon": "No URL plant",
                "island_id": "island_core",
                "doi": "10.1000/no-url",
                "source_api": "openalex",
                "is_open_access": "true",
                "candidate_pdf_url": "",
            },
        ]
    )


def test_public_pdf_candidates_requires_explicit_oa_and_url():
    candidates = locator.public_pdf_candidates(_leads())

    assert len(candidates) == 1
    assert candidates.iloc[0]["query_taxon"] == "Test plant"


def test_locator_records_pages_but_never_pdf_text_or_trait_value():
    candidates = locator.public_pdf_candidates(_leads())

    def fake_fetcher(url: str):
        assert url == "https://example.org/open.pdf"
        return 200, "application/pdf", b"%PDF-fake"

    def fake_scanner(payload: bytes):
        assert payload.startswith(b"%PDF")
        return 9, [{"page": 4, "matched_terms": ["flower", "pollination"], "score": 2}]

    result = locator.locate_candidates(
        (row for _, row in candidates.iterrows()),
        max_pdfs=5,
        fetcher=fake_fetcher,
        scanner=fake_scanner,
    )

    row = result.iloc[0]
    assert row["locator_status"] == "candidate_pages_located"
    assert row["candidate_pages_json"] == '[{"page": 4, "matched_terms": ["flower", "pollination"], "score": 2}]'
    assert "trait_value" not in result.columns
    assert "does not verify a trait value" in row["do_not_infer"]


def test_locator_keeps_one_failed_source_as_a_receipt():
    candidates = locator.public_pdf_candidates(_leads())

    def failing_fetcher(_url: str):
        raise RuntimeError("offline")

    result = locator.locate_candidates(
        (row for _, row in candidates.iterrows()),
        max_pdfs=5,
        fetcher=failing_fetcher,
    )

    row = result.iloc[0]
    assert row["locator_status"] == "access_or_parse_failed"
    assert "offline" in row["locator_note"]
