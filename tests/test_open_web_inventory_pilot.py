from __future__ import annotations

from pathlib import Path

from island_v2.open_web_evidence import DomainRecord, Page
from island_v2.open_web_inventory_pilot import (
    _accepted_species_from_url,
    _species_keys_from_url,
    discover_species_urls,
)


def test_species_key_supports_query_encoded_flora_treatments() -> None:
    url = (
        "https://plantnet.rbgsyd.nsw.gov.au/cgi-bin/NSWfl.pl?"
        "page=nswfl&lvl=sp&name=Acanthus~mollis"
    )
    assert _species_keys_from_url(url) == {"acanthus mollis"}


def test_inventory_backfill_never_runs_from_pr_synchronization() -> None:
    workflow = Path(".github/workflows/open-web-discovered-domain-backfill.yml").read_text(
        encoding="utf-8"
    )
    assert "backfill:\n    if: github.event_name == 'workflow_dispatch'" in workflow


def test_inventory_url_passes_exact_accepted_species_to_identity_gate() -> None:
    url = (
        "https://plantnet.rbgsyd.nsw.gov.au/cgi-bin/NSWfl.pl?"
        "page=nswfl&lvl=sp&name=Acanthus~mollis"
    )
    assert (
        _accepted_species_from_url(
            url,
            {
                "acanthus mollis": "Acanthus mollis",
                "acanthus spinosus": "Acanthus spinosus",
            },
        )
        == "Acanthus mollis"
    )


class FakeFetcher:
    def __init__(self, pages: dict[str, Page]) -> None:
        self.pages = pages
        self.calls: list[str] = []

    def fetch(self, url: str) -> Page:
        self.calls.append(url)
        return self.pages[url]


def _page(url: str, links: tuple[str, ...]) -> Page:
    return Page(
        requested_url=url,
        final_url=url,
        status_code=200,
        content_type="text/html",
        title="Index",
        text="index",
        language="en",
        retrieved_at_utc="2026-07-26T00:00:00Z",
        content_sha256="abc",
        links=links,
    )


def test_inventory_discovers_species_urls_without_fetching_treatments_twice() -> None:
    index = "https://flora.example/list"
    species_one = "https://flora.example/species/alpha/one"
    species_two = "https://flora.example/species/beta/two"
    outside = "https://elsewhere.example/species/gamma/three"
    record = DomainRecord(
        domain="flora.example",
        source_name="Example Flora",
        source_type="specialist_regional_flora",
        region="Test",
        language="en",
        organization_type="nonprofit",
        robots_access_notes="test",
        page_pattern=r"^https://flora\.example/species/[a-z-]+/[a-z-]+$",
        treatment_end_pattern="",
        inventory_urls=(index,),
        inventory_depth=1,
        scientific_name_displayed=True,
        cultivar_descriptions_present=False,
        provenance_quality="test",
        recommended_evidence_tier="B",
        review_status="approved",
        allowed_traits="all",
    )
    fetcher = FakeFetcher(
        {index: _page(index, (species_one, species_two, outside))}
    )
    urls, audit = discover_species_urls(
        record,
        fetcher=fetcher,  # type: ignore[arg-type]
        max_species_pages=10,
    )
    assert urls == [species_one, species_two]
    assert fetcher.calls == [index]
    assert len(audit) == 1


def test_inventory_respects_species_page_cap() -> None:
    index = "https://flora.example/list"
    children = tuple(
        f"https://flora.example/species/genus/species-{number}"
        for number in range(10)
    )
    record = DomainRecord(
        domain="flora.example",
        source_name="Example Flora",
        source_type="specialist_regional_flora",
        region="Test",
        language="en",
        organization_type="nonprofit",
        robots_access_notes="test",
        page_pattern=r"^https://flora\.example/species/[a-z-]+/[a-z0-9-]+$",
        treatment_end_pattern="",
        inventory_urls=(index,),
        inventory_depth=1,
        scientific_name_displayed=True,
        cultivar_descriptions_present=False,
        provenance_quality="test",
        recommended_evidence_tier="B",
        review_status="approved",
        allowed_traits="all",
    )
    urls, _ = discover_species_urls(
        record,
        fetcher=FakeFetcher({index: _page(index, children)}),  # type: ignore[arg-type]
        max_species_pages=3,
    )
    assert len(urls) == 3


def test_inventory_excludes_completed_species_pages() -> None:
    index = "https://flora.example/list"
    species_one = "https://flora.example/species/alpha/one"
    species_two = "https://flora.example/species/beta/two"
    species_three = "https://flora.example/species/gamma/three"
    record = DomainRecord(
        domain="flora.example",
        source_name="Example Flora",
        source_type="specialist_regional_flora",
        region="Test",
        language="en",
        organization_type="nonprofit",
        robots_access_notes="test",
        page_pattern=r"^https://flora\.example/species/[a-z-]+/[a-z-]+$",
        treatment_end_pattern="",
        inventory_urls=(index,),
        inventory_depth=1,
        scientific_name_displayed=True,
        cultivar_descriptions_present=False,
        provenance_quality="test",
        recommended_evidence_tier="B",
        review_status="approved",
        allowed_traits="all",
    )
    urls, _ = discover_species_urls(
        record,
        fetcher=FakeFetcher(  # type: ignore[arg-type]
            {index: _page(index, (species_one, species_two, species_three))}
        ),
        max_species_pages=2,
        completed_species_urls={species_one},
    )
    assert urls == [species_two, species_three]


def test_inventory_targets_only_unresolved_species_url_keys() -> None:
    index = "https://flora.example/list"
    species_one = "https://flora.example/species/alpha/one"
    species_two = "https://flora.example/species/beta/two"
    species_three = "https://flora.example/species/gamma/three"
    record = DomainRecord(
        domain="flora.example",
        source_name="Example Flora",
        source_type="specialist_regional_flora",
        region="Test",
        language="en",
        organization_type="nonprofit",
        robots_access_notes="test",
        page_pattern=r"^https://flora\.example/species/[a-z-]+/[a-z-]+$",
        treatment_end_pattern="",
        inventory_urls=(index,),
        inventory_depth=1,
        scientific_name_displayed=True,
        cultivar_descriptions_present=False,
        provenance_quality="test",
        recommended_evidence_tier="B",
        review_status="approved",
        allowed_traits="all",
    )
    urls, _ = discover_species_urls(
        record,
        fetcher=FakeFetcher(  # type: ignore[arg-type]
            {index: _page(index, (species_one, species_two, species_three))}
        ),
        max_species_pages=10,
        target_species_names={"Beta two", "Gamma three"},
    )
    assert urls == [species_two, species_three]


def test_inventory_targets_hyphenated_binomial_slug() -> None:
    index = "https://flora.example/sitemap.xml"
    species_one = "https://flora.example/plants/alpha-one"
    species_two = "https://flora.example/plants/beta-two"
    record = DomainRecord(
        domain="flora.example",
        source_name="Example Flora",
        source_type="university_extension",
        region="Test",
        language="en",
        organization_type="university",
        robots_access_notes="test",
        page_pattern=r"^https://flora\.example/plants/[a-z-]+$",
        treatment_end_pattern="",
        inventory_urls=(index,),
        inventory_depth=0,
        scientific_name_displayed=True,
        cultivar_descriptions_present=True,
        provenance_quality="test",
        recommended_evidence_tier="B",
        review_status="approved",
        allowed_traits="all",
    )
    urls, _ = discover_species_urls(
        record,
        fetcher=FakeFetcher(  # type: ignore[arg-type]
            {index: _page(index, (species_one, species_two))}
        ),
        max_species_pages=10,
        target_species_names={"Beta two"},
    )
    assert urls == [species_two]
