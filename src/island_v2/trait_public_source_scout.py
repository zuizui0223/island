"""Free source discovery for staged trait taxa.

Uses only public metadata/access endpoints (OpenAlex, Crossref, optional
Unpaywall) and standard-library HTTP. The resulting rows are discovery receipts,
never trait values or curation decisions. The design mirrors bita's separation
between seed harvesting, public-source discovery, and biological adjudication.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Iterable
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

OPENALEX_WORKS = "https://api.openalex.org/works"
CROSSREF_WORKS = "https://api.crossref.org/works/"
UNPAYWALL_WORKS = "https://api.unpaywall.org/v2/"
USER_AGENT = "island-floral-v2 public-trait-source-scout/0.1"
REQUIRED_TAXON_COLUMNS = {"accepted_species", "genus", "family"}

DOMAIN_TERMS = {
    "floral_signal": ("flower", "floral", "colour", "color", "corolla"),
    "floral_architecture": ("flower", "floral", "corolla", "inflorescence"),
    "pollen_vector": ("pollination", "pollinator", "pollen transfer"),
    "pollination": ("pollination", "pollinator", "floral visitor"),
    "reproductive_assurance": (
        "breeding system",
        "self incompatibility",
        "selfing",
        "mating system",
        "autonomous selfing",
    ),
}

SEED_FIELDS = [
    "seed_id", "accepted_species", "genus", "family", "island_id", "canonical_island_name",
    "trait_domain", "search_query", "openalex_work_id", "title", "authors", "publication_year",
    "doi", "landing_page_url", "open_access_url", "is_open_access", "cited_by_count",
    "abstract_available", "metadata_term_hits", "seed_status", "seed_warning",
    "trait_candidate_status", "release_gate",
]
RECEIPT_FIELDS = [
    "seed_id", "accepted_species", "doi", "provider", "source_kind", "resolution_status",
    "request_url", "source_identifier", "title", "landing_page_url", "content_url", "content_type",
    "license_label", "relation_to_article", "notes",
]


@dataclass(frozen=True)
class SourceReceipt:
    provider: str
    source_kind: str
    resolution_status: str
    request_url: str
    source_identifier: str
    title: str
    landing_page_url: str
    content_url: str
    content_type: str
    license_label: str
    relation_to_article: str
    notes: str


@dataclass
class LiteratureSeed:
    seed_id: str
    accepted_species: str
    genus: str
    family: str
    island_id: str
    canonical_island_name: str
    trait_domain: str
    search_query: str
    openalex_work_id: str
    title: str
    authors: str
    publication_year: str
    doi: str
    landing_page_url: str
    open_access_url: str
    is_open_access: bool
    cited_by_count: int
    abstract_available: bool
    metadata_term_hits: str
    trait_candidate_status: str
    release_gate: str
    seed_status: str = "M0_source_candidate_needs_content_screen"
    seed_warning: str = (
        "Metadata and public-access links are discovery leads only; they do not establish a trait value, "
        "taxon scope, island establishment, or analysis eligibility."
    )

    def row(self) -> dict[str, object]:
        payload = asdict(self)
        return {key: payload[key] for key in SEED_FIELDS}


def _text(value: object) -> str:
    return str(value or "").strip()


def _url(value: object) -> str:
    value = _text(value)
    return value if value.startswith(("https://", "http://")) else ""


def normalise_doi(value: object) -> str:
    doi = _text(value).lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if doi.startswith(prefix):
            return doi[len(prefix) :].strip()
    return doi


def request_json(url: str, params: dict[str, str] | None = None) -> tuple[int, Any]:
    request_url = f"{url}?{urlencode(params)}" if params else url
    request = Request(request_url, headers={"Accept": "application/json", "User-Agent": USER_AGENT})
    with urlopen(request, timeout=45) as response:  # nosec B310: fixed public metadata endpoints
        return int(getattr(response, "status", response.getcode())), json.loads(response.read().decode("utf-8"))


def _authors(work: dict[str, Any]) -> str:
    names: list[str] = []
    for authorship in work.get("authorships") or []:
        author = authorship.get("author") if isinstance(authorship, dict) else None
        name = author.get("display_name") if isinstance(author, dict) else None
        if isinstance(name, str) and name.strip():
            names.append(name.strip())
        if len(names) == 5:
            break
    return "; ".join(names)


def reconstruct_abstract(work: dict[str, Any]) -> str:
    inverted = work.get("abstract_inverted_index")
    if not isinstance(inverted, dict):
        return ""
    placed: list[tuple[int, str]] = []
    for word, positions in inverted.items():
        if not isinstance(word, str) or not isinstance(positions, list):
            continue
        for position in positions:
            if isinstance(position, int) and position >= 0:
                placed.append((position, word))
    return " ".join(word for _, word in sorted(placed))


def _locations(work: dict[str, Any]) -> tuple[str, str, bool]:
    primary = work.get("primary_location") if isinstance(work.get("primary_location"), dict) else {}
    best_oa = work.get("best_oa_location") if isinstance(work.get("best_oa_location"), dict) else {}
    open_access = work.get("open_access") if isinstance(work.get("open_access"), dict) else {}
    return (
        _url(primary.get("landing_page_url") or primary.get("pdf_url")),
        _url(best_oa.get("pdf_url") or best_oa.get("landing_page_url")),
        bool(open_access.get("is_oa")),
    )


def trait_domains(ontology_path: Path) -> list[str]:
    document = yaml.safe_load(ontology_path.read_text(encoding="utf-8")) or {}
    traits = document.get("traits")
    if not isinstance(traits, dict) or not traits:
        raise typer.BadParameter("Trait ontology must contain a traits mapping")
    return sorted({
        _text(spec.get("domain"))
        for spec in traits.values()
        if isinstance(spec, dict) and _text(spec.get("domain")) in DOMAIN_TERMS
    })


def build_queries(taxa: pd.DataFrame, domains: Iterable[str]) -> list[dict[str, str]]:
    missing = REQUIRED_TAXON_COLUMNS.difference(taxa.columns)
    if missing:
        raise typer.BadParameter(f"Taxon input missing columns: {', '.join(sorted(missing))}")
    queries: list[dict[str, str]] = []
    for taxon in taxa.fillna("").drop_duplicates("accepted_species").to_dict("records"):
        species = _text(taxon.get("accepted_species"))
        if not species:
            continue
        for domain in domains:
            terms = DOMAIN_TERMS[domain]
            queries.append({
                "accepted_species": species,
                "genus": _text(taxon.get("genus")),
                "family": _text(taxon.get("family")),
                "island_id": _text(taxon.get("island_id")),
                "canonical_island_name": _text(taxon.get("canonical_island_name")),
                "trait_domain": domain,
                "search_query": f'"{species}" ({" OR ".join(terms)})',
                "trait_candidate_status": _text(taxon.get("trait_candidate_status")) or "staged_not_curated",
                "release_gate": _text(taxon.get("release_gate")),
            })
    return queries


def seed_from_work(work: dict[str, Any], query: dict[str, str], index: int) -> LiteratureSeed:
    title = _text(work.get("display_name") or work.get("title"))
    abstract = reconstruct_abstract(work)
    hits = [term for term in DOMAIN_TERMS[query["trait_domain"]] if term in f"{title} {abstract}".lower()]
    landing, open_access_url, is_oa = _locations(work)
    openalex_id = _text(work.get("id"))
    suffix = openalex_id.rsplit("/", 1)[-1] if openalex_id else str(index)
    return LiteratureSeed(
        seed_id=f"{query['accepted_species'].replace(' ', '_')}:{query['trait_domain']}:{suffix}",
        accepted_species=query["accepted_species"], genus=query["genus"], family=query["family"],
        island_id=query["island_id"], canonical_island_name=query["canonical_island_name"],
        trait_domain=query["trait_domain"], search_query=query["search_query"], openalex_work_id=openalex_id,
        title=title, authors=_authors(work), publication_year=_text(work.get("publication_year")),
        doi=normalise_doi(work.get("doi")), landing_page_url=landing, open_access_url=open_access_url,
        is_open_access=is_oa, cited_by_count=int(work.get("cited_by_count") or 0),
        abstract_available=bool(abstract), metadata_term_hits=";".join(hits),
        trait_candidate_status=query["trait_candidate_status"], release_gate=query["release_gate"],
    )


def harvest_openalex(
    queries: Iterable[dict[str, str]], *, per_query: int,
    fetch_json: Callable[[str, dict[str, str] | None], tuple[int, Any]] = request_json,
) -> tuple[list[LiteratureSeed], list[dict[str, object]]]:
    seeds: dict[tuple[str, str, str], LiteratureSeed] = {}
    reports: list[dict[str, object]] = []
    for query in queries:
        params = {
            "search": query["search_query"],
            "per-page": str(per_query),
            "select": ",".join((
                "id", "display_name", "publication_year", "doi", "authorships", "cited_by_count",
                "open_access", "primary_location", "best_oa_location", "abstract_inverted_index",
            )),
        }
        try:
            status, payload = fetch_json(OPENALEX_WORKS, params)
            works = payload.get("results") if isinstance(payload, dict) else None
            if status >= 400 or not isinstance(works, list):
                raise RuntimeError(f"OpenAlex HTTP {status}")
            accepted = 0
            for index, work in enumerate(works, start=1):
                if not isinstance(work, dict):
                    continue
                seed = seed_from_work(work, query, index)
                key = (seed.accepted_species, seed.trait_domain, seed.openalex_work_id or seed.title)
                if key not in seeds:
                    seeds[key] = seed
                    accepted += 1
            reports.append({
                "accepted_species": query["accepted_species"], "trait_domain": query["trait_domain"],
                "search_query": query["search_query"], "status": "success",
                "works_returned": len(works), "usable_seeds": accepted,
            })
        except Exception as error:
            reports.append({
                "accepted_species": query["accepted_species"], "trait_domain": query["trait_domain"],
                "search_query": query["search_query"], "status": "query_failed",
                "message": f"{type(error).__name__}: {error}",
            })
    return list(seeds.values()), reports


def crossref_receipts(
    doi: str, *, fetch_json: Callable[[str, dict[str, str] | None], tuple[int, Any]] = request_json,
) -> list[SourceReceipt]:
    doi = normalise_doi(doi)
    if not doi:
        return []
    url = f"{CROSSREF_WORKS}{quote(doi, safe='')}"
    try:
        status, payload = fetch_json(url, None)
        message = payload.get("message") if isinstance(payload, dict) else None
        if status >= 400 or not isinstance(message, dict):
            raise RuntimeError(f"Crossref HTTP {status}")
        title = _text((message.get("title") or [""])[0]) if isinstance(message.get("title"), list) else _text(message.get("title"))
        landing = _url(message.get("URL")) or f"https://doi.org/{doi}"
        rows = [SourceReceipt(
            provider="Crossref", source_kind="article_metadata", resolution_status="metadata_recovered",
            request_url=url, source_identifier=_text(message.get("DOI")) or doi, title=title,
            landing_page_url=landing, content_url="", content_type="", license_label="",
            relation_to_article="exact_article_doi", notes="Metadata recovered; not full text or trait evidence.",
        )]
        for link in message.get("link") or []:
            if not isinstance(link, dict):
                continue
            content = _url(link.get("URL"))
            if content:
                rows.append(SourceReceipt(
                    provider="Crossref", source_kind="publisher_content_link", resolution_status="content_link_discovered",
                    request_url=url, source_identifier=_text(message.get("DOI")) or doi, title=title,
                    landing_page_url=landing, content_url=content, content_type=_text(link.get("content-type")),
                    license_label="", relation_to_article="exact_article_doi",
                    notes="Content link discovered; accessibility and trait content remain untested.",
                ))
        return rows
    except Exception as error:
        return [SourceReceipt(
            provider="Crossref", source_kind="article_metadata", resolution_status="access_failed",
            request_url=url, source_identifier="", title="", landing_page_url="", content_url="", content_type="",
            license_label="", relation_to_article="not_checked", notes=f"{type(error).__name__}: {error}",
        )]


def unpaywall_receipts(
    doi: str, *, email: str,
    fetch_json: Callable[[str, dict[str, str] | None], tuple[int, Any]] = request_json,
) -> list[SourceReceipt]:
    doi = normalise_doi(doi)
    if not doi or not email:
        return []
    url = f"{UNPAYWALL_WORKS}{quote(doi, safe='')}?email={quote(email, safe='@._+-')}"
    try:
        status, payload = fetch_json(url, None)
        if status >= 400 or not isinstance(payload, dict):
            raise RuntimeError(f"Unpaywall HTTP {status}")
        locations: list[dict[str, Any]] = []
        if isinstance(payload.get("best_oa_location"), dict):
            locations.append(payload["best_oa_location"])
        locations.extend(item for item in payload.get("oa_locations") or [] if isinstance(item, dict))
        rows: list[SourceReceipt] = []
        seen: set[str] = set()
        for location in locations:
            pdf, landing = _url(location.get("url_for_pdf")), _url(location.get("url"))
            chosen = pdf or landing
            if not chosen or chosen in seen:
                continue
            seen.add(chosen)
            rows.append(SourceReceipt(
                provider="Unpaywall", source_kind="open_access_location",
                resolution_status="public_fulltext_candidate" if pdf else "public_landing_candidate",
                request_url=url, source_identifier=_text(payload.get("doi")) or doi, title=_text(payload.get("title")),
                landing_page_url=landing, content_url=pdf, content_type="application/pdf" if pdf else "",
                license_label=_text(location.get("license")), relation_to_article="exact_article_doi",
                notes="Public-access candidate only; file contents remain untested.",
            ))
        return rows or [SourceReceipt(
            provider="Unpaywall", source_kind="open_access_location", resolution_status="not_found",
            request_url=url, source_identifier=_text(payload.get("doi")) or doi, title=_text(payload.get("title")),
            landing_page_url="", content_url="", content_type="", license_label="",
            relation_to_article="exact_article_doi", notes="No DOI-exact public-access candidate returned.",
        )]
    except Exception as error:
        return [SourceReceipt(
            provider="Unpaywall", source_kind="open_access_location", resolution_status="access_failed",
            request_url=url, source_identifier="", title="", landing_page_url="", content_url="", content_type="",
            license_label="", relation_to_article="not_checked", notes=f"{type(error).__name__}: {error}",
        )]


def source_receipts(
    seeds: Iterable[LiteratureSeed], *, unpaywall_email: str = "",
    fetch_json: Callable[[str, dict[str, str] | None], tuple[int, Any]] = request_json,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for seed in seeds:
        receipts = [*crossref_receipts(seed.doi, fetch_json=fetch_json), *unpaywall_receipts(seed.doi, email=unpaywall_email, fetch_json=fetch_json)]
        for receipt in receipts:
            rows.append({"seed_id": seed.seed_id, "accepted_species": seed.accepted_species, "doi": seed.doi, **asdict(receipt)})
    return rows


def write_outputs(directory: Path, queries: list[dict[str, str]], seeds: list[LiteratureSeed], reports: list[dict[str, object]], receipts: list[dict[str, object]]) -> dict[str, object]:
    directory.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(queries).to_csv(directory / "trait_source_query_registry.csv", index=False)
    pd.DataFrame([seed.row() for seed in seeds], columns=SEED_FIELDS).to_csv(directory / "trait_literature_seeds.csv", index=False)
    pd.DataFrame(reports).to_csv(directory / "trait_source_query_report.csv", index=False)
    pd.DataFrame(receipts, columns=RECEIPT_FIELDS).to_csv(directory / "trait_public_source_receipts.csv", index=False)
    result = {
        "n_query_rows": len(queries), "n_literature_seeds": len(seeds), "n_doi_exact_receipts": len(receipts),
        "n_open_access_seed_candidates": sum(bool(seed.open_access_url) for seed in seeds),
        "n_query_failures": sum(row.get("status") == "query_failed" for row in reports),
        "decision_boundary": "Seeds and receipts are literature/access leads only, never trait values or analysis inclusion decisions.",
    }
    (directory / "trait_source_scout_summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


@app.command("run")
def run(
    taxa_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    ontology_path: Path = typer.Option(Path("config/trait_ontology.yml")),
    per_query: int = typer.Option(10, min=1, max=50),
    unpaywall_email: str = typer.Option(""),
) -> None:
    """Harvest free OpenAlex/Crossref/Unpaywall literature and access leads."""
    taxa = pd.read_csv(taxa_csv, dtype=str).fillna("")
    queries = build_queries(taxa, trait_domains(ontology_path))
    seeds, reports = harvest_openalex(queries, per_query=per_query)
    receipts = source_receipts(seeds, unpaywall_email=unpaywall_email)
    typer.echo(json.dumps(write_outputs(output_dir, queries, seeds, reports, receipts), sort_keys=True))


if __name__ == "__main__":
    app()
