"""Free literature source discovery for staged Core-pilot taxa.

Given a staged-taxa CSV for an island that has **passed Core-pilot nomination**,
this collects only open bibliographic *leads* for each taxon so a human curator
has somewhere to start reading:

- literature seeds  — OpenAlex / Crossref works (title, DOI, year, venue, authors);
- open-access receipts — Unpaywall / OpenAlex OA status and a best open URL;
- candidate page locators — a best-guess open PDF URL plus an explicit
  "a human must open this and find the passage" note.

It NEVER emits a trait value, a native / establishment status, a Bombus
applicability decision, or an analysis-inclusion flag. Every row it writes is an
unreviewed M0 lead. It uses only free public APIs (OpenAlex, Crossref,
Unpaywall) — there are no paid LLM calls.

Operationally it refuses to run unless the target island appears in a Core-pilot
nomination report's ``eligible_island_ids``, and it is bounded to a small number
of staged taxa (default 5). This keeps it a post-nomination, human-in-the-loop
tool rather than an unbounded scraper.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

# A JSON getter: (url, params) -> parsed JSON dict. Injected so the pure
# discovery logic is testable without network access.
JsonGetter = Callable[[str, dict[str, Any]], dict[str, Any]]

OPENALEX_WORKS_URL = "https://api.openalex.org/works"
CROSSREF_WORKS_URL = "https://api.crossref.org/works"
UNPAYWALL_URL = "https://api.unpaywall.org/v2"

REQUIRED_TAXA_COLUMNS = {"accepted_species"}
# The only columns from the staged-taxa CSV we carry through as provenance.
PASSTHROUGH_TAXA_COLUMNS = ["island_id", "genus", "family", "trait_candidate_status"]

# Every lead is explicitly unreviewed and non-authoritative.
LEAD_STATUS = "unreviewed_literature_seed"
PAGE_LOCATOR_PLACEHOLDER = "HUMAN_REVIEW_REQUIRED: open source and locate the trait passage"

# Provisional source-tier HINT, using the study's own reliability vocabulary
# (see prompts/trait_evidence_extraction_v2.md and config/attrition_audit.yml).
# It grades the likely reliability of the SOURCE, never a trait value, and a
# human confirms the real grade on review. The analysis evidence tracks accept
# A/B (direct-conservative) or A/B/C (direct-broad-web); D is below every track,
# so leads are ranked to surface A/B/C first and push unvetted web last.
GRADE_A = "A_primary_or_monograph"
GRADE_B = "B_curated_database_or_institution"
GRADE_C = "C_curated_specialist_web"
GRADE_D = "D_unvetted_web"
GRADE_NONE = "none"
_GRADE_PRIORITY = {GRADE_A: 1, GRADE_B: 2, GRADE_C: 3, GRADE_D: 4, GRADE_NONE: 5}
# Curated / institutional hosts → likely B. Community / unvetted hosts → D.
_B_HOST_HINTS = (
    "powo.science.kew", "kew.org", "tropicos.org", "gbif.org", "ipni.org",
    "worldfloraonline.org", "catalogueoflife", "efloras.org", "biodiversitylibrary.org",
    "plants.usda.gov", "jstor.org/stable/community.", "globalplants", "ala.org.au",
    "herbari", "botanicgarden", "naturalis", ".gov/", "floraofaustralia",
)
_D_HOST_HINTS = (
    "wikipedia.org", "wikimedia.org", "wikispecies", "blogspot", "wordpress",
    "pinterest", "researchgate.net", "facebook.com", "reddit.com",
)


def grade_source_hint(source_api: str, venue: str, doi: str, *urls: str) -> tuple[str, int]:
    """Return a provisional (source_reliability_hint, priority_rank) for a lead.

    This is a navigation aid for reviewers, not a decision: it never inspects or
    emits a trait value, and the human reviewer sets the authoritative grade.
    """
    haystack = " ".join(u for u in (venue, *urls) if u).lower()
    if any(hint in haystack for hint in _D_HOST_HINTS):
        return GRADE_D, _GRADE_PRIORITY[GRADE_D]
    if any(hint in haystack for hint in _B_HOST_HINTS):
        return GRADE_B, _GRADE_PRIORITY[GRADE_B]
    scholarly = str(source_api).lower() in {"openalex", "crossref"}
    if scholarly and doi and venue.strip():
        return GRADE_A, _GRADE_PRIORITY[GRADE_A]
    if scholarly and doi:
        return GRADE_C, _GRADE_PRIORITY[GRADE_C]
    if haystack:
        return GRADE_C, _GRADE_PRIORITY[GRADE_C]
    return GRADE_NONE, _GRADE_PRIORITY[GRADE_NONE]

OUTPUT_COLUMNS = [
    "query_taxon",
    "source_api",
    "seed_type",
    "title",
    "doi",
    "publication_year",
    "venue",
    "authors",
    "is_open_access",
    "oa_status",
    "open_access_url",
    "candidate_pdf_url",
    "candidate_page_locator",
    "provisional_source_reliability_hint",
    "priority_rank",
    "lead_status",
    "provenance_note",
]
# Anything that would smuggle a finalized value into the output is forbidden.
PROHIBITED_OUTPUT_COLUMNS = {
    "flower_primary_color",
    "floral_symmetry",
    "floral_form",
    "trait",
    "trait_value",
    "trait_state",
    "native_status",
    "establishment_means",
    "is_native",
    "bombus_applicability",
    "applicability",
    "analysis_included",
    "accepted",
    "review_status_accepted",
}


@app.callback()
def main() -> None:
    """Discover free literature leads for staged Core-pilot taxa (no trait values)."""


@dataclass
class LiteratureSeed:
    """One unreviewed bibliographic lead for a taxon."""

    query_taxon: str
    source_api: str
    seed_type: str = "literature_seed"
    title: str = ""
    doi: str = ""
    publication_year: str = ""
    venue: str = ""
    authors: str = ""
    is_open_access: str = ""
    oa_status: str = ""
    open_access_url: str = ""
    candidate_pdf_url: str = ""
    candidate_page_locator: str = PAGE_LOCATOR_PLACEHOLDER
    provisional_source_reliability_hint: str = GRADE_NONE
    priority_rank: int = _GRADE_PRIORITY[GRADE_NONE]
    lead_status: str = LEAD_STATUS
    provenance_note: str = ""
    passthrough: dict[str, str] = field(default_factory=dict)

    def apply_source_grade_hint(self) -> None:
        """Set the provisional source-tier hint and priority from resolved fields."""
        self.provisional_source_reliability_hint, self.priority_rank = grade_source_hint(
            self.source_api, self.venue, self.doi, self.open_access_url, self.candidate_pdf_url
        )

    def to_row(self) -> dict[str, str]:
        row = {k: v for k, v in asdict(self).items() if k != "passthrough"}
        row.update(self.passthrough)
        return row


def _clean_doi(raw: Any) -> str:
    if not raw:
        return ""
    doi = str(raw).strip()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if doi.lower().startswith(prefix):
            doi = doi[len(prefix):]
    return doi.strip().lower()


def parse_openalex_works(payload: dict[str, Any], query_taxon: str, limit: int) -> list[LiteratureSeed]:
    """Parse an OpenAlex /works response into unreviewed literature seeds."""
    seeds: list[LiteratureSeed] = []
    for work in (payload.get("results") or [])[:limit]:
        oa = work.get("open_access") or {}
        location = work.get("primary_location") or work.get("best_oa_location") or {}
        source = location.get("source") or {}
        authorships = work.get("authorships") or []
        authors = "; ".join(
            (a.get("author") or {}).get("display_name", "") for a in authorships if a.get("author")
        )
        is_oa = oa.get("is_oa")
        seeds.append(
            LiteratureSeed(
                query_taxon=query_taxon,
                source_api="openalex",
                title=str(work.get("display_name") or work.get("title") or "").strip(),
                doi=_clean_doi(work.get("doi")),
                publication_year=str(work.get("publication_year") or "").strip(),
                venue=str(source.get("display_name") or "").strip(),
                authors=authors,
                is_open_access="" if is_oa is None else str(bool(is_oa)).lower(),
                oa_status=str(oa.get("oa_status") or "").strip(),
                open_access_url=str(oa.get("oa_url") or "").strip(),
                candidate_pdf_url=str(location.get("pdf_url") or "").strip(),
            )
        )
    return seeds


def parse_crossref_works(payload: dict[str, Any], query_taxon: str, limit: int) -> list[LiteratureSeed]:
    """Parse a Crossref /works response into unreviewed literature seeds."""
    items = ((payload.get("message") or {}).get("items")) or []
    seeds: list[LiteratureSeed] = []
    for item in items[:limit]:
        title_list = item.get("title") or []
        container = item.get("container-title") or []
        parts = (item.get("published") or item.get("issued") or {}).get("date-parts") or [[]]
        year = ""
        if parts and parts[0]:
            year = str(parts[0][0])
        authors = "; ".join(
            f"{a.get('given', '')} {a.get('family', '')}".strip() for a in (item.get("author") or [])
        )
        seeds.append(
            LiteratureSeed(
                query_taxon=query_taxon,
                source_api="crossref",
                title=str(title_list[0]).strip() if title_list else "",
                doi=_clean_doi(item.get("DOI")),
                publication_year=year,
                venue=str(container[0]).strip() if container else "",
                authors=authors,
            )
        )
    return seeds


def apply_unpaywall(seed: LiteratureSeed, payload: dict[str, Any]) -> LiteratureSeed:
    """Attach an Unpaywall open-access receipt to a seed (no trait inference)."""
    if not payload:
        return seed
    is_oa = payload.get("is_oa")
    seed.is_open_access = "" if is_oa is None else str(bool(is_oa)).lower()
    seed.oa_status = str(payload.get("oa_status") or seed.oa_status).strip()
    best = payload.get("best_oa_location") or {}
    if best:
        seed.open_access_url = str(best.get("url") or seed.open_access_url).strip()
        seed.candidate_pdf_url = str(best.get("url_for_pdf") or seed.candidate_pdf_url).strip()
    return seed


def discover_for_taxon(
    taxon: str,
    getter: JsonGetter,
    contact_email: str,
    max_seeds: int,
    passthrough: dict[str, str] | None = None,
) -> tuple[list[LiteratureSeed], list[str]]:
    """Collect open literature leads for one taxon. Returns (seeds, errors)."""
    passthrough = passthrough or {}
    seeds: list[LiteratureSeed] = []
    errors: list[str] = []

    try:
        oa_payload = getter(
            OPENALEX_WORKS_URL,
            {"search": taxon, "per_page": max_seeds, "mailto": contact_email},
        )
        seeds.extend(parse_openalex_works(oa_payload, taxon, max_seeds))
    except Exception as exc:  # noqa: BLE001 - one source failing must not abort the taxon
        errors.append(f"openalex:{taxon}:{exc}")

    if len(seeds) < max_seeds:
        try:
            cr_payload = getter(
                CROSSREF_WORKS_URL,
                {"query.bibliographic": taxon, "rows": max_seeds, "mailto": contact_email},
            )
            seeds.extend(parse_crossref_works(cr_payload, taxon, max_seeds - len(seeds)))
        except Exception as exc:  # noqa: BLE001
            errors.append(f"crossref:{taxon}:{exc}")

    for seed in seeds:
        seed.passthrough = dict(passthrough)
        if not seed.doi:
            seed.provenance_note = "No DOI resolved; open-access status not checked."
        else:
            try:
                up_payload = getter(f"{UNPAYWALL_URL}/{seed.doi}", {"email": contact_email})
                apply_unpaywall(seed, up_payload)
            except Exception as exc:  # noqa: BLE001
                errors.append(f"unpaywall:{seed.doi}:{exc}")
                seed.provenance_note = "Unpaywall lookup failed; open-access status unknown."
        # Grade the source tier only after any Unpaywall URLs are resolved.
        seed.apply_source_grade_hint()

    return seeds, errors


def _guard_output(frame: pd.DataFrame) -> None:
    leaked = PROHIBITED_OUTPUT_COLUMNS.intersection(frame.columns)
    if leaked:
        raise typer.BadParameter(
            "trait source discovery output must never contain finalized-value columns: "
            f"{', '.join(sorted(leaked))}"
        )


def discover_sources(
    staged_taxa: pd.DataFrame,
    getter: JsonGetter,
    contact_email: str,
    max_taxa: int,
    max_seeds_per_taxon: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build the literature-lead table and a run report for staged taxa."""
    missing = REQUIRED_TAXA_COLUMNS.difference(staged_taxa.columns)
    if missing:
        raise typer.BadParameter(f"staged taxa table missing columns: {', '.join(sorted(missing))}")
    _guard_output(staged_taxa)

    taxa = staged_taxa.copy()
    taxa["accepted_species"] = taxa["accepted_species"].fillna("").astype(str).str.strip()
    taxa = taxa.loc[taxa["accepted_species"] != ""].drop_duplicates("accepted_species")
    if len(taxa) > max_taxa:
        taxa = taxa.head(max_taxa)

    passthrough_cols = [c for c in PASSTHROUGH_TAXA_COLUMNS if c in taxa.columns]
    rows: list[dict[str, str]] = []
    all_errors: list[str] = []
    for record in taxa.to_dict("records"):
        taxon = str(record["accepted_species"]).strip()
        passthrough = {c: str(record.get(c, "") or "").strip() for c in passthrough_cols}
        seeds, errors = discover_for_taxon(
            taxon, getter, contact_email, max_seeds_per_taxon, passthrough
        )
        all_errors.extend(errors)
        rows.extend(seed.to_row() for seed in seeds)

    ordered = OUTPUT_COLUMNS + passthrough_cols
    frame = pd.DataFrame(rows, columns=ordered) if rows else pd.DataFrame(columns=ordered)
    _guard_output(frame)
    if len(frame):
        # Surface A/B/C sources before D (unvetted web) within each taxon.
        frame["priority_rank"] = pd.to_numeric(frame["priority_rank"], errors="coerce").fillna(5).astype(int)
        frame = frame.sort_values(["query_taxon", "priority_rank"]).reset_index(drop=True)

    def _n_hint(grade: str) -> int:
        return int((frame["provisional_source_reliability_hint"] == grade).sum()) if len(frame) else 0

    report = {
        "note": (
            "Unreviewed literature leads only. No trait value, native/establishment "
            "status, Bombus applicability, or analysis inclusion is decided here. The "
            "provisional_source_reliability_hint is a navigation aid; a human sets the "
            "authoritative grade. Analysis tracks accept A/B (direct-conservative) or "
            "A/B/C (direct-broad-web); D_unvetted_web is below every track."
        ),
        "n_taxa_queried": int(len(taxa)),
        "n_literature_seeds": int(len(frame)),
        "n_seeds_with_doi": int((frame["doi"] != "").sum()) if len(frame) else 0,
        "n_open_access_seeds": int((frame["is_open_access"] == "true").sum()) if len(frame) else 0,
        "n_hint_A_primary_or_monograph": _n_hint(GRADE_A),
        "n_hint_B_curated_or_institution": _n_hint(GRADE_B),
        "n_hint_C_specialist_web": _n_hint(GRADE_C),
        "n_hint_D_unvetted_web": _n_hint(GRADE_D),
        "n_hint_track_eligible_A_B_C": _n_hint(GRADE_A) + _n_hint(GRADE_B) + _n_hint(GRADE_C),
        "n_lookup_errors": len(all_errors),
        "lookup_errors": all_errors[:20],
    }
    return frame, report


def _require_nominated_island(nomination_report: Path, island_id: str) -> None:
    report = json.loads(nomination_report.read_text(encoding="utf-8"))
    eligible = report.get("eligible_island_ids") or []
    if island_id not in eligible:
        raise typer.BadParameter(
            f"island_id {island_id!r} is not in the Core-pilot nomination report's "
            f"eligible_island_ids ({eligible}); run source discovery only after the "
            "island passes nomination."
        )


def _httpx_getter() -> JsonGetter:
    import httpx

    client = httpx.Client(
        timeout=30.0,
        headers={"User-Agent": "island-floral-v2 trait-source-discovery (mailto provided per request)"},
        follow_redirects=True,
    )

    def getter(url: str, params: dict[str, Any]) -> dict[str, Any]:
        response = client.get(url, params=params)
        response.raise_for_status()
        return response.json()

    return getter


@app.command("discover")
def discover(
    staged_taxa_csv: Path = typer.Option(..., exists=True, help="Staged Core-pilot taxa CSV (needs accepted_species)."),
    nomination_report: Path = typer.Option(..., exists=True, help="core_pilot_nomination_report.json for gating."),
    island_id: str = typer.Option(..., help="Nominated island_id (must be eligible in the report)."),
    contact_email: str = typer.Option(..., help="Contact email for the OpenAlex/Unpaywall polite pool."),
    output_dir: Path = typer.Option(...),
    max_taxa: int = typer.Option(5, min=1, help="Bounded number of staged taxa to query."),
    max_seeds_per_taxon: int = typer.Option(5, min=1, help="Max literature seeds per taxon."),
) -> None:
    """Write unreviewed literature leads for a nominated island's staged taxa."""
    _require_nominated_island(nomination_report, island_id)
    staged_taxa = pd.read_csv(staged_taxa_csv, dtype=str).fillna("")
    frame, report = discover_sources(
        staged_taxa, _httpx_getter(), contact_email, max_taxa, max_seeds_per_taxon
    )
    report["island_id"] = island_id
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "trait_source_leads.csv", index=False)
    (output_dir / "trait_source_discovery_report.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"{report['n_literature_seeds']} unreviewed literature lead(s) for "
        f"{report['n_taxa_queried']} taxon(s) on {island_id}. "
        "No trait/native/applicability value was decided."
    )


if __name__ == "__main__":
    app()
