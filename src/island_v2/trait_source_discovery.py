"""Free, trait-group-targeted literature source discovery for staged taxa.

Given a staged-taxa CSV for an island that has **passed Core-pilot nomination**,
this collects open bibliographic *leads* - never trait values - so a human
curator has somewhere to start reading. It is a source-lead / literature-scout
stage, explicitly upstream of trait evidence review.

The retrieval is **species x trait-group**, not species-only, so a paper that is
high-quality but irrelevant to floral / pollination / reproductive traits
(micropropagation, tissue culture, nutrient content, nanoparticles, genomics,
pathology, ...) is demoted rather than surfaced. Three query groups are run:

- ``A_floral_morphology_colour``   - floral phenotype (M0);
- ``B_pollination_pollen_vector``  - pollination / pollen vector (M1);
- ``C_reproductive_assurance``     - self-incompatibility / selfing / mating.

Each lead carries the trait group it was retrieved for, the query template, a
title/abstract relevance score (keyword match against the group), the matched
keywords, a likely evidence type, and a provisional source-reliability hint.
Leads with ``title_relevance_score == 0`` are ranked below any relevant lead so
they never head a taxon's trait review.

It NEVER emits a trait value, a native / establishment status, a Bombus
applicability decision, or an analysis-inclusion flag. Every row is an
unreviewed lead. It uses only free public APIs (OpenAlex, Crossref, Unpaywall);
there are no paid LLM calls. It refuses to run unless the target island is in a
Core-pilot nomination report's ``eligible_island_ids`` and is bounded to a small
number of staged taxa.
"""

from __future__ import annotations

import json
import re
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

# --------------------------------------------------------------------------- #
# Trait-group retrieval. Each group targets one trait layer; a lead's relevance
# is scored by matching the group's keywords against the returned title/abstract.
# The `api_terms` bias the free-text query; the scoring is what actually demotes
# off-topic results, so it is robust to the APIs' loose boolean handling.
# --------------------------------------------------------------------------- #
TRAIT_GROUP_QUERIES: dict[str, dict[str, Any]] = {
    "A_floral_morphology_colour": {
        "query_template": (
            '"{species}" (flower OR floral OR corolla OR inflorescence OR morphology) '
            "(flora OR revision OR taxonomic treatment OR description)"
        ),
        "api_terms": "flower floral corolla inflorescence morphology flora revision description",
        "likely_evidence_type": "floral_morphology_or_colour",
        "keywords": (
            "flower", "floral", "corolla", "inflorescence", "perianth", "petal", "tepal",
            "morpholog", "colour", "color", "symmetr", "flora", "revision",
            "taxonomic treatment", "description", "monograph",
        ),
    },
    "B_pollination_pollen_vector": {
        "query_template": (
            '"{species}" (pollination OR pollinator OR floral visitor OR pollen vector '
            "OR breeding biology)"
        ),
        "api_terms": "pollination pollinator floral visitor pollen vector breeding biology",
        "likely_evidence_type": "pollination_or_pollen_vector",
        "keywords": (
            "pollinat", "pollen vector", "floral visitor", "flower visitor", "breeding biology",
            "pollen deposition", "pollen contact", "visitation", "anemophil", "entomophil",
            "melittophil", "ornithophil", "chiropterophil",
        ),
    },
    "C_reproductive_assurance": {
        "query_template": (
            '"{species}" ("self-incompatibility" OR "self-compatible" OR "autonomous selfing" '
            "OR bagging OR hand-pollination OR mating system OR herkogamy OR dichogamy "
            "OR cleistogamy)"
        ),
        "api_terms": (
            "self-incompatibility self-compatible autonomous selfing bagging hand-pollination "
            "mating system herkogamy dichogamy cleistogamy"
        ),
        "likely_evidence_type": "reproductive_assurance",
        "keywords": (
            "self-incompatib", "self incompatib", "self-compatib", "self compatib",
            "autonomous selfing", "autogam", "bagging experiment", "bagging", "hand-pollinat",
            "hand pollinat", "mating system", "herkogam", "dichogam", "cleistogam",
            "geitonogam", "xenogam", "apomix",
        ),
    },
}


def _match_keywords(text: str, keywords: tuple[str, ...]) -> list[str]:
    lowered = (text or "").lower()
    return [keyword for keyword in keywords if keyword in lowered]


def _taxon_tokens(query_taxon: str) -> tuple[str, str]:
    """Return the lowercased (genus, species_epithet) of a binomial, else ('','')."""
    tokens = [t for t in re.split(r"\s+", str(query_taxon or "").strip()) if t.isalpha()]
    if len(tokens) < 2:
        return "", ""
    return tokens[0].lower(), tokens[1].lower()


def _word_present(word: str, lowered_text: str) -> bool:
    """Whole-word presence test (avoids 'acer' matching 'aceraceae')."""
    return bool(word) and re.search(rf"\b{re.escape(word)}\b", lowered_text) is not None


# Taxon-relevance tiers (v2.1). S = the species is named with confidence (review
# target); A = only the genus is named / flora-monograph rescue (quarantined, not the
# main packet); B = the species is not named at all (background literature, retained
# but never reviewed).
TIER_SPECIES = "S"
TIER_GENUS_FLORA = "A"
TIER_BACKGROUND = "B"


def taxon_relevance_tier(taxon_relevance_score: int) -> str:
    """Map a taxon_relevance_score to a review tier (S/A/B)."""
    score = int(taxon_relevance_score or 0)
    if score >= 3:
        return TIER_SPECIES
    if score >= 1:
        return TIER_GENUS_FLORA
    return TIER_BACKGROUND


def score_taxon_relevance(query_taxon: str, title_text: str, abstract_text: str) -> tuple[int, str]:
    """Score whether a lead actually concerns the TARGET species (not the trait).

    Trait-keyword relevance says a paper is about flowers/pollination; it does not
    say the paper is about *this* taxon. A lead whose title/abstract never names the
    genus is almost never about the target species, so it scores 0 and is gated out
    of the review packet. Returns ``(score, match_kind)``; higher is stronger:

    - 4 binomial_title      full "Genus epithet" in the title
    - 3 genus_epithet_title genus and epithet both in the title (not adjacent)
    - 3 binomial_abstract   full binomial in the abstract
    - 2 genus_epithet_abstract  genus and epithet both in the abstract
    - 1 genus_title / genus_abstract  genus present, epithet absent
    - 0 none                genus absent (epithet alone is too ambiguous to count)
    """
    genus, epithet = _taxon_tokens(query_taxon)
    if not genus:
        return 0, "none"
    title = (title_text or "").lower()
    abstract = (abstract_text or "").lower()
    binomial = f"{genus} {epithet}"

    if epithet and binomial in title:
        return 4, "binomial_title"
    if epithet and _word_present(genus, title) and _word_present(epithet, title):
        return 3, "genus_epithet_title"
    if epithet and binomial in abstract:
        return 3, "binomial_abstract"
    if epithet and _word_present(genus, abstract) and _word_present(epithet, abstract):
        return 2, "genus_epithet_abstract"
    if _word_present(genus, title):
        return 1, "genus_title"
    if _word_present(genus, abstract):
        return 1, "genus_abstract"
    return 0, "none"


def openalex_abstract(work: dict[str, Any]) -> str:
    """Reconstruct an OpenAlex abstract from its inverted index (for scoring)."""
    inverted = work.get("abstract_inverted_index")
    if not isinstance(inverted, dict) or not inverted:
        return ""
    positions: dict[int, str] = {}
    for word, indices in inverted.items():
        for index in indices or []:
            positions[int(index)] = word
    return " ".join(positions[i] for i in sorted(positions))


def _strip_markup(text: str) -> str:
    return re.sub(r"<[^>]+>", " ", str(text or ""))


# Provisional source-tier HINT, using the study's own reliability vocabulary
# (see prompts/trait_evidence_extraction_v2.md and config/attrition_audit.yml).
# It grades the likely reliability of the SOURCE, never a trait value, and a
# human confirms the real grade on review.
GRADE_A = "A_primary_or_monograph"
GRADE_B = "B_curated_database_or_institution"
GRADE_C = "C_curated_specialist_web"
GRADE_D = "D_unvetted_web"
GRADE_NONE = "none"
_GRADE_PRIORITY = {GRADE_A: 1, GRADE_B: 2, GRADE_C: 3, GRADE_D: 4, GRADE_NONE: 5}
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
    "query_trait_group",
    "query_template",
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
    "taxon_relevance_score",
    "taxon_match_kind",
    "taxon_relevance_tier",
    "title_relevance_score",
    "abstract_relevance_score",
    "trait_keywords_matched",
    "likely_evidence_type",
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
    """Discover free, trait-group-targeted literature leads (no trait values)."""


@dataclass
class LiteratureSeed:
    """One unreviewed, trait-group-targeted bibliographic lead for a taxon."""

    query_taxon: str
    source_api: str
    query_trait_group: str = ""
    query_template: str = ""
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
    taxon_relevance_score: int = 0
    taxon_match_kind: str = "none"
    taxon_relevance_tier: str = TIER_BACKGROUND
    title_relevance_score: int = 0
    abstract_relevance_score: int = 0
    trait_keywords_matched: str = ""
    likely_evidence_type: str = ""
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

    def score_trait_relevance(self, title_text: str, abstract_text: str, group_cfg: dict[str, Any]) -> None:
        """Score how relevant this lead is to its trait group (never a trait value)."""
        keywords = tuple(group_cfg.get("keywords", ()))
        title_hits = _match_keywords(title_text, keywords)
        abstract_hits = _match_keywords(abstract_text, keywords)
        self.title_relevance_score = len(title_hits)
        self.abstract_relevance_score = len(abstract_hits)
        self.trait_keywords_matched = "|".join(sorted(set(title_hits) | set(abstract_hits)))
        self.likely_evidence_type = str(group_cfg.get("likely_evidence_type", ""))

    def score_taxon_relevance(self, title_text: str, abstract_text: str) -> None:
        """Score whether this lead concerns the target species (never a trait value)."""
        self.taxon_relevance_score, self.taxon_match_kind = score_taxon_relevance(
            self.query_taxon, title_text, abstract_text
        )
        self.taxon_relevance_tier = taxon_relevance_tier(self.taxon_relevance_score)

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


def parse_openalex_works(
    payload: dict[str, Any],
    query_taxon: str,
    limit: int,
    trait_group: str = "",
    group_cfg: dict[str, Any] | None = None,
) -> list[LiteratureSeed]:
    """Parse an OpenAlex /works response into trait-scored literature seeds."""
    group_cfg = group_cfg or {}
    template = str(group_cfg.get("query_template", "")).replace("{species}", query_taxon)
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
        title = str(work.get("display_name") or work.get("title") or "").strip()
        seed = LiteratureSeed(
            query_taxon=query_taxon,
            source_api="openalex",
            query_trait_group=trait_group,
            query_template=template,
            title=title,
            doi=_clean_doi(work.get("doi")),
            publication_year=str(work.get("publication_year") or "").strip(),
            venue=str(source.get("display_name") or "").strip(),
            authors=authors,
            is_open_access="" if is_oa is None else str(bool(is_oa)).lower(),
            oa_status=str(oa.get("oa_status") or "").strip(),
            open_access_url=str(oa.get("oa_url") or "").strip(),
            candidate_pdf_url=str(location.get("pdf_url") or "").strip(),
        )
        abstract_text = openalex_abstract(work)
        seed.score_trait_relevance(title, abstract_text, group_cfg)
        seed.score_taxon_relevance(title, abstract_text)
        seeds.append(seed)
    return seeds


def parse_crossref_works(
    payload: dict[str, Any],
    query_taxon: str,
    limit: int,
    trait_group: str = "",
    group_cfg: dict[str, Any] | None = None,
) -> list[LiteratureSeed]:
    """Parse a Crossref /works response into trait-scored literature seeds."""
    group_cfg = group_cfg or {}
    template = str(group_cfg.get("query_template", "")).replace("{species}", query_taxon)
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
        title = str(title_list[0]).strip() if title_list else ""
        seed = LiteratureSeed(
            query_taxon=query_taxon,
            source_api="crossref",
            query_trait_group=trait_group,
            query_template=template,
            title=title,
            doi=_clean_doi(item.get("DOI")),
            publication_year=year,
            venue=str(container[0]).strip() if container else "",
            authors=authors,
        )
        abstract_text = _strip_markup(item.get("abstract", ""))
        seed.score_trait_relevance(title, abstract_text, group_cfg)
        seed.score_taxon_relevance(title, abstract_text)
        seeds.append(seed)
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
    max_seeds_per_group: int,
    passthrough: dict[str, str] | None = None,
    trait_groups: dict[str, dict[str, Any]] | None = None,
) -> tuple[list[LiteratureSeed], list[str]]:
    """Collect trait-group-targeted leads for one taxon. Returns (seeds, errors)."""
    passthrough = passthrough or {}
    groups = trait_groups or TRAIT_GROUP_QUERIES
    seeds: list[LiteratureSeed] = []
    errors: list[str] = []

    for group_key, cfg in groups.items():
        # v2.1: quote the binomial so the API treats the scientific name as an exact
        # phrase, biasing retrieval toward on-species literature rather than any work
        # that merely shares the genus or the trait keywords.
        query = f'"{taxon}" {cfg.get("api_terms", "")}'.strip()
        try:
            oa_payload = getter(
                OPENALEX_WORKS_URL,
                {"search": query, "per_page": max_seeds_per_group, "mailto": contact_email},
            )
            seeds.extend(parse_openalex_works(oa_payload, taxon, max_seeds_per_group, group_key, cfg))
        except Exception as exc:  # noqa: BLE001 - one source failing must not abort the taxon
            errors.append(f"openalex:{group_key}:{taxon}:{exc}")
        try:
            cr_payload = getter(
                CROSSREF_WORKS_URL,
                {"query.bibliographic": query, "rows": max_seeds_per_group, "mailto": contact_email},
            )
            seeds.extend(parse_crossref_works(cr_payload, taxon, max_seeds_per_group, group_key, cfg))
        except Exception as exc:  # noqa: BLE001
            errors.append(f"crossref:{group_key}:{taxon}:{exc}")

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
    """Build the trait-group-targeted lead table and a run report for staged taxa.

    ``max_seeds_per_taxon`` bounds the seeds requested per trait group per API.
    """
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
        for column in (
            "taxon_relevance_score", "title_relevance_score", "abstract_relevance_score", "priority_rank"
        ):
            frame[column] = pd.to_numeric(frame[column], errors="coerce").fillna(0).astype(int)
        # Within a taxon and trait group, rank by TAXON relevance first (is the lead
        # actually about this species?), then trait relevance. A taxon_relevance_score
        # of 0 (the genus never appears) can never head the review, and a higher source
        # grade only breaks ties among otherwise-equal leads.
        frame = frame.sort_values(
            ["query_taxon", "query_trait_group", "taxon_relevance_score",
             "title_relevance_score", "abstract_relevance_score", "priority_rank"],
            ascending=[True, True, False, False, False, True],
        ).reset_index(drop=True)

    def _n_hint(grade: str) -> int:
        return int((frame["provisional_source_reliability_hint"] == grade).sum()) if len(frame) else 0

    def _n_group(group: str) -> int:
        return int((frame["query_trait_group"] == group).sum()) if len(frame) else 0

    n_title_relevant = int((frame["title_relevance_score"] > 0).sum()) if len(frame) else 0
    n_taxon_relevant = int((frame["taxon_relevance_score"] > 0).sum()) if len(frame) else 0
    if len(frame):
        matched_species = set(
            frame.loc[frame["taxon_relevance_score"] > 0, "query_taxon"].astype(str)
        )
    else:
        matched_species = set()
    n_taxon_matched_species = len(matched_species)
    n_zero_taxon_species = int(len(taxa)) - n_taxon_matched_species
    report = {
        "note": (
            "Unreviewed, trait-group-targeted literature leads only. A source lead is "
            "NOT trait evidence; no trait value, native/establishment status, Bombus "
            "applicability, or analysis inclusion is decided here. Retrieval is "
            "species x trait-group. Leads are ranked by taxon relevance first: a "
            "taxon_relevance_score==0 lead (the genus never appears in title/abstract) "
            "does not credibly concern the target species, is ranked last, and is gated "
            "out of the review packet. title_relevance_score==0 leads also never head review."
        ),
        "n_taxa_queried": int(len(taxa)),
        "n_literature_seeds": int(len(frame)),
        "n_taxon_relevant_seeds": n_taxon_relevant,
        "n_zero_taxon_relevance_seeds": int(len(frame)) - n_taxon_relevant,
        "n_taxon_matched_species": n_taxon_matched_species,
        "n_zero_taxon_species": n_zero_taxon_species,
        "n_seeds_tier_S_species_confident": int((frame["taxon_relevance_tier"] == TIER_SPECIES).sum()) if len(frame) else 0,
        "n_seeds_tier_A_genus_flora": int((frame["taxon_relevance_tier"] == TIER_GENUS_FLORA).sum()) if len(frame) else 0,
        "n_seeds_tier_B_background": int((frame["taxon_relevance_tier"] == TIER_BACKGROUND).sum()) if len(frame) else 0,
        "n_seeds_group_A_floral_morphology_colour": _n_group("A_floral_morphology_colour"),
        "n_seeds_group_B_pollination_pollen_vector": _n_group("B_pollination_pollen_vector"),
        "n_seeds_group_C_reproductive_assurance": _n_group("C_reproductive_assurance"),
        "n_title_relevant_seeds": n_title_relevant,
        "n_zero_title_relevance_seeds": int(len(frame)) - n_title_relevant,
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
    max_seeds_per_taxon: int = typer.Option(5, min=1, help="Max leads per trait group per API."),
) -> None:
    """Write unreviewed, trait-group-targeted literature leads for a nominated island."""
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
        f"{report['n_literature_seeds']} unreviewed lead(s) ({report['n_title_relevant_seeds']} "
        f"title-relevant) across 3 trait groups for {report['n_taxa_queried']} taxon(s) on "
        f"{island_id}. No trait/native/applicability value was decided."
    )


if __name__ == "__main__":
    app()
