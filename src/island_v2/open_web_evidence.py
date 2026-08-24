"""Fetch source pages and extract fail-closed species-trait evidence.

The public interface accepts discovery URLs plus a master taxonomy and returns
auditable candidates.  Search snippets are intentionally absent from that
interface.  Unknown domains, fuzzy names, family conflicts, cultivar-only
pages, unsupported ontology values, and incomplete provenance remain audit
rows rather than trait evidence.
"""

# Typer's documented command style evaluates Option metadata in defaults.
# ruff: noqa: B008

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import time
import urllib.parse
import urllib.robotparser
from collections import Counter, defaultdict, deque
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import httpx
import typer
import yaml
from bs4 import BeautifulSoup
from pypdf import PdfReader

CONTRACT = "open_web_species_trait_evidence_v1"
USER_AGENT = "island-floral-traits/0.1 (research crawler; contact via GitHub repository)"
MAX_PAGE_BYTES = 8 * 1024 * 1024
SCIENTIFIC_NAME = re.compile(
    r"(?<![A-Za-z])([A-Z][a-z]{2,}(?:\s+×)?\s+[a-z][a-z-]{2,})(?![A-Za-z-])"
)
CULTIVAR = re.compile(
    r"\b(cultivar|cultivars|variety|horticultural hybrid|hybrid cultivar|cv\.)\b|"
    r"園芸品種|栽培品種|交配種|品种|栽培品种|cultivar|hybride horticole",
    re.IGNORECASE,
)
GENERIC_HYBRID = re.compile(r"(?:\s|^)[×x]\s*[a-z]|hybrid species", re.IGNORECASE)
INFRASPECIFIC_RANK = re.compile(
    r"\b(?:subsp|ssp|var|f|form|forma|subvar)\.?\s+[a-z]", re.IGNORECASE
)
SCIENTIFIC_NAME_LABEL = re.compile(
    r"学\s*名|scientific\s+name|nombre\s+cient[ií]fico|nome\s+cient[ií]fico|"
    r"nom\s+scientifique|wissenschaftlicher\s+name|学名",
    re.IGNORECASE,
)
FAMILY_NAME = re.compile(r"\b([A-Z][a-z]+aceae)\b")
NEGATION = re.compile(
    r"\b(no|not|without|lacks?|never|non[- ]?)\b|"
    r"(?:ではない|を欠く|がない|しない|无|不|sin|sans|kein(?:e|en|er)?|nicht)",
    re.IGNORECASE,
)
NON_FLORAL_ORGANS = re.compile(
    r"\b(leaf|leaves|fruit|fruits|seed|seeds|stem|bark)\b|"
    r"(葉|果実|種子|茎|樹皮|叶|果|种子|hoja|fruto|feuille|fruit|blatt|frucht)",
    re.IGNORECASE,
)
FLORAL_CONTEXT = re.compile(
    r"\b(flower|flowers|floral|corolla|corollas|petal|petals|inflorescence|"
    r"pollinat|self[- ]|mating|breeding|cleistogam|nectar)\b|"
    r"(花|花冠|花弁|花序|受粉|自殖|交配|閉鎖花|花蜜|传粉|自交|花瓣|花序|"
    r"flor|corola|infloresc|pollinis|blüte|krone)",
    re.IGNORECASE,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _text(value: object) -> str:
    return " ".join(str(value or "").replace("\xa0", " ").split())


def _now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _sha_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _name_key(value: object) -> str:
    cleaned = _text(value).replace("× ", "×")
    tokens = cleaned.split()
    if len(tokens) < 2:
        return ""
    if tokens[1] == "×" and len(tokens) >= 3:
        return f"{tokens[0].casefold()} ×{tokens[2].casefold()}"
    return f"{tokens[0].casefold()} {tokens[1].casefold()}"


def canonical_url(value: object) -> str:
    parsed = urllib.parse.urlsplit(_text(value))
    if parsed.scheme.casefold() not in {"http", "https"} or not parsed.netloc:
        return ""
    kept = [
        (key, val)
        for key, val in urllib.parse.parse_qsl(parsed.query, keep_blank_values=True)
        if not key.casefold().startswith(("utm_", "ref", "source", "fbclid", "gclid"))
    ]
    path = re.sub(r"/+", "/", parsed.path).rstrip("/") or "/"
    return urllib.parse.urlunsplit(
        (
            parsed.scheme.casefold(),
            parsed.netloc.casefold().removeprefix("www."),
            path,
            urllib.parse.urlencode(kept),
            "",
        )
    )


def _normalised_lineage_text(value: object) -> str:
    return re.sub(
        r"\s+",
        " ",
        re.sub(r"[^\w\s]", " ", _text(value).casefold()),
    ).strip()


def source_lineage(
    url: str,
    excerpt: str,
    *,
    citation: str = "",
    content_fingerprint: str = "",
) -> tuple[str, str]:
    """Identify original-source lineage beyond URL-only equality."""

    doi = re.search(
        r"\b10\.\d{4,9}/[-._;()/:A-Za-z0-9]+\b",
        f"{url} {citation} {excerpt}",
    )
    if doi:
        return f"doi:{doi.group(0).rstrip('.,;').casefold()}", "doi"
    normalised_citation = _normalised_lineage_text(citation)
    normalised_excerpt = _normalised_lineage_text(excerpt)
    if normalised_citation and normalised_excerpt:
        fingerprint = _sha_text(f"{normalised_citation}|{normalised_excerpt}")
        return f"citation-excerpt:{fingerprint}", "citation_normalized_excerpt"
    if _text(content_fingerprint):
        return (
            f"content:{_text(content_fingerprint).casefold()}",
            "content_fingerprint",
        )
    canonical = canonical_url(url)
    if canonical:
        return f"url:{canonical}", "canonical_url"
    return f"excerpt:{_sha_text(_text(excerpt).casefold())}", "excerpt_sha256"


@dataclass(frozen=True)
class DomainRecord:
    domain: str
    source_name: str
    source_type: str
    region: str
    language: str
    organization_type: str
    robots_access_notes: str
    page_pattern: str
    treatment_end_pattern: str
    inventory_urls: tuple[str, ...]
    inventory_depth: int
    scientific_name_displayed: bool
    cultivar_descriptions_present: bool
    provenance_quality: str
    recommended_evidence_tier: str
    review_status: str
    allowed_traits: tuple[str, ...] | str

    def allows_trait(self, trait_name: str) -> bool:
        return self.allowed_traits == "all" or trait_name in self.allowed_traits


class DomainRegistry:
    """Registry seam: URL classification is centralized and fail-closed."""

    def __init__(self, records: Sequence[DomainRecord]) -> None:
        self._records = {record.domain.casefold(): record for record in records}

    @classmethod
    def load(cls, path: Path) -> DomainRegistry:
        payload = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        rows = payload.get("domains", []) if isinstance(payload, dict) else []
        records: list[DomainRecord] = []
        for row in rows:
            allowed = row.get("allowed_traits", [])
            records.append(
                DomainRecord(
                    domain=_text(row.get("domain")).casefold(),
                    source_name=_text(row.get("source_name")),
                    source_type=_text(row.get("source_type")),
                    region=_text(row.get("region")),
                    language=_text(row.get("language")),
                    organization_type=_text(row.get("organization_type")),
                    robots_access_notes=_text(row.get("robots_access_notes")),
                    page_pattern=_text(row.get("page_pattern")),
                    treatment_end_pattern=_text(row.get("treatment_end_pattern")),
                    inventory_urls=tuple(row.get("inventory_urls") or ()),
                    inventory_depth=int(row.get("inventory_depth") or 0),
                    scientific_name_displayed=bool(row.get("scientific_name_displayed")),
                    cultivar_descriptions_present=bool(row.get("cultivar_descriptions_present")),
                    provenance_quality=_text(row.get("provenance_quality")),
                    recommended_evidence_tier=_text(row.get("recommended_evidence_tier")),
                    review_status=_text(row.get("review_status")),
                    allowed_traits=(
                        "all"
                        if allowed == "all"
                        else tuple(_text(item) for item in (allowed or []))
                    ),
                )
            )
        return cls(records)

    def lookup(self, url: object) -> DomainRecord | None:
        host = urllib.parse.urlsplit(_text(url)).netloc.casefold().split(":", 1)[0]
        host = host.removeprefix("www.")
        record = self._records.get(host)
        if (
            record
            and record.review_status == "approved"
            and re.search(record.page_pattern, canonical_url(url))
        ):
            return record
        return None

    def all(self) -> tuple[DomainRecord, ...]:
        return tuple(self._records.values())


@dataclass(frozen=True)
class NameResolution:
    submitted_name: str
    accepted_species: str
    family: str
    method: str
    status: str
    backbone_count: int
    backbone_names: str


class StrictNameResolver:
    """Exact accepted names plus cross-backbone-agreed exact synonyms only."""

    def __init__(
        self,
        master_rows: Sequence[Mapping[str, object]],
        synonym_rows: Sequence[Mapping[str, object]] = (),
    ) -> None:
        self._accepted: dict[str, tuple[str, str]] = {}
        for row in master_rows:
            species = _text(row.get("accepted_species"))
            key = _name_key(species)
            if key:
                self._accepted[key] = (species, _text(row.get("family")))
        grouped: dict[str, list[Mapping[str, object]]] = defaultdict(list)
        for row in synonym_rows:
            grouped[
                _name_key(row.get("synonym") or row.get("basionym") or row.get("submitted_name"))
            ].append(row)
        self._synonyms: dict[str, NameResolution] = {}
        for key, rows in grouped.items():
            if not key:
                continue
            accepted = {
                _text(row.get("accepted_species") or row.get("accepted_name")) for row in rows
            }
            accepted.discard("")
            backbones = {_text(row.get("backbone") or row.get("source")).casefold() for row in rows}
            backbones.discard("")
            families = {_text(row.get("family")) for row in rows if _text(row.get("family"))}
            if len(accepted) != 1 or len(backbones) < 2 or len(families) > 1:
                continue
            accepted_name = next(iter(accepted))
            master = self._accepted.get(_name_key(accepted_name))
            if master is None or (families and master[1] and master[1] not in families):
                continue
            self._synonyms[key] = NameResolution(
                submitted_name=_text(
                    rows[0].get("synonym")
                    or rows[0].get("basionym")
                    or rows[0].get("submitted_name")
                ),
                accepted_species=master[0],
                family=master[1],
                method="strict_cross_backbone_synonym",
                status="accepted",
                backbone_count=len(backbones),
                backbone_names=";".join(sorted(backbones)),
            )

    def resolve(self, name: object) -> NameResolution:
        submitted = _text(name)
        key = _name_key(submitted)
        if key in self._accepted:
            accepted, family = self._accepted[key]
            return NameResolution(
                submitted_name=submitted,
                accepted_species=accepted,
                family=family,
                method="exact_accepted_name",
                status="accepted",
                backbone_count=1,
                backbone_names="master_denominator",
            )
        if key in self._synonyms:
            return self._synonyms[key]
        return NameResolution(
            submitted_name=submitted,
            accepted_species="",
            family="",
            method="none",
            status="rejected_no_exact_cross_backbone_match",
            backbone_count=0,
            backbone_names="",
        )


@dataclass(frozen=True)
class Page:
    requested_url: str
    final_url: str
    status_code: int
    content_type: str
    title: str
    text: str
    language: str
    retrieved_at_utc: str
    content_sha256: str
    error: str = ""
    links: tuple[str, ...] = ()


class PageFetcher:
    def __init__(
        self,
        *,
        client: httpx.Client | None = None,
        pause_seconds: float = 0.25,
        max_bytes: int = MAX_PAGE_BYTES,
    ) -> None:
        self._client = client or httpx.Client(
            timeout=40,
            follow_redirects=True,
            headers={"User-Agent": USER_AGENT, "Accept": "text/html,application/pdf;q=0.9"},
        )
        self._pause = max(0.0, pause_seconds)
        self._max_bytes = max_bytes
        self._last_request = 0.0
        self._robots: dict[str, urllib.robotparser.RobotFileParser | None] = {}
        self._cache: dict[str, Page] = {}

    def _wait(self) -> None:
        remaining = self._pause - (time.monotonic() - self._last_request)
        if remaining > 0:
            time.sleep(remaining)

    def _request(self, url: str) -> httpx.Response:
        self._wait()
        response = self._client.get(url)
        self._last_request = time.monotonic()
        return response

    def _robots_allowed(self, url: str) -> bool:
        parsed = urllib.parse.urlsplit(url)
        origin = f"{parsed.scheme}://{parsed.netloc}"
        if origin not in self._robots:
            robots_url = f"{origin}/robots.txt"
            try:
                response = self._request(robots_url)
                if response.status_code == 404:
                    self._robots[origin] = None
                elif response.status_code == 200:
                    parser = urllib.robotparser.RobotFileParser()
                    parser.set_url(robots_url)
                    parser.parse(response.text.splitlines())
                    self._robots[origin] = parser
                else:
                    self._robots[origin] = False  # type: ignore[assignment]
            except httpx.HTTPError:
                self._robots[origin] = False  # type: ignore[assignment]
        parser = self._robots[origin]
        return parser is None or bool(parser and parser.can_fetch(USER_AGENT, url))

    def fetch(self, url: str) -> Page:
        requested = canonical_url(url)
        if requested in self._cache:
            return self._cache[requested]
        retrieved = _now()
        if not requested:
            return Page(url, "", 0, "", "", "", "", retrieved, "", "invalid_url")
        if not self._robots_allowed(requested):
            return Page(requested, requested, 0, "", "", "", "", retrieved, "", "robots_disallowed")
        try:
            response = self._request(requested)
            response.raise_for_status()
        except httpx.HTTPError as exc:
            return Page(
                requested,
                requested,
                getattr(getattr(exc, "response", None), "status_code", 0) or 0,
                "",
                "",
                "",
                "",
                retrieved,
                "",
                f"{type(exc).__name__}:{_text(exc)}",
            )
        raw = response.content
        if len(raw) > self._max_bytes:
            return Page(
                requested,
                canonical_url(response.url),
                response.status_code,
                response.headers.get("content-type", ""),
                "",
                "",
                "",
                retrieved,
                _sha_text(raw[: self._max_bytes].hex()),
                "page_too_large",
            )
        content_type = response.headers.get("content-type", "").casefold()
        try:
            if "pdf" in content_type or requested.casefold().endswith(".pdf"):
                reader = PdfReader(io.BytesIO(raw))
                text = "\n".join(page.extract_text() or "" for page in reader.pages)
                title = _text(reader.metadata.title if reader.metadata else "")
                language = ""
            else:
                soup = BeautifulSoup(raw, "html.parser")
                links = tuple(
                    canonical_url(urllib.parse.urljoin(str(response.url), link.get("href", "")))
                    for link in soup.select("a[href]")
                )
                for element in soup(["script", "style", "nav", "noscript"]):
                    element.decompose()
                title = _text(soup.title.get_text(" ", strip=True) if soup.title else "")
                language = _text((soup.html or {}).get("lang") if soup.html else "")
                text = "\n".join(
                    line
                    for line in (_text(item) for item in soup.get_text("\n").splitlines())
                    if line
                )
                links = tuple(link for link in links if link)
        except (ValueError, UnicodeError, OSError) as exc:
            return Page(
                requested,
                canonical_url(response.url),
                response.status_code,
                content_type,
                "",
                "",
                "",
                retrieved,
                hashlib.sha256(raw).hexdigest(),
                f"parse_error:{type(exc).__name__}",
            )
        page = Page(
            requested_url=requested,
            final_url=canonical_url(response.url),
            status_code=response.status_code,
            content_type=content_type,
            title=title,
            text=text,
            language=language,
            retrieved_at_utc=retrieved,
            content_sha256=hashlib.sha256(raw).hexdigest(),
            links=links if "links" in locals() else (),
        )
        self._cache[requested] = page
        return page


def page_scientific_names(page: Page, *, title_only: bool = False) -> list[str]:
    haystack = page.title if title_only else f"{page.title}\n{page.text[:12000]}"
    candidates = SCIENTIFIC_NAME.findall(haystack)
    excluded = {"Flora of", "Plants of", "Journal of", "University of"}
    return list(
        dict.fromkeys(
            name for name in (_text(candidate) for candidate in candidates) if name not in excluded
        )
    )


def identity_gate(
    page: Page,
    *,
    searched_name: str,
    resolver: StrictNameResolver,
    expected_species: str = "",
) -> tuple[NameResolution, str, str]:
    """Return resolution, matched page name, and rejection reason."""

    if page.error or not page.text:
        return resolver.resolve(""), "", "page_fetch_or_parse_failure"
    if CULTIVAR.search(page.title) or (
        CULTIVAR.search(page.text[:1500]) and ("'" in page.title or "cv." in page.title.casefold())
    ):
        return resolver.resolve(""), "", "cultivar_or_horticultural_hybrid_page"
    if GENERIC_HYBRID.search(page.title):
        return resolver.resolve(""), "", "hybrid_page"
    names = page_scientific_names(page)
    title_names = page_scientific_names(page, title_only=True)
    searched_key = _name_key(searched_name)
    if not expected_species:
        # Inventory crawling has no species attached to the lead.  Treat the
        # page title as the subject boundary; names buried in references or
        # comparison prose must never become the page's focal species.
        if not title_names:
            return resolver.resolve(""), "", "subject_scientific_name_absent_from_title"
        names = title_names
    else:
        names = sorted(
            names,
            key=lambda name: (
                0 if name in title_names else 1,
                0 if _name_key(name) == searched_key else 1,
            ),
        )
    expected_key = _name_key(expected_species)
    for page_name in names:
        resolution = resolver.resolve(page_name)
        if resolution.status != "accepted":
            continue
        if expected_key and _name_key(resolution.accepted_species) != expected_key:
            continue
        target_key = expected_key or _name_key(resolution.accepted_species)
        for label in SCIENTIFIC_NAME_LABEL.finditer(page.text[:2500]):
            labeled_window = page.text[label.end() : label.end() + 180]
            labeled_match = re.search(
                r"([A-Z][a-z]{2,}\s+(?:[a-z][a-z-]{2,}|sp\.))([^\n]{0,100})",
                labeled_window,
            )
            if not labeled_match:
                continue
            labeled_name = _text(labeled_match.group(1))
            labeled_tail = labeled_match.group(2)
            if re.search(r"(?:\bcv\.?\s*)?['‘’\"][^'‘’\"]+['‘’\"]", labeled_tail):
                return resolution, page_name, "labeled_scientific_name_is_cultivar"
            labeled_resolution = resolver.resolve(labeled_name)
            if (
                labeled_resolution.status != "accepted"
                or _name_key(labeled_resolution.accepted_species) != target_key
            ):
                return resolution, page_name, "labeled_scientific_name_mismatch"
            break
        if GENERIC_HYBRID.search(resolution.accepted_species):
            return resolution, page_name, "accepted_species_is_hybrid"
        title_index = page.title.casefold().find(page_name.casefold())
        if title_index >= 0:
            title_tail = page.title[title_index + len(page_name) :][:100]
            if INFRASPECIFIC_RANK.search(title_tail):
                return resolution, page_name, "page_subject_rank_is_infraspecific"
        body_index = page.text.casefold().find(page_name.casefold())
        if body_index >= 0:
            subject_line_tail = page.text[body_index + len(page_name) :].split("\n", 1)[0]
            if re.search(r"(?:\bcv\.?\s*)?['‘’\"][^'‘’\"]+['‘’\"]", subject_line_tail[:100]):
                return resolution, page_name, "page_subject_is_cultivar"
        explicit_families = set(FAMILY_NAME.findall(f"{page.title} {page.text[:4000]}"))
        if explicit_families and resolution.family and resolution.family not in explicit_families:
            return resolution, page_name, "family_conflict"
        return resolution, page_name, ""
    if searched_key and searched_key not in {_name_key(name) for name in names}:
        return resolver.resolve(searched_name), "", "searched_name_absent_from_page"
    return resolver.resolve(searched_name), "", "name_match_not_exact_or_cross_backbone_synonym"


def _sentences(text: str) -> list[str]:
    return [
        _text(sentence)
        for sentence in re.split(
            r"(?<=[。！？；])\s*|(?<=[.!?;])\s+|\n+",
            text,
        )
        if _text(sentence)
    ]


def _paragraphs(text: str) -> list[str]:
    return [_text(paragraph) for paragraph in re.split(r"\n+", text) if _text(paragraph)]


def _term_match(sentence: str, term: str) -> bool:
    return _text(term).casefold() in sentence.casefold()


def _negated(sentence: str, term: str) -> bool:
    index = sentence.casefold().find(_text(term).casefold())
    if index < 0:
        return False
    before = sentence[max(0, index - 45) : index]
    after = sentence[index + len(_text(term)) : index + len(_text(term)) + 45]
    return bool(
        NEGATION.search(before)
        or re.search(
            r"(?:が|は|も)?(?:ない|無い)|を欠く|なし|"
            r"\b(?:absent|missing|lacking|not present)\b",
            after,
            re.IGNORECASE,
        )
    )


def _collapse_state_matches(
    matches: Sequence[tuple[str, str, str]],
    *,
    trait_name: str,
    source_text: str,
) -> list[tuple[str, str, str]]:
    deduped = list(dict.fromkeys(matches))
    if trait_name not in {"flower_primary_color", "floral_form", "floral_symmetry"}:
        return deduped
    states = sorted({normalized for _, normalized, _ in deduped})
    if len(states) <= 1:
        return deduped
    excerpts = [excerpt for _, _, excerpt in deduped]
    starts = [source_text.find(excerpt) for excerpt in excerpts]
    starts = [index for index in starts if index >= 0]
    if not starts:
        return deduped
    start = min(starts)
    end = max(
        source_text.find(excerpt) + len(excerpt)
        for excerpt in excerpts
        if source_text.find(excerpt) >= 0
    )
    variable = (
        "multicolored_variable" if trait_name == "flower_primary_color" else "multistate_variable"
    )
    return [(";".join(states), variable, source_text[start:end])]


def primary_treatment_text(page: Page, *, treatment_end_pattern: str = "") -> str:
    """Exclude configured appended catalogues from a focal species treatment."""

    pattern = treatment_end_pattern or (
        r"\bfamily\s*:?\s*[A-Z][a-z]+aceae\s*[－–—-]\s*genus\s+[A-Z][a-z]+"
    )
    marker = re.search(pattern, page.text, re.IGNORECASE)
    if marker and marker.start() > 0:
        return page.text[: marker.start()]
    return page.text


NON_FLOWER_SIZE_BRIDGE = re.compile(
    r"\b(?:flower|floral)\s+(?:stalk|pedicel|peduncle)|"
    r"\b(?:pedicel|peduncle|inflorescence|bract|lemma|glume|awn|ovary|anther|"
    r"stamen|filament|fruit|seed|shoot)\b|"
    r"小?花柄|花茎|花序|花穂|花托|シュート|苞|護[穎頴]|[穎頴]|芒|"
    r"葯|雄しべ|花糸|葉|茎|果|痩果|子房|種子",
    re.IGNORECASE,
)
TUBE_ONLY_BRIDGE = re.compile(
    r"\b(?:corolla|floral)\s+tube\b|筒部|花冠(?:の)?筒|花筒|萼筒|"
    r"tubo de la corola|tube de la corolle|Kronröhre|tubus corollae",
    re.IGNORECASE,
)


def _direct_flower_measurement(
    sentence: str,
    *,
    context: re.Match[str],
    measure: re.Match[str],
    trait_name: str,
) -> bool:
    """Reject dimensions of nearby non-floral organs and nested structures."""

    bridge = sentence[context.start() : context.start() + measure.start()]
    if re.search(r"似て|resembl|花に比べ", sentence, re.IGNORECASE):
        return False
    if trait_name == "tube_depth_class":
        if re.search(r"\bspathe\b|仏炎苞", sentence, re.IGNORECASE):
            return False
        tube_mentions = list(TUBE_ONLY_BRIDGE.finditer(bridge))
        if not tube_mentions:
            return False
        tube_tail = bridge[tube_mentions[-1].end() :]
        if len(tube_tail) > 40:
            return False
        if re.search(
            r"\b(?:petal|lobe|scale|limb|style|stigma|stamen|filament|ovary|anther)\b|"
            r"花被片|花弁|裂片|鱗片|拡大部|萼片|花柱|柱頭|雄しべ|花糸|子房|葯",
            tube_tail,
            re.IGNORECASE,
        ):
            return False
        after_measure = sentence[context.start() + measure.end() :][:30]
        if re.search(
            r"^(?:\D{0,18})(?:鱗片|花被片|花弁|裂片|拡大部|萼片|"
            r"花柱|柱頭|雄しべ|花糸|子房|葯|"
            r"\b(?:scales?|petals?|lobes?|limbs?|styles?|stigmas?|stamens?|"
            r"filaments?|ovaries|anthers?)\b)",
            after_measure,
            re.IGNORECASE,
        ):
            return False
        return bool(
            re.search(
                r"長さ|\blength\b|longitud|comprimento|longueur|Länge",
                tube_tail,
                re.IGNORECASE,
            )
        )
    if NON_FLOWER_SIZE_BRIDGE.search(bridge):
        return False
    if TUBE_ONLY_BRIDGE.search(bridge):
        return False
    if re.search(
        r"(?:花|flowers?)\D{0,24}(?:\d+\s*(?:～|〜|-|to)?\s*\d*|数)\s*個",
        bridge,
    ):
        return False
    if re.search(r"花(?:が|を)\s*(?:つき|密生し|少なく|多く)\s*[、,]", bridge):
        return False
    if re.search(r"花(?:は)?\s*(?:無く|なく)|\bflowers?\s+(?:absent|lacking)\b", bridge):
        return False
    if re.search(r"花被(?:に|で)\s*(?:包|覆)", bridge):
        return False
    if re.search(r"花の下部|基部", bridge):
        return False
    return bool(
        re.search(
            r"直径|長さ|幅|径|\b(?:diameter|length|wide|across)\b|"
            r"longitud|diámetro|comprimento|diamètre|longueur|Länge|Durchmesser",
            bridge,
            re.IGNORECASE,
        )
    )


def extract_rules(
    page: Page,
    *,
    trait_name: str,
    config: Mapping[str, Any],
    treatment_end_pattern: str = "",
) -> list[tuple[str, str, str]]:
    """Return raw value, normalized value, exact source sentence."""

    treatment_text = primary_treatment_text(
        page,
        treatment_end_pattern=treatment_end_pattern,
    )
    subject_keys = {_name_key(name) for name in page_scientific_names(page, title_only=True)}
    if trait_name in {"flower_size_class", "tube_depth_class"}:
        measurements: list[tuple[str, str, str]] = []
        seen_measurement = False
        stop_after_comparison = False
        for paragraph in _paragraphs(treatment_text):
            if trait_name == "tube_depth_class" and re.search(
                r"\bspathe\b|仏炎苞", paragraph, re.IGNORECASE
            ):
                continue
            paragraph_references_other_species = False
            for sentence in _sentences(paragraph):
                sentence_keys = {
                    _name_key(name) for name in SCIENTIFIC_NAME.findall(sentence) if _name_key(name)
                }
                if sentence_keys and any(key not in subject_keys for key in sentence_keys):
                    paragraph_references_other_species = True
                if paragraph_references_other_species:
                    if seen_measurement:
                        stop_after_comparison = True
                        break
                    continue
                if trait_name == "flower_size_class":
                    context_pattern = re.compile(
                        r"\bflower(?!ing| stalk| pedicel| peduncle)|corolla|petal|花冠|花弁|花瓣|"
                        r"花被|(?<!頭)花(?!序|穂|期|時|托|茎|粉|柱|糸|被|柄)|"
                        r"\bflor\b|\bcorola\b|\bblüte\b|\bkrone\b",
                        re.IGNORECASE,
                    )
                else:
                    context_pattern = re.compile(
                        r"corolla tube|floral tube|花冠(?:の)?筒|筒部|花筒|"
                        r"tubo de la corola|tube de la corolle|Kronröhre|tubus corollae",
                        re.IGNORECASE,
                    )
                context = None
                excerpt = ""
                measure = None
                for candidate_context in context_pattern.finditer(sentence):
                    candidate_excerpt = sentence[
                        candidate_context.start() : candidate_context.start() + 120
                    ]
                    candidate_measure = re.search(
                        r"(\d+(?:\.\d+)?)\s*(?:[-–—~〜～]\s*(\d+(?:\.\d+)?))?\s*"
                        r"(mm|㎜|cm|㎝)",
                        candidate_excerpt,
                        re.IGNORECASE,
                    )
                    if not candidate_measure:
                        continue
                    if not _direct_flower_measurement(
                        sentence,
                        context=candidate_context,
                        measure=candidate_measure,
                        trait_name=trait_name,
                    ):
                        continue
                    context = candidate_context
                    excerpt = candidate_excerpt
                    measure = candidate_measure
                    break
                if context is None or measure is None:
                    continue
                low = float(measure.group(1))
                high = float(measure.group(2) or measure.group(1))
                midpoint = (low + high) / 2
                if measure.group(3).casefold() in {"cm", "㎝"}:
                    midpoint *= 10
                if trait_name == "flower_size_class":
                    normalized = (
                        "very_small"
                        if midpoint <= 5
                        else "small"
                        if midpoint <= 15
                        else "medium"
                        if midpoint <= 30
                        else "large"
                        if midpoint <= 60
                        else "very_large"
                    )
                else:
                    normalized = (
                        "shallow" if midpoint <= 5 else "intermediate" if midpoint <= 20 else "deep"
                    )
                measurements.append((measure.group(0), normalized, excerpt))
                seen_measurement = True
            if stop_after_comparison:
                break
        return list(dict.fromkeys(measurements))

    trait_terms = config.get("extraction_terms", {}).get(trait_name, {})
    if not isinstance(trait_terms, Mapping):
        return []
    matches: list[tuple[str, str, str]] = []
    seen_match = False
    stop_after_comparison = False
    for paragraph in _paragraphs(treatment_text):
        paragraph_references_other_species = False
        paragraph_matches: list[tuple[str, str, str]] = []
        for sentence in _sentences(paragraph):
            sentence_keys = {
                _name_key(name) for name in SCIENTIFIC_NAME.findall(sentence) if _name_key(name)
            }
            if sentence_keys and any(key not in subject_keys for key in sentence_keys):
                paragraph_references_other_species = True
            if paragraph_references_other_species:
                if seen_match:
                    stop_after_comparison = True
                    break
                continue
            if not FLORAL_CONTEXT.search(sentence):
                continue
            sentence_without_terminal = sentence.rstrip("。.!！?？；;")
            has_terminal = len(sentence_without_terminal) != len(sentence)
            if len(sentence_without_terminal) < 5 or (
                not has_terminal
                and (
                    re.match(r"^(?:は|が|を|に|の|も)", sentence_without_terminal)
                    or re.search(r"(?:は|が|を|に|の)$", sentence_without_terminal)
                )
            ):
                continue
            if (
                trait_name == "flower_primary_color"
                and NON_FLORAL_ORGANS.search(sentence)
                and not re.search(
                    r"flower|corolla|petal|花|花冠|花弁|花瓣|flor|corola|blüte",
                    sentence,
                    re.IGNORECASE,
                )
            ):
                continue
            sentence_match_count = len(paragraph_matches)
            for normalized, terms in trait_terms.items():
                for term in terms or []:
                    if not _term_match(sentence, term):
                        continue
                    if trait_name == "floral_form":
                        index = sentence.casefold().find(_text(term).casefold())
                        before = sentence[max(0, index - 60) : index]
                        if (
                            NON_FLORAL_ORGANS.search(before)
                            or re.search(r"似て|resembl", sentence, re.IGNORECASE)
                            or re.search(
                                r"多くの種|多くの種類|\b(?:many|most)\s+species\b",
                                sentence,
                                re.IGNORECASE,
                            )
                        ):
                            continue
                    if trait_name == "flower_primary_color":
                        if re.search(
                            r"(?:花に見える|花のように見える).{0,30}(?:雄しべ|葯)|"
                            r"\bappears?\s+(?:to\s+be\s+)?(?:a\s+)?flower.{0,30}"
                            r"(?:stamen|anther)",
                            sentence,
                            re.IGNORECASE,
                        ):
                            continue
                        if re.search(
                            rf"(?:乾(?:く|燥)|\b(?:dry|dried|drying)\b).{{0,20}}"
                            rf"{re.escape(_text(term))}",
                            sentence,
                            re.IGNORECASE,
                        ):
                            continue
                        index = sentence.casefold().find(_text(term).casefold())
                        local = sentence[max(0, index - 20) : index + len(_text(term)) + 35]
                        before = sentence[max(0, index - 30) : index]
                        after = sentence[index + len(_text(term)) : index + len(_text(term)) + 15]
                        explicit_flower_subject = bool(
                            re.search(
                                r"(?:花冠(?:裂片)?|花弁|花被(?:片)?|(?<!頭)花)"
                                r"\s*(?:は|が|を|に|[:：,、])",
                                before,
                            )
                        )
                        explicit_flower_object = bool(
                            re.match(
                                r"\s*(?:の)?\s*(?:花(?=[をがはにの、。\s]|$)|花冠|花弁|花被)",
                                after,
                            )
                        )
                        if not re.search(
                            r"\b(?:flower|corolla|petal)\b|花冠|花弁|花被|"
                            r"(?<!頭)花(?!序|粉|柱|糸|柄)",
                            local,
                            re.IGNORECASE,
                        ):
                            continue
                        if not explicit_flower_object and (
                            NON_FLORAL_ORGANS.search(local)
                            or re.search(
                                r"\b(?:ovary|anther|stamen|calyx|style|stigma|gland|hair|"
                                r"pubescence|lemma|glume|spikelet|inflorescence|bract|awn)\b|"
                                r"子房|葯|雄しべ|萼|花柱|柱頭|腺|毛|護[穎頴]|[穎頴]|"
                                r"小穂|花穂|花序|鱗片|芒|苞",
                                local,
                                re.IGNORECASE,
                            )
                        ):
                            continue
                        if not (explicit_flower_subject or explicit_flower_object):
                            continue
                    if trait_name == "floral_symmetry":
                        index = sentence.casefold().find(_text(term).casefold())
                        context = sentence[max(0, index - 45) : index + len(_text(term)) + 45]
                        if not re.search(
                            r"\bflower(?!ing| stalk| pedicel| peduncle)|corolla|petal|"
                            r"花冠|花弁|花瓣|(?<!頭)花(?!序|茎|粉|柱|糸|被|柄)|"
                            r"flor|corola|blüte|krone",
                            context,
                            re.IGNORECASE,
                        ):
                            continue
                    negative_value = normalized in {"absent", "rewardless"}
                    if _negated(sentence, term) and not negative_value:
                        continue
                    paragraph_matches.append((_text(term), _text(normalized), sentence))
                    break
            if len(paragraph_matches) > sentence_match_count:
                seen_match = True
        matches.extend(
            _collapse_state_matches(
                paragraph_matches,
                trait_name=trait_name,
                source_text=paragraph,
            )
        )
        if stop_after_comparison:
            break
    return list(dict.fromkeys(matches))


def _quality(tier: str, source_type: str) -> str:
    if tier == "A" and source_type in {
        "original_research",
        "thesis",
        "public_flora",
        "public_botanical_database",
        "botanic_garden_database",
    }:
        return "high"
    return "medium"


def _candidate(
    *,
    resolution: NameResolution,
    page_name: str,
    searched_name: str,
    trait_name: str,
    raw_value: str,
    normalized_value: str,
    excerpt: str,
    page: Page,
    domain: DomainRecord,
    query: str,
    search_rank: object,
    search_provider: str,
    extraction_version: str,
) -> dict[str, object]:
    citation = f"{domain.source_name}. {page.title}".strip(" .")
    normalised_excerpt = _normalised_lineage_text(excerpt)
    lineage, lineage_method = source_lineage(
        page.final_url,
        excerpt,
        citation=citation,
        content_fingerprint=page.content_sha256,
    )
    return {
        "accepted_species": resolution.accepted_species,
        "searched_name": searched_name,
        "matched_page_name": page_name,
        "trait_name": trait_name,
        "raw_value": raw_value,
        "normalized_value": normalized_value,
        "state_set_json": (
            json.dumps(raw_value.split(";"), ensure_ascii=False)
            if normalized_value in {"multicolored_variable", "multistate_variable"}
            else ""
        ),
        "source_url": page.final_url,
        "page_title": page.title,
        "domain": domain.domain,
        "source_tier": domain.recommended_evidence_tier,
        "source_type": domain.source_type,
        "organization_type": domain.organization_type,
        "language": page.language or domain.language,
        "supporting_excerpt": excerpt,
        "source_citation": citation,
        "normalized_excerpt_sha256": _sha_text(normalised_excerpt),
        "content_fingerprint": page.content_sha256,
        "retrieval_date": page.retrieved_at_utc[:10],
        "query": query,
        "search_rank": search_rank,
        "search_provider": search_provider,
        "name_match_method": resolution.method,
        "name_match_backbones": resolution.backbone_names,
        "wild_cultivated_cultivar_status": "wild_or_unspecified_not_cultivar_limited",
        "source_lineage": lineage,
        "source_lineage_method": lineage_method,
        "page_content_sha256": page.content_sha256,
        "extraction_rule_or_model_version": extraction_version,
        "confidence": _quality(domain.recommended_evidence_tier, domain.source_type),
        "review_status": "pending_manual_review",
        "candidate_class": "reported",
        "evidence_scope": "species_direct",
        "family_inference": False,
        "global_fallback": False,
    }


def process_pages(
    leads: Sequence[Mapping[str, object]],
    *,
    registry: DomainRegistry,
    resolver: StrictNameResolver,
    config: Mapping[str, Any],
    fetcher: PageFetcher,
    default_traits: Sequence[str] = (),
) -> dict[str, list[dict[str, object]]]:
    candidates: list[dict[str, object]] = []
    identity_rows: list[dict[str, object]] = []
    fetch_rows: list[dict[str, object]] = []
    domain_rows: Counter[tuple[str, str]] = Counter()
    for lead in leads:
        url = _text(lead.get("result_url") or lead.get("source_url"))
        domain = registry.lookup(url)
        if domain is None:
            identity_rows.append(
                {
                    "source_url": url,
                    "accepted_species": _text(lead.get("accepted_species")),
                    "trait_name": _text(lead.get("trait_name")),
                    "status": "rejected",
                    "reason": "domain_not_approved_or_page_pattern_mismatch",
                }
            )
            host = urllib.parse.urlsplit(url).netloc.casefold().removeprefix("www.")
            domain_rows[(host, "unapproved")] += 1
            continue
        page = fetcher.fetch(url)
        fetch_rows.append(
            {
                "requested_url": page.requested_url,
                "final_url": page.final_url,
                "status_code": page.status_code,
                "content_type": page.content_type,
                "page_title": page.title,
                "language": page.language or domain.language,
                "retrieved_at_utc": page.retrieved_at_utc,
                "content_sha256": page.content_sha256,
                "text_length": len(page.text),
                "error": page.error,
            }
        )
        searched_name = _text(
            lead.get("searched_name")
            or lead.get("accepted_species")
            or lead.get("matched_page_name")
        )
        expected = _text(lead.get("accepted_species"))
        resolution, page_name, reason = identity_gate(
            page,
            searched_name=searched_name,
            resolver=resolver,
            expected_species=expected,
        )
        if reason:
            identity_rows.append(
                {
                    "source_url": page.final_url or url,
                    "accepted_species": expected,
                    "searched_name": searched_name,
                    "matched_page_name": page_name,
                    "trait_name": _text(lead.get("trait_name")),
                    "status": "rejected",
                    "reason": reason,
                    "name_match_method": resolution.method,
                }
            )
            domain_rows[(domain.domain, "identity_rejected")] += 1
            continue
        trait_names = (
            [_text(lead.get("trait_name"))]
            if _text(lead.get("trait_name"))
            else list(default_traits)
        )
        accepted_count = 0
        for trait_name in trait_names:
            if not trait_name or not domain.allows_trait(trait_name):
                identity_rows.append(
                    {
                        "source_url": page.final_url,
                        "accepted_species": resolution.accepted_species,
                        "searched_name": searched_name,
                        "matched_page_name": page_name,
                        "trait_name": trait_name,
                        "status": "rejected",
                        "reason": "source_tier_trait_not_allowed",
                        "name_match_method": resolution.method,
                    }
                )
                continue
            extracted = extract_rules(
                page,
                trait_name=trait_name,
                config=config,
                treatment_end_pattern=domain.treatment_end_pattern,
            )
            if not extracted:
                identity_rows.append(
                    {
                        "source_url": page.final_url,
                        "accepted_species": resolution.accepted_species,
                        "searched_name": searched_name,
                        "matched_page_name": page_name,
                        "trait_name": trait_name,
                        "status": "no_evidence",
                        "reason": "page_found_but_no_supported_trait_statement",
                        "name_match_method": resolution.method,
                    }
                )
                continue
            for raw_value, normalized, excerpt in extracted:
                candidates.append(
                    _candidate(
                        resolution=resolution,
                        page_name=page_name,
                        searched_name=searched_name,
                        trait_name=trait_name,
                        raw_value=raw_value,
                        normalized_value=normalized,
                        excerpt=excerpt,
                        page=page,
                        domain=domain,
                        query=_text(lead.get("query")),
                        search_rank=lead.get("rank") or lead.get("search_rank") or "",
                        search_provider=_text(lead.get("provider") or lead.get("search_provider")),
                        extraction_version="rules_only_multilingual_v1",
                    )
                )
                accepted_count += 1
        identity_rows.append(
            {
                "source_url": page.final_url,
                "accepted_species": resolution.accepted_species,
                "searched_name": searched_name,
                "matched_page_name": page_name,
                "trait_name": _text(lead.get("trait_name")),
                "status": "candidate_extracted" if accepted_count else "no_evidence",
                "reason": "" if accepted_count else "no_allowed_supported_statement",
                "name_match_method": resolution.method,
            }
        )
        domain_rows[(domain.domain, "candidate" if accepted_count else "no_evidence")] += 1
    deduped, duplicates = dedupe_candidates(candidates)
    domain_audit = [
        {"domain": domain, "outcome": outcome, "count": count}
        for (domain, outcome), count in sorted(domain_rows.items())
    ]
    return {
        "candidates": deduped,
        "lineage_duplicates": duplicates,
        "identity_audit": identity_rows,
        "fetch_audit": fetch_rows,
        "domain_audit": domain_audit,
    }


def dedupe_candidates(
    rows: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    selected: dict[tuple[str, str, str, str], dict[str, object]] = {}
    duplicates: list[dict[str, object]] = []
    excerpt_seen: dict[tuple[str, str, str, str], str] = {}
    citation_excerpt_seen: dict[tuple[str, str, str, str], str] = {}
    content_seen: dict[tuple[str, str, str, str], str] = {}
    for raw in rows:
        row = dict(raw)
        key = (
            _text(row.get("accepted_species")),
            _text(row.get("trait_name")),
            _text(row.get("normalized_value")),
            _text(row.get("source_lineage")),
        )
        excerpt_fingerprint = _text(row.get("normalized_excerpt_sha256")) or _sha_text(
            _normalised_lineage_text(row.get("supporting_excerpt"))
        )
        excerpt_key = (
            key[0],
            key[1],
            key[2],
            excerpt_fingerprint,
        )
        citation_text = _normalised_lineage_text(row.get("source_citation"))
        citation_fingerprint = _sha_text(citation_text)
        citation_excerpt_key = (
            key[0],
            key[1],
            key[2],
            _sha_text(f"{citation_fingerprint}|{excerpt_fingerprint}"),
        )
        content_key = (
            key[0],
            key[1],
            key[2],
            _text(row.get("content_fingerprint") or row.get("page_content_sha256")).casefold(),
        )
        if citation_text:
            mirror_lineage = citation_excerpt_seen.get(citation_excerpt_key)
            if mirror_lineage and mirror_lineage != key[3]:
                row["duplicate_reason"] = "same_citation_and_normalized_excerpt"
                row["duplicate_of_lineage"] = mirror_lineage
                duplicates.append(row)
                continue
            citation_excerpt_seen[citation_excerpt_key] = key[3]
        if content_key[3]:
            mirror_lineage = content_seen.get(content_key)
            if mirror_lineage and mirror_lineage != key[3]:
                row["duplicate_reason"] = "same_content_fingerprint"
                row["duplicate_of_lineage"] = mirror_lineage
                duplicates.append(row)
                continue
            content_seen[content_key] = key[3]
        mirror_lineage = excerpt_seen.get(excerpt_key)
        if mirror_lineage and mirror_lineage != key[3]:
            row["duplicate_reason"] = "same_excerpt_republished_or_mirrored"
            row["duplicate_of_lineage"] = mirror_lineage
            duplicates.append(row)
            continue
        excerpt_seen[excerpt_key] = key[3]
        if key in selected:
            row["duplicate_reason"] = "same_source_lineage"
            row["duplicate_of_lineage"] = key[3]
            duplicates.append(row)
            continue
        selected[key] = row
    return list(selected.values()), duplicates


def crawl_registry_inventory(
    record: DomainRecord,
    *,
    fetcher: PageFetcher,
    max_urls: int,
) -> tuple[list[str], list[dict[str, object]]]:
    """Generic bounded inventory crawl; never follows outside approved patterns."""

    pending = deque((url, 0) for url in record.inventory_urls)
    seen: set[str] = set()
    accepted: list[str] = []
    audit: list[dict[str, object]] = []
    while pending and len(accepted) < max_urls:
        url, depth = pending.popleft()
        canonical = canonical_url(url)
        if not canonical or canonical in seen:
            continue
        seen.add(canonical)
        page = fetcher.fetch(canonical)
        audit.append(
            {
                "url": canonical,
                "depth": depth,
                "status_code": page.status_code,
                "error": page.error,
                "text_length": len(page.text),
            }
        )
        if page.error:
            continue
        if re.search(record.page_pattern, canonical):
            accepted.append(canonical)
            continue
        if depth >= record.inventory_depth or "html" not in page.content_type:
            continue
        for child in page.links:
            if not child:
                continue
            child_host = urllib.parse.urlsplit(child).netloc.removeprefix("www.")
            if child_host.casefold() != record.domain.casefold():
                continue
            if child not in seen:
                pending.append((child, depth + 1))
    return sorted(dict.fromkeys(accepted)), audit


def read_csv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def write_csv(
    path: Path, rows: Sequence[Mapping[str, object]], columns: Sequence[str] = ()
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(columns or (list(rows[0]) if rows else []))
    with path.open("w", encoding="utf-8", newline="") as handle:
        if not fieldnames:
            return
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def load_yaml(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise TypeError(f"{path} must contain a mapping")
    return payload


def promote_reviewed(
    candidates: Sequence[Mapping[str, object]],
    reviews: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    decisions = {
        _text(row.get("candidate_id")): _text(row.get("decision")).casefold() for row in reviews
    }
    promoted: list[dict[str, object]] = []
    for raw in candidates:
        row = dict(raw)
        candidate_id = (
            _text(row.get("candidate_id"))
            or _sha_text(
                "|".join(
                    _text(row.get(column))
                    for column in (
                        "accepted_species",
                        "trait_name",
                        "normalized_value",
                        "source_lineage",
                        "supporting_excerpt",
                    )
                )
            )[:24]
        )
        row["candidate_id"] = candidate_id
        if decisions.get(candidate_id) != "accept":
            continue
        row["review_status"] = "accepted_manual_audit"
        promoted.append(row)
    return promoted


@app.command("extract")
def extract_command(
    leads_csv: Path = typer.Option(..., exists=True),
    master_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_yaml: Path = typer.Option(Path("config/open_web_trait_acquisition.yml"), exists=True),
    registry_yaml: Path = typer.Option(Path("config/open_web_domain_registry.yml"), exists=True),
    synonym_csv: Path | None = typer.Option(None),
    pause_seconds: float = typer.Option(0.25, min=0, max=30),
) -> None:
    master = read_csv(master_csv)
    synonyms = read_csv(synonym_csv) if synonym_csv else []
    result = process_pages(
        read_csv(leads_csv),
        registry=DomainRegistry.load(registry_yaml),
        resolver=StrictNameResolver(master, synonyms),
        config=load_yaml(config_yaml),
        fetcher=PageFetcher(pause_seconds=pause_seconds),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    filenames = {
        "candidates": "open_web_evidence_candidates.csv",
        "lineage_duplicates": "source_lineage_duplicates.csv",
        "identity_audit": "species_identity_audit.csv",
        "fetch_audit": "page_fetch_audit.csv",
        "domain_audit": "domain_outcome_audit.csv",
    }
    for key, filename in filenames.items():
        rows = result[key]
        if key == "candidates":
            for row in rows:
                row["candidate_id"] = _sha_text(
                    "|".join(
                        _text(row.get(column))
                        for column in (
                            "accepted_species",
                            "trait_name",
                            "normalized_value",
                            "source_lineage",
                            "supporting_excerpt",
                        )
                    )
                )[:24]
        write_csv(output_dir / filename, rows)
    report = {
        "contract": CONTRACT,
        "completed_at_utc": _now(),
        "leads": len(read_csv(leads_csv)),
        "pages_fetched": len(result["fetch_audit"]),
        "identity_audit_rows": len(result["identity_audit"]),
        "candidates_pending_review": len(result["candidates"]),
        "lineage_duplicates_removed": len(result["lineage_duplicates"]),
        "accepted_evidence_rows": 0,
        "policy": {
            "snippets_used_as_evidence": False,
            "unknown_domains_accepted": False,
            "fuzzy_names_accepted": False,
            "family_inference": False,
            "global_fallback": False,
            "manual_review_required_for_promotion": True,
        },
    }
    (output_dir / "extraction_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    app()
