"""Pilot broad-web Medium acquisition and audited genus-inference potential.

The search engine is only a URL discovery mechanism.  Accepted evidence always
comes from a fetched publisher/institution page or PDF, never from a result
snippet.  A second lane inspects exact GBIF checklist usages because useful
descriptions may be attached to a source checklist rather than the backbone.
"""

from __future__ import annotations

import hashlib
import io
import ipaddress
import json
import math
import re
import socket
import time
import urllib.robotparser
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Any
from urllib.parse import parse_qs, quote_plus, unquote, urlsplit, urlunsplit

import httpx
import pandas as pd
import typer
from bs4 import BeautifulSoup

from island_v2.direct_web_reacquisition import (
    DIRECT_EVIDENCE_COLUMNS,
    _file_sha256,
    _text,
    resolve_wfo_exact_synonyms,
)
from island_v2.integrated_trait_coverage import TRAIT_TO_AXIS
from island_v2.v1_category_search import (
    EVIDENCE_COLUMNS as V1_EVIDENCE_COLUMNS,
)
from island_v2.v1_category_search import (
    SOURCE_COLUMNS,
    V2_FIELD_VALUE_MAP,
    extract_evidence_from_sources,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT_VERSION = "broad_web_medium_pilot_v1"
USER_AGENT = "island-floral-traits/0.1 (broad-web-medium-pilot)"
GBIF_API = "https://api.gbif.org/v1"
DUCKDUCKGO_HTML = "https://html.duckduckgo.com/html/"
AXES = (
    "floral_structural_complexity",
    "flower_colour",
    "reproductive_assurance",
)
AXIS_QUERY_TERMS = {
    "floral_structural_complexity": (
        '"flower shape" corolla symmetry tubular campanulate actinomorphic '
        "zygomorphic flower diameter"
    ),
    "flower_colour": '"flower color" "flower colour" corolla petal colour',
    "reproductive_assurance": (
        '"self-compatible" "self-incompatible" autogamy "autonomous selfing" '
        '"mating system" outcrossing'
    ),
}
AXIS_MARKERS = {
    "floral_structural_complexity": re.compile(
        r"\b(?:flower|flowers|floral|floret|florets|corolla|corollas|petal|petals|"
        r"inflorescence|inflorescences)\b.{0,180}\b(?:actinomorphic|zygomorphic|"
        r"radially symmetric|bilaterally symmetric|tubular|tube[- ]shaped|"
        r"campanulate|bell[- ]shaped|funnelform|funnel[- ]shaped|"
        r"trumpet[- ]shaped|infundibuliform|rotate|cup[- ]shaped|bowl[- ]shaped|"
        r"papilionaceous|spurred|urceolate|raceme|spike|panicle|capitulum|"
        r"\d+(?:\.\d+)?\s*(?:mm|cm)\s*(?:long|wide|across|diameter))\b",
        re.IGNORECASE,
    ),
    "flower_colour": re.compile(
        r"\b(?:flower|flowers|floral|floret|florets|corolla|corollas|petal|petals|"
        r"tepal|tepals|perianth|inflorescence)\b.{0,140}\b(?:white|cream|yellow|"
        r"orange|red|pink|purple|violet|blue|green|brown|scarlet|crimson|"
        r"lavender|mauve|lilac|golden)\b",
        re.IGNORECASE,
    ),
    "reproductive_assurance": re.compile(
        r"\b(?:self[- ]compatible|self[- ]compatibility|self[- ]incompatible|"
        r"self[- ]incompatibility|autogam(?:y|ous)|autonomous selfing|"
        r"spontaneous selfing|mixed mating|facultative selfing|"
        r"predominantly selfing|predominantly outcrossing|obligate selfing|"
        r"cleistogam(?:y|ous))\b",
        re.IGNORECASE,
    ),
}
DOI_RE = re.compile(r"\b10\.\d{4,9}/[-._;()/:A-Z0-9]+\b", re.IGNORECASE)
MEASUREMENT_RE = re.compile(
    r"\b(?:flower|flowers|corolla)\b.{0,80}?"
    r"(?P<low>\d+(?:\.\d+)?)\s*(?:[-–—]\s*(?P<high>\d+(?:\.\d+)?))?\s*"
    r"(?P<unit>mm|cm)\s*(?:long|wide|across|in diameter|diameter)\b",
    re.IGNORECASE,
)
BLOCKED_HOST_FRAGMENTS = {
    "amazon.",
    "ebay.",
    "etsy.",
    "facebook.",
    "instagram.",
    "pinterest.",
    "reddit.",
    "tiktok.",
    "x.com",
    "youtube.",
}
SCIENTIFIC_HOST_FRAGMENTS = {
    "biomedcentral.com",
    "cambridge.org",
    "doi.org",
    "frontiersin.org",
    "jstor.org",
    "nature.com",
    "onlinelibrary.wiley.com",
    "oup.com",
    "plos.org",
    "pubmed.ncbi.nlm.nih.gov",
    "researchsquare.com",
    "sciencedirect.com",
    "springer.com",
    "tandfonline.com",
}
CURATED_HOST_FRAGMENTS = {
    "biodiversitylibrary.org",
    "botanicgardens.org",
    "efloras.org",
    "gbif.org",
    "gobotany.nativeplanttrust.org",
    "infoflora.ch",
    "kew.org",
    "missouribotanicalgarden.org",
    "natureserve.org",
    "nzpcn.org.nz",
    "pladias.cz",
    "plants.ces.ncsu.edu",
    "sanbi.org",
    "treatment.plazi.org",
    "worldfloraonline.org",
}
UNTRUSTED_SOURCE_TYPE = "general_web_unreviewed"
FETCHED_SOURCE_COLUMNS = SOURCE_COLUMNS + [
    "search_name",
    "name_role",
    "name_resolution_lineage",
    "target_axis",
    "source_host",
    "retrieved_at_epoch",
    "raw_content_sha256",
    "rights_status",
    "robots_status",
]
DISCOVERY_COLUMNS = [
    "accepted_species",
    "target_axis",
    "search_name",
    "name_role",
    "source_url",
    "source_host",
    "source_type",
    "fetch_status",
]
ERROR_COLUMNS = [
    "accepted_species",
    "target_axis",
    "search_name",
    "source",
    "source_url",
    "error_class",
    "message",
]
STRUCTURED_MATCH_COLUMNS = [
    "accepted_species",
    "provider",
    "source_url",
    "index_match_method",
    "fetch_status",
    "evidence_rows",
]
GBIF_USAGE_COLUMNS = [
    "accepted_species",
    "target_axis",
    "search_name",
    "usage_key",
    "dataset_key",
    "canonical_name",
    "taxonomic_status",
    "references",
    "description_records",
    "status",
]
PLADIAS_VALUE_MAP = {
    ("Flower symmetry", "actinomorphic"): ("floral_symmetry", "actinomorphic"),
    ("Flower symmetry", "zygomorphic"): ("floral_symmetry", "zygomorphic"),
    ("Generative reproduction type", "allogamy"): (
        "mating_system",
        "predominantly_outcrossing",
    ),
    ("Generative reproduction type", "allogamy self-incompatibility"): (
        "self_incompatibility",
        "SI",
    ),
    ("Generative reproduction type", "facultative allogamy"): (
        "mating_system",
        "mixed_mating",
    ),
    ("Generative reproduction type", "autogamy"): (
        "mating_system",
        "predominantly_selfing",
    ),
    ("Generative reproduction type", "facultative autogamy"): (
        "mating_system",
        "predominantly_selfing",
    ),
    ("Generative reproduction type", "mixed mating"): (
        "mating_system",
        "mixed_mating",
    ),
}
COLOUR_TOKEN_MAP = {
    "white": "white",
    "cream": "white",
    "yellow": "yellow_orange",
    "orange": "yellow_orange",
    "red": "red_pink",
    "pink": "red_pink",
    "purple": "blue_purple",
    "violet": "blue_purple",
    "blue": "blue_purple",
    "green": "green_brown_inconspicuous",
    "brown": "green_brown_inconspicuous",
}


@dataclass(frozen=True)
class FetchResult:
    url: str
    title: str
    text: str
    content_sha256: str
    content_type: str
    status: str
    error: str = ""


def _strict_binomial(value: object) -> str:
    name = _text(value)
    return name if re.fullmatch(r"[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+", name) else ""


def _canonical_url(value: object) -> str:
    raw = _text(value)
    parts = urlsplit(raw)
    if parts.scheme not in {"http", "https"} or not parts.hostname:
        return ""
    host = parts.hostname.casefold().rstrip(".")
    port = f":{parts.port}" if parts.port and parts.port not in {80, 443} else ""
    path = re.sub(r"/+", "/", unquote(parts.path or "/")).rstrip("/") or "/"
    return urlunsplit((parts.scheme.casefold(), host + port, path, "", ""))


def _unwrap_search_url(value: object) -> str:
    raw = _text(value)
    if raw.startswith("//"):
        raw = "https:" + raw
    parsed = urlsplit(raw)
    if parsed.hostname and parsed.hostname.endswith("duckduckgo.com"):
        target = parse_qs(parsed.query).get("uddg", [""])[0]
        if target:
            return _canonical_url(unquote(target))
    return _canonical_url(raw)


def parse_duckduckgo_urls(html: str, *, max_results: int) -> list[str]:
    """Return destination URLs only; snippets and ranks are never evidence."""
    soup = BeautifulSoup(html, "html.parser")
    urls: list[str] = []
    for result in soup.select(".result"):
        link = result.select_one(".result__a")
        url = _unwrap_search_url(link.get("href", "") if link else "")
        if not url or url in urls:
            continue
        urls.append(url)
        if len(urls) >= max_results:
            break
    return urls


def _host_class(url: str) -> tuple[str, str]:
    host = (urlsplit(url).hostname or "").casefold()
    if not host or any(fragment in host for fragment in BLOCKED_HOST_FRAGMENTS):
        return "", "blocked"
    if (
        host.endswith((".gov", ".edu", ".ac.uk"))
        or ".gov." in host
        or ".edu." in host
        or ".ac." in host
    ):
        return "general_web_institutional", "medium"
    if any(fragment in host for fragment in SCIENTIFIC_HOST_FRAGMENTS):
        return "general_web_scientific", "medium"
    if any(fragment in host for fragment in CURATED_HOST_FRAGMENTS):
        return "general_web_curated", "medium"
    return UNTRUSTED_SOURCE_TYPE, "pending"


def _safe_public_url(url: str) -> bool:
    parts = urlsplit(url)
    host = parts.hostname
    if parts.scheme not in {"http", "https"} or not host:
        return False
    if host.casefold() in {"localhost", "localhost.localdomain"}:
        return False
    try:
        addresses = {item[4][0] for item in socket.getaddrinfo(host, parts.port or 443)}
    except socket.gaierror:
        return False
    for address in addresses:
        try:
            ip = ipaddress.ip_address(address)
        except ValueError:
            return False
        if not ip.is_global:
            return False
    return True


def _robots_allowed(
    client: httpx.Client,
    url: str,
    cache: dict[str, tuple[bool, str]],
) -> tuple[bool, str]:
    parts = urlsplit(url)
    origin = f"{parts.scheme}://{parts.netloc}"
    if origin in cache:
        return cache[origin]
    robots_url = origin + "/robots.txt"
    try:
        response = client.get(robots_url)
        if response.status_code == 404:
            cache[origin] = (True, "not_found_allow")
            return cache[origin]
        response.raise_for_status()
        parser = urllib.robotparser.RobotFileParser()
        parser.set_url(robots_url)
        parser.parse(response.text.splitlines())
        allowed = parser.can_fetch(USER_AGENT, url)
        cache[origin] = (allowed, "allowed" if allowed else "disallowed")
    except Exception as exc:  # noqa: BLE001
        cache[origin] = (False, f"robots_unavailable:{type(exc).__name__}")
    return cache[origin]


def _visible_html(response: httpx.Response) -> tuple[str, str]:
    soup = BeautifulSoup(response.text, "html.parser")
    title = _text(soup.title.get_text(" ", strip=True) if soup.title else "")
    for element in soup(["script", "style", "noscript", "svg", "canvas"]):
        element.decompose()
    return title, _text(soup.get_text(" ", strip=True))


def _pdf_text(content: bytes, *, max_pages: int = 40) -> tuple[str, str]:
    from pypdf import PdfReader

    reader = PdfReader(io.BytesIO(content))
    pages = [_text(page.extract_text() or "") for page in reader.pages[:max_pages]]
    title = _text((reader.metadata or {}).get("/Title"))
    return title, _text(" ".join(pages))


def fetch_document(client: httpx.Client, url: str, *, max_bytes: int = 12_000_000) -> FetchResult:
    if not _safe_public_url(url):
        return FetchResult(url, "", "", "", "", "rejected", "non_public_or_invalid_url")
    try:
        response = client.get(url)
        response.raise_for_status()
        content = response.content
        if len(content) > max_bytes:
            return FetchResult(url, "", "", "", "", "rejected", "content_too_large")
        content_type = response.headers.get("content-type", "").casefold()
        digest = hashlib.sha256(content).hexdigest()
        if "pdf" in content_type or urlsplit(str(response.url)).path.casefold().endswith(".pdf"):
            title, text = _pdf_text(content)
        elif "html" in content_type or not content_type:
            title, text = _visible_html(response)
        elif content_type.startswith("text/"):
            title, text = "", _text(response.text)
        else:
            return FetchResult(
                str(response.url), "", "", digest, content_type, "rejected", "unsupported_type"
            )
        if not text:
            return FetchResult(
                str(response.url), title, "", digest, content_type, "rejected", "empty_text"
            )
        return FetchResult(str(response.url), title, text, digest, content_type, "fetched")
    except Exception as exc:  # noqa: BLE001
        return FetchResult(url, "", "", "", "", "error", f"{type(exc).__name__}:{exc}")


def _contains_exact_name(text: str, name: str) -> bool:
    return bool(
        re.search(
            rf"(?<![A-Za-z.-]){re.escape(name)}(?![A-Za-z.-])",
            text,
            flags=re.IGNORECASE,
        )
    )


def extract_direct_passages(
    *,
    accepted_species: str,
    search_name: str,
    target_axis: str,
    title: str,
    text: str,
    max_passages: int = 12,
) -> list[str]:
    """Keep explicit target-axis assertions from an exact taxon page/article."""
    marker = AXIS_MARKERS[target_axis]
    title_exact = _contains_exact_name(title, search_name)
    sentences = [
        _text(sentence)
        for sentence in re.split(r"(?<=[.!?;])\s+|\s*\n+\s*", text)
        if _text(sentence)
    ]
    passages: list[str] = []
    for sentence in sentences:
        if not marker.search(sentence):
            continue
        if not title_exact and not _contains_exact_name(sentence, search_name):
            continue
        if sentence not in passages:
            passages.append(sentence[:1600])
        if len(passages) >= max_passages:
            break
    if not passages and title_exact:
        # HTML flattening can merge headings and descriptions into long chunks.
        for match in marker.finditer(text):
            start = max(0, match.start() - 250)
            end = min(len(text), match.end() + 450)
            excerpt = _text(text[start:end])
            if excerpt and excerpt not in passages:
                passages.append(excerpt)
            if len(passages) >= max_passages:
                break
    if not passages:
        return []
    joined = " ".join(passages)
    if not (_contains_exact_name(title, search_name) or _contains_exact_name(joined, search_name)):
        return []
    # accepted_species may differ only when WFO supplied this exact synonym.
    if search_name.casefold() == accepted_species.casefold():
        return passages
    return passages


def _doi_lineage(*values: object) -> str:
    for value in values:
        match = DOI_RE.search(_text(value))
        if match:
            return f"doi:{match.group(0).rstrip('.,;)').casefold()}"
    return ""


def _url_lineage(url: str, content_sha256: str) -> str:
    canonical = _canonical_url(url)
    return f"url:{canonical}" if canonical else f"sha256:{content_sha256}"


def _colour_value(raw: object) -> str:
    text = _text(raw).casefold()
    states = {
        canonical
        for token, canonical in COLOUR_TOKEN_MAP.items()
        if re.search(rf"\b{re.escape(token)}\b", text)
    }
    if len(states) > 1:
        return "multicolored_variable"
    return next(iter(states), "")


def _measurement_class(raw: object) -> str:
    match = re.search(
        r"(?P<low>\d+(?:\.\d+)?)\s*(?:[-–—]\s*(?P<high>\d+(?:\.\d+)?))?\s*"
        r"(?P<unit>mm|cm)\b",
        _text(raw),
        flags=re.IGNORECASE,
    )
    if not match:
        return ""
    low = float(match.group("low"))
    high = float(match.group("high") or low)
    midpoint_mm = ((low + high) / 2) * (10 if match.group("unit").casefold() == "cm" else 1)
    return next(
        label
        for maximum, label in (
            (2, "very_small"),
            (10, "small"),
            (30, "medium"),
            (100, "large"),
            (float("inf"), "very_large"),
        )
        if midpoint_mm <= maximum
    )


def _structured_evidence_row(
    *,
    species: str,
    trait_name: str,
    value: str,
    provider: str,
    url: str,
    feature: str,
    raw_value: str,
    content_sha256: str,
) -> dict[str, str]:
    axis = TRAIT_TO_AXIS[trait_name]
    excerpt = f"{feature}: {raw_value}"
    return {
        "accepted_species": species,
        "axis": axis,
        "trait_name": trait_name,
        "normalized_value": value,
        "evidence_quality": "medium",
        "evidence_scope": "species_direct",
        "source_group": "structured_institutional_web",
        "source_provider": provider,
        "source_url": url,
        "source_citation": (
            f"{provider}; exact species page; feature={feature}; "
            f"publisher_content_sha256={content_sha256}"
        ),
        "source_excerpt": excerpt,
        "source_record_id": hashlib.sha256(
            f"{provider}|{species}|{feature}|{raw_value}".encode()
        ).hexdigest(),
        "source_lineage": f"url:{_canonical_url(url)}",
        "name_resolution_lineage": "",
        "name_match_method": "accepted_name_exact",
        "inference_rule": "",
    }


def parse_pladias_evidence(
    html: str,
    *,
    accepted_species: str,
    source_url: str,
) -> pd.DataFrame:
    soup = BeautifulSoup(html, "html.parser")
    title = _text(soup.title.get_text(" ", strip=True) if soup.title else "")
    if not _contains_exact_name(title, accepted_species):
        return pd.DataFrame(columns=DIRECT_EVIDENCE_COLUMNS)
    digest = hashlib.sha256(html.encode("utf-8")).hexdigest()
    rows: list[dict[str, str]] = []
    for item in soup.select("li.featureDetail"):
        labels = item.select("div.d-flex > div > span.font-italic")
        value_node = item.select_one("div.d-flex > div > b")
        if not labels or value_node is None:
            continue
        feature = " ".join(node.get_text(" ", strip=True) for node in labels).strip(" :")
        raw_value = _text(value_node.get_text(" ", strip=True))
        mapped: tuple[str, str] | None = None
        if feature == "Flower colour":
            colour = _colour_value(raw_value)
            mapped = ("flower_primary_color", colour) if colour else None
        else:
            mapped = PLADIAS_VALUE_MAP.get((feature, raw_value.casefold()))
        if mapped is None:
            continue
        trait_name, value = mapped
        rows.append(
            _structured_evidence_row(
                species=accepted_species,
                trait_name=trait_name,
                value=value,
                provider="Pladias",
                url=source_url,
                feature=feature,
                raw_value=raw_value,
                content_sha256=digest,
            )
        )
    return pd.DataFrame(rows, columns=DIRECT_EVIDENCE_COLUMNS).drop_duplicates()


def parse_gobotany_evidence(
    html: str,
    *,
    accepted_species: str,
    source_url: str,
) -> pd.DataFrame:
    soup = BeautifulSoup(html, "html.parser")
    title = _text(soup.title.get_text(" ", strip=True) if soup.title else "")
    if not _contains_exact_name(title, accepted_species):
        return pd.DataFrame(columns=DIRECT_EVIDENCE_COLUMNS)
    digest = hashlib.sha256(html.encode("utf-8")).hexdigest()
    rows: list[dict[str, str]] = []
    for item in soup.select("dl"):
        feature_node = item.find("dt")
        value_node = item.find("dd")
        if feature_node is None or value_node is None:
            continue
        feature = _text(feature_node.get_text(" ", strip=True))
        raw_value = _text(value_node.get_text(" ", strip=True))
        mapped: tuple[str, str] | None = None
        if feature == "Flower petal color":
            colour = _colour_value(raw_value)
            mapped = ("flower_primary_color", colour) if colour else None
        elif feature == "Flower symmetry":
            folded = raw_value.casefold()
            radial = "radially symmetrical" in folded or "two or more ways" in folded
            bilateral = "bilaterally symmetrical" in folded or "only one way" in folded
            if radial and not bilateral:
                mapped = ("floral_symmetry", "actinomorphic")
            elif bilateral and not radial:
                mapped = ("floral_symmetry", "zygomorphic")
        elif feature == "Corolla morphology":
            folded = raw_value.casefold()
            for marker, value in (
                ("campanulate", "bell_campanulate"),
                ("bell-shaped", "bell_campanulate"),
                ("funnelform", "funnel_trumpet"),
                ("trumpet", "funnel_trumpet"),
                ("tubular", "tubular"),
                ("rotate", "rotate"),
                ("urceolate", "urceolate"),
            ):
                if marker in folded:
                    mapped = ("floral_form", value)
                    break
        elif feature == "Flower diameter":
            size = _measurement_class(raw_value)
            mapped = ("flower_size_class", size) if size else None
        elif feature == "Cleistogamous flowers":
            folded = raw_value.casefold()
            if re.search(r"\b(?:no|none|absent|there are no)\b", folded):
                mapped = ("cleistogamy", "absent")
            elif re.search(r"\b(?:yes|present|there are)\b", folded):
                mapped = ("cleistogamy", "present")
        if mapped is None:
            continue
        trait_name, value = mapped
        rows.append(
            _structured_evidence_row(
                species=accepted_species,
                trait_name=trait_name,
                value=value,
                provider="Go Botany / Native Plant Trust",
                url=source_url,
                feature=feature,
                raw_value=raw_value,
                content_sha256=digest,
            )
        )
    return pd.DataFrame(rows, columns=DIRECT_EVIDENCE_COLUMNS).drop_duplicates()


def _gobotany_sitemap_index(client: httpx.Client) -> dict[str, str]:
    sitemap_url = "https://gobotany.nativeplanttrust.org/sitemap.txt"
    response = client.get(sitemap_url)
    response.raise_for_status()
    index: dict[str, str] = {}
    for raw_url in response.text.splitlines():
        url = _canonical_url(raw_url)
        parts = [unquote(part) for part in urlsplit(url).path.strip("/").split("/")]
        if len(parts) != 3 or parts[0].casefold() != "species":
            continue
        species = f"{parts[1].capitalize()} {parts[2].casefold()}"
        if _strict_binomial(species):
            index.setdefault(species.casefold(), url + "/")
    return index


def _structured_domain_evidence(
    targets: pd.DataFrame,
    *,
    workers: int,
    include_pladias: bool,
    include_gobotany: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    headers = {
        "User-Agent": USER_AGENT,
        "Accept": "text/html,application/xhtml+xml",
    }
    sitemap_errors: list[dict[str, str]] = []
    gobotany_index: dict[str, str] = {}
    if include_gobotany:
        try:
            with httpx.Client(
                timeout=httpx.Timeout(30, connect=10),
                follow_redirects=True,
                headers=headers,
            ) as client:
                gobotany_index = _gobotany_sitemap_index(client)
        except Exception as exc:  # noqa: BLE001
            sitemap_errors.append(
                {
                    "accepted_species": "",
                    "target_axis": "",
                    "search_name": "",
                    "source": "gobotany_sitemap",
                    "source_url": "https://gobotany.nativeplanttrust.org/sitemap.txt",
                    "error_class": type(exc).__name__,
                    "message": _text(exc),
                }
            )

    def fetch_species(
        species: str,
    ) -> tuple[list[pd.DataFrame], list[dict[str, str]], list[dict[str, str]]]:
        requests: list[tuple[str, str, Any, str]] = []
        if include_pladias:
            requests.append(
                (
                    "pladias",
                    f"https://pladias.cz/en/taxon/data/{quote_plus(species)}",
                    parse_pladias_evidence,
                    "exact_accepted_name_url",
                )
            )
        gobotany_url = gobotany_index.get(species.casefold(), "")
        if gobotany_url:
            requests.append(
                (
                    "gobotany",
                    gobotany_url,
                    parse_gobotany_evidence,
                    "official_sitemap_exact_binomial",
                )
            )
        frames: list[pd.DataFrame] = []
        errors: list[dict[str, str]] = []
        matches: list[dict[str, str]] = []
        with httpx.Client(
            timeout=httpx.Timeout(20, connect=10),
            follow_redirects=True,
            headers=headers,
        ) as client:
            for provider, url, parser, match_method in requests:
                match = {
                    "accepted_species": species,
                    "provider": provider,
                    "source_url": url,
                    "index_match_method": match_method,
                    "fetch_status": "",
                    "evidence_rows": "0",
                }
                try:
                    response = client.get(url)
                    if response.status_code in {302, 404}:
                        match["fetch_status"] = f"http_{response.status_code}"
                        matches.append(match)
                        continue
                    response.raise_for_status()
                    frame = parser(
                        response.text,
                        accepted_species=species,
                        source_url=str(response.url),
                    )
                    match["fetch_status"] = "fetched"
                    match["evidence_rows"] = str(len(frame))
                    matches.append(match)
                    if not frame.empty:
                        frames.append(frame)
                except Exception as exc:  # noqa: BLE001
                    match["fetch_status"] = "error"
                    matches.append(match)
                    errors.append(
                        {
                            "accepted_species": species,
                            "target_axis": "",
                            "search_name": species,
                            "source": provider,
                            "source_url": url,
                            "error_class": type(exc).__name__,
                            "message": _text(exc),
                        }
                    )
        return frames, errors, matches

    evidence_parts: list[pd.DataFrame] = []
    error_rows: list[dict[str, str]] = list(sitemap_errors)
    match_rows: list[dict[str, str]] = []
    target_species = set(targets["accepted_species"].map(_text))
    unique_species = sorted(
        target_species
        if include_pladias
        else {species for species in target_species if species.casefold() in gobotany_index}
    )
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(fetch_species, species): species for species in unique_species}
        for future in as_completed(futures):
            try:
                frames, errors, matches = future.result()
                evidence_parts.extend(frames)
                error_rows.extend(errors)
                match_rows.extend(matches)
            except Exception as exc:  # noqa: BLE001
                species = futures[future]
                error_rows.append(
                    {
                        "accepted_species": species,
                        "target_axis": "",
                        "search_name": species,
                        "source": "structured_domains",
                        "source_url": "",
                        "error_class": type(exc).__name__,
                        "message": _text(exc),
                    }
                )
    evidence = (
        pd.concat(evidence_parts, ignore_index=True).drop_duplicates()
        if evidence_parts
        else pd.DataFrame(columns=DIRECT_EVIDENCE_COLUMNS)
    )
    return (
        evidence,
        pd.DataFrame(error_rows, columns=ERROR_COLUMNS),
        pd.DataFrame(match_rows, columns=STRUCTURED_MATCH_COLUMNS),
    )


def select_pilot_targets(
    coverage: pd.DataFrame,
    master: pd.DataFrame,
    direct_evidence: pd.DataFrame,
    *,
    limit: int = 100,
) -> pd.DataFrame:
    """Select an axis-balanced, genus-diverse, high-yield 100-species pilot."""
    required_coverage = {"accepted_species", "axis", "after_quality"}
    if missing := required_coverage.difference(coverage.columns):
        raise ValueError(f"coverage is missing columns: {sorted(missing)}")
    if coverage.duplicated(["accepted_species", "axis"]).any():
        raise ValueError("coverage contains duplicate accepted_species x axis")
    if limit < 3:
        raise ValueError("limit must be at least three")

    cov = coverage[list(required_coverage)].copy().fillna("")
    cov["accepted_species"] = cov["accepted_species"].map(_text)
    cov["axis"] = cov["axis"].map(_text)
    cov["after_quality"] = cov["after_quality"].map(lambda value: _text(value).casefold())
    cov["genus"] = cov["accepted_species"].str.split().str[0]

    direct = direct_evidence.copy().fillna("")
    for column in ("accepted_species", "axis", "quality", "evidence_scope"):
        if column not in direct:
            direct[column] = ""
    direct = direct.loc[
        direct["quality"].isin({"high", "medium"})
        & direct["evidence_scope"].isin({"species_direct", "synonym_direct"})
    ].copy()
    direct["genus"] = direct["accepted_species"].map(_text).str.split().str[0]
    donors = (
        direct.drop_duplicates(["accepted_species", "axis"])
        .groupby(["genus", "axis"])["accepted_species"]
        .nunique()
        .rename("direct_donor_species")
    )
    unresolved_in_genus = (
        cov.loc[cov["after_quality"].eq("")]
        .groupby(["genus", "axis"])["accepted_species"]
        .nunique()
        .rename("unresolved_congeneric_cells")
    )

    meta = master.copy().fillna("")
    for column in ("accepted_species", "genus", "family", "n_records"):
        if column not in meta:
            meta[column] = ""
    meta = meta.drop_duplicates("accepted_species")
    meta["accepted_species"] = meta["accepted_species"].map(_text)
    meta["n_records"] = pd.to_numeric(meta["n_records"], errors="coerce").fillna(0)

    candidates = cov.loc[cov["after_quality"].eq("")].merge(
        meta[["accepted_species", "family", "n_records"]],
        on="accepted_species",
        how="left",
    )
    candidates = candidates.join(donors, on=["genus", "axis"]).join(
        unresolved_in_genus,
        on=["genus", "axis"],
    )
    candidates[["direct_donor_species", "unresolved_congeneric_cells"]] = candidates[
        ["direct_donor_species", "unresolved_congeneric_cells"]
    ].fillna(0)
    candidates = candidates.loc[candidates["accepted_species"].map(_strict_binomial).ne("")]
    candidates["score"] = (
        candidates["n_records"].map(lambda value: math.log1p(float(value)))
        + candidates["unresolved_congeneric_cells"].map(
            lambda value: 2.0 * math.log1p(float(value))
        )
        + candidates["direct_donor_species"].map(
            lambda value: 8.0 if int(value) == 2 else (3.0 if int(value) == 1 else 0.0)
        )
    )
    quotas = {axis: limit // len(AXES) for axis in AXES}
    for axis in AXES[: limit % len(AXES)]:
        quotas[axis] += 1

    selected_parts: list[pd.DataFrame] = []
    used_species: set[str] = set()
    for axis in AXES:
        pool = candidates.loc[
            candidates["axis"].eq(axis) & ~candidates["accepted_species"].isin(used_species)
        ].sort_values(
            ["score", "n_records", "accepted_species"],
            ascending=[False, False, True],
        )
        # First pass: one focal species per genus maximizes both source and
        # genus-rule diversity.  A second pass fills any quota shortfall.
        diverse = pool.drop_duplicates("genus", keep="first").head(quotas[axis])
        if len(diverse) < quotas[axis]:
            used = set(diverse["accepted_species"])
            extra = pool.loc[~pool["accepted_species"].isin(used)].head(quotas[axis] - len(diverse))
            diverse = pd.concat([diverse, extra], ignore_index=True)
        selected_parts.append(diverse)
        used_species.update(diverse["accepted_species"])
    selected = pd.concat(selected_parts, ignore_index=True)
    selected = selected.rename(columns={"axis": "target_axis"})
    selected["selection_reason"] = "axis_balanced_webability_plus_genus_leverage"
    return selected[
        [
            "accepted_species",
            "target_axis",
            "genus",
            "family",
            "n_records",
            "direct_donor_species",
            "unresolved_congeneric_cells",
            "score",
            "selection_reason",
        ]
    ].sort_values(["target_axis", "score"], ascending=[True, False])


def select_all_unresolved_targets(
    coverage: pd.DataFrame,
    master: pd.DataFrame,
    *,
    offset: int = 0,
    limit: int = 0,
) -> pd.DataFrame:
    """Select every species with at least one unresolved axis exactly once."""
    required_coverage = {"accepted_species", "axis", "after_quality"}
    if missing := required_coverage.difference(coverage.columns):
        raise ValueError(f"coverage is missing columns: {sorted(missing)}")
    if coverage.duplicated(["accepted_species", "axis"]).any():
        raise ValueError("coverage contains duplicate accepted_species x axis")
    if offset < 0 or limit < 0:
        raise ValueError("offset and limit must be non-negative")

    missing = coverage.copy().fillna("")
    missing["accepted_species"] = missing["accepted_species"].map(_text)
    missing["axis"] = missing["axis"].map(_text)
    missing["after_quality"] = missing["after_quality"].map(lambda value: _text(value).casefold())
    missing = missing.loc[
        missing["after_quality"].eq("") & missing["accepted_species"].map(_strict_binomial).ne("")
    ]
    grouped = (
        missing.groupby("accepted_species", sort=True)["axis"]
        .agg(lambda values: "|".join(sorted(set(values))))
        .rename("target_axis")
        .reset_index()
    )
    meta = master.copy().fillna("")
    for column in ("accepted_species", "genus", "family", "n_records"):
        if column not in meta:
            meta[column] = ""
    meta = meta.drop_duplicates("accepted_species")
    selected = grouped.merge(
        meta[["accepted_species", "genus", "family", "n_records"]],
        on="accepted_species",
        how="left",
    )
    selected["genus"] = selected["genus"].where(
        selected["genus"].map(_text).ne(""),
        selected["accepted_species"].str.split().str[0],
    )
    selected["n_records"] = pd.to_numeric(
        selected["n_records"],
        errors="coerce",
    ).fillna(0)
    selected["direct_donor_species"] = ""
    selected["unresolved_congeneric_cells"] = ""
    selected["score"] = selected["n_records"].map(lambda value: math.log1p(float(value)))
    selected["selection_reason"] = "all_species_with_any_unresolved_axis"
    selected = selected.sort_values("accepted_species").iloc[
        offset : (offset + limit) if limit else None
    ]
    return selected[
        [
            "accepted_species",
            "target_axis",
            "genus",
            "family",
            "n_records",
            "direct_donor_species",
            "unresolved_congeneric_cells",
            "score",
            "selection_reason",
        ]
    ]


def _search_names(
    targets: pd.DataFrame,
    *,
    resolve_synonyms: bool,
    max_aliases_per_species: int,
    wfo_pause_seconds: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    accepted = pd.DataFrame(
        {
            "accepted_species": targets["accepted_species"],
            "search_name": targets["accepted_species"],
            "name_role": "accepted",
            "name_resolution_lineage": "",
        }
    )
    if not resolve_synonyms:
        return accepted, pd.DataFrame(columns=["accepted_species", "source", "error"])
    aliases, errors = resolve_wfo_exact_synonyms(
        targets["accepted_species"].tolist(),
        min_interval_seconds=wfo_pause_seconds,
        max_aliases_per_species=max_aliases_per_species,
    )
    if aliases.empty:
        return accepted, errors
    aliases = aliases[
        [
            "accepted_species",
            "search_name",
            "name_role",
            "name_resolution_lineage",
        ]
    ]
    names = pd.concat([accepted, aliases], ignore_index=True).drop_duplicates(
        ["accepted_species", "search_name"]
    )
    return names, errors


def _gbif_exact_usages(
    client: httpx.Client,
    *,
    accepted_species: str,
    search_name: str,
    target_axis: str,
    max_usages: int,
) -> list[dict[str, str]]:
    response = client.get(
        f"{GBIF_API}/species/search",
        params={"q": search_name, "rank": "SPECIES", "limit": 100},
    )
    response.raise_for_status()
    rows: list[dict[str, str]] = []
    for record in response.json().get("results", []):
        canonical = _text(record.get("canonicalName"))
        if (
            canonical.casefold() != search_name.casefold()
            or _text(record.get("rank")).casefold() != "species"
            or _text(record.get("kingdom")).casefold() != "plantae"
            or not _text(record.get("key"))
        ):
            continue
        rows.append(
            {
                "accepted_species": accepted_species,
                "target_axis": target_axis,
                "search_name": search_name,
                "usage_key": _text(record.get("key")),
                "dataset_key": _text(record.get("datasetKey")),
                "canonical_name": canonical,
                "taxonomic_status": _text(record.get("taxonomicStatus") or record.get("status")),
                "references": _text(record.get("references")),
                "description_records": "0",
                "status": "matched",
            }
        )
        if len(rows) >= max_usages:
            break
    return rows


def _gbif_description_sources(
    targets: pd.DataFrame,
    *,
    max_usages_per_species: int,
    workers: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    usage_rows: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    with httpx.Client(timeout=45, follow_redirects=True, headers=headers) as client:
        for row in targets.itertuples(index=False):
            try:
                usage_rows.extend(
                    _gbif_exact_usages(
                        client,
                        accepted_species=row.accepted_species,
                        search_name=row.accepted_species,
                        target_axis=row.target_axis,
                        max_usages=max_usages_per_species,
                    )
                )
            except Exception as exc:  # noqa: BLE001
                errors.append(
                    {
                        "accepted_species": row.accepted_species,
                        "target_axis": row.target_axis,
                        "search_name": row.accepted_species,
                        "source": "gbif_species_search",
                        "source_url": f"{GBIF_API}/species/search",
                        "error_class": type(exc).__name__,
                        "message": _text(exc),
                    }
                )

    def fetch_usage(usage: dict[str, str]) -> tuple[dict[str, str], list[dict[str, str]]]:
        url = f"{GBIF_API}/species/{usage['usage_key']}/descriptions"
        with httpx.Client(timeout=45, follow_redirects=True, headers=headers) as client:
            response = client.get(url)
            response.raise_for_status()
            descriptions = response.json().get("results", [])
        return usage, descriptions

    source_rows: list[dict[str, str]] = []
    if usage_rows:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(fetch_usage, row): row for row in usage_rows}
            for future in as_completed(futures):
                usage = futures[future]
                try:
                    usage, descriptions = future.result()
                    usage["description_records"] = str(len(descriptions))
                    usage["status"] = "description_hit" if descriptions else "no_description"
                    for record in descriptions:
                        description = _text(record.get("description"))
                        if not description:
                            continue
                        source = _text(record.get("source"))
                        lineage = _doi_lineage(source)
                        source_rows.append(
                            {
                                "accepted_species": usage["accepted_species"],
                                "source_text": description,
                                "source_url": (
                                    f"https://www.gbif.org/species/{usage['usage_key']}"
                                ),
                                "source_citation": (
                                    f"{source} | GBIF usageKey={usage['usage_key']} | "
                                    f"datasetKey={usage['dataset_key']} | "
                                    f"descriptionKey={_text(record.get('key'))} | "
                                    f"original_lineage={lineage}"
                                ),
                                "source_type": "gbif_checklist_description",
                                "evidence_scope": "species_direct",
                                "search_name": usage["search_name"],
                                "name_role": "accepted",
                                "name_resolution_lineage": "",
                                "target_axis": usage["target_axis"],
                                "source_host": "gbif.org",
                                "retrieved_at_epoch": str(int(time.time())),
                                "raw_content_sha256": hashlib.sha256(
                                    description.encode("utf-8")
                                ).hexdigest(),
                                "rights_status": "dataset_specific",
                                "robots_status": "official_api",
                            }
                        )
                except Exception as exc:  # noqa: BLE001
                    usage["status"] = "error"
                    errors.append(
                        {
                            "accepted_species": usage["accepted_species"],
                            "target_axis": usage["target_axis"],
                            "search_name": usage["search_name"],
                            "source": "gbif_description",
                            "source_url": (f"{GBIF_API}/species/{usage['usage_key']}/descriptions"),
                            "error_class": type(exc).__name__,
                            "message": _text(exc),
                        }
                    )
    return (
        pd.DataFrame(source_rows, columns=FETCHED_SOURCE_COLUMNS),
        pd.DataFrame(usage_rows, columns=GBIF_USAGE_COLUMNS),
        pd.DataFrame(errors, columns=ERROR_COLUMNS),
    )


def _general_web_sources(
    targets: pd.DataFrame,
    search_names: pd.DataFrame,
    *,
    max_results_per_query: int,
    max_aliases_to_search: int,
    query_pause_seconds: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    target_axis = targets.set_index("accepted_species")["target_axis"].to_dict()
    names_by_species = {
        species: group.head(1 + max_aliases_to_search).to_dict("records")
        for species, group in search_names.groupby("accepted_species", sort=False)
    }
    source_rows: list[dict[str, str]] = []
    discovery_rows: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    robots_cache: dict[str, tuple[bool, str]] = {}
    headers = {
        "User-Agent": USER_AGENT,
        "Accept": "text/html,application/xhtml+xml,application/pdf,text/plain;q=0.8,*/*;q=0.1",
    }
    with httpx.Client(timeout=45, follow_redirects=True, headers=headers) as client:
        for accepted_species in targets["accepted_species"]:
            axis = target_axis[accepted_species]
            accepted_hit = False
            for name_record in names_by_species.get(accepted_species, []):
                if accepted_hit and name_record["name_role"] != "accepted":
                    break
                search_name = name_record["search_name"]
                query = f'"{search_name}" {AXIS_QUERY_TERMS[axis]}'
                try:
                    response = client.get(
                        DUCKDUCKGO_HTML,
                        params={"q": query},
                    )
                    response.raise_for_status()
                    urls = parse_duckduckgo_urls(
                        response.text,
                        max_results=max_results_per_query,
                    )
                except Exception as exc:  # noqa: BLE001
                    errors.append(
                        {
                            "accepted_species": accepted_species,
                            "target_axis": axis,
                            "search_name": search_name,
                            "source": "duckduckgo_discovery",
                            "source_url": DUCKDUCKGO_HTML,
                            "error_class": type(exc).__name__,
                            "message": _text(exc),
                        }
                    )
                    urls = []
                for url in urls:
                    source_type, quality = _host_class(url)
                    host = (urlsplit(url).hostname or "").casefold()
                    discovered = {
                        "accepted_species": accepted_species,
                        "target_axis": axis,
                        "search_name": search_name,
                        "name_role": name_record["name_role"],
                        "source_url": url,
                        "source_host": host,
                        "source_type": source_type,
                        "fetch_status": "",
                    }
                    if not source_type:
                        discovered["fetch_status"] = "blocked_domain"
                        discovery_rows.append(discovered)
                        continue
                    allowed, robots_status = _robots_allowed(client, url, robots_cache)
                    if not allowed:
                        discovered["fetch_status"] = robots_status
                        discovery_rows.append(discovered)
                        continue
                    fetched = fetch_document(client, url)
                    discovered["fetch_status"] = fetched.status
                    discovery_rows.append(discovered)
                    if fetched.status != "fetched":
                        errors.append(
                            {
                                "accepted_species": accepted_species,
                                "target_axis": axis,
                                "search_name": search_name,
                                "source": "publisher_fetch",
                                "source_url": url,
                                "error_class": fetched.status,
                                "message": fetched.error,
                            }
                        )
                        continue
                    passages = extract_direct_passages(
                        accepted_species=accepted_species,
                        search_name=search_name,
                        target_axis=axis,
                        title=fetched.title,
                        text=fetched.text,
                    )
                    if not passages:
                        continue
                    if quality == "medium":
                        accepted_hit = True
                    source_rows.append(
                        {
                            "accepted_species": accepted_species,
                            "source_text": " ".join(passages),
                            "source_url": fetched.url,
                            "source_citation": (
                                f"{fetched.title or search_name} | discovered_query_axis={axis} | "
                                f"publisher_content_sha256={fetched.content_sha256}"
                            ),
                            "source_type": source_type,
                            "evidence_scope": "species_direct",
                            "search_name": search_name,
                            "name_role": name_record["name_role"],
                            "name_resolution_lineage": name_record["name_resolution_lineage"],
                            "target_axis": axis,
                            "source_host": host,
                            "retrieved_at_epoch": str(int(time.time())),
                            "raw_content_sha256": fetched.content_sha256,
                            "rights_status": "unknown_minimal_excerpt_only",
                            "robots_status": robots_status,
                        }
                    )
                if query_pause_seconds:
                    time.sleep(query_pause_seconds)
    return (
        pd.DataFrame(source_rows, columns=FETCHED_SOURCE_COLUMNS),
        pd.DataFrame(discovery_rows, columns=DISCOVERY_COLUMNS),
        pd.DataFrame(errors, columns=ERROR_COLUMNS),
    )


def _measurement_evidence(sources: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for record in sources.fillna("").to_dict("records"):
        if record.get("target_axis") != "floral_structural_complexity":
            continue
        if record.get("source_type") == UNTRUSTED_SOURCE_TYPE:
            continue
        for match in MEASUREMENT_RE.finditer(record.get("source_text", "")):
            low = float(match.group("low"))
            high = float(match.group("high") or low)
            midpoint_mm = ((low + high) / 2.0) * (
                10 if match.group("unit").casefold() == "cm" else 1
            )
            value = next(
                label
                for maximum, label in (
                    (2.0, "very_small"),
                    (10.0, "small"),
                    (30.0, "medium"),
                    (100.0, "large"),
                    (float("inf"), "very_large"),
                )
                if midpoint_mm <= maximum
            )
            rows.append(
                {
                    "species": record["accepted_species"],
                    "field": "flower_size_class",
                    "value": value,
                    "matched_term": match.group(0),
                    "source_type": record["source_type"],
                    "source_url": record["source_url"],
                    "source_citation": record["source_citation"],
                    "source_excerpt": match.group(0),
                    "evidence_type": "flora",
                    "confidence": "medium",
                    "source_tier": "A_flora",
                    "rule_id": "explicit_flower_measurement",
                }
            )
    return pd.DataFrame(rows, columns=V1_EVIDENCE_COLUMNS)


def _direct_evidence(
    sources: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    extraction_sources = sources[SOURCE_COLUMNS].copy()
    v1 = extract_evidence_from_sources(extraction_sources)
    measurements = _measurement_evidence(sources)
    v1 = pd.concat([v1, measurements], ignore_index=True).drop_duplicates()
    meta = sources.drop_duplicates(["accepted_species", "source_url"]).set_index(
        ["accepted_species", "source_url"]
    )
    rows: list[dict[str, str]] = []
    for record in v1.fillna("").to_dict("records"):
        if record["confidence"] != "medium":
            continue
        mapped = V2_FIELD_VALUE_MAP.get((record["field"], record["value"]))
        if mapped is None and record["field"] == "flower_size_class":
            mapped = (
                "M0_floral_phenotype",
                "flower_size_class",
                record["value"],
            )
        if mapped is None:
            continue
        _, trait_name, value = mapped
        axis = TRAIT_TO_AXIS.get(trait_name, "")
        if not axis:
            continue
        source_meta = meta.loc[(record["species"], record["source_url"])]
        if isinstance(source_meta, pd.DataFrame):
            source_meta = source_meta.iloc[0]
        lineage = _doi_lineage(
            record["source_citation"],
            record["source_url"],
            record["source_excerpt"],
        ) or _url_lineage(
            record["source_url"],
            source_meta["raw_content_sha256"],
        )
        role = source_meta["name_role"]
        rows.append(
            {
                "accepted_species": record["species"],
                "axis": axis,
                "trait_name": trait_name,
                "normalized_value": value,
                "evidence_quality": "medium",
                "evidence_scope": "species_direct" if role == "accepted" else "synonym_direct",
                "source_group": "broad_general_web",
                "source_provider": source_meta["source_host"] or record["source_type"],
                "source_url": record["source_url"],
                "source_citation": record["source_citation"],
                "source_excerpt": record["source_excerpt"],
                "source_record_id": hashlib.sha256(
                    (
                        record["species"]
                        + "|"
                        + trait_name
                        + "|"
                        + record["source_url"]
                        + "|"
                        + record["source_excerpt"]
                    ).encode("utf-8")
                ).hexdigest(),
                "source_lineage": lineage,
                "name_resolution_lineage": source_meta["name_resolution_lineage"],
                "name_match_method": (
                    "accepted_name_exact" if role == "accepted" else "exact_synonym_to_accepted"
                ),
                "inference_rule": "",
            }
        )
    direct = pd.DataFrame(rows, columns=DIRECT_EVIDENCE_COLUMNS)
    if not direct.empty:
        direct = direct.drop_duplicates(
            ["accepted_species", "axis", "trait_name", "normalized_value", "source_lineage"]
        )
    return v1, direct


def run_pilot(
    *,
    coverage_csv: Path,
    master_csv: Path,
    integrated_evidence_csv: Path,
    output_dir: Path,
    limit: int = 100,
    resolve_synonyms: bool = True,
    max_aliases_per_species: int = 3,
    max_aliases_to_search: int = 2,
    wfo_pause_seconds: float = 0.25,
    query_pause_seconds: float = 0.25,
    max_results_per_query: int = 5,
    max_gbif_usages_per_species: int = 25,
    workers: int = 12,
    include_structured_domains: bool = True,
    include_pladias: bool = False,
    include_gobotany: bool = True,
    include_gbif: bool = True,
    include_general_web: bool = True,
    selection_mode: str = "pilot",
    offset: int = 0,
) -> dict[str, Any]:
    coverage = pd.read_csv(coverage_csv, dtype=str).fillna("")
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    integrated = pd.read_csv(integrated_evidence_csv, dtype=str).fillna("")
    if selection_mode == "pilot":
        targets = select_pilot_targets(coverage, master, integrated, limit=limit)
    elif selection_mode == "all_unresolved":
        targets = select_all_unresolved_targets(
            coverage,
            master,
            offset=offset,
            limit=limit,
        )
    else:
        raise ValueError("selection_mode must be pilot or all_unresolved")
    search_names, wfo_errors = _search_names(
        targets,
        resolve_synonyms=resolve_synonyms and include_general_web,
        max_aliases_per_species=max_aliases_per_species,
        wfo_pause_seconds=wfo_pause_seconds,
    )
    if include_gbif:
        gbif_sources, gbif_usages, gbif_errors = _gbif_description_sources(
            targets,
            max_usages_per_species=max_gbif_usages_per_species,
            workers=workers,
        )
    else:
        gbif_sources = pd.DataFrame(columns=FETCHED_SOURCE_COLUMNS)
        gbif_usages = pd.DataFrame(columns=GBIF_USAGE_COLUMNS)
        gbif_errors = pd.DataFrame(columns=ERROR_COLUMNS)
    if include_general_web:
        web_sources, discoveries, web_errors = _general_web_sources(
            targets,
            search_names,
            max_results_per_query=max_results_per_query,
            max_aliases_to_search=max_aliases_to_search,
            query_pause_seconds=query_pause_seconds,
        )
    else:
        web_sources = pd.DataFrame(columns=FETCHED_SOURCE_COLUMNS)
        discoveries = pd.DataFrame(columns=DISCOVERY_COLUMNS)
        web_errors = pd.DataFrame(columns=ERROR_COLUMNS)
    if include_structured_domains:
        structured_direct, structured_errors, structured_matches = _structured_domain_evidence(
            targets,
            workers=workers,
            include_pladias=include_pladias,
            include_gobotany=include_gobotany,
        )
    else:
        structured_direct = pd.DataFrame(columns=DIRECT_EVIDENCE_COLUMNS)
        structured_errors = pd.DataFrame(columns=ERROR_COLUMNS)
        structured_matches = pd.DataFrame(columns=STRUCTURED_MATCH_COLUMNS)
    sources = pd.concat([gbif_sources, web_sources], ignore_index=True).drop_duplicates(
        ["accepted_species", "source_url", "raw_content_sha256"]
    )
    v1, extracted_direct = _direct_evidence(sources)
    direct = pd.concat(
        [structured_direct, extracted_direct],
        ignore_index=True,
    ).drop_duplicates()
    errors = pd.concat(
        [gbif_errors, web_errors, structured_errors],
        ignore_index=True,
    )
    if not wfo_errors.empty:
        converted = pd.DataFrame(
            {
                "accepted_species": wfo_errors["accepted_species"],
                "target_axis": "",
                "search_name": wfo_errors["accepted_species"],
                "source": wfo_errors["source"],
                "source_url": "",
                "error_class": "wfo_resolution",
                "message": wfo_errors["error"],
            }
        )
        errors = pd.concat([errors, converted], ignore_index=True)

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "targets": output_dir / "pilot_targets.csv",
        "search_names": output_dir / "exact_search_names.csv.gz",
        "gbif_usages": output_dir / "gbif_exact_usages.csv.gz",
        "discoveries": output_dir / "publisher_url_discovery.csv.gz",
        "sources": output_dir / "fetched_source_passages.csv.gz",
        "structured_matches": output_dir / "structured_domain_matches.csv.gz",
        "structured_evidence": output_dir / "structured_domain_medium_evidence.csv.gz",
        "v1_evidence": output_dir / "broad_web_v1_evidence.csv.gz",
        "direct_evidence": output_dir / "broad_web_medium_evidence.csv.gz",
        "errors": output_dir / "lookup_errors.csv.gz",
    }
    targets.to_csv(paths["targets"], index=False)
    search_names.to_csv(paths["search_names"], index=False)
    gbif_usages.to_csv(paths["gbif_usages"], index=False)
    discoveries.to_csv(paths["discoveries"], index=False)
    sources.to_csv(paths["sources"], index=False)
    structured_matches.to_csv(paths["structured_matches"], index=False)
    structured_direct.to_csv(paths["structured_evidence"], index=False)
    v1.to_csv(paths["v1_evidence"], index=False)
    direct.to_csv(paths["direct_evidence"], index=False)
    errors.to_csv(paths["errors"], index=False)

    accepted_keys = set(zip(direct["accepted_species"], direct["axis"])) if len(direct) else set()
    baseline = coverage.set_index(["accepted_species", "axis"])["after_quality"]
    gained_keys = {
        key for key in accepted_keys if key in baseline.index and not _text(baseline.loc[key])
    }
    upgraded_keys = {
        key
        for key in accepted_keys
        if key in baseline.index and _text(baseline.loc[key]).casefold() == "low"
    }
    by_axis = {
        axis: len({species for species, candidate_axis in gained_keys if candidate_axis == axis})
        for axis in AXES
    }
    target_axis_counts = {
        axis: int(
            targets["target_axis"]
            .map(lambda value, target_axis=axis: target_axis in _text(value).split("|"))
            .sum()
        )
        for axis in AXES
    }
    summary = {
        "contract_version": CONTRACT_VERSION,
        "policy": {
            "search_results_are_discovery_only": True,
            "publisher_page_or_pdf_required": True,
            "exact_accepted_or_wfo_synonym_required": True,
            "family_inference_used": False,
            "global_fallback_used": False,
        },
        "inputs": {
            "coverage_csv": str(coverage_csv),
            "coverage_sha256": _file_sha256(coverage_csv),
            "master_csv": str(master_csv),
            "master_sha256": _file_sha256(master_csv),
            "integrated_evidence_csv": str(integrated_evidence_csv),
            "integrated_evidence_sha256": _file_sha256(integrated_evidence_csv),
        },
        "targets": {
            "selection_mode": selection_mode,
            "offset": offset,
            "limit": limit,
            "species": int(targets["accepted_species"].nunique()),
            "by_axis": target_axis_counts,
            "genera": int(targets["genus"].nunique()),
        },
        "retrieval": {
            "structured_domains_enabled": include_structured_domains,
            "pladias_enabled": include_pladias,
            "gobotany_enabled": include_gobotany,
            "gbif_enabled": include_gbif,
            "general_web_enabled": include_general_web,
            "search_names": len(search_names),
            "wfo_synonyms": int(search_names["name_role"].eq("synonym").sum()),
            "gbif_exact_usages": len(gbif_usages),
            "gbif_description_hits": int(gbif_usages["status"].eq("description_hit").sum()),
            "publisher_urls_discovered": len(discoveries),
            "publisher_documents_fetched": int(discoveries["fetch_status"].eq("fetched").sum()),
            "species_direct_source_passages": len(sources),
            "structured_domain_matches": len(structured_matches),
            "structured_domain_pages_fetched": int(
                structured_matches["fetch_status"].eq("fetched").sum()
            ),
            "structured_domain_evidence_rows": len(structured_direct),
            "errors": len(errors),
        },
        "accepted": {
            "quality": "medium",
            "evidence_rows": len(direct),
            "pure_gain_species_axis": len(gained_keys),
            "validated_low_to_medium_upgrades": len(upgraded_keys),
            "by_axis": by_axis,
        },
        "files": {
            name: {
                "path": str(path),
                "sha256": _file_sha256(path),
            }
            for name, path in paths.items()
        },
    }
    summary_path = output_dir / "broad_web_medium_pilot_summary.json"
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


@app.command("select")
def select_command(
    coverage_csv: Annotated[
        Path,
        typer.Option(exists=True, dir_okay=False),
    ],
    master_csv: Annotated[
        Path,
        typer.Option(exists=True, dir_okay=False),
    ],
    integrated_evidence_csv: Annotated[
        Path,
        typer.Option(exists=True, dir_okay=False),
    ],
    output_csv: Annotated[Path, typer.Option(dir_okay=False)],
    limit: Annotated[int, typer.Option(min=3, max=1000)] = 100,
) -> None:
    targets = select_pilot_targets(
        pd.read_csv(coverage_csv, dtype=str).fillna(""),
        pd.read_csv(master_csv, dtype=str).fillna(""),
        pd.read_csv(integrated_evidence_csv, dtype=str).fillna(""),
        limit=limit,
    )
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    targets.to_csv(output_csv, index=False)
    typer.echo(json.dumps({"targets": len(targets)}, ensure_ascii=False))


@app.command("run")
def run_command(
    coverage_csv: Annotated[
        Path,
        typer.Option(exists=True, dir_okay=False),
    ],
    master_csv: Annotated[
        Path,
        typer.Option(exists=True, dir_okay=False),
    ],
    integrated_evidence_csv: Annotated[
        Path,
        typer.Option(exists=True, dir_okay=False),
    ],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
    limit: Annotated[int, typer.Option(min=0, max=110_000)] = 100,
    resolve_synonyms: Annotated[
        bool,
        typer.Option("--resolve-synonyms/--no-resolve-synonyms"),
    ] = True,
    max_aliases_per_species: Annotated[
        int,
        typer.Option(min=0, max=20),
    ] = 3,
    max_aliases_to_search: Annotated[
        int,
        typer.Option(min=0, max=10),
    ] = 2,
    wfo_pause_seconds: Annotated[
        float,
        typer.Option(min=0.0, max=10),
    ] = 0.25,
    query_pause_seconds: Annotated[
        float,
        typer.Option(min=0.0, max=10),
    ] = 0.25,
    max_results_per_query: Annotated[int, typer.Option(min=1, max=10)] = 5,
    max_gbif_usages_per_species: Annotated[
        int,
        typer.Option(min=1, max=100),
    ] = 25,
    workers: Annotated[int, typer.Option(min=1, max=32)] = 12,
    include_structured_domains: Annotated[
        bool,
        typer.Option("--structured-domains/--no-structured-domains"),
    ] = True,
    include_pladias: Annotated[
        bool,
        typer.Option("--pladias/--no-pladias"),
    ] = False,
    include_gobotany: Annotated[
        bool,
        typer.Option("--gobotany/--no-gobotany"),
    ] = True,
    include_gbif: Annotated[
        bool,
        typer.Option("--gbif/--no-gbif"),
    ] = True,
    include_general_web: Annotated[
        bool,
        typer.Option("--general-web/--no-general-web"),
    ] = True,
    selection_mode: Annotated[str, typer.Option("--selection-mode")] = "pilot",
    offset: Annotated[int, typer.Option(min=0, max=110_000)] = 0,
) -> None:
    summary = run_pilot(
        coverage_csv=coverage_csv,
        master_csv=master_csv,
        integrated_evidence_csv=integrated_evidence_csv,
        output_dir=output_dir,
        limit=limit,
        resolve_synonyms=resolve_synonyms,
        max_aliases_per_species=max_aliases_per_species,
        max_aliases_to_search=max_aliases_to_search,
        wfo_pause_seconds=wfo_pause_seconds,
        query_pause_seconds=query_pause_seconds,
        max_results_per_query=max_results_per_query,
        max_gbif_usages_per_species=max_gbif_usages_per_species,
        workers=workers,
        include_structured_domains=include_structured_domains,
        include_pladias=include_pladias,
        include_gobotany=include_gobotany,
        include_gbif=include_gbif,
        include_general_web=include_general_web,
        selection_mode=selection_mode,
        offset=offset,
    )
    typer.echo(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    app()
