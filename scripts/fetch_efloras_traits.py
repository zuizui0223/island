#!/usr/bin/env python3
"""Acquire eFloras species treatments and reviewable trait candidates.

The crawler follows eFloras' taxonomic hierarchy, inventories names before
matching them to the Island v2 master, and validates every selected leaf against
the species treatment itself.  Rule-extracted values remain explicitly
unreviewed and always retain the source URL and verbatim evidence excerpt.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import random
import re
import sys
import time
import unicodedata
from collections import Counter, defaultdict, deque
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import parse_qs, urlencode, urljoin, urlparse, urlunparse

import requests
from bs4 import BeautifulSoup

BASE_URL = "http://www.efloras.org/"
HOME_URL = urljoin(BASE_URL, "index.aspx")
USER_AGENT = "island-v2-efloras-harvester/2.0 (+https://github.com/zuizui0223/island)"

RANK_MARKERS = {"subsp", "subsp.", "ssp", "ssp.", "var", "var.", "f", "f.", "forma"}
BOTANICAL_TOKEN = re.compile(r"^[A-Za-z][A-Za-z-]*$")
FLORAL_ORGAN = re.compile(
    r"\b(flowers?|corollas?|petals?|perianths?|tepals?|caly(?:x|ces)|sepals?|"
    r"lips?|labella|labellums?)\b",
    re.I,
)
DESCRIPTION_SIGNAL = re.compile(
    r"\b(herbs?|shrubs?|trees?|plants?|climbers?|stems?|roots?|rhizomes?|caudex|"
    r"leaves|leaf|flowers?|inflorescences?|corollas?|petals?|perianths?|tepals?|"
    r"caly(?:x|ces)|sepals?|stamens?|pistils?|ovaries|capsules?)\b",
    re.I,
)
NONFLORAL_COLOR_CONTEXT = re.compile(
    r"\b(leaves|leaflets?|fruits?|seeds?|stems?|branches?|ovaries|ovary|capsules?)\b",
    re.I,
)
INDUMENTUM_AFTER_COLOR = re.compile(
    r"^\s+(?:(?:densely|sparsely|appressed|spreading|minute|short|long)\s+)*"
    r"(?:hairs?|scales?|trichomes?|pubescence|puberulent|pubescent|tomentose|villous)\b",
    re.I,
)
FRAGMENT_SPLIT = re.compile(r"(?<=[;!?])\s+|(?<=[a-z\]\)])\.\s+(?=[A-Z])|\n+")

COLOR_TERMS = {
    "white": "white",
    "whitish": "white",
    "cream": "cream",
    "creamy": "cream",
    "yellow": "yellow",
    "yellowish": "yellow",
    "orange": "orange",
    "red": "red",
    "reddish": "red",
    "pink": "pink",
    "pinkish": "pink",
    "purple": "purple",
    "purplish": "purple",
    "violet": "violet",
    "blue": "blue",
    "bluish": "blue",
    "green": "green",
    "greenish": "green",
    "brown": "brown",
    "brownish": "brown",
    "black": "black",
    "blackish": "black",
}
FORM_TERMS = {
    "campanulate": "campanulate",
    "bell-shaped": "campanulate",
    "tubular": "tubular",
    "funnelform": "funnelform",
    "funnel-shaped": "funnelform",
    "rotate": "rotate",
    "salverform": "salverform",
    "urceolate": "urceolate",
    "bilabiate": "bilabiate",
    "two-lipped": "bilabiate",
    "2-lipped": "bilabiate",
    "papilionaceous": "papilionaceous",
    "spurred": "spurred",
}
SYMMETRY_TERMS = {
    "actinomorphic": "radial",
    "radially symmetric": "radial",
    "radially symmetrical": "radial",
    "zygomorphic": "bilateral",
    "bilaterally symmetric": "bilateral",
    "bilaterally symmetrical": "bilateral",
}
SEX_TERMS = (
    "bisexual",
    "hermaphroditic",
    "unisexual",
    "monoecious",
    "dioecious",
    "andromonoecious",
    "gynomonoecious",
)

DISCOVERED_FIELDS = [
    "flora_id",
    "taxon_id",
    "original_taxon_name",
    "listing_rank",
    "family",
    "genus",
    "parent_taxon_id",
    "treatment_url",
    "browse_url",
    "discovery_url",
]
MATCH_FIELDS = [
    "flora_id",
    "taxon_id",
    "original_taxon_name",
    "listing_rank",
    "accepted_species",
    "match_status",
    "match_candidates",
    "selection_status",
]
TREATMENT_FIELDS = [
    "accepted_species",
    "original_taxon_name",
    "taxon_heading",
    "taxonomic_rank",
    "match_status",
    "flora_id",
    "flora_name",
    "taxon_id",
    "source_url",
    "description",
    "description_chars",
    "trait_candidate_rows",
    "treatment_confirmed",
]
CANDIDATE_FIELDS = [
    "accepted_species",
    "original_taxon_name",
    "trait_name",
    "candidate_value",
    "flora_id",
    "flora_name",
    "taxon_id",
    "source_url",
    "source_excerpt",
    "source_record_id",
    "evidence_scope",
    "evidence_status",
    "extraction_method",
    "name_match_method",
]
MANIFEST_FIELDS = [
    "flora_id",
    "flora_name",
    "landing_url",
    "status",
    "has_volume_structure",
    "has_browse_structure",
]


class AcquisitionError(RuntimeError):
    """Raised after diagnostic artifacts have been written."""


def clean(value: Any) -> str:
    text = html.unescape(str(value or ""))
    text = unicodedata.normalize("NFKC", text).strip().strip('"')
    return re.sub(r"\s+", " ", text)


def exact_key(value: Any) -> str:
    return clean(value).casefold()


def binomial_key(value: Any) -> str | None:
    """Return a conservative genus+epithet key, never collapsing infraspecies."""
    text = clean(value).replace("×", " x ")
    tokens = [token.strip(".,;:()[]{}") for token in text.split()]
    lowered = [token.casefold() for token in tokens]
    if "x" in lowered or any(token in RANK_MARKERS for token in lowered):
        return None
    botanical = [token for token in tokens if BOTANICAL_TOKEN.fullmatch(token)]
    if len(botanical) < 2:
        return None
    genus, epithet = botanical[:2]
    if epithet.casefold() in {"sp", "spp", "cf", "aff"}:
        return None
    return f"{genus.casefold()} {epithet.casefold()}"


def taxonomic_rank(name: str) -> str:
    text = clean(name).replace("×", " x ")
    tokens = [token.strip(".,;:()[]{}") for token in text.split()]
    lowered = [token.casefold() for token in tokens]
    if "x" in lowered:
        return "hybrid"
    if any(token in RANK_MARKERS for token in lowered):
        return "infraspecies"
    botanical = [token for token in tokens if BOTANICAL_TOKEN.fullmatch(token)]
    if len(botanical) >= 2 and botanical[0][:1].isupper() and botanical[1][:1].islower():
        return "species"
    if len(botanical) == 1 and botanical[0][:1].isupper():
        return "family" if botanical[0].casefold().endswith("aceae") else "genus"
    return "unknown"


def query_value(url: str, key: str) -> str:
    return parse_qs(urlparse(url).query).get(key, [""])[0]


def with_query(url: str, **updates: Any) -> str:
    parsed = urlparse(url)
    query = parse_qs(parsed.query, keep_blank_values=True)
    for key, value in updates.items():
        if value is None:
            query.pop(key, None)
        else:
            query[key] = [str(value)]
    encoded = urlencode(query, doseq=True)
    return urlunparse(parsed._replace(query=encoded))


def canonical_browse_url(url: str) -> str:
    return with_query(url, page=None)


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    rows: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            rows.append(json.loads(line))
    return rows


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    content = json.dumps(payload, ensure_ascii=False, indent=2) + "\n"
    temporary.write_text(content, encoding="utf-8")
    for attempt in range(6):
        try:
            temporary.replace(path)
            return
        except PermissionError:
            if attempt < 5:
                time.sleep(0.1 * (attempt + 1))
    # OneDrive/antivirus can hold the temporary file open on Windows.  A direct
    # final write is preferable to losing a resumable checkpoint altogether.
    path.write_text(content, encoding="utf-8")
    temporary.unlink(missing_ok=True)


def write_errors(path: Path, errors: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for error in errors:
            handle.write(json.dumps(error, ensure_ascii=False) + "\n")


class Client:
    def __init__(
        self,
        delay_min: float = 1.5,
        delay_max: float = 3.5,
        timeout: float = 120,
        attempts: int = 6,
    ):
        self.delay_min = delay_min
        self.delay_max = delay_max
        self.timeout = timeout
        self.attempts = attempts
        self.session = requests.Session()
        self.session.headers["User-Agent"] = USER_AGENT
        self.http_errors = 0
        self.requests_made = 0
        self.errors: list[dict[str, Any]] = []

    def get(self, url: str, stage: str = "request") -> requests.Response:
        last: Exception | None = None
        for attempt in range(1, self.attempts + 1):
            self.requests_made += 1
            try:
                response = self.session.get(url, timeout=self.timeout)
                if response.status_code == 200:
                    if self.delay_max > 0:
                        time.sleep(random.uniform(self.delay_min, self.delay_max))
                    return response
                self.http_errors += 1
                last = RuntimeError(f"HTTP {response.status_code}")
                retry_after = response.headers.get("Retry-After", "")
                wait = float(retry_after) if retry_after.isdigit() else min(60.0, 2.0 ** (attempt - 1))
            except Exception as exc:  # requests exposes several transport exceptions
                self.http_errors += 1
                last = exc
                wait = min(60.0, 2.0 ** (attempt - 1))
            if attempt < self.attempts:
                time.sleep(wait)
        error = {"stage": stage, "url": url, "error": str(last), "attempts": self.attempts}
        self.errors.append(error)
        raise RuntimeError(f"GET failed: {url}: {last}")


class MasterMatcher:
    def __init__(self, accepted_names: Iterable[str], synonyms: Iterable[tuple[str, str]] = ()):
        names = sorted({clean(name) for name in accepted_names if clean(name)})
        self.accepted = set(names)
        self.exact: dict[str, list[str]] = defaultdict(list)
        self.binomial: dict[str, list[str]] = defaultdict(list)
        self.synonyms: dict[str, list[str]] = defaultdict(list)
        for name in names:
            self.exact[exact_key(name)].append(name)
            key = binomial_key(name)
            if key:
                self.binomial[key].append(name)
        for alias, accepted in synonyms:
            alias, accepted = clean(alias), clean(accepted)
            if alias and accepted in self.accepted:
                self.synonyms[exact_key(alias)].append(accepted)

    @classmethod
    def from_csv(
        cls,
        master_csv: Path | None,
        name_column: str,
        synonyms_csv: Path | None = None,
    ) -> MasterMatcher | None:
        if master_csv is None:
            return None
        rows = read_csv(master_csv)
        if not rows or name_column not in rows[0]:
            raise ValueError(f"Missing master column {name_column!r} in {master_csv}")
        synonyms: list[tuple[str, str]] = []
        if synonyms_csv:
            synonym_rows = read_csv(synonyms_csv)
            if synonym_rows:
                alias_column = next(
                    (name for name in ("original_taxon_name", "synonym", "scientific_name", "name") if name in synonym_rows[0]),
                    None,
                )
                accepted_column = next(
                    (name for name in ("accepted_species", "accepted_name") if name in synonym_rows[0]),
                    None,
                )
                if not alias_column or not accepted_column:
                    raise ValueError("Synonym CSV needs an alias column and accepted_species/accepted_name")
                synonyms = [(row[alias_column], row[accepted_column]) for row in synonym_rows]
        return cls((row[name_column] for row in rows), synonyms)

    @staticmethod
    def _result(candidates: list[str], status: str) -> dict[str, str]:
        unique = sorted(set(candidates))
        if len(unique) == 1:
            return {"accepted_species": unique[0], "match_status": status, "match_candidates": unique[0]}
        if len(unique) > 1:
            return {
                "accepted_species": "",
                "match_status": "ambiguous",
                "match_candidates": "|".join(unique),
            }
        return {"accepted_species": "", "match_status": "unmatched", "match_candidates": ""}

    def match(self, source_name: str) -> dict[str, str]:
        exact = self.exact.get(exact_key(source_name), [])
        if exact:
            return self._result(exact, "exact_accepted_name_match")
        key = binomial_key(source_name)
        normalized = self.binomial.get(key, []) if key else []
        if normalized:
            return self._result(normalized, "normalized_binomial_match")
        synonym = self.synonyms.get(exact_key(source_name), [])
        if synonym:
            return self._result(synonym, "synonym_match")
        return self._result([], "unmatched")


@dataclass
class BrowseTask:
    url: str
    parent_rank: str
    parent_taxon_id: str = ""
    family: str = ""
    genus: str = ""
    depth: int = 0


def infer_child_rank(parent_rank: str, label: str) -> str:
    observed = taxonomic_rank(label)
    if parent_rank == "flora":
        return "family"
    if parent_rank == "family":
        return "genus"
    if parent_rank == "genus":
        return observed if observed in {"species", "infraspecies", "hybrid"} else "unknown"
    if parent_rank == "species":
        return observed if observed != "unknown" else "infraspecies"
    return observed


def parse_listing_page(
    page_html: str,
    page_url: str,
    flora_id: int,
    task: BrowseTask,
) -> list[dict[str, str]]:
    soup = BeautifulSoup(page_html, "html.parser")
    rows = soup.select("tr.underline") or soup.find_all("tr")
    parsed: list[dict[str, str]] = []
    for row in rows:
        taxon_link = None
        for link in row.find_all("a", href=True):
            absolute = urljoin(page_url, link["href"])
            if "florataxon.aspx" not in absolute:
                continue
            linked_flora = query_value(absolute, "flora_id")
            taxon_id = query_value(absolute, "taxon_id")
            if taxon_id.isdigit() and (not linked_flora or linked_flora == str(flora_id)):
                taxon_link = (link, absolute, taxon_id)
                break
        if not taxon_link:
            continue
        link, treatment_url, taxon_id = taxon_link
        label = clean(link.get_text(" ", strip=True))
        if not label:
            continue
        browse_url = ""
        for child in row.find_all("a", href=True):
            absolute = urljoin(page_url, child["href"])
            if "browse.aspx" in absolute and query_value(absolute, "start_taxon_id") == taxon_id:
                browse_url = absolute
                break
        rank = infer_child_rank(task.parent_rank, label)
        family = label if rank == "family" else task.family
        genus = label if rank == "genus" else task.genus
        parsed.append(
            {
                "flora_id": str(flora_id),
                "taxon_id": taxon_id,
                "original_taxon_name": label,
                "listing_rank": rank,
                "family": family,
                "genus": genus,
                "parent_taxon_id": task.parent_taxon_id,
                "treatment_url": treatment_url,
                "browse_url": browse_url,
                "discovery_url": page_url,
            }
        )
    return parsed


def listing_entries(
    client: Client,
    task: BrowseTask,
    flora_id: int,
    max_pages: int,
) -> tuple[list[dict[str, str]], int]:
    found: dict[str, dict[str, str]] = {}
    signatures: set[tuple[str, ...]] = set()
    pages_visited = 0
    for page_number in range(1, max_pages + 1):
        url = with_query(task.url, page=page_number)
        response = client.get(url, stage="browse")
        pages_visited += 1
        entries = parse_listing_page(response.text, response.url, flora_id, task)
        signature = tuple(entry["taxon_id"] for entry in entries)
        if not signature or signature in signatures:
            break
        signatures.add(signature)
        for entry in entries:
            found.setdefault(entry["taxon_id"], entry)
    return list(found.values()), pages_visited


def _matched_species_count(taxa: Iterable[dict[str, str]], matcher: MasterMatcher | None) -> int:
    species = [row for row in taxa if row["listing_rank"] == "species"]
    if matcher is None:
        return len(species)
    return sum(matcher.match(row["original_taxon_name"])["match_status"] not in {"unmatched", "ambiguous"} for row in species)


def crawl_taxa(
    *,
    client: Client,
    flora_id: int,
    matcher: MasterMatcher | None,
    max_pages: int,
    max_depth: int,
    match_target: int,
    checkpoint_path: Path,
    checkpoint_every: int,
    resume: bool,
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    queue: deque[BrowseTask]
    taxa: dict[str, dict[str, str]] = {}
    visited: set[str] = set()
    pages_visited = 0
    resumed = False
    if resume and checkpoint_path.exists():
        saved = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        if saved.get("complete"):
            return saved.get("taxa", []), saved.get("metrics", {})
        queue = deque(BrowseTask(**row) for row in saved.get("queue", []))
        taxa = {row["taxon_id"]: row for row in saved.get("taxa", [])}
        visited = set(saved.get("visited", []))
        pages_visited = int(saved.get("pages_visited", 0))
        resumed = True
    else:
        queue = deque([BrowseTask(url=f"{BASE_URL}browse.aspx?flora_id={flora_id}", parent_rank="flora")])

    tasks_processed = 0
    stop_reason = "queue_exhausted"
    while queue:
        # Depth-first traversal reaches a confirmed species quickly in bounded
        # smoke/pilot runs.  An unbounded production run still visits the same
        # complete hierarchy.
        task = queue.pop()
        canonical = canonical_browse_url(task.url)
        if canonical in visited or task.depth > max_depth:
            continue
        try:
            entries, page_count = listing_entries(client, task, flora_id, max_pages)
        except Exception:
            # Preserve the failed node for a later resumed run after the remote
            # service or network recovers; do not spin on it in this process.
            queue.appendleft(task)
            stop_reason = "browse_request_error"
            break
        visited.add(canonical)
        pages_visited += page_count
        for entry in entries:
            taxa.setdefault(entry["taxon_id"], entry)
            if entry["browse_url"] and task.depth < max_depth and entry["listing_rank"] not in {"infraspecies", "hybrid"}:
                queue.append(
                    BrowseTask(
                        url=entry["browse_url"],
                        parent_rank=entry["listing_rank"],
                        parent_taxon_id=entry["taxon_id"],
                        family=entry["family"],
                        genus=entry["genus"],
                        depth=task.depth + 1,
                    )
                )
        tasks_processed += 1
        if match_target > 0 and _matched_species_count(taxa.values(), matcher) >= match_target:
            stop_reason = "max_species_match_target"
            break
        if tasks_processed % max(1, checkpoint_every) == 0:
            write_json(
                checkpoint_path,
                {
                    "complete": False,
                    "queue": [asdict(item) for item in queue],
                    "visited": sorted(visited),
                    "taxa": list(taxa.values()),
                    "pages_visited": pages_visited,
                },
            )

    complete = not queue and stop_reason == "queue_exhausted"
    metrics = {
        "browse_pages_visited": pages_visited,
        "browse_nodes_visited": len(visited),
        "inventory_complete": complete,
        "inventory_stop_reason": stop_reason,
        "crawl_resumed": resumed,
    }
    ordered = sorted(taxa.values(), key=lambda row: (row["listing_rank"], row["original_taxon_name"], row["taxon_id"]))
    write_json(
        checkpoint_path,
        {
            "complete": complete,
            "queue": [asdict(item) for item in queue],
            "visited": sorted(visited),
            "taxa": ordered,
            "pages_visited": pages_visited,
            "metrics": metrics,
        },
    )
    return ordered, metrics


def parse_treatment(page_html: str, source_url: str, min_description_chars: int = 20) -> dict[str, Any]:
    soup = BeautifulSoup(page_html, "html.parser")
    node = soup.find(id="lblTaxonDesc")
    description = clean(node.get_text(" ", strip=True)) if node else ""
    heading_node = node.find("b") if node else None
    heading = clean(heading_node.get_text(" ", strip=True)) if heading_node else ""
    rank = taxonomic_rank(heading)
    narrative_chars = max(0, len(description) - len(heading))
    confirmed = bool(
        node
        and heading
        and rank in {"species", "infraspecies", "hybrid"}
        and narrative_chars >= min_description_chars
        and DESCRIPTION_SIGNAL.search(description)
    )
    return {
        "taxon_heading": heading,
        "taxonomic_rank": rank,
        "description": description,
        "description_chars": len(description),
        "narrative_chars": narrative_chars,
        "source_url": source_url,
        "treatment_confirmed": confirmed,
    }


def _term_hits(text: str, mapping: dict[str, str]) -> list[str]:
    low = text.casefold()
    return sorted(
        {
            normalized
            for phrase, normalized in mapping.items()
            if re.search(rf"(?<![A-Za-z]){re.escape(phrase)}(?![A-Za-z])", low)
        }
    )


def _color_hits(text: str) -> list[str]:
    """Extract floral-organ colours while rejecting nearby non-floral material."""
    low = text.casefold()
    organ_positions = [match.start() for match in FLORAL_ORGAN.finditer(text)]
    nonfloral_positions = [match.start() for match in NONFLORAL_COLOR_CONTEXT.finditer(text)]
    hits: set[str] = set()
    for phrase, normalized in COLOR_TERMS.items():
        for match in re.finditer(rf"(?<![A-Za-z]){re.escape(phrase)}(?![A-Za-z])", low):
            previous_organs = [position for position in organ_positions if position < match.start()]
            if not previous_organs:
                continue
            previous_nonfloral = [position for position in nonfloral_positions if position < match.start()]
            if previous_nonfloral and max(previous_nonfloral) > max(previous_organs):
                continue
            if INDUMENTUM_AFTER_COLOR.search(low[match.end() : match.end() + 50]):
                continue
            hits.add(normalized)
    return sorted(hits)


def trait_candidates(
    *,
    flora_id: int,
    flora_name: str,
    accepted_species: str,
    original_taxon_name: str,
    taxon_id: str,
    description: str,
    source_url: str,
    name_match_method: str,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for fragment in FRAGMENT_SPLIT.split(description):
        excerpt = clean(fragment)
        if not 8 <= len(excerpt) <= 2000:
            continue
        low = excerpt.casefold()
        hits: list[tuple[str, str]] = []
        if FLORAL_ORGAN.search(excerpt):
            colors = _color_hits(excerpt)
            if colors:
                hits.append(("flower_primary_color", "|".join(colors)))
            forms = _term_hits(excerpt, FORM_TERMS)
            if forms:
                hits.append(("floral_form", "|".join(forms)))
            symmetry = _term_hits(excerpt, SYMMETRY_TERMS)
            if symmetry:
                hits.append(("floral_symmetry", "|".join(symmetry)))
        sex = sorted({term for term in SEX_TERMS if re.search(rf"\b{re.escape(term)}\b", low)})
        if sex:
            hits.append(("sex_system", "|".join(sex)))
        if re.search(r"\bself[- ]incompat(?:ible|ibility)?\b", low):
            hits.append(("self_incompatibility", "self_incompatible"))
        elif re.search(r"\bself[- ]compat(?:ible|ibility)?\b", low):
            hits.append(("self_incompatibility", "self_compatible"))

        for trait_name, candidate_value in hits:
            key = (trait_name, candidate_value, excerpt)
            if key in seen:
                continue
            seen.add(key)
            rows.append(
                {
                    "accepted_species": accepted_species,
                    "original_taxon_name": original_taxon_name,
                    "trait_name": trait_name,
                    "candidate_value": candidate_value,
                    "flora_id": str(flora_id),
                    "flora_name": flora_name,
                    "taxon_id": taxon_id,
                    "source_url": source_url,
                    "source_excerpt": excerpt,
                    "source_record_id": f"efloras:{flora_id}:{taxon_id}",
                    "evidence_scope": "species_direct",
                    "evidence_status": "source_backed_rule_extracted_unreviewed",
                    "extraction_method": "rule_based_floral_context_v2",
                    "name_match_method": name_match_method,
                }
            )
    return rows


def flora_name_from_landing(page_html: str, fallback: str) -> str:
    soup = BeautifulSoup(page_html, "html.parser")
    title = clean(soup.title.get_text(" ", strip=True) if soup.title else "")
    title = re.sub(r"\s*@\s*efloras\.org\s*$", "", title, flags=re.I)
    return title or fallback


def discover_floras(
    *,
    client: Client,
    output_dir: Path,
    requested_ids: set[int] | None = None,
) -> list[dict[str, Any]]:
    home = client.get(HOME_URL, stage="flora_manifest")
    soup = BeautifulSoup(home.text, "html.parser")
    listed: dict[int, tuple[str, str]] = {}
    for link in soup.find_all("a", href=True):
        absolute = urljoin(home.url, link["href"])
        flora_id = query_value(absolute, "flora_id")
        if "flora_page.aspx" in absolute and flora_id.isdigit():
            listed.setdefault(int(flora_id), (clean(link.get_text(" ", strip=True)), absolute))
    ids = sorted(requested_ids if requested_ids is not None else listed)
    rows: list[dict[str, Any]] = []
    for flora_id in ids:
        listed_name, landing_url = listed.get(
            flora_id,
            (f"eFloras flora {flora_id}", f"{BASE_URL}flora_page.aspx?flora_id={flora_id}"),
        )
        status = "ok"
        has_volume = False
        has_browse = False
        flora_name = listed_name
        try:
            landing = client.get(landing_url, stage="flora_manifest_landing")
            landing_soup = BeautifulSoup(landing.text, "html.parser")
            flora_name = flora_name_from_landing(landing.text, listed_name)
            has_volume = any("volume_page.aspx" in urljoin(landing.url, link["href"]) for link in landing_soup.find_all("a", href=True))
            browse = client.get(f"{BASE_URL}browse.aspx?flora_id={flora_id}", stage="flora_manifest_browse")
            browse_soup = BeautifulSoup(browse.text, "html.parser")
            has_browse = bool(browse_soup.select("tr.underline") and browse_soup.find("a", href=re.compile(r"taxon_id=")))
            if flora_id not in listed:
                status = "not_listed" if not has_browse else "ok_unlisted"
            elif not has_browse:
                status = "unsupported_no_browse"
        except Exception as exc:
            status = "http_error"
            client.errors.append({"stage": "flora_manifest", "flora_id": flora_id, "url": landing_url, "error": str(exc)})
        rows.append(
            {
                "flora_id": flora_id,
                "flora_name": flora_name,
                "landing_url": landing_url,
                "status": status,
                "has_volume_structure": has_volume,
                "has_browse_structure": has_browse,
            }
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "flora_manifest.csv", rows, MANIFEST_FIELDS)
    write_json(output_dir / "flora_manifest.json", rows)
    write_errors(output_dir / "flora_manifest_errors.jsonl", client.errors)
    return rows


def _match_inventory(
    taxa: list[dict[str, str]],
    flora_id: int,
    matcher: MasterMatcher | None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    report: list[dict[str, str]] = []
    selected: list[dict[str, str]] = []
    for taxon in taxa:
        if taxon["listing_rank"] not in {"species", "infraspecies", "hybrid"}:
            continue
        match = (
            matcher.match(taxon["original_taxon_name"])
            if matcher
            else {
                "accepted_species": taxon["original_taxon_name"],
                "match_status": "no_master_filter",
                "match_candidates": taxon["original_taxon_name"],
            }
        )
        selection_status = "eligible"
        if taxon["listing_rank"] != "species":
            selection_status = f"excluded_{taxon['listing_rank']}"
        elif match["match_status"] in {"unmatched", "ambiguous"}:
            selection_status = match["match_status"]
        row = {
            "flora_id": str(flora_id),
            "taxon_id": taxon["taxon_id"],
            "original_taxon_name": taxon["original_taxon_name"],
            "listing_rank": taxon["listing_rank"],
            **match,
            "selection_status": selection_status,
        }
        report.append(row)
        if selection_status == "eligible":
            selected.append({**taxon, **match})

    priority = {"exact_accepted_name_match": 0, "normalized_binomial_match": 1, "synonym_match": 2, "no_master_filter": 3}
    selected.sort(key=lambda row: (row["accepted_species"], priority.get(row["match_status"], 9), row["taxon_id"]))
    unique: list[dict[str, str]] = []
    seen_species: set[str] = set()
    for row in selected:
        if row["accepted_species"] in seen_species:
            for audit in report:
                if audit["taxon_id"] == row["taxon_id"]:
                    audit["selection_status"] = "duplicate_flora_match"
                    break
            continue
        seen_species.add(row["accepted_species"])
        unique.append(row)
    return report, unique


def _persist_treatment_state(
    *,
    output_dir: Path,
    treatments: list[dict[str, Any]],
    candidates: list[dict[str, str]],
    errors: list[dict[str, Any]],
    processed_taxon_ids: set[str],
    request_total: int,
) -> None:
    treatments.sort(key=lambda row: (str(row["accepted_species"]), str(row["taxon_id"])))
    candidates.sort(key=lambda row: (row["accepted_species"], row["trait_name"], row["taxon_id"], row["source_excerpt"]))
    write_csv(output_dir / "species_treatments.csv", treatments, TREATMENT_FIELDS)
    write_csv(output_dir / "trait_candidates.csv", candidates, CANDIDATE_FIELDS)
    write_errors(output_dir / "errors.jsonl", errors)
    write_json(
        output_dir / "treatment_checkpoint.json",
        {
            "processed_taxon_ids": sorted(processed_taxon_ids),
            "species_treatment_requests": request_total,
            "species_treatments_saved": len(treatments),
            "trait_candidate_rows": len(candidates),
        },
    )


def run_flora(args: argparse.Namespace) -> dict[str, Any]:
    flora_id = int(args.flora_id)
    shard_count = int(args.shard_count)
    shard_index = int(args.shard_index)
    if shard_count < 1 or not 0 <= shard_index < shard_count:
        raise ValueError("shard-index must satisfy 0 <= shard-index < shard-count")

    output_dir = args.output_dir / f"flora-{flora_id}"
    if shard_count > 1:
        output_dir = output_dir / f"shard-{shard_index:04d}-of-{shard_count:04d}"
    output_dir.mkdir(parents=True, exist_ok=True)
    html_dir = output_dir / "html"
    if args.save_html:
        html_dir.mkdir(parents=True, exist_ok=True)

    client = Client(args.delay_min, args.delay_max, args.timeout, args.attempts)
    matcher = MasterMatcher.from_csv(args.master_csv, args.master_name_column, args.synonyms_csv)
    landing_url = f"{BASE_URL}flora_page.aspx?flora_id={flora_id}"
    try:
        landing = client.get(landing_url, stage="landing")
        flora_name = flora_name_from_landing(landing.text, f"eFloras flora {flora_id}")
    except Exception as exc:
        report = {"status": "failed_landing", "flora_id": flora_id, "flora_name": "", "error": str(exc), "HTTP_errors": client.http_errors}
        write_json(output_dir / "coverage_report.json", report)
        write_errors(output_dir / "errors.jsonl", client.errors)
        raise AcquisitionError(str(exc)) from exc

    match_target = args.max_species * shard_count if args.max_species > 0 else 0
    taxa, crawl_metrics = crawl_taxa(
        client=client,
        flora_id=flora_id,
        matcher=matcher,
        max_pages=args.max_pages,
        max_depth=args.max_depth,
        match_target=match_target,
        checkpoint_path=output_dir / "crawl_checkpoint.json",
        checkpoint_every=args.checkpoint_every,
        resume=args.resume,
    )
    write_csv(output_dir / "discovered_taxa.csv", taxa, DISCOVERED_FIELDS)
    match_report, eligible = _match_inventory(taxa, flora_id, matcher)
    sharded = [row for index, row in enumerate(eligible) if index % shard_count == shard_index]
    if args.max_species > 0:
        sharded = sharded[: args.max_species]
    selected_ids = {row["taxon_id"] for row in sharded}
    for row in match_report:
        if row["selection_status"] == "eligible":
            row["selection_status"] = "selected" if row["taxon_id"] in selected_ids else "not_in_current_shard_or_limit"
    write_csv(output_dir / "master_match_report.csv", match_report, MATCH_FIELDS)

    treatment_path = output_dir / "species_treatments.csv"
    candidate_path = output_dir / "trait_candidates.csv"
    treatments: list[dict[str, Any]] = read_csv(treatment_path) if args.resume else []
    candidates: list[dict[str, str]] = read_csv(candidate_path) if args.resume else []
    processed_taxon_ids: set[str] = set()
    request_total = 0
    checkpoint_path = output_dir / "treatment_checkpoint.json"
    if args.resume and checkpoint_path.exists():
        checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        processed_taxon_ids = set(checkpoint.get("processed_taxon_ids", []))
        request_total = int(checkpoint.get("species_treatment_requests", 0))
    errors = read_jsonl(output_dir / "errors.jsonl") if args.resume else []
    errors.extend(error for error in client.errors if error not in errors)

    for index, selected in enumerate(sharded, start=1):
        taxon_id = selected["taxon_id"]
        if taxon_id in processed_taxon_ids:
            continue
        request_total += 1
        try:
            response = client.get(selected["treatment_url"], stage="species_treatment")
            parsed = parse_treatment(response.text, response.url, args.min_description_chars)
        except Exception as exc:
            errors.append(
                {
                    "stage": "species_treatment",
                    "flora_id": flora_id,
                    "taxon_id": taxon_id,
                    "original_taxon_name": selected["original_taxon_name"],
                    "url": selected["treatment_url"],
                    "error": str(exc),
                }
            )
            continue

        if not parsed["treatment_confirmed"]:
            errors.append(
                {
                    "stage": "treatment_validation",
                    "flora_id": flora_id,
                    "taxon_id": taxon_id,
                    "original_taxon_name": selected["original_taxon_name"],
                    "taxon_heading": parsed["taxon_heading"],
                    "taxonomic_rank": parsed["taxonomic_rank"],
                    "description_chars": parsed["description_chars"],
                    "narrative_chars": parsed["narrative_chars"],
                    "error": "not_a_confirmed_species_treatment",
                }
            )
            processed_taxon_ids.add(taxon_id)
            continue

        heading_match = matcher.match(parsed["taxon_heading"]) if matcher else {
            "accepted_species": parsed["taxon_heading"],
            "match_status": "no_master_filter",
            "match_candidates": parsed["taxon_heading"],
        }
        if parsed["taxonomic_rank"] != "species" or heading_match["match_status"] in {"unmatched", "ambiguous"}:
            errors.append(
                {
                    "stage": "treatment_validation",
                    "flora_id": flora_id,
                    "taxon_id": taxon_id,
                    "original_taxon_name": selected["original_taxon_name"],
                    "taxon_heading": parsed["taxon_heading"],
                    "taxonomic_rank": parsed["taxonomic_rank"],
                    "error": "rank_or_heading_master_mismatch",
                }
            )
            processed_taxon_ids.add(taxon_id)
            continue
        accepted_species = heading_match["accepted_species"]
        extracted = trait_candidates(
            flora_id=flora_id,
            flora_name=flora_name,
            accepted_species=accepted_species,
            original_taxon_name=selected["original_taxon_name"],
            taxon_id=taxon_id,
            description=parsed["description"],
            source_url=parsed["source_url"],
            name_match_method=heading_match["match_status"],
        )
        candidates.extend(extracted)
        treatments.append(
            {
                "accepted_species": accepted_species,
                "original_taxon_name": selected["original_taxon_name"],
                "taxon_heading": parsed["taxon_heading"],
                "taxonomic_rank": parsed["taxonomic_rank"],
                "match_status": heading_match["match_status"],
                "flora_id": flora_id,
                "flora_name": flora_name,
                "taxon_id": taxon_id,
                "source_url": parsed["source_url"],
                "description": parsed["description"],
                "description_chars": parsed["description_chars"],
                "trait_candidate_rows": len(extracted),
                "treatment_confirmed": True,
            }
        )
        if args.save_html:
            (html_dir / f"{taxon_id}.html").write_text(response.text, encoding="utf-8")
        processed_taxon_ids.add(taxon_id)
        if index % max(1, args.checkpoint_every) == 0:
            _persist_treatment_state(
                output_dir=output_dir,
                treatments=treatments,
                candidates=candidates,
                errors=errors,
                processed_taxon_ids=processed_taxon_ids,
                request_total=request_total,
            )
            print(f"{index}/{len(sharded)} selected treatments; saved={len(treatments)} candidates={len(candidates)}", flush=True)

    # De-duplicate deterministically when resuming from a previous partial run.
    treatments = list({(row["flora_id"], row["taxon_id"]): row for row in treatments}.values())
    candidates = list(
        {
            (row["flora_id"], row["taxon_id"], row["trait_name"], row["candidate_value"], row["source_excerpt"]): row
            for row in candidates
        }.values()
    )
    errors.extend(error for error in client.errors if error not in errors)
    _persist_treatment_state(
        output_dir=output_dir,
        treatments=treatments,
        candidates=candidates,
        errors=errors,
        processed_taxon_ids=processed_taxon_ids,
        request_total=request_total,
    )

    status_counts = Counter(row["match_status"] for row in match_report if row["listing_rank"] == "species")
    species_names = {row["original_taxon_name"] for row in taxa if row["listing_rank"] == "species"}
    species_with_traits = {row["accepted_species"] for row in candidates}
    traits = sorted({row["trait_name"] for row in candidates})
    pending_selected = sorted(selected_ids.difference(processed_taxon_ids))
    if crawl_metrics.get("inventory_stop_reason") == "browse_request_error":
        status = "failed_browse_request"
    elif pending_selected:
        status = "partial_treatment_errors"
    else:
        status = "ok" if treatments else ("no_master_matches" if not sharded else "failed_no_treatments")
    report = {
        "status": status,
        "flora_id": flora_id,
        "flora_name": flora_name,
        "landing_url": landing_url,
        **crawl_metrics,
        "taxa_links_discovered": len(taxa),
        "species_names_discovered": len(species_names),
        "exact_master_matches": status_counts.get("exact_accepted_name_match", 0),
        "normalized_matches": status_counts.get("normalized_binomial_match", 0),
        "synonym_matches": status_counts.get("synonym_match", 0),
        "unmatched_names": status_counts.get("unmatched", 0),
        "ambiguous_names": status_counts.get("ambiguous", 0),
        "eligible_unique_master_species": len(eligible),
        "selected_for_shard": len(sharded),
        "shard_index": shard_index,
        "shard_count": shard_count,
        "species_treatment_requests": request_total,
        "species_treatments_saved": len(treatments),
        "pending_selected_treatments": len(pending_selected),
        "empty_descriptions": sum(
            error.get("error") == "not_a_confirmed_species_treatment"
            and int(error.get("narrative_chars", 0)) < args.min_description_chars
            for error in errors
        ),
        "HTTP_errors": client.http_errors,
        "species_with_any_trait": len(species_with_traits),
        "trait_candidate_rows": len(candidates),
        "unique_species_trait_cells": len({(row["accepted_species"], row["trait_name"]) for row in candidates}),
        "species_by_trait": {
            trait: len({row["accepted_species"] for row in candidates if row["trait_name"] == trait})
            for trait in traits
        },
        "error_records": len(errors),
        "source_record_contract": "eFloras lblTaxonDesc with bold accepted heading",
    }
    write_json(output_dir / "coverage_report.json", report)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if status != "ok":
        raise AcquisitionError(
            f"flora {flora_id}: status={status}; selected={len(sharded)} treatments_saved={len(treatments)}"
        )
    return report


def combine_results(root: Path, output_dir: Path) -> dict[str, Any]:
    treatment_paths = sorted(root.rglob("species_treatments.csv"))
    candidate_paths = sorted(root.rglob("trait_candidates.csv"))
    report_paths = sorted(root.rglob("coverage_report.json"))
    treatments = [row for path in treatment_paths for row in read_csv(path)]
    candidates = [row for path in candidate_paths for row in read_csv(path)]
    reports = []
    for path in report_paths:
        try:
            reports.append(json.loads(path.read_text(encoding="utf-8")))
        except Exception as exc:
            reports.append({"status": "report_read_error", "path": str(path), "error": str(exc)})
    treatments = list({(row.get("flora_id", ""), row.get("taxon_id", "")): row for row in treatments}.values())
    candidates = list(
        {
            (
                row.get("flora_id", ""),
                row.get("taxon_id", ""),
                row.get("trait_name", ""),
                row.get("candidate_value", ""),
                row.get("source_excerpt", ""),
            ): row
            for row in candidates
        }.values()
    )
    flora_by_species: dict[str, set[str]] = defaultdict(set)
    for row in treatments:
        flora_by_species[row.get("accepted_species", "")].add(row.get("flora_id", ""))
    traits = sorted({row.get("trait_name", "") for row in candidates if row.get("trait_name")})
    report = {
        "flora_reports": len(reports),
        "flora_statuses": dict(Counter(str(row.get("status", "unknown")) for row in reports)),
        "species_treatment_rows": len(treatments),
        "unique_species_excluding_cross_flora_duplicates": len({row.get("accepted_species", "") for row in treatments if row.get("accepted_species")}),
        "species_repeated_across_floras": sum(len(floras) > 1 for floras in flora_by_species.values()),
        "trait_candidate_rows": len(candidates),
        "unique_species_trait_cells": len({(row.get("accepted_species", ""), row.get("trait_name", "")) for row in candidates}),
        "species_by_trait": {
            trait: len({row.get("accepted_species", "") for row in candidates if row.get("trait_name") == trait})
            for trait in traits
        },
        "per_flora": reports,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "all_species_treatments.csv", treatments, TREATMENT_FIELDS)
    write_csv(output_dir / "all_trait_candidates.csv", candidates, CANDIDATE_FIELDS)
    write_json(output_dir / "combined_coverage.json", report)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return report


def comma_ids(value: str) -> set[int] | None:
    cleaned = clean(value)
    if not cleaned:
        return None
    return {int(item.strip()) for item in cleaned.split(",") if item.strip()}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--discover-floras", action="store_true")
    mode.add_argument("--combine-root", type=Path)
    mode.add_argument("--flora-id", type=int)
    parser.add_argument("--flora-ids", default="", help="Comma-separated manifest filter; missing IDs remain explicit failures")
    parser.add_argument("--output-dir", type=Path, default=Path("data/v2/staging/traits/bulk/efloras"))
    parser.add_argument("--master-csv", type=Path)
    parser.add_argument("--master-name-column", default="accepted_species")
    parser.add_argument("--synonyms-csv", type=Path)
    parser.add_argument("--max-pages", type=int, default=50)
    parser.add_argument("--max-depth", type=int, default=5)
    parser.add_argument("--max-species", type=int, default=0)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--checkpoint-every", type=int, default=10)
    parser.add_argument("--delay-min", type=float, default=1.5)
    parser.add_argument("--delay-max", type=float, default=3.5)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--attempts", type=int, default=6)
    parser.add_argument("--min-description-chars", type=int, default=20)
    parser.add_argument("--save-html", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--resume", action=argparse.BooleanOptionalAction, default=True)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        if args.discover_floras:
            client = Client(args.delay_min, args.delay_max, args.timeout, args.attempts)
            rows = discover_floras(client=client, output_dir=args.output_dir, requested_ids=comma_ids(args.flora_ids))
            print(json.dumps(rows, ensure_ascii=False, indent=2))
        elif args.combine_root:
            combine_results(args.combine_root, args.output_dir)
        else:
            run_flora(args)
    except (AcquisitionError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
