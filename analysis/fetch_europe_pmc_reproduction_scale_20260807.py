"""Fetch and screen Europe PMC full text for species-direct reproduction traits.

The acquisition unit is a high-information genus x trait queue entry, not a
fixed species batch. Search results are discovery only. Candidate evidence is
emitted only from fetched full-text XML when an exact master species name (or
an article-resolved abbreviation) and an explicit trait phrase occur in the
same sentence. Contrastive clauses are separated before attribution.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

REST = "https://www.ebi.ac.uk/europepmc/webservices/rest"
USER_AGENT = "island-trait-acquisition/0.1 (Europe PMC official REST API)"
REPRO_TRAITS = {
    "self_incompatibility",
    "autonomous_selfing_capacity",
    "mating_system",
    "cleistogamy",
}
QUERY_TERMS = {
    "self_incompatibility": (
        '"self-incompatibility" OR "self incompatibility" OR '
        '"self-incompatible" OR "self-compatible" OR "self compatible" OR '
        '"self-fertile"'
    ),
    "autonomous_selfing_capacity": (
        '"autonomous selfing" OR "autonomous self-pollination" OR '
        '"spontaneous selfing" OR "automatic self-pollination" OR '
        '"delayed selfing" OR "bagged flowers"'
    ),
    "mating_system": (
        '"mating system" OR "mixed mating" OR "selfing rate" OR '
        '"outcrossing rate" OR "predominantly selfing" OR '
        '"predominantly outcrossing"'
    ),
    "cleistogamy": ('"cleistogamy" OR "cleistogamous" OR "chasmogamous"'),
}
VALUE_PATTERNS = {
    "self_incompatibility": (
        ("SI", re.compile(r"\bself[- ]incompatible\b", re.IGNORECASE)),
        ("SC", re.compile(r"\bself[- ]compatible\b|\bself[- ]fertile\b", re.IGNORECASE)),
    ),
    "autonomous_selfing_capacity": (
        (
            "absent",
            re.compile(
                r"\b(no|not|without|unable to)\b.{0,35}"
                r"\bautonomous self(?:ing|-pollination)\b",
                re.IGNORECASE,
            ),
        ),
        (
            "delayed",
            re.compile(r"\bdelayed self(?:ing|-pollination)\b", re.IGNORECASE),
        ),
        (
            "facilitated",
            re.compile(r"\bfacilitated selfing\b|\bgeitonogamous selfing\b", re.IGNORECASE),
        ),
        (
            "autonomous",
            re.compile(
                r"\bautonomous self(?:ing|-pollination)\b|"
                r"\bautomatic self-pollination\b|\bspontaneous selfing\b",
                re.IGNORECASE,
            ),
        ),
    ),
    "mating_system": (
        (
            "predominantly_outcrossing",
            re.compile(r"\b(predominantly|mainly|mostly) outcross(?:ing|ed)\b", re.IGNORECASE),
        ),
        (
            "predominantly_selfing",
            re.compile(r"\b(predominantly|mainly|mostly) self(?:ing|ed)\b", re.IGNORECASE),
        ),
        (
            "obligate_selfing",
            re.compile(r"\bobligate(?:ly)? self(?:ing|er)\b", re.IGNORECASE),
        ),
        (
            "mixed_mating",
            re.compile(r"\bmixed[- ]mating(?: system)?\b", re.IGNORECASE),
        ),
    ),
    "cleistogamy": (
        (
            "facultative",
            re.compile(
                r"\b(chasmogamous.{0,80}cleistogamous|"
                r"cleistogamous.{0,80}chasmogamous|facultative cleistogam)\b",
                re.IGNORECASE,
            ),
        ),
        (
            "obligate",
            re.compile(r"\bobligate(?:ly)? cleistogam(?:ous|y)\b", re.IGNORECASE),
        ),
        (
            "absent",
            re.compile(r"\b(no|not|without)\b.{0,30}\bcleistogam(?:ous|y)\b", re.IGNORECASE),
        ),
    ),
}


def text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def sha(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def request_bytes(url: str, *, attempts: int = 4) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT, "Accept": "*/*"},
    )
    last: Exception | None = None
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError) as exc:
            last = exc
            if attempt + 1 == attempts:
                break
            time.sleep(2**attempt)
    raise RuntimeError(f"request failed after {attempts} attempts: {url}: {last}")


def search(query: str) -> dict[str, Any]:
    url = (
        REST
        + "/search?"
        + urllib.parse.urlencode(
            {
                "query": query,
                "format": "json",
                "pageSize": 1000,
                "resultType": "core",
            }
        )
    )
    return json.loads(request_bytes(url).decode("utf-8"))


def build_tasks(
    queue_path: Path,
    top_n: int,
    completed_query_logs: list[Path] | None = None,
    max_tasks: int | None = None,
    current_direct_path: Path | None = None,
    require_agreement: bool = False,
    support_levels: tuple[int, ...] = (2,),
) -> pd.DataFrame:
    queue = pd.read_csv(queue_path, dtype=str).fillna("").head(top_n)
    # Historical queues search one genus x trait.  Later acquisition waves may
    # nominate an exact unresolved congener so the same literature backend can
    # look for genuinely new direct evidence instead of rediscovering the two
    # species that already support the rule.  Empty values preserve the frozen
    # Wave37 query and checkpoint contract byte-for-byte at the task level.
    if "target_species" not in queue:
        queue["target_species"] = ""
    if "query_scope" not in queue:
        queue["query_scope"] = ""
    if current_direct_path is not None:
        direct = pd.read_csv(current_direct_path, dtype=str).fillna("")
        direct = direct.loc[direct["resolution_status"].eq("resolved")].drop_duplicates(
            ["accepted_species", "genus", "axis", "trait_name"]
        )
        current = direct.groupby(["genus", "axis", "trait_name"], as_index=False).agg(
            latest_support=("accepted_species", "nunique"),
            latest_value_count=("state_set", "nunique"),
        )
        queue["_queue_order"] = range(len(queue))
        queue = queue.merge(
            current,
            on=["genus", "axis", "trait_name"],
            how="left",
            validate="many_to_one",
        ).sort_values("_queue_order", kind="stable")
        queue["current_support"] = (
            pd.to_numeric(queue["latest_support"], errors="coerce")
            .fillna(0)
            .astype(int)
            .astype(str)
        )
        if require_agreement:
            queue = queue.loc[queue["latest_value_count"].eq(1)].copy()
    allowed_support = {str(value) for value in support_levels}
    if not allowed_support or not allowed_support <= {"1", "2"}:
        raise ValueError("support_levels must contain only 1 and/or 2")
    queue = queue.loc[
        queue["trait_name"].isin(REPRO_TRAITS) & queue["current_support"].isin(allowed_support)
    ].copy()
    queue["query_name"] = queue["genus"]
    exact = queue["query_scope"].eq("exact_species") & queue["target_species"].ne("")
    queue.loc[exact, "query_name"] = queue.loc[exact, "target_species"]
    queue["query"] = [
        f'("{name}") AND ({QUERY_TERMS[trait]}) AND OPEN_ACCESS:Y AND IN_EPMC:Y'
        for name, trait in zip(queue["query_name"], queue["trait_name"], strict=True)
    ]
    task_key = ["genus", "trait_name", "query_scope", "target_species"]
    queue = queue.drop_duplicates(task_key)
    completed: set[tuple[str, str, str, str]] = set()
    for path in completed_query_logs or []:
        columns = set(pd.read_csv(path, nrows=0).columns)
        usecols = [column for column in task_key if column in columns]
        log = pd.read_csv(path, usecols=usecols, dtype=str).fillna("")
        for column in task_key:
            if column not in log:
                log[column] = ""
        completed.update(log[task_key].itertuples(index=False, name=None))
    if completed:
        keys = queue[task_key].apply(tuple, axis=1)
        queue = queue.loc[~keys.isin(completed)].copy()
    if max_tasks is not None:
        queue = queue.head(max_tasks)
    return queue.reset_index(drop=True)


def load_species(
    coverage_path: Path,
    genus_path: Path,
) -> tuple[dict[str, str], dict[str, list[str]], dict[str, str]]:
    coverage = (
        pd.read_csv(
            coverage_path,
            usecols=["accepted_species"],
            dtype=str,
        )
        .fillna("")
        .drop_duplicates()
    )
    genus_rows = (
        pd.read_csv(
            genus_path,
            usecols=["accepted_species", "genus"],
            dtype=str,
        )
        .fillna("")
        .drop_duplicates("accepted_species")
    )
    species = coverage.merge(genus_rows, on="accepted_species", how="left")
    species["genus"] = species["genus"].where(
        species["genus"].ne(""),
        species["accepted_species"].str.split().str[0],
    )
    canonical: dict[str, str] = {}
    by_genus: dict[str, list[str]] = defaultdict(list)
    genus_map: dict[str, str] = {}
    for row in species.itertuples(index=False):
        name = text(row.accepted_species).strip('" ')
        tokens = name.replace("× ", "×").split()
        if len(tokens) < 2 or tokens[1] == "×":
            continue
        binomial = f"{tokens[0]} {tokens[1]}"
        canonical[binomial.casefold()] = text(row.accepted_species)
        by_genus[text(row.genus)].append(binomial)
        genus_map[text(row.accepted_species)] = text(row.genus)
    return canonical, by_genus, genus_map


def xml_metadata(raw: bytes) -> tuple[str, str, str]:
    root = ET.fromstring(raw)
    title = text(" ".join(root.findtext(".//article-title", default="").split()))
    doi = ""
    pmcid = ""
    for node in root.findall(".//article-id"):
        kind = text(node.attrib.get("pub-id-type")).casefold()
        if kind == "doi":
            doi = text(node.text)
        elif kind in {"pmc", "pmcid"}:
            pmcid = text(node.text)
    return title, doi, pmcid


def paragraphs(raw: bytes) -> list[str]:
    root = ET.fromstring(raw)
    values: list[str] = []
    for node in root.findall(".//p") + root.findall(".//tr"):
        value = text(" ".join(node.itertext()))
        if value and value not in values:
            values.append(value)
    return values


def sentences(value: str) -> list[str]:
    return [
        text(part) for part in re.split(r"(?<=[.!?])\s+(?=[A-Z])", value) if len(text(part)) >= 20
    ]


def attribution_clauses(value: str) -> list[str]:
    """Split contrastive clauses that commonly assign different taxa values."""

    return [
        text(part)
        for part in re.split(
            r"\s*;\s*|\s*,?\s+\b(?:whereas|while|but)\b\s+",
            value,
            flags=re.IGNORECASE,
        )
        if text(part)
    ]


def document_genus_binomials(value: str, genus: str) -> list[str]:
    """Return explicit source binomials for a queued genus, including non-target taxa."""

    matches = re.finditer(
        rf"(?<![A-Za-z]){re.escape(genus)}\s+([a-z][a-z-]{{2,}})\b",
        value,
        flags=re.IGNORECASE,
    )
    excluded = {"species", "spec", "spp", "unknown"}
    return sorted(
        {
            f"{genus} {match.group(1)}"
            for match in matches
            if match.group(1).casefold() not in excluded
        }
    )


def article_candidates(
    raw: bytes,
    *,
    pmcid: str,
    title: str,
    doi: str,
    task_rows: pd.DataFrame,
    canonical: dict[str, str],
    by_genus: dict[str, list[str]],
    article_license: str,
    retrieved_at: str,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str, str, str]] = set()
    article_paragraphs = paragraphs(raw)
    article_text = text(" ".join(article_paragraphs))
    task_contexts: list[tuple[Any, str, tuple[tuple[str, re.Pattern[str]], ...], list[str]]] = []
    for task in task_rows.itertuples(index=False):
        trait = text(task.trait_name)
        task_genus = text(task.genus)
        target_species = text(getattr(task, "target_species", ""))
        query_scope = text(getattr(task, "query_scope", ""))
        task_contexts.append(
            (
                task,
                trait,
                VALUE_PATTERNS[trait],
                (
                    [target_species]
                    if query_scope == "exact_species" and target_species
                    else document_genus_binomials(article_text, task_genus)
                ),
            )
        )
    for paragraph_index, paragraph in enumerate(article_paragraphs):
        for task, trait, patterns, species_names in task_contexts:
            if not any(pattern.search(paragraph) for _, pattern in patterns):
                continue
            present: dict[str, tuple[str, str]] = {}
            for name in species_names:
                fixed_master_name = canonical.get(name.casefold(), "")
                if match := re.search(
                    rf"\b{re.escape(name)}\b",
                    paragraph,
                    re.IGNORECASE,
                ):
                    present[name] = (
                        match.group(0),
                        (
                            "accepted_binomial_exact_in_full_text"
                            if fixed_master_name
                            else "source_binomial_exact_pending_two_backbone"
                        ),
                    )
                    continue
                genus, epithet = name.split(maxsplit=1)
                abbreviation = re.search(
                    rf"\b{re.escape(genus[0])}\.\s*{re.escape(epithet)}\b",
                    paragraph,
                    re.IGNORECASE,
                )
                if not abbreviation:
                    continue
                full_name_established = bool(
                    re.search(rf"\b{re.escape(name)}\b", article_text, re.IGNORECASE)
                )
                genus_established_here = bool(
                    re.search(rf"\b{re.escape(genus)}\b", paragraph, re.IGNORECASE)
                )
                if full_name_established or genus_established_here:
                    present[name] = (
                        abbreviation.group(0),
                        (
                            (
                                "accepted_binomial_established_then_unambiguous_abbreviation"
                                if fixed_master_name
                                else "source_binomial_established_then_abbreviation_pending_two_backbone"
                            )
                            if full_name_established
                            else "full_genus_and_unambiguous_abbreviation_in_same_paragraph"
                        ),
                    )
            if not present:
                continue
            paragraph_sentences = sentences(paragraph)
            for value, pattern in patterns:
                for sentence_index, sentence in enumerate(paragraph_sentences):
                    sentence_clauses = attribution_clauses(sentence)
                    for clause in sentence_clauses:
                        match = pattern.search(clause)
                        if not match:
                            continue
                        same_clause = [
                            name
                            for name in present
                            if re.search(
                                rf"\b{re.escape(present[name][0])}\b",
                                clause,
                                re.IGNORECASE,
                            )
                        ]
                        same_sentence = (
                            [
                                name
                                for name in present
                                if re.search(
                                    rf"\b{re.escape(present[name][0])}\b",
                                    sentence,
                                    re.IGNORECASE,
                                )
                            ]
                            if len(sentence_clauses) == 1
                            else []
                        )
                        attributed = same_clause or same_sentence
                        if not attributed:
                            continue
                        excerpt = sentence
                        for matched_name in attributed:
                            accepted = canonical.get(matched_name.casefold(), matched_name)
                            in_fixed_master = matched_name.casefold() in canonical
                            lineage = f"doi:{doi}" if doi else f"pmcid:{pmcid}"
                            key = (accepted, trait, value, lineage)
                            if key in seen:
                                continue
                            seen.add(key)
                            rows.append(
                                {
                                    "candidate_id": (
                                        "EPMC-"
                                        + hashlib.sha256("|".join(key).encode("utf-8")).hexdigest()[
                                            :20
                                        ]
                                    ),
                                    "accepted_species": accepted,
                                    "source_species_name": matched_name,
                                    "in_fixed_master": str(in_fixed_master).lower(),
                                    "rule_support_only": str(not in_fixed_master).lower(),
                                    "genus": text(task.genus),
                                    "trait_name": trait,
                                    "axis": "reproductive_assurance",
                                    "raw_value": value,
                                    "normalized_value": value,
                                    "pmcid": pmcid,
                                    "doi": doi,
                                    "article_title": title,
                                    "source_url": f"https://europepmc.org/articles/{pmcid}",
                                    "fulltext_xml_url": f"{REST}/{pmcid}/fullTextXML",
                                    "exact_supporting_quote": excerpt,
                                    "matched_name": present[matched_name][0],
                                    "name_match_method": present[matched_name][1],
                                    "source_lineage": lineage,
                                    "article_license": article_license,
                                    "content_sha256": sha(raw),
                                    "paragraph_index": str(paragraph_index),
                                    "sentence_index": str(sentence_index),
                                    "attribution_method": (
                                        "same_clause_exact_species_and_trait"
                                        if same_clause
                                        else "same_sentence_exact_species_and_trait"
                                    ),
                                    "retrieved_at_utc": retrieved_at,
                                    "review_status": "source_backed_fulltext_review_pending",
                                    "decision_reason": (
                                        "Explicit source binomial and trait phrase in the same "
                                        "fetched Europe PMC sentence; experimental context, "
                                        "conflict review, and two-backbone identity review are "
                                        "required before non-master congeners support a rule."
                                    ),
                                }
                            )
    return rows


def write_manifest(output: Path) -> None:
    rows = []
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "file_manifest.sha256":
            relative = path.relative_to(output).as_posix()
            rows.append(f"{sha(path.read_bytes())}  {relative}")
    (output / "file_manifest.sha256").write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", type=Path, required=True)
    parser.add_argument("--coverage", type=Path, required=True)
    parser.add_argument("--genus-map", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument(
        "--completed-query-log",
        type=Path,
        action="append",
        default=[],
        help="Exclude terminal genus x trait tasks recorded by earlier runs.",
    )
    parser.add_argument("--max-tasks", type=int)
    parser.add_argument(
        "--current-direct",
        type=Path,
        help="Recalculate current support and state agreement from this direct ledger.",
    )
    parser.add_argument(
        "--require-agreement",
        action="store_true",
        help="Keep only genus x trait tasks whose current species values agree.",
    )
    parser.add_argument(
        "--support-level",
        type=int,
        action="append",
        choices=(1, 2),
        dest="support_levels",
        help=(
            "Current direct-species support to search. Repeat for both 1 and 2; "
            "the default remains support 2 for backward compatibility."
        ),
    )
    parser.add_argument("--delay", type=float, default=0.15)
    parser.add_argument(
        "--search-workers",
        type=int,
        default=1,
        help="Concurrent official Europe PMC search requests; cached tasks are not resent.",
    )
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    raw_dir = args.output / "raw_xml"
    raw_dir.mkdir(exist_ok=True)
    search_dir = args.output / "search_results"
    search_dir.mkdir(exist_ok=True)

    tasks = build_tasks(
        args.queue,
        args.top_n,
        args.completed_query_log,
        args.max_tasks,
        args.current_direct,
        args.require_agreement,
        tuple(args.support_levels or (2,)),
    )
    canonical, by_genus, _ = load_species(args.coverage, args.genus_map)
    query_rows: list[dict[str, Any]] = []
    articles: dict[str, dict[str, Any]] = {}
    task_key = ["genus", "trait_name", "query_scope", "target_species"]
    article_tasks: dict[str, set[tuple[str, str, str, str]]] = defaultdict(set)

    def process_search(task: Any) -> tuple[dict[str, Any], list[dict[str, Any]]]:
        search_key = hashlib.sha256(
            (
                f"{task.genus}|{task.trait_name}|{task.query_scope}|"
                f"{task.target_species}|{task.query}"
            ).encode()
        ).hexdigest()[:24]
        search_path = search_dir / f"{search_key}.json.gz"
        try:
            if search_path.exists():
                result = json.loads(gzip.decompress(search_path.read_bytes()).decode("utf-8"))
                search_status = "cached"
            else:
                result = search(text(task.query))
                search_path.write_bytes(
                    gzip.compress(
                        json.dumps(
                            result,
                            ensure_ascii=False,
                            sort_keys=True,
                        ).encode("utf-8"),
                        mtime=0,
                    )
                )
                search_status = "fetched"
            hits = int(result.get("hitCount") or 0)
            records = result.get("resultList", {}).get("result", [])
            status = "success"
            error = ""
        except Exception as exc:  # noqa: BLE001 - recorded, never silently skipped
            hits, records, status, error, search_status = (
                0,
                [],
                "error",
                str(exc),
                "error",
            )
        if search_status == "fetched":
            time.sleep(args.delay)
        return (
            {
                "queue_order": int(task.Index),
                "genus": task.genus,
                "trait_name": task.trait_name,
                "query_scope": task.query_scope,
                "target_species": task.target_species,
                "query": task.query,
                "status": status,
                "search_cache_status": search_status,
                "search_cache_path": str(search_path),
                "hit_count": hits,
                "returned_count": len(records),
                "cost_usd": 0.0,
                "error": error,
            },
            records,
        )

    search_workers = max(1, args.search_workers)
    search_futures = {}
    indexed_tasks = tasks.reset_index(drop=True).itertuples(index=True)
    with ThreadPoolExecutor(max_workers=search_workers) as executor:
        for task in indexed_tasks:
            search_futures[executor.submit(process_search, task)] = task.Index
        for completed, future in enumerate(as_completed(search_futures), start=1):
            query_row, records = future.result()
            query_rows.append(query_row)
            for record in records:
                pmcid = text(record.get("pmcid"))
                if not pmcid:
                    continue
                articles.setdefault(pmcid, record)
                article_tasks[pmcid].add(
                    tuple(text(query_row[column]) for column in task_key)
                )
            if completed % 50 == 0:
                print(f"[europe-pmc-search] {completed}/{len(tasks)} queries", flush=True)

    query_rows.sort(key=lambda row: int(row["queue_order"]))
    for row in query_rows:
        row.pop("queue_order", None)

    manifest_rows: list[dict[str, Any]] = []
    candidate_rows: list[dict[str, str]] = []
    task_lookup = tasks.set_index(task_key, drop=False)

    def process_article(
        pmcid: str,
        record: dict[str, Any],
    ) -> tuple[dict[str, Any], list[dict[str, str]]]:
        raw_path = raw_dir / f"{pmcid}.xml.gz"
        retrieved_at = now()
        status = "cached"
        error = ""
        try:
            if raw_path.exists():
                try:
                    raw = gzip.decompress(raw_path.read_bytes())
                except (EOFError, gzip.BadGzipFile):
                    raw = request_bytes(f"{REST}/{pmcid}/fullTextXML")
                    temporary = raw_path.with_suffix(".xml.gz.tmp")
                    temporary.write_bytes(gzip.compress(raw, mtime=0))
                    temporary.replace(raw_path)
                    status = "refetched_corrupt_cache"
                    time.sleep(args.delay)
            else:
                raw = request_bytes(f"{REST}/{pmcid}/fullTextXML")
                temporary = raw_path.with_suffix(".xml.gz.tmp")
                temporary.write_bytes(gzip.compress(raw, mtime=0))
                temporary.replace(raw_path)
                status = "fetched"
                time.sleep(args.delay)
            title, doi_xml, pmcid_xml = xml_metadata(raw)
            doi = doi_xml or text(record.get("doi"))
            license_text = text(record.get("license"))
            selected_tasks = pd.concat(
                [
                    task_lookup.loc[[pair]]
                    for pair in sorted(article_tasks[pmcid])
                    if pair in task_lookup.index
                ],
                ignore_index=True,
            )
            article_rows = article_candidates(
                raw,
                pmcid=pmcid_xml or pmcid,
                title=title or text(record.get("title")),
                doi=doi,
                task_rows=selected_tasks,
                canonical=canonical,
                by_genus=by_genus,
                article_license=license_text,
                retrieved_at=retrieved_at,
            )
        except Exception as exc:  # noqa: BLE001 - per-article failure is manifested
            raw = b""
            article_rows = []
            title, doi, license_text = "", text(record.get("doi")), text(record.get("license"))
            status, error = "error", str(exc)
        manifest_row = {
            "pmcid": pmcid,
            "doi": doi,
            "title": title or text(record.get("title")),
            "source_url": f"https://europepmc.org/articles/{pmcid}",
            "fulltext_xml_url": f"{REST}/{pmcid}/fullTextXML",
            "fetch_status": status,
            "article_license": license_text,
            "content_sha256": sha(raw) if raw else "",
            "raw_path": str(raw_path),
            "retrieved_at_utc": retrieved_at,
            "task_count": len(article_tasks[pmcid]),
            "error": error,
        }
        return manifest_row, article_rows

    futures = {}
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
        for pmcid, record in sorted(articles.items()):
            futures[executor.submit(process_article, pmcid, record)] = pmcid
        for index, future in enumerate(as_completed(futures), start=1):
            manifest_row, article_rows = future.result()
            manifest_rows.append(manifest_row)
            candidate_rows.extend(article_rows)
            if index % 50 == 0:
                print(f"[europe-pmc] {index}/{len(articles)} articles", flush=True)

    manifest_rows.sort(key=lambda row: text(row["pmcid"]))
    candidate_rows.sort(
        key=lambda row: (
            text(row["accepted_species"]),
            text(row["trait_name"]),
            text(row["normalized_value"]),
            text(row["source_lineage"]),
        )
    )

    query_frame = pd.DataFrame(query_rows)
    manifest = pd.DataFrame(manifest_rows)
    candidates = (
        pd.DataFrame(candidate_rows).drop_duplicates(
            ["accepted_species", "trait_name", "normalized_value", "source_lineage"]
        )
        if candidate_rows
        else pd.DataFrame()
    )
    query_frame.to_csv(args.output / "query_log.csv", index=False)
    manifest.to_csv(args.output / "article_manifest.csv", index=False)
    candidates.to_csv(
        args.output / "fulltext_reproductive_candidates.csv.gz",
        index=False,
        compression="gzip",
    )
    summary = {
        "contract": "europe_pmc_reproduction_source_scale_v1",
        "created_at_utc": now(),
        "queue_entries": len(tasks),
        "queries": len(query_frame),
        "search_hits": int(
            pd.to_numeric(query_frame["hit_count"], errors="coerce").fillna(0).sum()
        ),
        "unique_pmcids": len(articles),
        "xml_fetched_or_cached": int(
            manifest["fetch_status"].isin({"fetched", "cached", "refetched_corrupt_cache"}).sum()
        ),
        "xml_errors": int(manifest["fetch_status"].eq("error").sum()),
        "candidate_rows": len(candidates),
        "candidate_species": (
            int(candidates["accepted_species"].nunique()) if not candidates.empty else 0
        ),
        "candidate_species_trait": (
            int(candidates[["accepted_species", "trait_name"]].drop_duplicates().shape[0])
            if not candidates.empty
            else 0
        ),
        "candidate_by_trait": (
            candidates.groupby("trait_name").size().sort_index().to_dict()
            if not candidates.empty
            else {}
        ),
        "query_cost_usd": 0.0,
        "search_workers": search_workers,
        "workers": max(1, args.workers),
        "formal_promotion_status": "review_pending",
        "safeguards": {
            "search_snippets_as_evidence": False,
            "fulltext_required": True,
            "exact_master_binomial_required": True,
            "self_fertile_to_autonomous_selfing": False,
            "trait_substitution": False,
            "family_inference": False,
            "global_fallback": False,
        },
    }
    (args.output / "run_manifest.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    write_manifest(args.output)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
