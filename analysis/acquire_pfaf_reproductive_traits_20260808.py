"""Acquire exact species-level self-fertility statements from public PFAF pages.

The crawler is deliberately narrow: it starts from a precomputed strict-unresolved
queue, accepts only exact binomial pages present in PFAF's public sitemap, and
extracts the generated physical-characteristics statement.  A positive
``self-fertile`` statement maps to self-compatibility (SC), never to autonomous
selfing.  Negative, cultivar-qualified, or identity-ambiguous text is retained in
the fetch inventory but is not emitted as evidence.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import html
import json
import re
import time
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

import pandas as pd
import requests

USER_AGENT = (
    "island-trait-research/0.1 (scientific evidence audit; https://github.com/zuizui0223/island)"
)
SELF_FERTILE_RE = re.compile(
    r"\b(?:the plant|the species|it) is (?:generally )?self-fertile\b", re.IGNORECASE
)
PHYSICAL_RE = re.compile(
    r'<span[^>]+id="ContentPlaceHolder1_lblPhystatment"[^>]*>(.*?)</span>',
    re.IGNORECASE | re.DOTALL,
)
H1_RE = re.compile(
    r'<span[^>]+id="ContentPlaceHolder1_lbldisplatinname"[^>]*>(.*?)</span>',
    re.IGNORECASE | re.DOTALL,
)
CULTIVAR_RE = re.compile(r"\b(?:cultivar|cv\.|hybrid|variety only)\b", re.IGNORECASE)
NONSEXUAL_RE = re.compile(r"\b(?:apomict\w*|apomix\w*|without sexual)\b", re.IGNORECASE)
DIOECIOUS_RE = re.compile(r"\bdioecious\b", re.IGNORECASE)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _visible(fragment: str) -> str:
    fragment = re.sub(r"<(?:br|/p)\s*/?>", " ", fragment, flags=re.IGNORECASE)
    return " ".join(html.unescape(re.sub(r"<[^>]+>", " ", fragment)).split())


def sitemap_species(sitemap: Path) -> dict[str, str]:
    raw = sitemap.read_text(encoding="utf-8", errors="replace")
    urls = re.findall(r"<loc>(.*?)</loc>", raw, flags=re.IGNORECASE)
    result: dict[str, str] = {}
    for url in urls:
        if not re.search(r"/user/plant\.aspx\?LatinName=", url, re.IGNORECASE):
            continue
        name = parse_qs(urlsplit(url).query).get("LatinName", [""])[0].replace("+", " ").strip()
        if len(name.split()) != 2:
            continue
        result.setdefault(name, url)
    return result


def parse_page(species: str, url: str, payload: bytes) -> dict[str, object]:
    raw = payload.decode("utf-8", errors="replace")
    h1_match = H1_RE.search(raw)
    h1 = _visible(h1_match.group(1)) if h1_match else ""
    matched_page_name = " ".join(h1.split()[:2])
    identity_exact = matched_page_name.casefold() == species.casefold()
    physical_match = PHYSICAL_RE.search(raw)
    statement = _visible(physical_match.group(1)) if physical_match else ""
    hit = bool(SELF_FERTILE_RE.search(statement))
    cultivar_qualified = bool(CULTIVAR_RE.search(statement))
    nonsexual_qualified = bool(NONSEXUAL_RE.search(statement))
    dioecious_conflict = bool(DIOECIOUS_RE.search(statement))
    if not identity_exact:
        status = "identity_mismatch"
    elif not statement:
        status = "physical_statement_missing"
    elif cultivar_qualified:
        status = "cultivar_or_hybrid_qualified"
    elif hit and nonsexual_qualified:
        status = "nonsexual_reproduction_category_error"
    elif hit and dioecious_conflict:
        status = "sexual_system_internal_conflict"
    elif hit:
        status = "accepted_self_fertile_statement"
    elif re.search(r"\bnot self-fertile\b", statement, re.IGNORECASE):
        status = "negative_self_fertile_not_mapped_to_si"
    else:
        status = "no_positive_self_fertile_statement"
    quote = SELF_FERTILE_RE.search(statement)
    return {
        "accepted_species": species,
        "source_url": url,
        "matched_page_name": matched_page_name,
        "identity_exact": identity_exact,
        "supporting_excerpt": quote.group(0) if quote else "",
        "physical_statement": statement,
        "normalized_value": "SC" if status == "accepted_self_fertile_statement" else "",
        "trait_name": "self_incompatibility" if status == "accepted_self_fertile_statement" else "",
        "axis": "reproductive_assurance" if status == "accepted_self_fertile_statement" else "",
        "status": status,
        "cultivar_qualified": cultivar_qualified,
        "nonsexual_qualified": nonsexual_qualified,
        "dioecious_conflict": dioecious_conflict,
        "page_sha256": _sha256_bytes(payload),
        "content_fingerprint": _sha256_bytes(statement.casefold().encode("utf-8")),
    }


def fetch_one(species: str, url: str, cache_dir: Path, timeout: float) -> dict[str, object]:
    cache_path = cache_dir / f"{hashlib.sha256(url.encode('utf-8')).hexdigest()}.html"
    from_cache = cache_path.exists()
    try:
        if from_cache:
            payload = cache_path.read_bytes()
            http_status = 200
        else:
            response = requests.get(url, headers={"User-Agent": USER_AGENT}, timeout=timeout)
            http_status = response.status_code
            response.raise_for_status()
            payload = response.content
            cache_path.write_bytes(payload)
            time.sleep(0.15)
        row = parse_page(species, url, payload)
        row.update(
            {
                "http_status": http_status,
                "retrieved_at_utc": datetime.now(UTC).isoformat(),
                "from_cache": from_cache,
                "cache_file": cache_path.name,
                "error": "",
            }
        )
        return row
    except Exception as exc:  # noqa: BLE001 - every task must remain accounted for
        return {
            "accepted_species": species,
            "source_url": url,
            "matched_page_name": "",
            "identity_exact": False,
            "supporting_excerpt": "",
            "physical_statement": "",
            "normalized_value": "",
            "trait_name": "",
            "axis": "",
            "status": "fetch_error",
            "cultivar_qualified": False,
            "nonsexual_qualified": False,
            "dioecious_conflict": False,
            "page_sha256": "",
            "content_fingerprint": "",
            "http_status": "",
            "retrieved_at_utc": datetime.now(UTC).isoformat(),
            "from_cache": from_cache,
            "cache_file": cache_path.name,
            "error": f"{type(exc).__name__}: {exc}",
        }


def build(
    *,
    queue_csv: Path,
    sitemap_xml: Path,
    robots_txt: Path,
    cache_dir: Path,
    output_dir: Path,
    min_islands: int,
    max_workers: int,
    timeout: float,
) -> dict[str, object]:
    queue = pd.read_csv(queue_csv, dtype={"accepted_species": str}).fillna("")
    queue["n_islands"] = pd.to_numeric(queue["n_islands"], errors="raise").astype(int)
    queue = queue.loc[queue["n_islands"].ge(min_islands)].copy()
    sitemap = sitemap_species(sitemap_xml)
    queue["source_url"] = queue["accepted_species"].map(sitemap).fillna("")
    queue = queue.loc[queue["source_url"].ne("")].drop_duplicates("accepted_species")
    queue = queue.sort_values(["n_islands", "accepted_species"], ascending=[False, True])

    cache_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    tasks = [(row.accepted_species, row.source_url) for row in queue.itertuples(index=False)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        rows = list(
            executor.map(
                lambda task: fetch_one(task[0], task[1], cache_dir, timeout),
                tasks,
            )
        )
    inventory = queue.drop(columns="source_url").merge(
        pd.DataFrame(rows), on="accepted_species", how="left", validate="one_to_one"
    )
    inventory.to_csv(output_dir / "pfaf_fetch_inventory.csv.gz", index=False)
    evidence = inventory.loc[inventory["status"].eq("accepted_self_fertile_statement")].copy()
    evidence.to_csv(output_dir / "pfaf_self_fertile_candidates.csv.gz", index=False)
    summary: dict[str, object] = {
        "contract": "pfaf_exact_species_self_fertility_acquisition_v1",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "queue_rows": len(queue),
        "fetched_or_cached": int(inventory["http_status"].astype(str).eq("200").sum()),
        "fetch_errors": int(inventory["status"].eq("fetch_error").sum()),
        "identity_exact": int(inventory["identity_exact"].astype(bool).sum()),
        "accepted_self_fertile_statements": len(evidence),
        "accepted_species": int(evidence["accepted_species"].nunique()),
        "accepted_n_islands_sum": int(evidence["n_islands"].sum()),
        "status_counts": inventory["status"].value_counts().sort_index().to_dict(),
        "min_islands": min_islands,
        "mapping": "positive self-fertile -> self_incompatibility=SC; never autonomous selfing",
        "sitemap_sha256": _sha256_bytes(sitemap_xml.read_bytes()),
        "robots_sha256": _sha256_bytes(robots_txt.read_bytes()),
        "robots_url": "https://pfaf.org/robots.txt",
        "sitemap_url": "https://pfaf.org/sitemap.xml",
        "source_access": "public species pages allowed by robots.txt",
    }
    (output_dir / "pfaf_acquisition_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue-csv", type=Path, required=True)
    parser.add_argument("--sitemap-xml", type=Path, required=True)
    parser.add_argument("--robots-txt", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-islands", type=int, default=11)
    parser.add_argument("--max-workers", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=45.0)
    args = parser.parse_args()
    print(json.dumps(build(**vars(args)), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
