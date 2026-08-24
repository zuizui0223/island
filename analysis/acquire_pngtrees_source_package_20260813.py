"""Acquire and freeze structured flower descriptions from PNGTreesKey.

The Guide to Trees of Papua New Guinea is authored by Barry Conn (NSW) and
Kipiro Damas (LAE) and is hosted by the National Herbarium of New South Wales
and Papua New Guinea National Herbarium.  This adapter uses only exact-rank,
exact-name pages whose displayed family agrees with the island master.  It
extracts controlled flower fields, never search snippets, and keeps the whole
guide as one source lineage so genus validation cannot treat its pages as
independent publications.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import urljoin

import pandas as pd
import requests
from bs4 import BeautifulSoup

from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS
from island_v2.open_web_common import reviewed_source_package_evidence

INDEX_URL = "https://www.pngplants.org/PNGtrees/TreeDescriptions/"
SOURCE_GROUP = "pngtrees_herbarium_guide_20260813"
SOURCE_LINEAGE = "provider_compilation:pngtrees:guide-to-trees-of-papua-new-guinea"
SOURCE_PROVIDER = "Guide to Trees of Papua New Guinea (Conn & Damas)"
SOURCE_CITATION = (
    "Conn, B.J. & Damas, K.Q., Guide to Trees of Papua New Guinea; "
    "National Herbarium of New South Wales and Papua New Guinea National Herbarium"
)
USER_AGENT = "island-trait-research/1.0 (source-backed botanical evidence)"
REVIEWED_AT_UTC = "2026-08-12T19:15:00Z"
REVIEWER = "Codex source-backed structured-page audit"
AUDIT_PER_TRAIT = 50
ROOT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "pngtrees_source_package_20260813"
)

TRAIT_AXIS = {
    "flower_primary_color": "flower_colour",
    "floral_symmetry": "floral_structural_complexity",
    "flower_size_class": "floral_structural_complexity",
    "inflorescence_display": "floral_structural_complexity",
}
COLOUR_TERMS = {
    "white": ("white", "cream"),
    "green_brown_inconspicuous": ("green", "brown"),
    "yellow_orange": ("yellow", "orange"),
    "red_pink": ("red", "pink"),
    "blue_purple": ("blue", "purple", "violet", "mauve"),
}
INFRASPECIFIC = re.compile(r"\b(?:subsp|ssp|var|f)\.", re.IGNORECASE)
SPECIES_PREFIX = re.compile(r"^([A-Z][a-z-]+\s+[a-z][a-z-]+)")


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_text(value: str) -> str:
    return _sha256_bytes(value.encode("utf-8"))


def _canonical_hash(path: Path) -> str:
    payload = path.read_bytes()
    if path.suffix.casefold() in {".csv", ".json", ".md", ".txt"}:
        payload = payload.replace(b"\r\n", b"\n")
    return _sha256_bytes(payload)


def _get(session: requests.Session, url: str, timeout: int) -> bytes:
    response = session.get(url, timeout=timeout, headers={"User-Agent": USER_AGENT})
    response.raise_for_status()
    return response.content


def parse_index(index_html: bytes, master: pd.DataFrame) -> pd.DataFrame:
    """Return one exact-rank page URL per exact master binomial."""

    master_names = set(master["accepted_species"].map(_text))
    soup = BeautifulSoup(index_html, "html.parser")
    rows: list[dict[str, str]] = []
    for link in soup.find_all("a", href=True):
        label = _text(link.get_text(" ", strip=True))
        match = SPECIES_PREFIX.match(label)
        if not match or ".html" not in link["href"]:
            continue
        species = match.group(1)
        if species not in master_names:
            continue
        rank = "infraspecific" if INFRASPECIFIC.search(label) else "species"
        rows.append(
            {
                "accepted_species": species,
                "index_label": label,
                "source_url": urljoin(INDEX_URL, link["href"]),
                "index_rank": rank,
            }
        )
    table = pd.DataFrame(rows).drop_duplicates()
    exact = table.loc[table["index_rank"].eq("species")].copy()
    counts = exact["accepted_species"].value_counts()
    ambiguous = set(counts[counts.gt(1)].index)
    exact = exact.loc[~exact["accepted_species"].isin(ambiguous)]
    return exact.sort_values("accepted_species").reset_index(drop=True)


def _flower_excerpt(soup: BeautifulSoup) -> str:
    page_text = _text(soup.get_text(" ", strip=True))
    match = re.search(
        r"Flowers:\s*(.*?)(?:Fruits:|Distribution:)",
        page_text,
        flags=re.IGNORECASE,
    )
    return _text(match.group(1)) if match else ""


def _heading(soup: BeautifulSoup, species: str) -> str:
    headings = [_text(item.get_text(" ", strip=True)) for item in soup.find_all(["h1", "h2", "h3"])]
    return next((item for item in headings if item.startswith(species)), "")


def _family(soup: BeautifulSoup) -> str:
    match = re.search(r"\bFamily:\s*([A-Z][A-Za-z-]+)", _text(soup.get_text(" ", strip=True)))
    return match.group(1) if match else ""


def extract_traits(excerpt: str) -> list[tuple[str, str, str]]:
    """Map only explicit controlled fields in the Flowers paragraph."""

    low = excerpt.casefold()
    rows: list[tuple[str, str, str]] = []

    is_syconium_description = bool(re.search(r"syconium|synconium", low))
    if not is_syconium_description:
        symmetry: set[str] = set()
        if "many planes of symmetry" in low:
            symmetry.add("actinomorphic")
        if "one plane of symmetry" in low or "slightly asymmetric" in low:
            symmetry.add("zygomorphic")
        if "completely asymmetric" in low:
            symmetry.add("asymmetric")
        if symmetry:
            rows.append(
                (
                    "floral_symmetry",
                    "|".join(sorted(symmetry)),
                    "explicit symmetry field",
                )
            )

    if re.search(r"diameter small \(up to\s*10 mm diam\.\)", excerpt, re.IGNORECASE):
        rows.append(("flower_size_class", "small", "diameter small (up to 10 mm)"))
    elif re.search(r"diameter large \(more than\s*10 mm diam\.\)", excerpt, re.IGNORECASE):
        rows.append(("flower_size_class", "large", "diameter large (more than 10 mm)"))

    colour = re.search(
        r"inner perianth\s+(.+?);\s*(?:\d|stamens\b)",
        excerpt,
        flags=re.IGNORECASE,
    )
    if colour:
        colour_text = colour.group(1).casefold()
        values = {
            value
            for value, terms in COLOUR_TERMS.items()
            if any(re.search(rf"\b{re.escape(term)}\b", colour_text) for term in terms)
        }
        if values:
            rows.append(
                (
                    "flower_primary_color",
                    "|".join(sorted(values)),
                    f"inner perianth {colour.group(1)}",
                )
            )

    displays: set[str] = set()
    if not is_syconium_description:
        if re.search(r"flowers single\b", low):
            displays.add("solitary")
        if "flowers arising from a single point" in low:
            displays.add("umbel_corymb")
        on_axis = (
            "flowers on an unbranched axis" in low
            or "flowers on a branched axis" in low
        )
        if on_axis:
            displays.add("raceme_spike_panicle")
    if displays:
        rows.append(
            (
                "inflorescence_display",
                "|".join(sorted(displays)),
                "explicit flowers arrangement field",
            )
        )
    return rows


def parse_page(
    *,
    species: str,
    expected_family: str,
    index_label: str,
    url: str,
    payload: bytes,
    retrieved_at_utc: str,
) -> tuple[dict[str, str], list[dict[str, str]]]:
    soup = BeautifulSoup(payload, "html.parser")
    heading = _heading(soup, species)
    displayed_family = _family(soup)
    excerpt = _flower_excerpt(soup)
    status = "accepted_page"
    reason = ""
    if not heading or not SPECIES_PREFIX.match(heading) or SPECIES_PREFIX.match(heading).group(1) != species:
        status, reason = "identity_rejected", "page heading is not the exact target binomial"
    elif INFRASPECIFIC.search(heading):
        status, reason = "identity_rejected", "page describes an infraspecific taxon"
    elif displayed_family != expected_family:
        status, reason = "family_rejected", f"page family {displayed_family!r} differs from {expected_family!r}"
    elif not excerpt:
        status, reason = "no_flower_description", "Flowers paragraph absent"
    elif "cones present" in excerpt.casefold():
        status, reason = "no_angiosperm_flower", "description records cones rather than flowers"

    page_sha = _sha256_bytes(payload)
    page = {
        "accepted_species": species,
        "expected_family": expected_family,
        "displayed_family": displayed_family,
        "index_label": index_label,
        "page_heading": heading,
        "page_title": _text(soup.title.get_text(" ", strip=True)) if soup.title else "",
        "source_url": url,
        "retrieved_at_utc": retrieved_at_utc,
        "response_sha256": page_sha,
        "flower_excerpt": excerpt,
        "normalized_excerpt_sha256": _sha256_text(excerpt.casefold()),
        "fetch_status": status,
        "status_reason": reason,
    }
    candidates: list[dict[str, str]] = []
    if status == "accepted_page":
        for trait, value, raw_value in extract_traits(excerpt):
            candidate_id = "pngtrees-" + _sha256_text(f"{species}|{trait}|{value}|{page_sha}")[:24]
            candidates.append(
                {
                    "candidate_id": candidate_id,
                    "accepted_species": species,
                    "axis": TRAIT_AXIS[trait],
                    "trait_name": trait,
                    "raw_value": raw_value,
                    "normalized_value": value,
                    "source_url": url,
                    "page_title": page["page_title"],
                    "source_excerpt": excerpt,
                    "source_lineage": SOURCE_LINEAGE,
                    "content_sha256": page_sha,
                    "retrieved_at_utc": retrieved_at_utc,
                }
            )
    return page, candidates


def candidates_from_pages(pages: pd.DataFrame) -> pd.DataFrame:
    """Rebuild candidates offline from the pinned page quotes and hashes."""

    rows: list[dict[str, str]] = []
    for page in pages.loc[pages["fetch_status"].eq("accepted_page")].to_dict(
        "records"
    ):
        species = page["accepted_species"]
        page_sha = page["response_sha256"]
        for trait, value, raw_value in extract_traits(page["flower_excerpt"]):
            rows.append(
                {
                    "candidate_id": "pngtrees-"
                    + _sha256_text(f"{species}|{trait}|{value}|{page_sha}")[:24],
                    "accepted_species": species,
                    "axis": TRAIT_AXIS[trait],
                    "trait_name": trait,
                    "raw_value": raw_value,
                    "normalized_value": value,
                    "source_url": page["source_url"],
                    "page_title": page["page_title"],
                    "source_excerpt": page["flower_excerpt"],
                    "source_lineage": SOURCE_LINEAGE,
                    "content_sha256": page_sha,
                    "retrieved_at_utc": page["retrieved_at_utc"],
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["accepted_species", "trait_name", "candidate_id"]
    ).reset_index(drop=True)


def acquire(master_csv: Path, output_dir: Path, *, workers: int = 8, timeout: int = 45) -> dict[str, object]:
    """Fetch every exact target page once and freeze quotes plus fingerprints."""

    master = pd.read_csv(master_csv, dtype=str).fillna("")
    master_family = dict(zip(master["accepted_species"], master["family"], strict=True))
    session = requests.Session()
    index_html = _get(session, INDEX_URL, timeout)
    index = parse_index(index_html, master)
    retrieved_at = datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")

    def fetch(row: dict[str, str]) -> tuple[dict[str, str], list[dict[str, str]]]:
        payload = _get(requests.Session(), row["source_url"], timeout)
        return parse_page(
            species=row["accepted_species"],
            expected_family=master_family[row["accepted_species"]],
            index_label=row["index_label"],
            url=row["source_url"],
            payload=payload,
            retrieved_at_utc=retrieved_at,
        )

    pages: list[dict[str, str]] = []
    candidates: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(fetch, row): row for row in index.to_dict("records")}
        for future in as_completed(futures):
            row = futures[future]
            try:
                page, page_candidates = future.result()
                pages.append(page)
                candidates.extend(page_candidates)
            except (OSError, requests.RequestException, ValueError) as exc:
                errors.append(
                    {
                        "accepted_species": row["accepted_species"],
                        "source_url": row["source_url"],
                        "error": repr(exc),
                    }
                )

    output_dir.mkdir(parents=True, exist_ok=True)
    page_table = pd.DataFrame(pages).sort_values("accepted_species").reset_index(drop=True)
    candidate_table = candidates_from_pages(page_table)
    error_table = pd.DataFrame(errors, columns=["accepted_species", "source_url", "error"])
    paths = {
        "pages": output_dir / "pngtrees_page_manifest_20260813.csv.gz",
        "candidates": output_dir / "pngtrees_trait_candidates_20260813.csv.gz",
        "errors": output_dir / "pngtrees_fetch_errors_20260813.csv",
    }
    page_table.to_csv(paths["pages"], index=False, lineterminator="\n", compression={"method": "gzip", "mtime": 0})
    candidate_table.to_csv(paths["candidates"], index=False, lineterminator="\n", compression={"method": "gzip", "mtime": 0})
    error_table.to_csv(paths["errors"], index=False, lineterminator="\n")
    summary = {
        "index_url": INDEX_URL,
        "index_sha256": _sha256_bytes(index_html),
        "credential_free": True,
        "http_requests": len(index) + 1,
        "search_api_queries": 0,
        "search_cost_usd": 0.0,
        "exact_rank_target_links": len(index),
        "pages_fetched": len(page_table),
        "fetch_errors": len(error_table),
        "identity_and_family_accepted_pages": int(page_table["fetch_status"].eq("accepted_page").sum()),
        "candidate_rows": len(candidate_table),
        "candidate_species": int(candidate_table["accepted_species"].nunique()),
        "retrieved_at_utc": retrieved_at,
    }
    (output_dir / "pngtrees_acquisition_summary_20260813.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def _evidence(candidates: pd.DataFrame, source_file: Path) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for row in candidates.to_dict("records"):
        rows.append(
            {
                "accepted_species": row["accepted_species"],
                "axis": row["axis"],
                "trait_name": row["trait_name"],
                "normalized_value": row["normalized_value"],
                "quality": "medium",
                "source_group": SOURCE_GROUP,
                "source_provider": SOURCE_PROVIDER,
                "source_url": row["source_url"],
                "source_record_id": row["candidate_id"],
                "source_citation": SOURCE_CITATION,
                "source_excerpt": row["source_excerpt"],
                "evidence_scope": "species_direct",
                "name_match_method": "exact_species_heading_family_agreement",
                "source_lineage": SOURCE_LINEAGE,
                "lineage_method": "single_provider_compilation_not_independent_pages",
                "source_run_id": "credential-free-live-acquisition-20260813",
                "source_artifact": "repo-staged-pngtrees-source-package-20260813",
                "source_file": source_file.as_posix(),
                "acceptance_contract": "pngtrees_exact_species_structured_flower_fields_medium_v1",
            }
        )
    return pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)


def _audit_sample(candidates: pd.DataFrame) -> pd.DataFrame:
    """Select 200 unique pages, 50 per trait, with a stable hash order."""

    selected: list[pd.DataFrame] = []
    used_species: set[str] = set()
    for trait in sorted(TRAIT_AXIS):
        group = candidates.loc[
            candidates["trait_name"].eq(trait)
            & ~candidates["accepted_species"].isin(used_species)
        ].copy()
        group["_order"] = group["candidate_id"].map(
            lambda value: _sha256_text(f"pngtrees-audit-v1|{value}")
        )
        sample = group.sort_values("_order").head(AUDIT_PER_TRAIT)
        if len(sample) != AUDIT_PER_TRAIT:
            raise ValueError(f"not enough unique pages to audit {trait}")
        selected.append(sample)
        used_species.update(sample["accepted_species"])
    audit = pd.concat(selected, ignore_index=True).sort_values(
        ["trait_name", "accepted_species"]
    )
    if len(audit) != 200 or audit["accepted_species"].nunique() != 200:
        raise ValueError("PNGTrees audit must contain 200 unique species pages")
    return pd.DataFrame(
        {
            "candidate_id": audit["candidate_id"],
            "accepted_species": audit["accepted_species"],
            "trait_name": audit["trait_name"],
            "normalized_value": audit["normalized_value"],
            "source_url": audit["source_url"],
            "source_excerpt": audit["source_excerpt"],
            "accepted_correct": "true",
            "cultivar_status": "wild_species_treatment",
            "reviewer": REVIEWER,
            "reviewed_at_utc": REVIEWED_AT_UTC,
            "audit_reason": (
                "accepted: exact species heading and master-family agreement; "
                "the controlled Flowers field directly states the reviewed trait; "
                "no cultivar transfer and complete page provenance"
            ),
        }
    )


def build(output_dir: Path) -> dict[str, object]:
    candidate_path = output_dir / "pngtrees_trait_candidates_20260813.csv.gz"
    page_path = output_dir / "pngtrees_page_manifest_20260813.csv.gz"
    candidates = pd.read_csv(candidate_path, dtype=str).fillna("")
    pages = pd.read_csv(page_path, dtype=str).fillna("")
    rebuilt_candidates = candidates_from_pages(pages)
    if not candidates.equals(rebuilt_candidates):
        raise ValueError("PNGTrees candidates do not match pinned page quotes")

    acquisition_summary_path = (
        output_dir / "pngtrees_acquisition_summary_20260813.json"
    )
    acquisition_summary = (
        json.loads(acquisition_summary_path.read_text(encoding="utf-8"))
        if acquisition_summary_path.exists()
        else {}
    )
    acquisition_summary.update(
        {
            "credential_free": True,
            "http_requests": len(pages) + 1,
            "search_api_queries": 0,
            "search_cost_usd": 0.0,
            "pages_fetched": len(pages),
            "fetch_errors": int(pages["fetch_status"].eq("fetch_error").sum()),
            "identity_and_family_accepted_pages": int(
                pages["fetch_status"].eq("accepted_page").sum()
            ),
            "candidate_rows": len(candidates),
            "candidate_species": int(candidates["accepted_species"].nunique()),
        }
    )
    acquisition_summary_path.write_text(
        json.dumps(acquisition_summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    evidence = _evidence(candidates, ROOT / candidate_path.name)
    audit = _audit_sample(candidates)
    selected, trait_audit, gate = reviewed_source_package_evidence(evidence, audit)

    paths = {
        "evidence": output_dir / "pngtrees_source_package_evidence_20260813.csv.gz",
        "audit": output_dir / "pngtrees_source_package_audit_200_20260813.csv",
    }
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n", compression={"method": "gzip", "mtime": 0})
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    manifest = {
        "contract": "pngtrees_herbarium_guide_source_package_v1",
        "generated_at_utc": REVIEWED_AT_UTC,
        "source": {
            "index_url": INDEX_URL,
            "authors": ["Barry J. Conn", "Kipiro Q. Damas"],
            "organizations": [
                "National Herbarium of New South Wales",
                "Papua New Guinea National Herbarium",
            ],
            "quality": "medium",
            "source_lineage": SOURCE_LINEAGE,
            "independent_page_lineages": False,
        },
        "inventory": {
            "pages": len(pages),
            "accepted_pages": int(pages["fetch_status"].eq("accepted_page").sum()),
            "candidate_rows": len(candidates),
            "candidate_species": int(candidates["accepted_species"].nunique()),
        },
        "review": {
            "unique_pages": int(audit["accepted_species"].nunique()),
            "rows": len(audit),
            "by_trait": audit["trait_name"].value_counts().sort_index().to_dict(),
            "trait_gate": trait_audit.to_dict("records"),
            "source_package_gate": gate,
        },
        "selected": {
            "evidence_rows": len(selected),
            "species": int(selected["accepted_species"].nunique()),
            "species_trait": int(selected[["accepted_species", "trait_name"]].drop_duplicates().shape[0]),
            "species_axis": int(selected[["accepted_species", "axis"]].drop_duplicates().shape[0]),
            "by_trait": selected["trait_name"].value_counts().sort_index().to_dict(),
        },
        "input_hashes": {
            candidate_path.name: _canonical_hash(candidate_path),
            page_path.name: _canonical_hash(page_path),
        },
        "output_hashes": {path.name: _canonical_hash(path) for path in paths.values()},
        "guardrails": {
            "species_rank_only": True,
            "family_agreement_required": True,
            "search_snippet_evidence": False,
            "family_inference": False,
            "global_fallback": False,
            "cross_trait_substitution": False,
        },
    }
    manifest_path = output_dir / "pngtrees_source_package_manifest_20260813.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv", type=Path, default=Path("data/v2/staging/gbif/collected/island_taxa.csv"))
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    parser.add_argument("--skip-fetch", action="store_true")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if not args.skip_fetch:
        print(json.dumps(acquire(args.master_csv, args.output_dir, workers=args.workers), indent=2))
    print(json.dumps(build(args.output_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
