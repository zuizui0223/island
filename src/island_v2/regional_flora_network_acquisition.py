"""Acquire unresolved floral traits from a reusable regional-flora network.

The Flora of Zimbabwe, Zambia, Mozambique, Malawi and Botswana sites expose
the same family-checklist and species-treatment interface.  This module uses
that interface as a provider family, not as site-specific scraping code.  It
first inventories family indexes, intersects exact binomials with the fixed
island master, and only then retrieves pages that can fill a currently missing
strict axis or provide a third species for a support=2 genus x trait rule.

Outputs are review candidates.  They are never promoted to the formal ledger
by this acquisition command alone.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
import urllib.parse
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
import requests
from bs4 import BeautifulSoup

from island_v2.open_web_evidence import Page
from island_v2.open_web_structured_traits import extract_botanical_description

USER_AGENT = "island-floral-trait-research/1.0 (reproducible academic acquisition)"
TARGET_TRAITS = (
    "flower_primary_color",
    "floral_form",
    "floral_symmetry",
    "flower_size_class",
    "inflorescence_display",
    "tube_depth_class",
)
TRAIT_AXIS = {
    "flower_primary_color": "flower_colour",
    "floral_form": "floral_structural_complexity",
    "floral_symmetry": "floral_structural_complexity",
    "flower_size_class": "floral_structural_complexity",
    "inflorescence_display": "floral_structural_complexity",
    "tube_depth_class": "floral_structural_complexity",
}


@dataclass(frozen=True)
class FloraSite:
    key: str
    provider: str
    origin: str
    domain: str
    country: str

    def family_url(self, family_id: int) -> str:
        return (
            f"{self.origin}/speciesdata/utilities/utility-display-checklist.php?"
            f"family_id={family_id}"
        )

    def species_url(self, species_id: str) -> str:
        return f"{self.origin}/speciesdata/species.php?species_id={species_id}"


SITES = {
    item.key: item
    for item in (
        FloraSite(
            "zambia",
            "Flora of Zambia",
            "https://www.zambiaflora.com",
            "zambiaflora.com",
            "Zambia",
        ),
        FloraSite(
            "mozambique",
            "Flora of Mozambique",
            "https://www.mozambiqueflora.com",
            "mozambiqueflora.com",
            "Mozambique",
        ),
        FloraSite(
            "malawi",
            "Flora of Malawi",
            "https://www.malawiflora.com",
            "malawiflora.com",
            "Malawi",
        ),
        FloraSite(
            "botswana",
            "Flora of Botswana",
            "https://www.botswanaflora.com",
            "botswanaflora.com",
            "Botswana",
        ),
    )
}


@dataclass(frozen=True)
class FetchResult:
    url: str
    status: int
    payload: bytes
    retrieved_at_utc: str
    error: str = ""


def _now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _text(value: object) -> str:
    return " ".join(str(value or "").replace("\xa0", " ").split())


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _fetch(url: str, *, timeout: float = 75, attempts: int = 3) -> FetchResult:
    last_error = ""
    for attempt in range(attempts):
        try:
            response = requests.get(
                url,
                headers={"User-Agent": USER_AGENT},
                timeout=timeout,
            )
            if response.status_code == 200:
                return FetchResult(url, 200, response.content, _now())
            last_error = f"http_status_{response.status_code}"
        except requests.RequestException as exc:
            last_error = f"{type(exc).__name__}:{_text(exc)}"
        if attempt + 1 < attempts:
            time.sleep(1.5 * (attempt + 1))
    return FetchResult(url, 0, b"", _now(), last_error or "fetch_failed")


def parse_family_index(
    payload: bytes,
    *,
    site: FloraSite,
    family_id: int,
) -> tuple[str, list[dict[str, str]]]:
    """Return the displayed family and exact species links from one index."""

    soup = BeautifulSoup(payload, "html.parser")
    heading = soup.find("h1")
    family = _text(heading.get_text(" ", strip=True) if heading else "")
    family = re.sub(r"^Checklist:\s*", "", family, flags=re.IGNORECASE)
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for link in soup.select('a[href*="species.php?species_id="]'):
        name = _text(link.get_text(" ", strip=True)).lstrip("*").strip()
        match = re.search(r"(?:\?|&)species_id=(\d+)", str(link.get("href", "")))
        if not name or not match or len(name.split()) != 2:
            continue
        species_id = match.group(1)
        key = (name, species_id)
        if key in seen:
            continue
        seen.add(key)
        rows.append(
            {
                "site_key": site.key,
                "source_provider": site.provider,
                "domain": site.domain,
                "country": site.country,
                "family_id": str(family_id),
                "source_family": family,
                "source_scientific_name": name,
                "species_id": species_id,
                "family_url": site.family_url(family_id),
                "source_url": site.species_url(species_id),
            }
        )
    return family, rows


def _species_subject(soup: BeautifulSoup) -> str:
    heading = soup.find("h1")
    if not heading:
        return ""
    clone = BeautifulSoup(str(heading), "html.parser")
    for author in clone.select(".author"):
        author.decompose()
    return _text(clone.get_text(" ", strip=True))


def _species_family(soup: BeautifulSoup) -> str:
    for link in soup.select('a[href*="family.php?family_id="]'):
        value = _text(link.get_text(" ", strip=True))
        if value:
            return value
    return ""


def _description(soup: BeautifulSoup) -> str:
    marker = soup.select_one('a[href="about.php#descr"]')
    row = marker.find_parent("tr") if marker else None
    if not row:
        return ""
    cells = row.find_all("td")
    return _text(cells[-1].get_text(" ", strip=True)) if len(cells) >= 2 else ""


def parse_species_page(
    payload: bytes,
    *,
    site: FloraSite,
    expected_species: str,
    expected_family: str,
    source_url: str,
    missing_traits: Iterable[str],
    retrieved_at_utc: str,
) -> tuple[dict[str, str], list[dict[str, str]]]:
    """Apply fail-closed identity gates and extract organ-scoped evidence."""

    page_hash = _sha256(payload)
    soup = BeautifulSoup(payload, "html.parser")
    title = _text(soup.title.get_text(" ", strip=True) if soup.title else "")
    subject = _species_subject(soup)
    source_family = _species_family(soup)
    description = _description(soup)
    audit = {
        "site_key": site.key,
        "source_provider": site.provider,
        "accepted_species": expected_species,
        "expected_family": expected_family,
        "matched_page_name": subject,
        "matched_page_family": source_family,
        "source_url": source_url,
        "page_title": title,
        "content_sha256": page_hash,
        "retrieved_at_utc": retrieved_at_utc,
        "status": "",
        "error": "",
    }
    if subject != expected_species:
        audit["status"] = "rejected_name_mismatch"
        return audit, []
    if not source_family or source_family != expected_family:
        audit["status"] = "rejected_family_mismatch"
        return audit, []
    if not description:
        audit["status"] = "no_description"
        return audit, []
    if re.search(
        r"\b(?:cultivar|cultivated hybrid|horticultural hybrid|cv\.)\b",
        description,
        re.IGNORECASE,
    ):
        audit["status"] = "rejected_cultivar_or_hybrid_description"
        return audit, []

    page = Page(
        requested_url=source_url,
        final_url=source_url,
        status_code=200,
        content_type="text/html",
        title=title,
        text=f"Description:\n{description}",
        language="en",
        retrieved_at_utc=retrieved_at_utc,
        content_sha256=page_hash,
    )
    rows: list[dict[str, str]] = []
    species_id = urllib.parse.parse_qs(urllib.parse.urlsplit(source_url).query).get(
        "species_id", [expected_species.replace(" ", "_")]
    )[0]
    for trait in sorted(set(missing_traits).intersection(TARGET_TRAITS)):
        extracted = extract_botanical_description(page, trait_name=trait)
        if not extracted:
            continue
        values = sorted({normalized for _, normalized, _ in extracted})
        excerpts = list(dict.fromkeys(excerpt for _, _, excerpt in extracted))
        raw_values = list(dict.fromkeys(raw for raw, _, _ in extracted))
        rows.append(
            {
                "accepted_species": expected_species,
                "family": expected_family,
                "axis": TRAIT_AXIS[trait],
                "trait_name": trait,
                "normalized_value": "|".join(values),
                "raw_value": " || ".join(raw_values),
                "evidence_quality": "medium",
                "source_provider": site.provider,
                "source_url": source_url,
                "page_title": title,
                "source_citation": (
                    f"{site.provider} curated regional-flora species treatment"
                ),
                "source_excerpt": " || ".join(excerpts),
                "source_record_id": f"{site.key}-flora:{species_id}:{trait}",
                "source_lineage": f"provider_treatment:{site.key}_flora:{species_id}",
                "lineage_method": "curated_regional_flora_species_treatment",
                "source_tier": "B",
                "source_type": "curated_regional_flora_species_treatment",
                "domain": site.domain,
                "language": "en",
                "wild_cultivated_cultivar_status": (
                    "regional_flora_species_treatment_not_cultivar_limited"
                ),
                "content_sha256": page_hash,
                "content_sha256_basis": "retrieved_species_treatment_html_bytes",
                "retrieved_at_utc": retrieved_at_utc,
                "query": (
                    "credential_free_regional_flora_family_index_exact_species_page"
                ),
                "review_status": "pending_individual_review",
            }
        )
    audit["status"] = "candidate_extracted" if rows else "description_no_safe_trait_statement"
    return audit, rows


def _read_optional_csv(path: Path, **kwargs: object) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, **kwargs).fillna("") if path.exists() else pd.DataFrame()


def _resume_tables(resume_dir: Path | None) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if resume_dir is None:
        return (pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame())
    return (
        _read_optional_csv(resume_dir / "family_task_audit.csv"),
        _read_optional_csv(resume_dir / "regional_flora_inventory.csv.gz"),
        _read_optional_csv(resume_dir / "page_task_audit.csv.gz"),
        _read_optional_csv(resume_dir / "regional_flora_candidates.csv.gz"),
    )


def acquire(
    *,
    master_csv: Path,
    strict_coverage_csv: Path,
    direct_ledger_csv: Path,
    rule_queue_csv: Path,
    output_dir: Path,
    site_keys: Iterable[str],
    family_id_start: int,
    family_id_end: int,
    max_species_pages: int,
    workers: int,
    resume_dir: Path | None = None,
) -> dict[str, object]:
    if family_id_start < 1 or family_id_end <= family_id_start:
        raise ValueError("family ID range must be a non-empty half-open interval starting at 1+")
    selected_sites = [SITES[key] for key in site_keys]
    if not selected_sites:
        raise ValueError("at least one regional flora site is required")
    if max_species_pages < 1:
        raise ValueError("max_species_pages must be positive")

    output_dir.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    if master["accepted_species"].duplicated().any():
        raise ValueError("fixed master accepted species must be unique")
    family_by_species = master.set_index("accepted_species")["family"].to_dict()
    master_names = set(family_by_species)

    strict = pd.read_csv(strict_coverage_csv, dtype=str).fillna("")
    if len(strict) != 318_885 or strict["accepted_species"].nunique() != 106_295:
        raise ValueError("strict artifact denominator is not 106,295 x 3")
    strict_names = set(strict["accepted_species"])
    missing_from_master = strict_names.difference(master_names)
    if missing_from_master:
        raise ValueError(
            f"strict denominator has {len(missing_from_master)} species absent from fixed master"
        )
    axis_quality = strict.pivot(index="accepted_species", columns="axis", values="quality")
    axis_quality = axis_quality.reindex(sorted(strict_names)).fillna("")

    direct = pd.read_csv(
        direct_ledger_csv,
        usecols=["accepted_species", "trait_name", "resolution_status"],
        dtype=str,
    ).fillna("")
    completed_pairs = set(
        direct.loc[
            direct["resolution_status"].eq("resolved"),
            ["accepted_species", "trait_name"],
        ].itertuples(index=False, name=None)
    )
    queue = pd.read_csv(rule_queue_csv, dtype=str).fillna("")
    support2 = queue.loc[
        queue["current_support"].eq("2")
        & queue["trait_name"].isin(TARGET_TRAITS)
        & pd.to_numeric(queue["current_dominance"], errors="coerce").eq(1.0)
    ].copy()
    support2_potential = {
        (row.genus, row.trait_name): int(float(row.potential_cells_unlocked))
        for row in support2.itertuples(index=False)
    }

    old_family_audit, old_inventory, old_page_audit, old_candidates = _resume_tables(resume_dir)
    completed_family_tasks = set()
    if not old_family_audit.empty:
        completed_family_tasks = set(
            old_family_audit.loc[
                old_family_audit["status"].eq("success"), ["site_key", "family_id"]
            ].itertuples(index=False, name=None)
        )
    family_tasks = [
        (site, family_id)
        for site in selected_sites
        for family_id in range(family_id_start, family_id_end)
        if (site.key, str(family_id)) not in completed_family_tasks
    ]

    family_audit_rows: list[dict[str, str]] = []
    inventory_rows: list[dict[str, str]] = []

    def family_task(site: FloraSite, family_id: int) -> tuple[dict[str, str], list[dict[str, str]]]:
        result = _fetch(site.family_url(family_id))
        audit = {
            "site_key": site.key,
            "family_id": str(family_id),
            "source_url": result.url,
            "http_status": str(result.status),
            "status": "",
            "displayed_family": "",
            "listed_species": "0",
            "content_sha256": _sha256(result.payload) if result.payload else "",
            "retrieved_at_utc": result.retrieved_at_utc,
            "error": result.error,
        }
        if result.status != 200:
            audit["status"] = "fetch_failed"
            return audit, []
        family, rows = parse_family_index(result.payload, site=site, family_id=family_id)
        audit["displayed_family"] = family
        audit["listed_species"] = str(len(rows))
        audit["status"] = "success"
        return audit, rows

    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        futures = {
            pool.submit(family_task, site, family_id): (site.key, family_id)
            for site, family_id in family_tasks
        }
        for future in as_completed(futures):
            audit, rows = future.result()
            family_audit_rows.append(audit)
            inventory_rows.extend(rows)

    family_audit = pd.concat(
        [old_family_audit, pd.DataFrame(family_audit_rows)], ignore_index=True
    ).fillna("")
    if not family_audit.empty:
        family_audit = family_audit.drop_duplicates(["site_key", "family_id"], keep="last")
    inventory = pd.concat(
        [old_inventory, pd.DataFrame(inventory_rows)], ignore_index=True
    ).fillna("")
    if inventory.empty:
        inventory = pd.DataFrame(
            columns=[
                "site_key", "source_provider", "domain", "country", "family_id",
                "source_family", "source_scientific_name", "species_id", "family_url",
                "source_url",
            ]
        )
    inventory = inventory.drop_duplicates(["site_key", "species_id"], keep="last")
    inventory["accepted_species"] = inventory["source_scientific_name"]
    inventory["identity_status"] = "name_not_in_fixed_master"
    exact = inventory["accepted_species"].isin(strict_names)
    inventory.loc[exact, "expected_family"] = inventory.loc[exact, "accepted_species"].map(
        family_by_species
    )
    inventory.loc[exact, "identity_status"] = "accepted_name_exact_pending_page_family"

    completed_urls = set()
    if not old_page_audit.empty:
        completed_urls = set(
            old_page_audit.loc[
                ~old_page_audit["status"].isin({"fetch_failed", ""}), "source_url"
            ]
        )
    candidates = inventory.loc[exact & ~inventory["source_url"].isin(completed_urls)].copy()

    def missing_traits(species: str) -> list[str]:
        return [trait for trait in TARGET_TRAITS if (species, trait) not in completed_pairs]

    candidates["missing_traits"] = candidates["accepted_species"].map(
        lambda species: "|".join(missing_traits(species))
    )
    candidates = candidates.loc[candidates["missing_traits"].ne("")].copy()
    candidates["genus"] = candidates["accepted_species"].str.split().str[0]
    candidates["zero_axis"] = candidates["accepted_species"].map(
        lambda species: int((axis_quality.loc[species] == "").all())
    )
    candidates["missing_strict_axes"] = candidates["accepted_species"].map(
        lambda species: int(axis_quality.at[species, "flower_colour"] == "")
        + int(axis_quality.at[species, "floral_structural_complexity"] == "")
    )
    candidates["rule_unlock_potential"] = candidates.apply(
        lambda row: max(
            [
                support2_potential.get((row["genus"], trait), 0)
                for trait in row["missing_traits"].split("|")
            ]
            or [0]
        ),
        axis=1,
    )
    island_counts = pd.to_numeric(master.set_index("accepted_species")["n_islands"], errors="coerce")
    record_counts = pd.to_numeric(master.set_index("accepted_species")["n_records"], errors="coerce")
    candidates["n_islands"] = candidates["accepted_species"].map(island_counts).fillna(0)
    candidates["n_records"] = candidates["accepted_species"].map(record_counts).fillna(0)
    candidates = candidates.loc[
        candidates["missing_strict_axes"].gt(0) | candidates["rule_unlock_potential"].gt(0)
    ]
    candidates = candidates.sort_values(
        [
            "zero_axis", "missing_strict_axes", "rule_unlock_potential", "n_islands",
            "n_records", "accepted_species", "site_key",
        ],
        ascending=[False, False, False, False, False, True, True],
        kind="stable",
    )
    # One provider page per species is the efficient first pass.  A later
    # recovery pass can deliberately seek an independent lineage for no-hits.
    selected = candidates.drop_duplicates("accepted_species", keep="first").head(
        max_species_pages
    )

    page_audit_rows: list[dict[str, str]] = []
    candidate_rows: list[dict[str, str]] = []

    def page_task(row: object) -> tuple[dict[str, str], list[dict[str, str]]]:
        item = row._asdict()
        site = SITES[str(item["site_key"])]
        result = _fetch(str(item["source_url"]))
        if result.status != 200:
            return (
                {
                    "site_key": site.key,
                    "source_provider": site.provider,
                    "accepted_species": str(item["accepted_species"]),
                    "expected_family": str(item["expected_family"]),
                    "matched_page_name": "",
                    "matched_page_family": "",
                    "source_url": str(item["source_url"]),
                    "page_title": "",
                    "content_sha256": "",
                    "retrieved_at_utc": result.retrieved_at_utc,
                    "status": "fetch_failed",
                    "error": result.error,
                },
                [],
            )
        return parse_species_page(
            result.payload,
            site=site,
            expected_species=str(item["accepted_species"]),
            expected_family=str(item["expected_family"]),
            source_url=str(item["source_url"]),
            missing_traits=str(item["missing_traits"]).split("|"),
            retrieved_at_utc=result.retrieved_at_utc,
        )

    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        futures = {pool.submit(page_task, row): row.source_url for row in selected.itertuples()}
        for future in as_completed(futures):
            audit, rows = future.result()
            page_audit_rows.append(audit)
            candidate_rows.extend(rows)

    page_audit = pd.concat(
        [old_page_audit, pd.DataFrame(page_audit_rows)], ignore_index=True
    ).fillna("")
    if not page_audit.empty:
        page_audit = page_audit.drop_duplicates("source_url", keep="last")
    evidence = pd.concat(
        [old_candidates, pd.DataFrame(candidate_rows)], ignore_index=True
    ).fillna("")
    if not evidence.empty:
        evidence = evidence.drop_duplicates(
            ["accepted_species", "trait_name", "source_lineage"], keep="last"
        ).sort_values(["source_provider", "accepted_species", "trait_name"])

    paths = {
        "family_audit": output_dir / "family_task_audit.csv",
        "inventory": output_dir / "regional_flora_inventory.csv.gz",
        "page_audit": output_dir / "page_task_audit.csv.gz",
        "candidates": output_dir / "regional_flora_candidates.csv.gz",
    }
    family_audit.to_csv(paths["family_audit"], index=False, lineterminator="\n")
    inventory.to_csv(paths["inventory"], index=False, lineterminator="\n", compression="gzip")
    page_audit.to_csv(paths["page_audit"], index=False, lineterminator="\n", compression="gzip")
    evidence.to_csv(paths["candidates"], index=False, lineterminator="\n", compression="gzip")

    manifest: dict[str, object] = {
        "contract": "regional_flora_network_acquisition_v1",
        "created_at": _now(),
        "sites": [asdict(site) for site in selected_sites],
        "family_id_range": {"start": family_id_start, "end_exclusive": family_id_end},
        "fixed_denominator": {"species": 106_295, "species_axis": 318_885},
        "resume_dir": str(resume_dir or ""),
        "family_tasks": family_audit["status"].value_counts().sort_index().to_dict(),
        "inventory_rows": len(inventory),
        "fixed_master_exact_inventory_rows": int(exact.sum()),
        "selected_species_pages": len(selected),
        "page_tasks": page_audit["status"].value_counts().sort_index().to_dict(),
        "candidate_rows": len(evidence),
        "candidate_species": int(evidence["accepted_species"].nunique()) if not evidence.empty else 0,
        "by_provider": (
            evidence["source_provider"].value_counts().sort_index().to_dict()
            if not evidence.empty else {}
        ),
        "by_trait": (
            evidence["trait_name"].value_counts().sort_index().to_dict()
            if not evidence.empty else {}
        ),
        "search_api_queries": 0,
        "search_api_cost": 0,
        "promotion_status": "not_promoted_pending_individual_review",
        "guardrails": {
            "accepted_name_exact": True,
            "page_family_exact": True,
            "original_page_required": True,
            "exact_excerpt_required": True,
            "species_direct_only": True,
            "family_inference": False,
            "global_fallback": False,
            "cross_trait_substitution": False,
        },
        "inputs": {
            "master_csv": str(master_csv),
            "strict_coverage_csv": str(strict_coverage_csv),
            "direct_ledger_csv": str(direct_ledger_csv),
            "rule_queue_csv": str(rule_queue_csv),
        },
        "outputs": {
            label: {"path": str(path), "sha256": _sha256(path.read_bytes())}
            for label, path in paths.items()
        },
    }
    (output_dir / "regional_flora_acquisition_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--strict-coverage-csv", type=Path, required=True)
    parser.add_argument("--direct-ledger-csv", type=Path, required=True)
    parser.add_argument("--rule-queue-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sites", default=",".join(SITES))
    parser.add_argument("--family-id-start", type=int, default=1)
    parser.add_argument("--family-id-end", type=int, default=320)
    parser.add_argument("--max-species-pages", type=int, default=5000)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--resume-dir", type=Path)
    args = parser.parse_args()
    site_keys = [item.strip() for item in args.sites.split(",") if item.strip()]
    unknown = sorted(set(site_keys).difference(SITES))
    if unknown:
        raise ValueError(f"unknown regional flora sites: {unknown}")
    report = acquire(
        master_csv=args.master_csv,
        strict_coverage_csv=args.strict_coverage_csv,
        direct_ledger_csv=args.direct_ledger_csv,
        rule_queue_csv=args.rule_queue_csv,
        output_dir=args.output_dir,
        site_keys=site_keys,
        family_id_start=args.family_id_start,
        family_id_end=args.family_id_end,
        max_species_pages=args.max_species_pages,
        workers=args.workers,
        resume_dir=args.resume_dir,
    )
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
