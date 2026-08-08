"""Backfill an approved botanical site's species inventory.

The general-search lane discovers useful regional floras and specialist plant
sites.  Once a domain has been manually registered with an inventory URL and a
species-page pattern, this module performs a bounded, generic site inventory
backfill.  It is not specific to Go Botany or Flora of Mikawa.

All output remains review-only.  No candidate is added to the evidence ledger
until a separate manual audit and promotion step accepts it.
"""

from __future__ import annotations

import json
import re
import urllib.parse
from collections import deque
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.open_web_evidence import (
    DomainRecord,
    DomainRegistry,
    PageFetcher,
    canonical_url,
    load_yaml,
)
from island_v2.open_web_network_pilot import _load_baseline, _now, audit_sample
from island_v2.open_web_network_structured_pilot import _fetch_all_strict_traits

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT = "open_web_discovered_domain_inventory_v1"


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).replace("\xa0", " ").split())


def _species_keys_from_url(url: str) -> set[str]:
    """Return conservative binomial candidates encoded in a treatment URL."""

    parsed = urllib.parse.urlsplit(canonical_url(url))
    keys: set[str] = set()
    query = urllib.parse.parse_qs(parsed.query)
    for parameter in ("name", "species", "taxon"):
        for raw_value in query.get(parameter, []):
            value = _text(urllib.parse.unquote_plus(raw_value).replace("~", " "))
            tokens = value.split()
            if len(tokens) == 2:
                keys.add(value.casefold())
    if keys:
        return keys

    parts = [
        urllib.parse.unquote(part).replace("_", " ")
        for part in parsed.path.strip("/").split("/")
        if part
    ]
    if not parts:
        return keys
    slug_tokens = [token for token in parts[-1].split("-") if token]
    if len(slug_tokens) == 2:
        keys.add(" ".join(slug_tokens).casefold())
    if len(parts) < 2:
        return keys
    genus = _text(parts[-2])
    epithet = _text(parts[-1])
    if " " not in genus and " " not in epithet:
        keys.add(f"{genus} {epithet}".casefold())
    return keys


def discover_species_urls(
    record: DomainRecord,
    *,
    fetcher: PageFetcher,
    max_species_pages: int,
    completed_species_urls: set[str] | None = None,
    target_species_names: set[str] | None = None,
) -> tuple[list[str], list[dict[str, object]]]:
    """Discover species-page URLs without fetching each treatment twice."""

    completed = {
        canonical_url(url)
        for url in (completed_species_urls or set())
        if canonical_url(url)
    }
    targets = {
        _text(name).casefold() for name in (target_species_names or set()) if _text(name)
    }
    pending = deque((canonical_url(url), 0) for url in record.inventory_urls)
    seen_indexes: set[str] = set()
    species_urls: list[str] = []
    species_seen: set[str] = set()
    audit: list[dict[str, object]] = []

    while pending and len(species_urls) < max_species_pages:
        url, depth = pending.popleft()
        if not url or url in seen_indexes:
            continue
        seen_indexes.add(url)
        page = fetcher.fetch(url)
        audit.append(
            {
                "url": url,
                "depth": depth,
                "status_code": page.status_code,
                "error": page.error,
                "text_length": len(page.text),
                "links": len(page.links),
            }
        )
        if page.error:
            continue
        page_links = list(page.links)
        if "xml" in page.content_type.casefold() or url.casefold().endswith(".xml"):
            page_links.extend(re.findall(r"https?://[^\s<>]+", page.text))
        for child in page_links:
            child = canonical_url(child)
            if not child:
                continue
            child_host = urllib.parse.urlsplit(child).netloc.casefold().removeprefix("www.")
            if child_host != record.domain.casefold():
                continue
            if re.search(record.page_pattern, child):
                if targets and not (_species_keys_from_url(child) & targets):
                    continue
                if child not in species_seen:
                    species_seen.add(child)
                    if child in completed:
                        continue
                    species_urls.append(child)
                    if len(species_urls) >= max_species_pages:
                        break
                continue
            if depth < record.inventory_depth and child not in seen_indexes:
                pending.append((child, depth + 1))
    return species_urls, audit


def run_inventory(
    *,
    baseline_dir: Path,
    master_csv: Path,
    config_yaml: Path,
    registry_yaml: Path,
    output_dir: Path,
    domain: str,
    max_species_pages: int,
    inventory_pause_seconds: float,
    page_pause_seconds: float,
    completed_species_pages_csv: Path | None = None,
    target_coverage_csv: Path | None = None,
) -> dict[str, Any]:
    coverage, evidence = _load_baseline(baseline_dir)
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    universe = set(coverage["accepted_species"])
    master = master.loc[master["accepted_species"].isin(universe)].drop_duplicates(
        "accepted_species"
    )

    registry = DomainRegistry.load(registry_yaml)
    record = next((item for item in registry.all() if item.domain == domain), None)
    if record is None:
        raise ValueError(f"domain is not registered: {domain}")
    if record.review_status != "approved":
        raise ValueError(f"domain is not approved: {domain}")
    if not record.inventory_urls:
        raise ValueError(f"domain has no inventory URL: {domain}")

    inventory_fetcher = PageFetcher(pause_seconds=inventory_pause_seconds)
    completed_species_urls: set[str] = set()
    if completed_species_pages_csv is not None:
        completed_pages = pd.read_csv(
            completed_species_pages_csv,
            usecols=["source_url"],
            dtype=str,
        ).fillna("")
        completed_species_urls = {
            canonical_url(url)
            for url in completed_pages["source_url"]
            if canonical_url(url)
        }
    target_species_names: set[str] = set()
    if target_coverage_csv is not None:
        target_coverage = pd.read_csv(
            target_coverage_csv,
            usecols=["accepted_species", "quality"],
            dtype=str,
        ).fillna("")
        target_species_names = set(
            target_coverage.loc[
                target_coverage["quality"].eq(""), "accepted_species"
            ]
        )
        if not target_species_names:
            raise ValueError("target coverage contains no unresolved species")
    species_urls, inventory_audit = discover_species_urls(
        record,
        fetcher=inventory_fetcher,
        max_species_pages=max_species_pages,
        completed_species_urls=completed_species_urls,
        target_species_names=target_species_names,
    )
    leads = [
        {
            "accepted_species": "",
            "trait_name": "",
            "query": f"inventory:{domain}",
            "provider": "approved_discovered_domain_inventory",
            "rank": index + 1,
            "result_url": url,
            "result_domain": domain,
            "title": "",
            "snippet_discovery_only": "",
        }
        for index, url in enumerate(species_urls)
    ]
    candidates, page_audit, domain_audit = _fetch_all_strict_traits(
        leads,
        master=master,
        extraction_config=load_yaml(config_yaml),
        registry=registry,
        max_pages=max_species_pages,
        page_pause_seconds=page_pause_seconds,
    )

    candidate_frame = pd.DataFrame(candidates)
    baseline_direct = evidence.loc[evidence["quality"].isin({"high", "medium"})]
    baseline_pairs = set(
        zip(baseline_direct["accepted_species"], baseline_direct["trait_name"])
    )
    if not candidate_frame.empty:
        candidate_frame["already_direct_in_baseline"] = [
            (species, trait) in baseline_pairs
            for species, trait in zip(
                candidate_frame["accepted_species"], candidate_frame["trait_name"]
            )
        ]
        novel = candidate_frame.loc[~candidate_frame["already_direct_in_baseline"]].copy()
    else:
        novel = candidate_frame.copy()

    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(inventory_audit).to_csv(
        output_dir / "inventory_index_audit.csv", index=False
    )
    pd.DataFrame(page_audit).to_csv(output_dir / "species_page_audit.csv", index=False)
    pd.DataFrame(domain_audit).to_csv(output_dir / "domain_outcome_audit.csv", index=False)
    candidate_frame.to_csv(output_dir / "all_inventory_candidates.csv", index=False)
    novel.to_csv(output_dir / "novel_inventory_candidates.csv", index=False)
    audit_sample(novel, limit=100).to_csv(
        output_dir / "manual_audit_sample.csv", index=False
    )
    pd.DataFrame({"source_url": species_urls}).to_csv(
        output_dir / "discovered_species_pages.csv", index=False
    )

    report: dict[str, Any] = {
        "contract": CONTRACT,
        "completed_at_utc": _now(),
        "domain": domain,
        "source_name": record.source_name,
        "source_type": record.source_type,
        "source_tier": record.recommended_evidence_tier,
        "inventory_urls": list(record.inventory_urls),
        "inventory_index_pages_fetched": len(inventory_audit),
        "species_pages_discovered": len(species_urls),
        "species_pages_attempted": len(page_audit),
        "completed_species_pages_excluded": len(completed_species_urls),
        "completed_species_pages_csv": (
            str(completed_species_pages_csv)
            if completed_species_pages_csv is not None
            else ""
        ),
        "target_coverage_csv": (
            str(target_coverage_csv) if target_coverage_csv is not None else ""
        ),
        "target_unresolved_species": len(target_species_names),
        "inventory_targeting": (
            "unresolved_species_url_key" if target_species_names else "unfiltered"
        ),
        "identity_passed_pages": sum(
            row.get("page_status") == "identity_passed" for row in page_audit
        ),
        "all_candidate_rows": len(candidate_frame),
        "all_candidate_species_trait": (
            len(candidate_frame.drop_duplicates(["accepted_species", "trait_name"]))
            if not candidate_frame.empty
            else 0
        ),
        "novel_candidate_rows": len(novel),
        "novel_candidate_species": (
            int(novel["accepted_species"].nunique()) if not novel.empty else 0
        ),
        "novel_candidate_species_trait": (
            len(novel.drop_duplicates(["accepted_species", "trait_name"]))
            if not novel.empty
            else 0
        ),
        "novel_candidate_species_axis": (
            len(novel.drop_duplicates(["accepted_species", "axis"]))
            if not novel.empty
            else 0
        ),
        "novel_by_trait": (
            {
                _text(trait): int(count)
                for trait, count in novel.drop_duplicates(
                    ["accepted_species", "trait_name"]
                )["trait_name"].value_counts().items()
            }
            if not novel.empty
            else {}
        ),
        "accepted_evidence_rows": 0,
        "manual_review_required": True,
        "automatic_promotion": False,
        "family_inference": False,
        "global_fallback": False,
    }
    (output_dir / "inventory_pilot_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False, sort_keys=True))
    return report


@app.command()
def main(
    baseline_dir: Path = typer.Option(..., exists=True, file_okay=False),
    master_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(..., file_okay=False),
    domain: str = typer.Option(...),
    config_yaml: Path = typer.Option(
        Path("config/open_web_trait_acquisition.yml"), exists=True, dir_okay=False
    ),
    registry_yaml: Path = typer.Option(
        Path("config/open_web_domain_registry.yml"), exists=True, dir_okay=False
    ),
    max_species_pages: int = typer.Option(250, min=1, max=5000),
    inventory_pause_seconds: float = typer.Option(0.5, min=0.1, max=10),
    page_pause_seconds: float = typer.Option(0.35, min=0.1, max=10),
    completed_species_pages_csv: Path | None = typer.Option(
        None,
        exists=True,
        dir_okay=False,
        help="Prior discovered_species_pages.csv; listed URLs are never fetched again.",
    ),
    target_coverage_csv: Path | None = typer.Option(
        None,
        exists=True,
        dir_okay=False,
        help="Latest coverage ledger; only treatment URLs for unresolved species are fetched.",
    ),
) -> None:
    run_inventory(
        baseline_dir=baseline_dir,
        master_csv=master_csv,
        config_yaml=config_yaml,
        registry_yaml=registry_yaml,
        output_dir=output_dir,
        domain=domain,
        max_species_pages=max_species_pages,
        inventory_pause_seconds=inventory_pause_seconds,
        page_pause_seconds=page_pause_seconds,
        completed_species_pages_csv=completed_species_pages_csv,
        target_coverage_csv=target_coverage_csv,
    )


if __name__ == "__main__":
    app()
