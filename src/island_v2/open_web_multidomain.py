"""Credential-free multi-domain pilot for open-Web plant trait evidence.

The official search adapters remain the production discovery path. This module
adds a reproducible way to crawl several independently reviewed botanical-site
inventories in one run, so the pilot is not tied to any one example domain.
Every source page still passes the existing registry, identity, cultivar,
provenance, ontology, and source-lineage gates before it becomes a candidate.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
import typer

from island_v2.open_web_evidence import (
    DomainRecord,
    DomainRegistry,
    Page,
    PageFetcher,
    StrictNameResolver,
    crawl_registry_inventory,
    load_yaml,
    page_scientific_names,
    process_pages,
)
from island_v2.open_web_pilot import (
    DIRECT_QUALITY,
    TRAIT_AXIS,
    _read_baseline,
    _text,
    _validate_denominator,
)

CONTRACT = "open_web_multidomain_pilot_v1"
app = typer.Typer(add_completion=False, no_args_is_help=True)


def _now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _write_rows(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    columns = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _resolved_subject(
    page: Page,
    resolver: StrictNameResolver,
) -> tuple[str, str]:
    """Resolve one focal accepted species for an inventory-discovered page.

    Prefer a scientific binomial in the title. If the title uses only a local
    name, accept the first unique exact master name in the opening treatment.
    The downstream identity gate rechecks the page against this expected taxon.
    """

    title_names = page_scientific_names(page, title_only=True)
    body_names = page_scientific_names(page)
    ordered = list(dict.fromkeys([*title_names, *body_names]))
    accepted: list[tuple[str, str]] = []
    for name in ordered:
        resolution = resolver.resolve(name)
        if resolution.status != "accepted":
            continue
        pair = (resolution.accepted_species, name)
        if pair not in accepted:
            accepted.append(pair)
    if not accepted:
        return "", ""
    if title_names:
        title_keys = {name.casefold() for name in title_names}
        for species, matched_name in accepted:
            if matched_name.casefold() in title_keys:
                return species, matched_name
    unique_species = list(dict.fromkeys(species for species, _ in accepted))
    if len(unique_species) == 1:
        species = unique_species[0]
        matched = next(name for value, name in accepted if value == species)
        return species, matched
    return "", ""


def _hash_sample(
    candidates: pd.DataFrame,
    *,
    sample_pages: int,
    seed: str,
) -> pd.DataFrame:
    """Select a deterministic trait-stratified sample of unique pages."""

    if candidates.empty:
        return pd.DataFrame()
    table = candidates.copy()
    stable = (
        table["source_url"].map(_text)
        + "|"
        + table["candidate_id"].map(_text)
        + "|"
        + seed
    )
    table["_audit_hash"] = stable.map(
        lambda value: hashlib.sha256(value.encode("utf-8")).hexdigest()
    )
    table = table.sort_values(
        ["trait_name", "_audit_hash", "source_url", "candidate_id"]
    )
    pools = {
        trait: group.to_dict("records")
        for trait, group in table.groupby("trait_name", sort=True)
    }
    positions = {trait: 0 for trait in pools}
    selected: list[dict[str, object]] = []
    seen_urls: set[str] = set()
    while len(selected) < sample_pages:
        progressed = False
        for trait in sorted(pools):
            rows = pools[trait]
            while positions[trait] < len(rows):
                row = rows[positions[trait]]
                positions[trait] += 1
                url = _text(row.get("source_url"))
                if not url or url in seen_urls:
                    continue
                selected.append(row)
                seen_urls.add(url)
                progressed = True
                break
            if len(selected) >= sample_pages:
                break
        if not progressed:
            break
    return pd.DataFrame(selected).drop(columns=["_audit_hash"], errors="ignore")


def _domain_rows(
    record: DomainRecord,
    *,
    inventory_audit: Sequence[Mapping[str, object]],
    fetched_pages: int,
    candidate_rows: int,
    candidate_pages: int,
) -> dict[str, object]:
    return {
        "domain": record.domain,
        "source_name": record.source_name,
        "source_type": record.source_type,
        "region": record.region,
        "language": record.language,
        "organization_type": record.organization_type,
        "source_tier": record.recommended_evidence_tier,
        "inventory_audit_rows": len(inventory_audit),
        "fetched_pages": fetched_pages,
        "candidate_rows": candidate_rows,
        "candidate_pages": candidate_pages,
        "candidate_page_rate": (
            candidate_pages / fetched_pages if fetched_pages else None
        ),
    }


def run_multidomain_pilot(
    *,
    baseline_dir: Path,
    master_csv: Path,
    registry_yaml: Path,
    config_yaml: Path,
    output_dir: Path,
    domains: Sequence[str],
    max_pages_per_domain: int,
    audit_pages: int,
    pause_seconds: float,
    audit_seed: str,
) -> dict[str, object]:
    summary, coverage, evidence, _ = _read_baseline(baseline_dir)
    _validate_denominator(summary, coverage)
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    universe = set(coverage["accepted_species"])
    master = master.loc[master["accepted_species"].isin(universe)].drop_duplicates(
        "accepted_species"
    )
    resolver = StrictNameResolver(master.to_dict("records"))
    registry = DomainRegistry.load(registry_yaml)
    records = {record.domain: record for record in registry.all()}
    requested = list(dict.fromkeys(domain.strip() for domain in domains if domain.strip()))
    missing = [domain for domain in requested if domain not in records]
    if missing:
        raise ValueError(f"unregistered domains: {missing}")
    selected_records = [records[domain] for domain in requested]
    if any(not record.inventory_urls for record in selected_records):
        absent = [record.domain for record in selected_records if not record.inventory_urls]
        raise ValueError(f"domains without inventory URLs: {absent}")

    fetcher = PageFetcher(pause_seconds=pause_seconds)
    all_candidates: list[dict[str, object]] = []
    all_duplicates: list[dict[str, object]] = []
    all_identity: list[dict[str, object]] = []
    all_fetch: list[dict[str, object]] = []
    all_domain_audit: list[dict[str, object]] = []
    all_inventory: list[dict[str, object]] = []
    domain_summary: list[dict[str, object]] = []
    config = load_yaml(config_yaml)

    for record in selected_records:
        urls, inventory_audit = crawl_registry_inventory(
            record,
            fetcher=fetcher,
            max_urls=max_pages_per_domain,
        )
        leads: list[dict[str, object]] = []
        prefetch_failures = 0
        unresolved_subjects = 0
        for rank, url in enumerate(urls, start=1):
            page = fetcher.fetch(url)
            if page.error or not page.text:
                prefetch_failures += 1
                continue
            accepted_species, searched_name = _resolved_subject(page, resolver)
            if not accepted_species:
                unresolved_subjects += 1
                continue
            leads.append(
                {
                    "result_url": url,
                    "provider": "approved_multidomain_inventory",
                    "query": f"registry:{record.domain}",
                    "rank": rank,
                    "searched_name": searched_name,
                    "accepted_species": accepted_species,
                    "trait_name": "",
                }
            )
        result = process_pages(
            leads,
            registry=registry,
            resolver=resolver,
            config=config,
            fetcher=fetcher,
            default_traits=tuple(TRAIT_AXIS),
        )
        candidates = list(result["candidates"])
        for row in candidates:
            row["inventory_domain"] = record.domain
        for row in result["lineage_duplicates"]:
            row["inventory_domain"] = record.domain
        for row in result["identity_audit"]:
            row["inventory_domain"] = record.domain
        for row in result["fetch_audit"]:
            row["inventory_domain"] = record.domain
        for row in result["domain_audit"]:
            row["inventory_domain"] = record.domain
        for row in inventory_audit:
            row = dict(row)
            row["inventory_domain"] = record.domain
            all_inventory.append(row)
        all_candidates.extend(candidates)
        all_duplicates.extend(result["lineage_duplicates"])
        all_identity.extend(result["identity_audit"])
        all_fetch.extend(result["fetch_audit"])
        all_domain_audit.extend(result["domain_audit"])
        candidate_pages = len(
            {_text(row.get("source_url")) for row in candidates if row.get("source_url")}
        )
        summary_row = _domain_rows(
            record,
            inventory_audit=inventory_audit,
            fetched_pages=len(result["fetch_audit"]),
            candidate_rows=len(candidates),
            candidate_pages=candidate_pages,
        )
        summary_row["inventory_urls_discovered"] = len(urls)
        summary_row["prefetch_failures"] = prefetch_failures
        summary_row["unresolved_page_subjects"] = unresolved_subjects
        domain_summary.append(summary_row)

    candidates = pd.DataFrame(all_candidates)
    direct_pairs = set(
        zip(
            evidence.loc[evidence["quality"].isin(DIRECT_QUALITY), "accepted_species"],
            evidence.loc[evidence["quality"].isin(DIRECT_QUALITY), "trait_name"],
        )
    )
    if not candidates.empty:
        candidates = candidates.loc[
            [
                (species, trait) not in direct_pairs
                for species, trait in zip(
                    candidates["accepted_species"], candidates["trait_name"]
                )
            ]
        ].copy()
        candidates["candidate_id"] = candidates.apply(
            lambda row: hashlib.sha256(
                "|".join(
                    [
                        _text(row.get("accepted_species")),
                        _text(row.get("trait_name")),
                        _text(row.get("normalized_value")),
                        _text(row.get("source_lineage")),
                        _text(row.get("supporting_excerpt")),
                    ]
                ).encode("utf-8")
            ).hexdigest()[:24],
            axis=1,
        )
        candidates = candidates.drop_duplicates(
            [
                "accepted_species",
                "trait_name",
                "normalized_value",
                "source_lineage",
            ]
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    candidate_rows = candidates.to_dict("records") if not candidates.empty else []
    _write_rows(output_dir / "open_web_multidomain_candidates.csv", candidate_rows)
    _write_rows(output_dir / "source_lineage_duplicates.csv", all_duplicates)
    _write_rows(output_dir / "species_identity_audit.csv", all_identity)
    _write_rows(output_dir / "page_fetch_audit.csv", all_fetch)
    _write_rows(output_dir / "domain_outcome_audit.csv", all_domain_audit)
    _write_rows(output_dir / "registry_inventory_audit.csv", all_inventory)
    _write_rows(output_dir / "domain_summary.csv", domain_summary)

    audit = _hash_sample(
        candidates,
        sample_pages=audit_pages,
        seed=audit_seed,
    )
    if not audit.empty:
        columns = [
            "candidate_id",
            "accepted_species",
            "trait_name",
            "normalized_value",
            "source_url",
            "page_title",
            "source_tier",
            "source_type",
            "domain",
            "language",
            "supporting_excerpt",
            "name_match_method",
            "wild_cultivated_cultivar_status",
        ]
        audit = audit[[column for column in columns if column in audit.columns]]
        for column in (
            "decision",
            "species_identity_correct",
            "value_correct",
            "provenance_complete",
            "cultivar_contamination",
            "false_positive_reason",
            "reviewer",
            "reviewed_at_utc",
        ):
            audit[column] = ""
    audit.to_csv(output_dir / "manual_audit_multidomain.csv", index=False)

    report: dict[str, object] = {
        "contract": CONTRACT,
        "completed_at_utc": _now(),
        "baseline_run_contract": summary.get("contract", ""),
        "domains_requested": requested,
        "domain_count": len(requested),
        "max_pages_per_domain": max_pages_per_domain,
        "candidate_rows_pending_review": len(candidate_rows),
        "candidate_species": (
            int(candidates["accepted_species"].nunique()) if not candidates.empty else 0
        ),
        "candidate_species_trait": (
            int(candidates.drop_duplicates(["accepted_species", "trait_name"]).shape[0])
            if not candidates.empty
            else 0
        ),
        "candidate_species_axis": (
            int(
                candidates.assign(axis=candidates["trait_name"].map(TRAIT_AXIS))
                .dropna(subset=["axis"])
                .drop_duplicates(["accepted_species", "axis"])
                .shape[0]
            )
            if not candidates.empty
            else 0
        ),
        "audit_pages_sampled": len(audit),
        "audit_sampling": "deterministic_hash_trait_stratified_unique_page",
        "accepted_evidence_rows": 0,
        "manual_review_required": True,
        "official_search_queries": 0,
        "family_inference": False,
        "global_fallback": False,
        "domain_results": domain_summary,
    }
    (output_dir / "multidomain_pilot_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


@app.command("run")
def run_command(
    baseline_dir: Path = typer.Option(..., exists=True, file_okay=False),
    master_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    registry_yaml: Path = typer.Option(
        Path("config/open_web_domain_registry_v2.yml"),
        exists=True,
    ),
    config_yaml: Path = typer.Option(
        Path("config/open_web_trait_acquisition.yml"),
        exists=True,
    ),
    domain: list[str] = typer.Option([], "--domain"),
    max_pages_per_domain: int = typer.Option(150, min=1, max=2000),
    audit_pages: int = typer.Option(120, min=20, max=1000),
    pause_seconds: float = typer.Option(0.25, min=0, max=30),
    audit_seed: str = typer.Option("open-web-multidomain-v1"),
) -> None:
    registry = DomainRegistry.load(registry_yaml)
    domains = domain or [
        record.domain
        for record in registry.all()
        if record.review_status == "approved" and record.inventory_urls
    ]
    report = run_multidomain_pilot(
        baseline_dir=baseline_dir,
        master_csv=master_csv,
        registry_yaml=registry_yaml,
        config_yaml=config_yaml,
        output_dir=output_dir,
        domains=domains,
        max_pages_per_domain=max_pages_per_domain,
        audit_pages=audit_pages,
        pause_seconds=pause_seconds,
        audit_seed=audit_seed,
    )
    typer.echo(json.dumps(report, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    app()
