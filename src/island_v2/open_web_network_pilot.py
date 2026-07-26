"""Credential-free broad Web discovery pilot for unresolved plant traits.

This module complements the official Brave/Google adapters with a deliberately
small, rate-limited Marginalia Search pilot.  Marginalia is useful here because
its index emphasizes non-commercial and independent sites, including regional
floras and specialist personal botanical pages.  Search results are discovery
leads only.  Every candidate still requires a fetched source page, exact species
identity, an exact supporting quote, and manual domain/candidate review.

Unknown domains never become accepted evidence automatically.  They are written
as review-only candidates and domain-registry proposals.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import time
import urllib.parse
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import httpx
import pandas as pd
import typer

from island_v2.open_web_evidence import (
    DomainRegistry,
    PageFetcher,
    StrictNameResolver,
    canonical_url,
    extract_rules,
    identity_gate,
    load_yaml,
    source_lineage,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT = "open_web_network_discovery_pilot_v1"
MARGINALIA_ENDPOINT = "https://api.marginalia.nu/public/search"
USER_AGENT = "island-floral-traits/0.2 (small non-commercial discovery pilot)"

# These are the manuscript-facing three-axis mappings.  Reward and pollen-vector
# mode remain useful separate traits, but neither is allowed to resolve the
# floral-structural-complexity axis.
AXIS_TRAITS: dict[str, tuple[str, ...]] = {
    "flower_colour": ("flower_primary_color",),
    "floral_structural_complexity": (
        "floral_form",
        "floral_symmetry",
        "tube_depth_class",
        "flower_size_class",
        "inflorescence_display",
    ),
    "reproductive_assurance": (
        "self_incompatibility",
        "autonomous_selfing_capacity",
        "mating_system",
        "cleistogamy",
    ),
}
TRAIT_AXIS = {trait: axis for axis, traits in AXIS_TRAITS.items() for trait in traits}

SEARCH_TERMS: dict[str, tuple[str, ...]] = {
    "flower_primary_color": ("flower color", "corolla colour", "floral description"),
    "floral_form": ("flower shape", "corolla shape", "floral morphology"),
    "floral_symmetry": ("flower symmetry", "actinomorphic", "zygomorphic"),
    "tube_depth_class": ("corolla tube length", "floral tube", "flower description"),
    "flower_size_class": ("flower size", "corolla length", "floral description"),
    "inflorescence_display": ("inflorescence", "flower arrangement", "floral description"),
    "self_incompatibility": ("self compatible", "self incompatible", "breeding system"),
    "autonomous_selfing_capacity": ("autonomous selfing", "bagging experiment", "autofertility"),
    "mating_system": ("mating system", "breeding system", "selfing rate"),
    "cleistogamy": ("cleistogamy", "cleistogamous", "closed flowers"),
}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).replace("\xa0", " ").split())


def _now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _domain(url: object) -> str:
    return urllib.parse.urlsplit(_text(url)).netloc.casefold().removeprefix("www.")


def _write_csv(path: Path, rows: Sequence[Mapping[str, object]], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _load_baseline(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    coverage_path = root / "integrated_species_axis_coverage.csv.gz"
    evidence_path = root / "integrated_evidence_lineage.csv.gz"
    summary_path = root / "integrated_trait_coverage_summary.json"
    for path in (coverage_path, evidence_path, summary_path):
        if not path.exists():
            raise ValueError(f"baseline artifact missing {path.name}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    denominator = summary.get("denominator", {})
    if (
        int(denominator.get("species", 0)) != 106_295
        or int(denominator.get("axes", 0)) != 3
        or int(denominator.get("species_axis", 0)) != 318_885
    ):
        raise ValueError("baseline denominator is not 106,295 x 3")
    coverage = pd.read_csv(coverage_path, dtype=str).fillna("")
    evidence = pd.read_csv(evidence_path, dtype=str).fillna("")
    if len(coverage) != 318_885 or coverage.duplicated(["accepted_species", "axis"]).any():
        raise ValueError("baseline coverage ledger violates the fixed denominator")
    return coverage, evidence


def build_priority_queue(
    coverage: pd.DataFrame,
    evidence: pd.DataFrame,
    master: pd.DataFrame,
    *,
    species_limit: int,
    traits_per_species: int,
) -> pd.DataFrame:
    """Select high-information unresolved species-trait tasks.

    Priority is highest when two direct congeners agree and one more observation
    could unlock a supported genus-trait rule.  The queue remains species-trait
    specific and never treats reward or pollen vector as structural resolution.
    """

    direct = evidence.loc[
        evidence["quality"].isin({"high", "medium"})
        & evidence["trait_name"].isin(TRAIT_AXIS)
        & evidence["evidence_scope"].isin({"species_direct", "synonym_direct"})
    ].copy()
    direct["genus"] = direct["accepted_species"].str.split().str[0]
    votes = direct.drop_duplicates(
        ["accepted_species", "trait_name", "normalized_value", "source_lineage"]
    )
    species_values = (
        votes.groupby(["genus", "trait_name", "accepted_species"])["normalized_value"]
        .agg(lambda values: sorted({_text(value) for value in values if _text(value)}))
        .reset_index()
    )
    species_values = species_values.loc[species_values["normalized_value"].map(len).eq(1)].copy()
    species_values["value"] = species_values["normalized_value"].str[0]

    genus_stats: dict[tuple[str, str], dict[str, object]] = {}
    for (genus, trait), group in species_values.groupby(["genus", "trait_name"]):
        counts = group["value"].value_counts()
        genus_stats[(genus, trait)] = {
            "support": int(group["accepted_species"].nunique()),
            "dominance": float(counts.iloc[0] / counts.sum()),
            "distribution": json.dumps(counts.to_dict(), ensure_ascii=False, sort_keys=True),
        }

    master = master.drop_duplicates("accepted_species").copy()
    if "genus" not in master.columns:
        master["genus"] = ""
    master["genus"] = master["genus"].where(
        master["genus"].ne(""), master["accepted_species"].str.split().str[0]
    )
    unresolved = coverage.loc[coverage["after_quality"].eq("")].merge(
        master[[column for column in ("accepted_species", "genus", "family", "n_records") if column in master]],
        on="accepted_species",
        how="left",
    )
    unresolved["genus"] = unresolved["genus"].where(
        unresolved["genus"].ne(""), unresolved["accepted_species"].str.split().str[0]
    )
    congener_counts = unresolved.groupby(["genus", "axis"])["accepted_species"].nunique().to_dict()

    rows: list[dict[str, object]] = []
    for target in unresolved.itertuples(index=False):
        axis = _text(target.axis)
        for trait in AXIS_TRAITS[axis]:
            stats = genus_stats.get((_text(target.genus), trait), {})
            support = int(stats.get("support", 0))
            dominance = float(stats.get("dominance", 0.0))
            congeners = int(congener_counts.get((_text(target.genus), axis), 1))
            near_rule = support == 2 and dominance == 1.0
            records = float(_text(getattr(target, "n_records", 0)) or 0)
            endpoint_weight = 1.4 if axis == "reproductive_assurance" else 1.0
            score = (
                15.0 * int(near_rule)
                + 2.5 * math.log1p(congeners)
                + min(3.0, support * 0.6)
                + 0.1 * math.log1p(records)
                + endpoint_weight
            )
            rows.append(
                {
                    "accepted_species": target.accepted_species,
                    "genus": target.genus,
                    "family": _text(getattr(target, "family", "")),
                    "axis": axis,
                    "trait_name": trait,
                    "current_support": support,
                    "current_dominance": dominance,
                    "current_value_distribution": stats.get("distribution", "{}"),
                    "unresolved_congener_count": congeners,
                    "potential_species_trait_unlocked": congeners if near_rule else 1,
                    "priority_reason": (
                        "two_agree_one_more_unlocks_rule"
                        if near_rule
                        else "no_congener_direct_evidence"
                        if support == 0
                        else "one_supporting_congener"
                        if support == 1
                        else "conflicted_or_under_threshold"
                    ),
                    "priority_score": round(score, 6),
                }
            )
    queue = pd.DataFrame(rows).sort_values(
        ["priority_score", "potential_species_trait_unlocked", "accepted_species", "trait_name"],
        ascending=[False, False, True, True],
    )
    species_order = queue.drop_duplicates("accepted_species").head(species_limit)["accepted_species"]
    selected = queue.loc[queue["accepted_species"].isin(set(species_order))].copy()
    selected["species_rank"] = selected["accepted_species"].map(
        {species: index + 1 for index, species in enumerate(species_order)}
    )
    selected = (
        selected.sort_values(["species_rank", "priority_score"], ascending=[True, False])
        .groupby("accepted_species", sort=False)
        .head(traits_per_species)
        .reset_index(drop=True)
    )
    return selected


def marginalia_search(
    queue: pd.DataFrame,
    *,
    max_queries: int,
    results_per_query: int,
    pause_seconds: float,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    """Query the rate-limited public Marginalia endpoint.

    The shared public key may return 503 when exhausted.  Such failures are
    recorded and never converted into biological absence.
    """

    hits: list[dict[str, object]] = []
    errors: list[dict[str, object]] = []
    licenses: set[str] = set()
    client = httpx.Client(timeout=35, follow_redirects=True, headers={"User-Agent": USER_AGENT})
    used = 0
    for row in queue.itertuples(index=False):
        if used >= max_queries:
            break
        terms = SEARCH_TERMS[_text(row.trait_name)]
        query = f'"{_text(row.accepted_species)}" ({" ".join(terms)})'
        encoded = urllib.parse.quote(query, safe="")
        url = f"{MARGINALIA_ENDPOINT}/{encoded}"
        used += 1
        try:
            response = client.get(url, params={"index": 0, "count": results_per_query})
            if response.status_code == 503:
                errors.append(
                    {
                        "accepted_species": row.accepted_species,
                        "trait_name": row.trait_name,
                        "query": query,
                        "status": 503,
                        "error": "shared_public_rate_limit_or_unavailable",
                    }
                )
                time.sleep(max(2.0, pause_seconds))
                continue
            response.raise_for_status()
            payload = response.json()
            if _text(payload.get("license")):
                licenses.add(_text(payload.get("license")))
            results = payload.get("results", []) if isinstance(payload, Mapping) else []
            for rank, item in enumerate(results[:results_per_query], start=1):
                if not isinstance(item, Mapping):
                    continue
                result_url = canonical_url(item.get("url"))
                if not result_url:
                    continue
                hits.append(
                    {
                        "accepted_species": row.accepted_species,
                        "trait_name": row.trait_name,
                        "query": query,
                        "provider": "marginalia_public",
                        "rank": rank,
                        "result_url": result_url,
                        "result_domain": _domain(result_url),
                        "title": _text(item.get("title")),
                        "snippet_discovery_only": _text(item.get("description")),
                    }
                )
        except (httpx.HTTPError, ValueError, json.JSONDecodeError) as exc:
            errors.append(
                {
                    "accepted_species": row.accepted_species,
                    "trait_name": row.trait_name,
                    "query": query,
                    "status": 0,
                    "error": f"{type(exc).__name__}:{_text(exc)}",
                }
            )
        time.sleep(pause_seconds)

    deduped: dict[tuple[str, str, str], dict[str, object]] = {}
    for hit in hits:
        key = (_text(hit["accepted_species"]), _text(hit["trait_name"]), _text(hit["result_url"]))
        deduped.setdefault(key, hit)
    report = {
        "provider": "Marginalia Search legacy public API",
        "queries_used": used,
        "hits": len(deduped),
        "errors": len(errors),
        "licenses_reported": sorted(licenses),
        "snippets_are_evidence": False,
        "public_shared_rate_limit": True,
    }
    return list(deduped.values()), errors, report


def fetch_candidates(
    hits: Sequence[Mapping[str, object]],
    *,
    master: pd.DataFrame,
    extraction_config: Mapping[str, Any],
    registry: DomainRegistry,
    max_pages: int,
    page_pause_seconds: float,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    resolver = StrictNameResolver(master.drop_duplicates("accepted_species").to_dict("records"))
    fetcher = PageFetcher(pause_seconds=page_pause_seconds)
    candidates: list[dict[str, object]] = []
    page_audit: list[dict[str, object]] = []
    domain_stats: dict[str, Counter[str]] = defaultdict(Counter)
    seen_urls: set[str] = set()

    for hit in hits:
        if len(seen_urls) >= max_pages:
            break
        url = canonical_url(hit.get("result_url"))
        if not url or url in seen_urls:
            continue
        seen_urls.add(url)
        domain = _domain(url)
        domain_stats[domain]["search_hits"] += 1
        page = fetcher.fetch(url)
        if page.error:
            page_audit.append(
                {
                    **dict(hit),
                    "page_status": "fetch_failed",
                    "identity_reason": page.error,
                    "content_sha256": page.content_sha256,
                }
            )
            domain_stats[domain]["fetch_failed"] += 1
            continue
        resolution, page_name, reason = identity_gate(
            page,
            searched_name=_text(hit.get("accepted_species")),
            expected_species=_text(hit.get("accepted_species")),
            resolver=resolver,
        )
        page_audit.append(
            {
                **dict(hit),
                "page_status": "identity_rejected" if reason else "identity_passed",
                "identity_reason": reason,
                "matched_page_name": page_name,
                "name_match_method": resolution.method,
                "content_sha256": page.content_sha256,
            }
        )
        if reason:
            domain_stats[domain]["identity_rejected"] += 1
            continue
        domain_stats[domain]["identity_passed"] += 1
        trait = _text(hit.get("trait_name"))
        registry_record = registry.lookup(url)
        extracted = extract_rules(page, trait_name=trait, config=extraction_config)
        if not extracted:
            domain_stats[domain]["identity_pass_no_trait"] += 1
            continue
        for raw_value, normalized_value, excerpt in extracted:
            lineage, lineage_method = source_lineage(url, excerpt)
            candidate_id = _sha(
                "|".join(
                    [
                        resolution.accepted_species,
                        trait,
                        normalized_value,
                        lineage,
                        excerpt,
                    ]
                )
            )[:24]
            approved = registry_record is not None
            candidates.append(
                {
                    "candidate_id": candidate_id,
                    "accepted_species": resolution.accepted_species,
                    "trait_name": trait,
                    "axis": TRAIT_AXIS[trait],
                    "raw_value": raw_value,
                    "normalized_value": normalized_value,
                    "source_url": url,
                    "page_title": page.title,
                    "domain": domain,
                    "language": page.language or "unknown_not_recorded",
                    "supporting_excerpt": excerpt,
                    "source_lineage": lineage,
                    "lineage_method": lineage_method,
                    "name_match_method": resolution.method,
                    "matched_page_name": page_name,
                    "source_tier": (
                        registry_record.recommended_evidence_tier if approved else "D_pending_domain_review"
                    ),
                    "source_type": (
                        registry_record.source_type if approved else "unclassified_open_web"
                    ),
                    "domain_review_status": "approved_registry" if approved else "pending_registry_review",
                    "evidence_status": "review_only_not_promoted",
                    "content_sha256": page.content_sha256,
                    "retrieved_at_utc": page.retrieved_at_utc,
                    "query": _text(hit.get("query")),
                    "search_rank": int(hit.get("rank") or 0),
                }
            )
            domain_stats[domain]["candidate_rows"] += 1

    domain_rows: list[dict[str, object]] = []
    for domain, counts in sorted(domain_stats.items()):
        identity_passed = counts["identity_passed"]
        candidate_rows = counts["candidate_rows"]
        domain_rows.append(
            {
                "domain": domain,
                "search_hits": counts["search_hits"],
                "fetch_failed": counts["fetch_failed"],
                "identity_passed": identity_passed,
                "identity_rejected": counts["identity_rejected"],
                "identity_pass_no_trait": counts["identity_pass_no_trait"],
                "candidate_rows": candidate_rows,
                "registry_proposal": (
                    "priority_manual_review"
                    if identity_passed >= 2 and candidate_rows >= 2
                    else "retain_discovery_only"
                ),
                "automatic_acceptance": False,
            }
        )
    return candidates, page_audit, domain_rows


def audit_sample(candidates: pd.DataFrame, limit: int = 100) -> pd.DataFrame:
    """Return a deterministic trait-stratified sample with unique pages."""

    if candidates.empty:
        return candidates.copy()
    table = candidates.copy()
    table["_sample_hash"] = table.apply(
        lambda row: _sha(f"{row['source_url']}|{row['candidate_id']}"), axis=1
    )
    pools = {
        trait: group.sort_values("_sample_hash").to_dict("records")
        for trait, group in table.groupby("trait_name", sort=True)
    }
    positions = {trait: 0 for trait in pools}
    selected: list[dict[str, object]] = []
    seen_urls: set[str] = set()
    while len(selected) < limit:
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
            if len(selected) >= limit:
                break
        if not progressed:
            break
    result = pd.DataFrame(selected).drop(columns=["_sample_hash"], errors="ignore")
    for column in (
        "decision",
        "species_identity_correct",
        "value_correct",
        "provenance_complete",
        "cultivar_contamination",
        "domain_quality_decision",
        "false_positive_reason",
        "reviewer",
        "reviewed_at_utc",
    ):
        result[column] = ""
    return result


def run(
    *,
    baseline_dir: Path,
    master_csv: Path,
    extraction_config_yaml: Path,
    registry_yaml: Path,
    output_dir: Path,
    species_limit: int,
    traits_per_species: int,
    max_queries: int,
    results_per_query: int,
    max_pages: int,
    query_pause_seconds: float,
    page_pause_seconds: float,
) -> dict[str, object]:
    coverage, evidence = _load_baseline(baseline_dir)
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    queue = build_priority_queue(
        coverage,
        evidence,
        master,
        species_limit=species_limit,
        traits_per_species=traits_per_species,
    )
    hits, search_errors, search_report = marginalia_search(
        queue,
        max_queries=max_queries,
        results_per_query=results_per_query,
        pause_seconds=query_pause_seconds,
    )
    registry = DomainRegistry.load(registry_yaml)
    candidates, page_audit, domains = fetch_candidates(
        hits,
        master=master,
        extraction_config=load_yaml(extraction_config_yaml),
        registry=registry,
        max_pages=max_pages,
        page_pause_seconds=page_pause_seconds,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    queue.to_csv(output_dir / "priority_queue.csv", index=False)
    _write_csv(
        output_dir / "search_hits.csv",
        hits,
        (
            "accepted_species",
            "trait_name",
            "query",
            "provider",
            "rank",
            "result_url",
            "result_domain",
            "title",
            "snippet_discovery_only",
        ),
    )
    _write_csv(
        output_dir / "search_errors.csv",
        search_errors,
        ("accepted_species", "trait_name", "query", "status", "error"),
    )
    candidate_columns = (
        "candidate_id",
        "accepted_species",
        "trait_name",
        "axis",
        "raw_value",
        "normalized_value",
        "source_url",
        "page_title",
        "domain",
        "language",
        "supporting_excerpt",
        "source_lineage",
        "lineage_method",
        "name_match_method",
        "matched_page_name",
        "source_tier",
        "source_type",
        "domain_review_status",
        "evidence_status",
        "content_sha256",
        "retrieved_at_utc",
        "query",
        "search_rank",
    )
    _write_csv(output_dir / "provisional_open_web_candidates.csv", candidates, candidate_columns)
    _write_csv(
        output_dir / "page_identity_audit.csv",
        page_audit,
        tuple(page_audit[0]) if page_audit else (
            "accepted_species",
            "trait_name",
            "result_url",
            "page_status",
            "identity_reason",
        ),
    )
    _write_csv(
        output_dir / "discovered_domain_registry_candidates.csv",
        domains,
        (
            "domain",
            "search_hits",
            "fetch_failed",
            "identity_passed",
            "identity_rejected",
            "identity_pass_no_trait",
            "candidate_rows",
            "registry_proposal",
            "automatic_acceptance",
        ),
    )
    audit_sample(pd.DataFrame(candidates), limit=100).to_csv(
        output_dir / "manual_audit_sample.csv", index=False
    )

    candidate_frame = pd.DataFrame(candidates)
    report: dict[str, object] = {
        "contract": CONTRACT,
        "completed_at_utc": _now(),
        "baseline_filled_species_axis": int(coverage["after_quality"].ne("").sum()),
        "baseline_unresolved_species_axis": int(coverage["after_quality"].eq("").sum()),
        "queue_species": int(queue["accepted_species"].nunique()),
        "queue_species_trait": len(queue),
        "search": search_report,
        "unique_domains_discovered": len({_text(hit.get("result_domain")) for hit in hits}),
        "source_pages_attempted": len(page_audit),
        "identity_passed_pages": sum(row.get("page_status") == "identity_passed" for row in page_audit),
        "provisional_candidate_rows": len(candidates),
        "provisional_species_trait": (
            len(candidate_frame.drop_duplicates(["accepted_species", "trait_name"]))
            if not candidate_frame.empty
            else 0
        ),
        "provisional_species_axis": (
            len(candidate_frame.drop_duplicates(["accepted_species", "axis"]))
            if not candidate_frame.empty
            else 0
        ),
        "priority_domains_for_manual_review": sum(
            row["registry_proposal"] == "priority_manual_review" for row in domains
        ),
        "accepted_evidence_rows": 0,
        "automatic_domain_approval": False,
        "manual_review_required": True,
        "family_inference": False,
        "global_fallback": False,
        "strict_axis_traits": AXIS_TRAITS,
        "separate_non_axis_traits": ["reward_type", "pollen_vector_mode"],
    }
    (output_dir / "network_pilot_report.json").write_text(
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
    extraction_config_yaml: Path = typer.Option(
        Path("config/open_web_trait_acquisition.yml"), exists=True, dir_okay=False
    ),
    registry_yaml: Path = typer.Option(
        Path("config/open_web_domain_registry.yml"), exists=True, dir_okay=False
    ),
    species_limit: int = typer.Option(80, min=1, max=500),
    traits_per_species: int = typer.Option(1, min=1, max=3),
    max_queries: int = typer.Option(60, min=1, max=500),
    results_per_query: int = typer.Option(5, min=1, max=10),
    max_pages: int = typer.Option(180, min=1, max=1000),
    query_pause_seconds: float = typer.Option(1.25, min=0.5, max=10),
    page_pause_seconds: float = typer.Option(0.35, min=0.1, max=10),
) -> None:
    run(
        baseline_dir=baseline_dir,
        master_csv=master_csv,
        extraction_config_yaml=extraction_config_yaml,
        registry_yaml=registry_yaml,
        output_dir=output_dir,
        species_limit=species_limit,
        traits_per_species=traits_per_species,
        max_queries=max_queries,
        results_per_query=results_per_query,
        max_pages=max_pages,
        query_pause_seconds=query_pause_seconds,
        page_pause_seconds=page_pause_seconds,
    )


if __name__ == "__main__":
    app()
