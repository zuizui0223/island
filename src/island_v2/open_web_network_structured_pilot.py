"""Broad-Web CLI that extracts every strict trait from each fetched page."""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence

from island_v2 import open_web_network_pilot as network
from island_v2.open_web_evidence import (
    PageFetcher,
    StrictNameResolver,
    canonical_url,
    identity_gate,
    source_lineage,
)
from island_v2.open_web_structured_traits import extract_structured_characteristics

_RULE_EXTRACTOR = network.extract_rules
_QUEUE_BUILDER = network.build_priority_queue


def _combined_extract_rules(
    page,
    *,
    trait_name: str,
    config,
    treatment_end_pattern: str = "",
):
    # Explicit key-value blocks have a clearer subject boundary than sentences
    # elsewhere on the page.  When present, they replace rather than supplement
    # free-text matches, preventing comparison keys or alternative options from
    # creating contradictory pseudo-records.
    structured_rows = extract_structured_characteristics(page, trait_name=trait_name)
    if structured_rows:
        return structured_rows
    return _RULE_EXTRACTOR(
        page,
        trait_name=trait_name,
        config=config,
        treatment_end_pattern=treatment_end_pattern,
    )


def _diversified_priority_queue(
    coverage,
    evidence,
    master,
    *,
    species_limit: int,
    traits_per_species: int,
):
    """Retain information value while preventing one large genus dominating discovery."""

    expanded_limit = min(500, max(species_limit * 20, species_limit))
    expanded = _QUEUE_BUILDER(
        coverage,
        evidence,
        master,
        species_limit=expanded_limit,
        traits_per_species=traits_per_species,
    ).sort_values(
        ["priority_score", "potential_species_trait_unlocked", "accepted_species"],
        ascending=[False, False, True],
    )
    species_rows = expanded.drop_duplicates("accepted_species")
    selected_species: list[str] = []
    genus_counts: Counter[str] = Counter()
    for row in species_rows.itertuples(index=False):
        genus = network._text(row.genus)
        if genus_counts[genus] >= 2:
            continue
        selected_species.append(row.accepted_species)
        genus_counts[genus] += 1
        if len(selected_species) >= species_limit:
            break
    if len(selected_species) < species_limit:
        chosen = set(selected_species)
        for species in species_rows["accepted_species"]:
            if species in chosen:
                continue
            selected_species.append(species)
            chosen.add(species)
            if len(selected_species) >= species_limit:
                break
    rank = {species: index + 1 for index, species in enumerate(selected_species)}
    selected = expanded.loc[expanded["accepted_species"].isin(rank)].copy()
    selected["species_rank"] = selected["accepted_species"].map(rank)
    return selected.sort_values(
        ["species_rank", "priority_score"], ascending=[True, False]
    ).reset_index(drop=True)


def _fetch_all_strict_traits(
    hits: Sequence[Mapping[str, object]],
    *,
    master,
    extraction_config,
    registry,
    max_pages: int,
    page_pause_seconds: float,
):
    """Fetch each URL once and extract all strict traits after identity passes."""

    resolver = StrictNameResolver(master.drop_duplicates("accepted_species").to_dict("records"))
    fetcher = PageFetcher(pause_seconds=page_pause_seconds)
    candidates: list[dict[str, object]] = []
    page_audit: list[dict[str, object]] = []
    domain_stats: dict[str, Counter[str]] = defaultdict(Counter)
    seen_urls: set[str] = set()
    seen_candidates: set[tuple[str, str, str, str]] = set()

    for hit in hits:
        if len(seen_urls) >= max_pages:
            break
        url = canonical_url(hit.get("result_url"))
        if not url or url in seen_urls:
            continue
        seen_urls.add(url)
        domain = network._domain(url)
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

        expected_species = network._text(hit.get("accepted_species"))
        resolution, page_name, reason = identity_gate(
            page,
            searched_name=expected_species,
            expected_species=expected_species,
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

        registry_record = registry.lookup(url)
        page_candidate_count = 0
        for trait in network.TRAIT_AXIS:
            extracted = _combined_extract_rules(
                page,
                trait_name=trait,
                config=extraction_config,
                treatment_end_pattern=(
                    registry_record.treatment_end_pattern if registry_record else ""
                ),
            )
            for raw_value, normalized_value, excerpt in extracted:
                lineage, lineage_method = source_lineage(url, excerpt)
                candidate_key = (
                    resolution.accepted_species,
                    trait,
                    normalized_value,
                    lineage,
                )
                if candidate_key in seen_candidates:
                    continue
                seen_candidates.add(candidate_key)
                candidate_id = network._sha("|".join([*candidate_key, excerpt]))[:24]
                approved = registry_record is not None and registry_record.allows_trait(trait)
                candidates.append(
                    {
                        "candidate_id": candidate_id,
                        "accepted_species": resolution.accepted_species,
                        "trait_name": trait,
                        "axis": network.TRAIT_AXIS[trait],
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
                            registry_record.recommended_evidence_tier
                            if approved
                            else "D_pending_domain_or_trait_review"
                        ),
                        "source_type": (
                            registry_record.source_type
                            if approved
                            else "unclassified_open_web"
                        ),
                        "domain_review_status": (
                            "approved_registry_trait"
                            if approved
                            else "pending_registry_or_trait_review"
                        ),
                        "evidence_status": "review_only_not_promoted",
                        "content_sha256": page.content_sha256,
                        "retrieved_at_utc": page.retrieved_at_utc,
                        "query": network._text(hit.get("query")),
                        "search_rank": int(hit.get("rank") or 0),
                    }
                )
                page_candidate_count += 1
                domain_stats[domain]["candidate_rows"] += 1
        if page_candidate_count == 0:
            domain_stats[domain]["identity_pass_no_trait"] += 1

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


network.extract_rules = _combined_extract_rules
network.build_priority_queue = _diversified_priority_queue
network.fetch_candidates = _fetch_all_strict_traits


if __name__ == "__main__":
    network.app()
