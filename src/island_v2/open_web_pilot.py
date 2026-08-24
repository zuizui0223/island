"""Prioritize, run, and audit an open-Web trait acquisition pilot.

This module owns the reproducible pilot interface.  It consumes the formal
integrated-coverage artifact, keeps the 106,295 x 3 denominator fixed, builds a
species x trait queue, and reports gains only from manually accepted
species-direct evidence.  Registry crawling is a credential-free pilot adapter;
production discovery remains an official-search adapter.
"""

# Typer's documented command style evaluates Option metadata in defaults.
# ruff: noqa: B008

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.open_web_common import (
    STRICT_AXIS_TRAITS,
    STRICT_TRAIT_AXIS,
    rebuild_with_common_all_evidence,
)
from island_v2.open_web_evidence import (
    DomainRegistry,
    PageFetcher,
    StrictNameResolver,
    crawl_registry_inventory,
    load_yaml,
    process_pages,
    promote_reviewed,
)

CONTRACT = "open_web_trait_pilot_v1"
AXIS_TRAITS = STRICT_AXIS_TRAITS
TRAIT_AXIS = STRICT_TRAIT_AXIS
DIRECT_QUALITY = {"high", "medium"}
MIN_AUDITED_PER_PRODUCTION_TRAIT = 10
BASELINE_FILES = {
    "summary": "integrated_trait_coverage_summary.json",
    "coverage": "integrated_species_axis_coverage.csv.gz",
    "evidence": "integrated_evidence_lineage.csv.gz",
    "bulk_audit": "bulk_candidate_acceptance_audit.csv.gz",
    "duplicates": "source_lineage_duplicates.csv.gz",
}

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _domain(value: object) -> str:
    from urllib.parse import urlsplit

    return urlsplit(_text(value)).netloc.casefold().removeprefix("www.") or "unknown"


def _language_from_evidence(row: pd.Series) -> str:
    explicit = _text(row.get("language"))
    if explicit:
        return explicit
    domain = _domain(row.get("source_url"))
    if domain.endswith("wikipedia.org"):
        return domain.split(".", 1)[0]
    if domain in {"efloras.org", "plants.ces.ncsu.edu", "missouribotanicalgarden.org"}:
        return "en"
    return "unknown_not_recorded"


def _source_type(row: pd.Series) -> str:
    explicit = _text(row.get("source_type"))
    if explicit:
        return explicit
    group = _text(row.get("source_group"))
    provider = _text(row.get("source_provider"))
    mapping = {
        "efloras": "public_flora",
        "eol_traitbank": "public_trait_database",
        "austraits": "research_trait_database",
        "gift": "research_trait_database",
        "usda_plants": "public_database",
        "meyer_2026": "original_research_dataset",
        "razanajatovo_2016": "original_research_dataset",
        "three_pass_medium": "mixed_public_web",
        "three_pass_high": "original_or_authoritative_source",
        "promoted_public_web": provider or "public_web",
        "validated_low": "validated_genus_inference",
    }
    return mapping.get(group, provider or "unknown_not_recorded")


def _read_baseline(
    directory: Path,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    root = Path(directory)
    missing = [name for name in BASELINE_FILES.values() if not (root / name).exists()]
    if missing:
        raise ValueError(f"baseline artifact is missing required files: {missing}")
    summary = json.loads((root / BASELINE_FILES["summary"]).read_text(encoding="utf-8"))
    coverage = pd.read_csv(root / BASELINE_FILES["coverage"], dtype=str).fillna("")
    evidence = pd.read_csv(root / BASELINE_FILES["evidence"], dtype=str).fillna("")
    bulk_audit = pd.read_csv(root / BASELINE_FILES["bulk_audit"], dtype=str).fillna("")
    return summary, coverage, evidence, bulk_audit


def _validate_denominator(summary: Mapping[str, Any], coverage: pd.DataFrame) -> None:
    denominator = summary.get("denominator", {})
    if (
        int(denominator.get("species", 0)) != 106_295
        or int(denominator.get("axes", 0)) != 3
        or int(denominator.get("species_axis", 0)) != 318_885
    ):
        raise ValueError("baseline denominator is not the fixed 106,295 x 3 contract")
    if len(coverage) != 318_885:
        raise ValueError(f"coverage ledger has {len(coverage)} rows, expected 318,885")
    if coverage.duplicated(["accepted_species", "axis"]).any():
        raise ValueError("coverage ledger has duplicate species-axis rows")


def baseline_audit(
    baseline_dir: Path,
    *,
    output_dir: Path,
) -> dict[str, Any]:
    summary, coverage, evidence, bulk_audit = _read_baseline(baseline_dir)
    _validate_denominator(summary, coverage)
    direct = evidence.loc[evidence["quality"].isin(DIRECT_QUALITY)].copy()
    trait_pairs = evidence[["accepted_species", "trait_name"]].drop_duplicates()
    direct_trait_pairs = direct[["accepted_species", "trait_name"]].drop_duplicates()
    n_unresolved_trait = 106_295 * len(TRAIT_AXIS) - len(trait_pairs)

    evidence["domain"] = evidence["source_url"].map(_domain)
    evidence["language_audit"] = evidence.apply(_language_from_evidence, axis=1)
    evidence["source_type_audit"] = evidence.apply(_source_type, axis=1)

    def grouped(column: str) -> list[dict[str, object]]:
        rows: list[dict[str, object]] = []
        for value, group in evidence.groupby(column, dropna=False):
            rows.append(
                {
                    column: _text(value) or "unknown_not_recorded",
                    "evidence_rows": len(group),
                    "species": int(group["accepted_species"].nunique()),
                    "species_trait": len(group.drop_duplicates(["accepted_species", "trait_name"])),
                }
            )
        return sorted(rows, key=lambda row: (-int(row["evidence_rows"]), str(row[column])))

    rejection_counts = bulk_audit.loc[
        bulk_audit["accepted"].astype(str).str.casefold().ne("true"), "reason"
    ].value_counts()
    name_failures = int(
        rejection_counts[
            rejection_counts.index.to_series().str.contains("name|species|synonym", case=False)
        ].sum()
    )
    provenance_failures = int(
        rejection_counts[
            rejection_counts.index.to_series().str.contains(
                "provenance|citation|excerpt|url|source_backed", case=False
            )
        ].sum()
    )
    ontology_failures = int(
        rejection_counts[
            rejection_counts.index.to_series().str.contains(
                "trait_outside|unmapped|value|axis", case=False
            )
        ].sum()
    )
    calculated_quality = (
        coverage.loc[coverage["after_quality"].ne(""), "after_quality"].value_counts().to_dict()
    )
    by_axis: dict[str, Any] = {}
    for axis, group in coverage.groupby("axis"):
        filled = int(group["after_quality"].ne("").sum())
        by_axis[axis] = {
            "filled_species": filled,
            "fill_rate": filled / 106_295,
            "unresolved_species": 106_295 - filled,
            "quality_counts": {
                quality: int(group["after_quality"].eq(quality).sum())
                for quality in ("high", "medium", "low")
            },
        }
    computed = {
        "contract": "open_web_baseline_audit_v1",
        "computed_at_utc": _now(),
        "source_artifact_directory": str(Path(baseline_dir).resolve()),
        "denominator": {"species": 106_295, "axes": 3, "species_axis": 318_885},
        "quality_counts": {
            "high": int(calculated_quality.get("high", 0)),
            "medium": int(calculated_quality.get("medium", 0)),
            "low": int(calculated_quality.get("low", 0)),
        },
        "filled_species_axis": int(coverage["after_quality"].ne("").sum()),
        "strict_fill_rate": float(coverage["after_quality"].ne("").mean()),
        "unresolved_species_axis": int(coverage["after_quality"].eq("").sum()),
        "by_axis": by_axis,
        "direct_species_trait": len(direct_trait_pairs),
        "all_tier_species_trait": len(trait_pairs),
        "unresolved_species_trait": int(n_unresolved_trait),
        "trait_denominator": int(106_295 * len(TRAIT_AXIS)),
        "lineage_deduplication": {
            "duplicate_rows": int(summary["input_audit"]["lineage_duplicate_rows"]),
            "duplicate_groups": int(summary["input_audit"]["lineage_duplicate_groups"]),
            "accepted_lineages_after": len(evidence),
            "rows_before_approx": int(
                len(evidence) + summary["input_audit"]["lineage_duplicate_rows"]
            ),
        },
        "bulk_rejection_scope": "bulk_candidate_acceptance_audit only",
        "name_match_failures": name_failures,
        "provenance_insufficient_exclusions": provenance_failures,
        "ontology_or_value_exclusions": ontology_failures,
        "provider_counts": grouped("source_provider"),
        "language_counts": grouped("language_audit"),
        "domain_counts": grouped("domain"),
        "source_type_counts": grouped("source_type_audit"),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "baseline_coverage_audit.json").write_text(
        json.dumps(computed, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    for name, rows in (
        ("provider_coverage.csv", computed["provider_counts"]),
        ("language_coverage.csv", computed["language_counts"]),
        ("domain_coverage.csv", computed["domain_counts"]),
        ("source_type_coverage.csv", computed["source_type_counts"]),
    ):
        pd.DataFrame(rows).to_csv(output_dir / name, index=False)
    return computed


def build_priority_queue(
    baseline_dir: Path,
    master_csv: Path,
    *,
    output_csv: Path,
    pilot_species: int = 300,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    summary, coverage, evidence, _ = _read_baseline(baseline_dir)
    _validate_denominator(summary, coverage)
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    universe = coverage[["accepted_species"]].drop_duplicates()
    master = universe.merge(master, on="accepted_species", how="left")
    master["genus"] = master["genus"].where(
        master["genus"].ne(""), master["accepted_species"].str.split().str[0]
    )
    evidence = evidence.loc[
        evidence["quality"].isin(DIRECT_QUALITY) & evidence["trait_name"].isin(TRAIT_AXIS)
    ].copy()
    evidence["genus"] = evidence["accepted_species"].str.split().str[0]
    votes = evidence.drop_duplicates(
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
        distribution = group["value"].value_counts()
        loo_outcomes: list[bool] = []
        for held_index, held in group.iterrows():
            training = group.drop(index=held_index)
            if training.empty:
                continue
            predicted = _text(training["value"].value_counts().index[0])
            loo_outcomes.append(predicted == _text(held["value"]))
        genus_stats[(genus, trait)] = {
            "support": int(group["accepted_species"].nunique()),
            "distribution": json.dumps(distribution.to_dict(), ensure_ascii=False, sort_keys=True),
            "dominance": float(distribution.iloc[0] / distribution.sum()),
            "agreement_at_two": bool(len(group) == 2 and len(distribution) == 1),
            "masked_accuracy": (sum(loo_outcomes) / len(loo_outcomes) if loo_outcomes else 0.0),
        }
    unresolved = coverage.loc[coverage["after_quality"].eq("")].merge(master, on="accepted_species")
    unresolved_by_genus_axis = (
        unresolved.groupby(["genus", "axis"])["accepted_species"].nunique().to_dict()
    )
    rows: list[dict[str, object]] = []
    for target in unresolved.itertuples(index=False):
        for trait in AXIS_TRAITS[_text(target.axis)]:
            stats = genus_stats.get((_text(target.genus), trait), {})
            support = int(stats.get("support", 0))
            dominance = float(stats.get("dominance", 0))
            congener_count = int(
                unresolved_by_genus_axis.get((_text(target.genus), _text(target.axis)), 1)
            )
            unlock = congener_count if support == 2 and dominance == 1 else 1
            endpoint_weight = 1.3 if target.axis == "reproductive_assurance" else 1.0
            occurrence_weight = math.log1p(float(_text(getattr(target, "n_records", 0)) or 0))
            near_rule = 1 if support == 2 and dominance == 1 else 0
            score = (
                12 * near_rule
                + 2.5 * math.log1p(congener_count)
                + endpoint_weight
                + 0.15 * occurrence_weight
                + min(3.0, support * 0.5)
            )
            rows.append(
                {
                    "accepted_species": target.accepted_species,
                    "genus": target.genus,
                    "family": getattr(target, "family", ""),
                    "axis": target.axis,
                    "trait_name": trait,
                    "current_support": support,
                    "current_value_distribution": stats.get("distribution", "{}"),
                    "current_dominance": dominance,
                    "masked_accuracy": float(stats.get("masked_accuracy", 0)),
                    "unresolved_congener_count": congener_count,
                    "potential_species_trait_unlocked": unlock,
                    "expected_information_gain": round(
                        near_rule * math.log2(1 + congener_count)
                        + (1 / (1 + support))
                        + endpoint_weight,
                        6,
                    ),
                    "query_cost_units": 1,
                    "domain_availability": "official_search_required_or_registry_match",
                    "priority_score": round(score, 6),
                    "priority_reason": (
                        "two_agree_one_more_unlocks_rule"
                        if near_rule
                        else {
                            0: "no_congener_direct_evidence",
                            1: "one_supporting_congener",
                            2: "two_congeners_not_unanimous",
                        }.get(support, "three_plus_but_not_currently_resolved")
                    ),
                    "recommended_query_source": (
                        "regional flora / institutional species page / exact synonym search"
                    ),
                }
            )
    queue = pd.DataFrame(rows).sort_values(
        [
            "priority_score",
            "potential_species_trait_unlocked",
            "expected_information_gain",
            "accepted_species",
            "trait_name",
        ],
        ascending=[False, False, False, True, True],
    )
    species_order = queue.drop_duplicates("accepted_species").head(pilot_species)[
        "accepted_species"
    ]
    selected = queue.loc[queue["accepted_species"].isin(set(species_order))].copy()
    selected["pilot_species_rank"] = selected["accepted_species"].map(
        {species: index + 1 for index, species in enumerate(species_order)}
    )
    selected = selected.sort_values(
        ["pilot_species_rank", "priority_score"], ascending=[True, False]
    )
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(output_csv, index=False)
    reason_rows = (
        queue.groupby(["axis", "priority_reason"], dropna=False)
        .agg(
            unresolved_species_trait=("accepted_species", "size"),
            unresolved_species=("accepted_species", "nunique"),
            genera=("genus", "nunique"),
        )
        .reset_index()
    )
    reason_rows.to_csv(output_csv.parent / "unresolved_reason_report.csv", index=False)

    settings = {
        "current": {
            "min_species": 3,
            "floral_structural_complexity": 0.80,
            "flower_colour": 0.90,
            "reproductive_assurance": 0.95,
            "floral_structural_complexity_masked": 0.75,
            "flower_colour_masked": 0.80,
            "reproductive_assurance_masked": 0.85,
        },
        "mild_relaxation": {
            "min_species": 3,
            "floral_structural_complexity": 0.75,
            "flower_colour": 0.80,
            "reproductive_assurance": 0.90,
            "floral_structural_complexity_masked": 0.75,
            "flower_colour_masked": 0.80,
            "reproductive_assurance_masked": 0.85,
        },
        "min_species_2_diagnostic_only": {
            "min_species": 2,
            "floral_structural_complexity": 0.75,
            "flower_colour": 0.80,
            "reproductive_assurance": 0.90,
            "floral_structural_complexity_masked": 0.75,
            "flower_colour_masked": 0.80,
            "reproductive_assurance_masked": 0.85,
        },
    }
    sensitivity_rows: list[dict[str, object]] = []
    for setting, thresholds in settings.items():
        eligible: set[tuple[str, str]] = set()
        for key, stats in genus_stats.items():
            axis = TRAIT_AXIS.get(key[1], "")
            if (
                int(stats["support"]) >= int(thresholds["min_species"])
                and float(stats["dominance"]) >= float(thresholds[axis])
                and float(stats["masked_accuracy"]) >= float(thresholds[f"{axis}_masked"])
            ):
                eligible.add(key)
        unlocked_trait = queue.loc[
            [
                (genus, trait) in eligible
                for genus, trait in zip(queue["genus"], queue["trait_name"])
            ]
        ]
        sensitivity_rows.append(
            {
                "setting": setting,
                "min_species": thresholds["min_species"],
                "structure_dominance": thresholds["floral_structural_complexity"],
                "colour_dominance": thresholds["flower_colour"],
                "reproductive_dominance": thresholds["reproductive_assurance"],
                "structure_masked_accuracy": thresholds["floral_structural_complexity_masked"],
                "colour_masked_accuracy": thresholds["flower_colour_masked"],
                "reproductive_masked_accuracy": thresholds["reproductive_assurance_masked"],
                "eligible_genus_trait_rules": len(eligible),
                "potential_resolved_species_trait": len(unlocked_trait),
                "potential_resolved_species_axis": len(
                    unlocked_trait.drop_duplicates(["accepted_species", "axis"])
                ),
                "analysis_role": (
                    "diagnostic_only"
                    if setting == "min_species_2_diagnostic_only"
                    else "candidate_secondary_robustness"
                ),
            }
        )
    pd.DataFrame(sensitivity_rows).to_csv(
        output_csv.parent / "threshold_sensitivity_report.csv", index=False
    )
    report = {
        "contract": "prioritized_open_web_queue_v1",
        "generated_at_utc": _now(),
        "pilot_species_requested": pilot_species,
        "pilot_species": int(selected["accepted_species"].nunique()),
        "species_trait_tasks": len(selected),
        "near_rule_species_trait_tasks": int(
            selected["priority_reason"].eq("two_agree_one_more_unlocks_rule").sum()
        ),
        "potential_species_trait_unlocked_sum": int(
            selected["potential_species_trait_unlocked"].sum()
        ),
        "unresolved_reason_rows": len(reason_rows),
        "threshold_settings": sensitivity_rows,
        "family_inference": False,
        "global_fallback": False,
    }
    (output_csv.parent / "priority_queue_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return selected, report


def _write_rows(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _stratified_unique_page_audit(
    candidates: pd.DataFrame,
    audit_pages: int,
    audit_seed: str = "open-web-audit-v2",
) -> pd.DataFrame:
    """Reproducibly sample unique pages within domain x trait strata."""

    ordered = candidates.copy()
    ordered["domain"] = ordered.get("domain", "unknown_not_recorded")
    ordered["_audit_hash"] = (
        ordered["candidate_id"].map(_text)
        + "|"
        + ordered["source_url"].map(_text)
        + "|"
        + audit_seed
    ).map(lambda value: hashlib.sha256(value.encode("utf-8")).hexdigest())
    ordered = ordered.sort_values(
        [
            "domain",
            "trait_name",
            "_audit_hash",
            "source_url",
            "candidate_id",
        ]
    )
    pools = {
        key: group.to_dict("records")
        for key, group in ordered.groupby(["domain", "trait_name"], sort=True)
    }
    positions = {key: 0 for key in pools}
    selected: list[dict[str, object]] = []
    seen_urls: set[str] = set()
    while len(selected) < audit_pages:
        progressed = False
        for key in sorted(pools):
            rows = pools[key]
            while positions[key] < len(rows):
                row = rows[positions[key]]
                positions[key] += 1
                url = _text(row.get("source_url"))
                if not url or url in seen_urls:
                    continue
                selected.append(row)
                seen_urls.add(url)
                progressed = True
                break
            if len(selected) >= audit_pages:
                break
        if not progressed:
            break
    return pd.DataFrame(selected).drop(columns=["_audit_hash"], errors="ignore")


def _production_approved_traits(by_trait: Mapping[str, Mapping[str, object]]) -> list[str]:
    return sorted(
        trait
        for trait, metrics in by_trait.items()
        if int(metrics["audited"]) >= MIN_AUDITED_PER_PRODUCTION_TRAIT
        and float(metrics["precision"]) >= 0.90
    )


def run_registry_pilot(
    *,
    baseline_dir: Path,
    master_csv: Path,
    config_yaml: Path,
    registry_yaml: Path,
    output_dir: Path,
    domain: str,
    pilot_species: int,
    max_pages: int,
    audit_pages: int,
    pause_seconds: float,
) -> dict[str, Any]:
    summary, coverage, evidence, _ = _read_baseline(baseline_dir)
    _validate_denominator(summary, coverage)
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    universe = set(coverage["accepted_species"])
    master = master.loc[master["accepted_species"].isin(universe)].drop_duplicates(
        "accepted_species"
    )
    registry = DomainRegistry.load(registry_yaml)
    record = next((item for item in registry.all() if item.domain == domain), None)
    if record is None or not record.inventory_urls:
        raise ValueError(f"{domain} has no approved registry inventory")
    fetcher = PageFetcher(pause_seconds=pause_seconds)
    urls, inventory_audit = crawl_registry_inventory(record, fetcher=fetcher, max_urls=max_pages)
    leads = [
        {
            "result_url": url,
            "provider": "approved_domain_registry_inventory",
            "query": f"registry:{domain}",
            "rank": index + 1,
            "searched_name": "",
            "accepted_species": "",
            "trait_name": "",
        }
        for index, url in enumerate(urls)
    ]
    result = process_pages(
        leads,
        registry=registry,
        resolver=StrictNameResolver(master.to_dict("records")),
        config=load_yaml(config_yaml),
        fetcher=fetcher,
        default_traits=tuple(TRAIT_AXIS),
    )
    candidates = pd.DataFrame(result["candidates"])
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
                for species, trait in zip(candidates["accepted_species"], candidates["trait_name"])
            ]
        ].copy()
        species_hash = candidates["accepted_species"].map(
            lambda value: hashlib.sha256(_text(value).encode("utf-8")).hexdigest()
        )
        candidates["_species_hash"] = species_hash
        selected_species = (
            candidates[["accepted_species", "_species_hash"]]
            .drop_duplicates()
            .sort_values("_species_hash")
            .head(pilot_species)["accepted_species"]
        )
        candidates = candidates.loc[candidates["accepted_species"].isin(set(selected_species))]
        candidates = candidates.drop(columns="_species_hash")
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
    output_dir.mkdir(parents=True, exist_ok=True)
    candidate_rows = candidates.to_dict("records") if not candidates.empty else []
    _write_rows(output_dir / "open_web_evidence_candidates.csv", candidate_rows)
    _write_rows(output_dir / "source_lineage_duplicates.csv", result["lineage_duplicates"])
    _write_rows(output_dir / "species_identity_audit.csv", result["identity_audit"])
    _write_rows(output_dir / "page_fetch_audit.csv", result["fetch_audit"])
    _write_rows(output_dir / "domain_outcome_audit.csv", result["domain_audit"])
    _write_rows(output_dir / "registry_inventory_audit.csv", inventory_audit)
    processed_pages = len(result["fetch_audit"])
    candidate_pages = len(
        {
            _text(row.get("source_url"))
            for row in result["candidates"]
            if _text(row.get("source_url"))
        }
    )
    registry_observed = {
        **{
            field: getattr(record, field)
            for field in (
                "domain",
                "source_name",
                "source_type",
                "region",
                "language",
                "organization_type",
                "robots_access_notes",
                "page_pattern",
                "treatment_end_pattern",
                "scientific_name_displayed",
                "cultivar_descriptions_present",
                "provenance_quality",
                "recommended_evidence_tier",
                "review_status",
            )
        },
        "pages_processed": processed_pages,
        "pages_with_candidate": candidate_pages,
        "success_rate": candidate_pages / processed_pages if processed_pages else None,
        "false_positive_audit": "pending_manual_review",
    }
    _write_rows(output_dir / "domain_registry_observed.csv", [registry_observed])

    if candidates.empty:
        audit = pd.DataFrame()
    else:
        audit = _stratified_unique_page_audit(candidates, audit_pages)
        audit = audit[
            [
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
        ]
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
    audit.to_csv(output_dir / "manual_audit_100_pages.csv", index=False)
    report = {
        "contract": CONTRACT,
        "completed_at_utc": _now(),
        "pilot_adapter": "approved_domain_registry_inventory",
        "domain": domain,
        "pilot_species_requested": pilot_species,
        "pilot_species_with_candidate": int(
            candidates["accepted_species"].nunique() if not candidates.empty else 0
        ),
        "inventory_pages_fetched": len(inventory_audit),
        "species_pages_discovered": len(urls),
        "source_pages_processed": len(result["fetch_audit"]),
        "candidate_rows_pending_review": len(candidate_rows),
        "candidate_species_trait_pending_review": int(
            len(candidates.drop_duplicates(["accepted_species", "trait_name"]))
            if not candidates.empty
            else 0
        ),
        "manual_audit_pages_required": audit_pages,
        "manual_audit_pages_sampled": len(audit),
        "discovered_domains": 1,
        "tier_counts_pending_review": {
            _text(tier): int(count)
            for tier, count in (
                candidates["source_tier"].value_counts().items() if not candidates.empty else []
            )
        },
        "accepted_evidence_rows": 0,
        "queries_used": 0,
        "query_cost_usd": 0.0,
        "official_search_credentials_used": False,
        "production_approved": False,
        "policy": {
            "registry_seed_is_not_general_search": True,
            "snippets_as_evidence": False,
            "manual_review_required": True,
            "family_inference": False,
            "global_fallback": False,
        },
    }
    (output_dir / "pilot_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


def _bool(value: object) -> bool:
    return _text(value).casefold() in {"1", "true", "yes", "y"}


def apply_manual_audit(
    *,
    baseline_dir: Path,
    candidates_csv: Path,
    audit_csv: Path,
    output_dir: Path,
) -> dict[str, Any]:
    summary, coverage, evidence, _ = _read_baseline(baseline_dir)
    _validate_denominator(summary, coverage)
    candidates = pd.read_csv(candidates_csv, dtype=str).fillna("")
    audits = pd.read_csv(audit_csv, dtype=str).fillna("")
    required = {
        "candidate_id",
        "decision",
        "species_identity_correct",
        "value_correct",
        "provenance_complete",
        "cultivar_contamination",
    }
    if missing := required.difference(audits.columns):
        raise ValueError(f"audit is missing columns: {sorted(missing)}")
    reviewed = audits.loc[audits["decision"].str.casefold().isin({"accept", "reject"})]
    accepted_audit = reviewed.loc[
        reviewed["decision"].str.casefold().eq("accept")
        & reviewed["species_identity_correct"].map(_bool)
        & reviewed["value_correct"].map(_bool)
        & reviewed["provenance_complete"].map(_bool)
        & ~reviewed["cultivar_contamination"].map(_bool)
    ]
    review_rows = accepted_audit.assign(decision="accept").to_dict("records")
    promoted = pd.DataFrame(promote_reviewed(candidates.to_dict("records"), review_rows))
    direct_baseline = evidence.loc[evidence["quality"].isin(DIRECT_QUALITY)].copy()
    baseline_trait = set(zip(direct_baseline["accepted_species"], direct_baseline["trait_name"]))
    new_trait = (
        {
            (row.accepted_species, row.trait_name)
            for row in promoted.itertuples(index=False)
            if (row.accepted_species, row.trait_name) not in baseline_trait
        }
        if not promoted.empty
        else set()
    )
    coverage_lookup = coverage.set_index(["accepted_species", "axis"])
    new_axes: set[tuple[str, str]] = set()
    low_upgrades: set[tuple[str, str]] = set()
    low_conflicts: set[tuple[str, str]] = set()
    for row in promoted.itertuples(index=False):
        axis = TRAIT_AXIS.get(row.trait_name, "")
        if not axis:
            continue
        key = (row.accepted_species, axis)
        baseline = coverage_lookup.loc[key]
        quality = _text(baseline["after_quality"])
        if not quality:
            new_axes.add(key)
        if quality == "low":
            low_upgrades.add(key)
            if _text(baseline["after_value"]) not in {"", row.normalized_value}:
                low_conflicts.add(key)
    precision_denominator = len(reviewed)
    precision = len(accepted_audit) / precision_denominator if precision_denominator else None
    cultivar_rate = (
        float(reviewed["cultivar_contamination"].map(_bool).mean()) if len(reviewed) else None
    )
    by_trait: dict[str, Any] = {}
    for trait, group in reviewed.groupby("trait_name"):
        good = group.loc[
            group["decision"].str.casefold().eq("accept")
            & group["species_identity_correct"].map(_bool)
            & group["value_correct"].map(_bool)
            & group["provenance_complete"].map(_bool)
            & ~group["cultivar_contamination"].map(_bool)
        ]
        by_trait[trait] = {
            "audited": len(group),
            "accepted_correct": len(good),
            "precision": len(good) / len(group) if len(group) else None,
        }
    production_approved_traits = _production_approved_traits(by_trait)
    grouped_precision: dict[str, list[dict[str, object]]] = {}
    for column in ("source_tier", "source_type", "domain", "language"):
        rows: list[dict[str, object]] = []
        if column not in reviewed.columns:
            grouped_precision[column] = rows
            continue
        for value, group in reviewed.groupby(column):
            good = group.loc[
                group["decision"].str.casefold().eq("accept")
                & group["species_identity_correct"].map(_bool)
                & group["value_correct"].map(_bool)
                & group["provenance_complete"].map(_bool)
                & ~group["cultivar_contamination"].map(_bool)
            ]
            rows.append(
                {
                    column: _text(value) or "unknown_not_recorded",
                    "audited": len(group),
                    "accepted_correct": len(good),
                    "precision": len(good) / len(group) if len(group) else None,
                    "cultivar_contamination_rate": float(
                        group["cultivar_contamination"].map(_bool).mean()
                    ),
                }
            )
        grouped_precision[column] = rows
    output_dir.mkdir(parents=True, exist_ok=True)
    promoted.to_csv(output_dir / "accepted_open_web_evidence.csv", index=False)
    shared_report, shared_tables = rebuild_with_common_all_evidence(
        coverage=coverage,
        evidence=evidence,
        promoted=promoted,
    )
    rebuilt_rules = shared_tables["rules"]
    rebuilt_low = shared_tables["rebuilt_low"]
    direct_conflicts = shared_tables["direct_conflicts"]
    low_rebuild = shared_report["validated_low"]
    rebuilt_rules.to_csv(output_dir / "rebuilt_validated_genus_rules.csv.gz", index=False)
    rebuilt_low.to_csv(output_dir / "rebuilt_validated_low_evidence.csv.gz", index=False)
    direct_conflicts.to_csv(
        output_dir / "direct_species_trait_conflicts_excluded_from_rules.csv.gz",
        index=False,
    )
    final_strict_filled = int(shared_report["coverage_after"]["filled_species_axis"])
    report = {
        "contract": "open_web_manual_audit_and_gain_v1",
        "completed_at_utc": _now(),
        "audited_pages_or_candidates": len(reviewed),
        "accepted_correct": len(accepted_audit),
        "precision": precision,
        "species_identity_accuracy": (
            float(reviewed["species_identity_correct"].map(_bool).mean()) if len(reviewed) else None
        ),
        "value_accuracy": (
            float(reviewed["value_correct"].map(_bool).mean()) if len(reviewed) else None
        ),
        "provenance_completeness": (
            float(reviewed["provenance_complete"].map(_bool).mean()) if len(reviewed) else None
        ),
        "cultivar_contamination_rate": cultivar_rate,
        "by_trait": by_trait,
        "minimum_audited_per_production_trait": MIN_AUDITED_PER_PRODUCTION_TRAIT,
        "production_approved_traits": production_approved_traits,
        "by_source_tier": grouped_precision["source_tier"],
        "by_source_type": grouped_precision["source_type"],
        "by_domain": grouped_precision["domain"],
        "by_language": grouped_precision["language"],
        "new_direct_species_trait": len(new_trait),
        "new_direct_species_axis": len(new_axes),
        "validated_low_upgraded_to_direct": len(low_upgrades),
        "old_low_conflicting_with_new_direct": len(low_conflicts),
        "validated_low_rebuild": low_rebuild,
        "shared_all_evidence": shared_report,
        "strict_filled_species_axis_before": int(coverage["after_quality"].ne("").sum()),
        "strict_filled_species_axis_after_direct_only": int(
            coverage["after_quality"].ne("").sum() + len(new_axes)
        ),
        "strict_filled_species_axis_after": final_strict_filled,
        "strict_coverage_before": float(coverage["after_quality"].ne("").mean()),
        "strict_coverage_after_direct_only": float(
            (coverage["after_quality"].ne("").sum() + len(new_axes)) / 318_885
        ),
        "strict_coverage_after": float(final_strict_filled / 318_885),
        "strict_coverage_increment": float(
            (final_strict_filled - coverage["after_quality"].ne("").sum()) / 318_885
        ),
        "production_approved": bool(
            len(reviewed) >= 100
            and precision is not None
            and precision >= 0.90
            and (cultivar_rate or 0.0) <= 0.02
            and bool(production_approved_traits)
        ),
    }
    (output_dir / "manual_audit_gain_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


@app.command("audit-baseline")
def audit_baseline_command(
    baseline_dir: Path = typer.Option(..., exists=True, file_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    typer.echo(json.dumps(baseline_audit(baseline_dir, output_dir=output_dir), sort_keys=True))


@app.command("build-queue")
def build_queue_command(
    baseline_dir: Path = typer.Option(..., exists=True, file_okay=False),
    master_csv: Path = typer.Option(..., exists=True),
    output_csv: Path = typer.Option(...),
    pilot_species: int = typer.Option(300, min=1, max=1000),
) -> None:
    _, report = build_priority_queue(
        baseline_dir, master_csv, output_csv=output_csv, pilot_species=pilot_species
    )
    typer.echo(json.dumps(report, ensure_ascii=False, sort_keys=True))


@app.command("registry-pilot")
def registry_pilot_command(
    baseline_dir: Path = typer.Option(..., exists=True, file_okay=False),
    master_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_yaml: Path = typer.Option(Path("config/open_web_trait_acquisition.yml"), exists=True),
    registry_yaml: Path = typer.Option(Path("config/open_web_domain_registry.yml"), exists=True),
    domain: str = typer.Option("mikawanoyasou.org"),
    pilot_species: int = typer.Option(300, min=1, max=1000),
    max_pages: int = typer.Option(600, min=1, max=5000),
    audit_pages: int = typer.Option(100, min=100, max=500),
    pause_seconds: float = typer.Option(0.25, min=0, max=30),
) -> None:
    report = run_registry_pilot(
        baseline_dir=baseline_dir,
        master_csv=master_csv,
        config_yaml=config_yaml,
        registry_yaml=registry_yaml,
        output_dir=output_dir,
        domain=domain,
        pilot_species=pilot_species,
        max_pages=max_pages,
        audit_pages=audit_pages,
        pause_seconds=pause_seconds,
    )
    typer.echo(json.dumps(report, ensure_ascii=False, sort_keys=True))


@app.command("apply-audit")
def apply_audit_command(
    baseline_dir: Path = typer.Option(..., exists=True, file_okay=False),
    candidates_csv: Path = typer.Option(..., exists=True),
    audit_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    report = apply_manual_audit(
        baseline_dir=baseline_dir,
        candidates_csv=candidates_csv,
        audit_csv=audit_csv,
        output_dir=output_dir,
    )
    typer.echo(json.dumps(report, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    app()
