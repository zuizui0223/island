"""Extract conservative, review-pending traits from acquired WFO flora text."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.open_web_common import SEARCH_TRAITS, STRICT_TRAIT_AXIS
from island_v2.open_web_evidence import Page, extract_rules, load_yaml


LANGUAGE_PRIORITY = {"en": 0, "EN": 0, "fre": 1, "es": 2, "pt": 3, "tr": 4}


def text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def sha(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def now_utc() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def load_baseline(
    coverage_path: Path,
    direct_traits_path: Path,
) -> tuple[dict[tuple[str, str], str], set[tuple[str, str]]]:
    coverage = pd.read_csv(
        coverage_path,
        usecols=["accepted_species", "axis", "quality"],
        dtype=str,
    ).fillna("")
    axis_quality = {
        (row.accepted_species, row.axis): row.quality
        for row in coverage.itertuples(index=False)
    }
    direct = pd.read_csv(
        direct_traits_path,
        usecols=["accepted_species", "trait_name", "resolution_status"],
        dtype=str,
    ).fillna("")
    resolved_traits = set(
        direct.loc[
            direct["resolution_status"].eq("resolved"),
            ["accepted_species", "trait_name"],
        ].itertuples(index=False, name=None)
    )
    return axis_quality, resolved_traits


def candidate_rows(
    descriptions: pd.DataFrame,
    config: dict[str, Any],
    *,
    axis_quality: dict[tuple[str, str], str],
    resolved_traits: set[tuple[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    no_hit_counter: Counter[tuple[str, str]] = Counter()
    for source in descriptions.itertuples(index=False):
        full_text = text(source.description)
        page = Page(
            requested_url=text(source.source_url),
            final_url=text(source.source_url),
            status_code=200,
            content_type="text/plain; source=official_flora_archive",
            title=f"{source.accepted_species} species treatment",
            text=full_text,
            language=text(source.language),
            retrieved_at_utc=text(source.retrieved_at_utc),
            content_sha256=text(source.description_sha256),
        )
        for trait_name in SEARCH_TRAITS:
            extracted = extract_rules(
                page,
                trait_name=trait_name,
                config=config,
            )
            if not extracted:
                no_hit_counter[(text(source.provider), trait_name)] += 1
                continue
            for raw_value, normalized_value, excerpt in extracted:
                axis = STRICT_TRAIT_AXIS.get(trait_name, "")
                baseline_quality = (
                    axis_quality.get((text(source.accepted_species), axis), "")
                    if axis
                    else ""
                )
                baseline_trait_resolved = (
                    text(source.accepted_species),
                    trait_name,
                ) in resolved_traits
                normalized_excerpt = " ".join(
                    "".join(
                        character.casefold()
                        if character.isalnum() or character.isspace()
                        else " "
                        for character in excerpt
                    ).split()
                )
                source_lineage = text(source.source_lineage)
                candidate_key = "|".join(
                    (
                        text(source.accepted_species),
                        trait_name,
                        text(normalized_value),
                        source_lineage,
                        sha(normalized_excerpt),
                    )
                )
                rows.append(
                    {
                        "candidate_id": f"WFOB-{sha(candidate_key)[:20]}",
                        "accepted_species": text(source.accepted_species),
                        "trait_name": trait_name,
                        "axis": axis,
                        "raw_value": text(raw_value),
                        "normalized_value": text(normalized_value),
                        "source_id": (
                            f"WFO:{text(source.provider)}:{text(source.wfo_id)}"
                        ),
                        "source_url": text(source.source_url),
                        "source_title": (
                            f"{text(source.accepted_species)} species treatment "
                            f"({text(source.provider).replace('_', ' ')})"
                        ),
                        "domain": "files.worldfloraonline.org",
                        "source_type": "official_flora_archive",
                        "source_tier": "A",
                        "proposed_quality": "high",
                        "source_lineage": source_lineage,
                        "exact_supporting_quote": text(excerpt),
                        "name_match_method": text(source.name_match_method),
                        "matched_wfo_name": text(source.matched_wfo_name),
                        "wfo_name_status": text(source.wfo_name_status),
                        "cultivar_status": "wild_or_flora_treatment",
                        "language": text(source.language),
                        "provider": text(source.provider),
                        "archive_url": text(source.archive_url),
                        "license": text(source.license),
                        "description_sha256": text(source.description_sha256),
                        "normalized_excerpt_sha256": sha(normalized_excerpt),
                        "retrieval_date": text(source.retrieved_at_utc)[:10],
                        "extraction_rule_or_model_version": (
                            "rules_only_multilingual_v1"
                        ),
                        "baseline_trait_status": (
                            "resolved_direct"
                            if baseline_trait_resolved
                            else "absent_direct"
                        ),
                        "baseline_axis_quality": baseline_quality or "unresolved",
                        "candidate_effect": (
                            "new_species_trait_and_axis"
                            if not baseline_trait_resolved
                            and axis
                            and not baseline_quality
                            else "new_species_trait"
                            if not baseline_trait_resolved
                            else "additional_direct_lineage"
                        ),
                        "model_decision": (
                            "source_backed_rule_extracted_unreviewed"
                        ),
                        "review_status": "unreviewed",
                        "decision_reason": (
                            "Official flora treatment, exact WFO species/synonym "
                            "identity, and exact source quote; manual context and "
                            "conflict review still required."
                        ),
                    }
                )

    raw = pd.DataFrame(rows)
    if raw.empty:
        return raw, pd.DataFrame()
    raw["_language_priority"] = raw["language"].map(
        lambda value: LANGUAGE_PRIORITY.get(text(value), 99)
    )
    raw = raw.sort_values(
        [
            "accepted_species",
            "trait_name",
            "normalized_value",
            "source_lineage",
            "_language_priority",
            "normalized_excerpt_sha256",
        ]
    )
    dedupe_columns = [
        "accepted_species",
        "trait_name",
        "normalized_value",
        "source_lineage",
    ]
    duplicate_mask = raw.duplicated(dedupe_columns, keep="first")
    duplicates = raw.loc[duplicate_mask].copy()
    if not duplicates.empty:
        duplicates["dedupe_reason"] = (
            "same_species_trait_value_and_provider_treatment_lineage"
        )
    deduped = raw.loc[~duplicate_mask].drop(
        columns=["_language_priority"],
        errors="ignore",
    )
    duplicates = duplicates.drop(columns=["_language_priority"], errors="ignore")
    return deduped.reset_index(drop=True), duplicates.reset_index(drop=True)


def summarize(
    descriptions: pd.DataFrame,
    candidates: pd.DataFrame,
    duplicates: pd.DataFrame,
) -> dict[str, Any]:
    strict = candidates.loc[candidates["axis"].ne("")].copy()
    novel_trait = candidates.loc[
        candidates["baseline_trait_status"].eq("absent_direct")
    ]
    novel_axis = strict.loc[strict["baseline_axis_quality"].eq("unresolved")]
    per_trait = (
        candidates.groupby("trait_name", sort=True)
        .agg(
            candidate_rows=("candidate_id", "size"),
            species=("accepted_species", "nunique"),
            providers=("provider", "nunique"),
        )
        .reset_index()
        .to_dict("records")
        if not candidates.empty
        else []
    )
    per_provider = (
        candidates.groupby("provider", sort=True)
        .agg(
            candidate_rows=("candidate_id", "size"),
            species=("accepted_species", "nunique"),
            species_traits=(
                "candidate_id",
                lambda series: candidates.loc[
                    series.index, ["accepted_species", "trait_name"]
                ].drop_duplicates().shape[0],
            ),
        )
        .reset_index()
        .to_dict("records")
        if not candidates.empty
        else []
    )
    return {
        "contract": "wfo_official_bulk_trait_candidate_extraction_v1",
        "created_at_utc": now_utc(),
        "formal_promotion_status": "not_promoted_review_pending",
        "acquired_descriptions": {
            "rows": len(descriptions),
            "species": int(descriptions["accepted_species"].nunique()),
            "providers": int(descriptions["provider"].nunique()),
        },
        "candidates": {
            "rows_after_lineage_deduplication": len(candidates),
            "lineage_translation_duplicates_removed": len(duplicates),
            "species": int(candidates["accepted_species"].nunique()),
            "species_traits": int(
                candidates[["accepted_species", "trait_name"]]
                .drop_duplicates()
                .shape[0]
            ),
            "strict_species_axes": int(
                strict[["accepted_species", "axis"]].drop_duplicates().shape[0]
            ),
            "formal_novel_species_traits_pending_review": int(
                novel_trait[["accepted_species", "trait_name"]]
                .drop_duplicates()
                .shape[0]
            ),
            "formal_novel_species_axes_pending_review": int(
                novel_axis[["accepted_species", "axis"]]
                .drop_duplicates()
                .shape[0]
            ),
            "non_axis_rows": int(candidates["axis"].eq("").sum()),
        },
        "per_trait": per_trait,
        "per_provider": per_provider,
        "query_accounting": {
            "search_api_queries": 0,
            "page_requests": 0,
            "official_archive_downloads": 30,
            "official_consolidated_exports": 1,
            "cost_usd": 0.0,
        },
        "policy": {
            "snippets_used": False,
            "family_inference": False,
            "global_fallback": False,
            "pollen_vector_or_reward_mixed_into_strict_axes": False,
            "manual_review_required_for_promotion": True,
        },
    }


def write_file_manifest(output_dir: Path) -> None:
    rows = []
    for path in sorted(output_dir.iterdir()):
        if not path.is_file() or path.name == "file_manifest.sha256":
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append(f"{digest}  {path.name}")
    (output_dir / "file_manifest.sha256").write_text(
        "\n".join(rows) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--descriptions", type=Path, nargs="+", required=True)
    parser.add_argument("--coverage", type=Path, required=True)
    parser.add_argument("--direct-traits", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()

    descriptions = pd.concat(
        [
            pd.read_csv(description_path, dtype=str).fillna("")
            for description_path in args.descriptions
        ],
        ignore_index=True,
    )
    axis_quality, resolved_traits = load_baseline(
        args.coverage,
        args.direct_traits,
    )
    config = load_yaml(args.config)
    workers = max(1, args.workers)
    if workers == 1:
        candidates, duplicates = candidate_rows(
            descriptions,
            config,
            axis_quality=axis_quality,
            resolved_traits=resolved_traits,
        )
    else:
        descriptions["_partition"] = descriptions["source_lineage"].map(
            lambda value: int(sha(text(value))[:12], 16) % workers
        )
        futures = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
            for partition in range(workers):
                chunk = descriptions.loc[
                    descriptions["_partition"].eq(partition)
                ].drop(columns=["_partition"])
                futures.append(
                    executor.submit(
                        candidate_rows,
                        chunk,
                        config,
                        axis_quality=axis_quality,
                        resolved_traits=resolved_traits,
                    )
                )
            results = [future.result() for future in futures]
        candidates = pd.concat(
            [result[0] for result in results if not result[0].empty],
            ignore_index=True,
        )
        duplicates = pd.concat(
            [result[1] for result in results if not result[1].empty],
            ignore_index=True,
        )
        descriptions = descriptions.drop(columns=["_partition"])
    args.output.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(
        args.output / "llm_evidence_candidates.csv.gz",
        index=False,
        compression="gzip",
    )
    duplicates.to_csv(
        args.output / "source_lineage_translation_duplicates.csv.gz",
        index=False,
        compression="gzip",
    )
    summary = summarize(descriptions, candidates, duplicates)
    (args.output / "run_manifest.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    write_file_manifest(args.output)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
