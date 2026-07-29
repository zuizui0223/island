"""Open-Web promotion through the PR #131 all-evidence implementation.

This module is intentionally orchestration-only.  It does not implement a
second genus-inference algorithm.  Reviewed Web evidence is converted to the
common evidence schema, then the source-lineage, direct-conflict, trait-rule,
Validated Low, coverage, and acquisition-queue functions from
``all_evidence_trait_audit`` are called directly.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.all_evidence_trait_audit import (
    acquisition_queue,
    apply_genus_rules,
    build_rule_audit,
    canonical_trait_name,
    compare_old_low,
    coverage_snapshot,
    dedupe_direct_lineages,
    load_ontology,
    normalise_state_set,
    resolve_direct_cells,
    species_axis_coverage,
    trait_axis,
)
from island_v2.integrated_trait_coverage import (
    AXES,
    AXIS_TRAITS,
    DIRECT_SCOPES,
    EVIDENCE_COLUMNS,
    QUALITY_RANK,
)

STRICT_AXIS_TRAITS = {
    axis: tuple(sorted(trait for trait in traits if trait not in {"flower_color", "flower_shape"}))
    for axis, traits in AXIS_TRAITS.items()
}
STRICT_TRAIT_AXIS = {trait: axis for axis, traits in STRICT_AXIS_TRAITS.items() for trait in traits}
NON_AXIS_TRAITS = ("pollen_vector_mode", "reward_type")
SEARCH_TRAITS = (*STRICT_TRAIT_AXIS, *NON_AXIS_TRAITS)

MIN_REVIEWED_PAGES = 200
MIN_REVIEWED_DOMAINS = 5
MIN_REVIEWED_PER_SCOPE = 10
MIN_PRECISION = 0.90
MAX_CULTIVAR_CONTAMINATION = 0.02


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _bool(value: object) -> bool:
    return _text(value).casefold() in {"1", "true", "yes", "y"}


def reviewed_audit_metrics(audit: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Calculate fail-closed precision and production approval by domain x trait."""

    required = {
        "candidate_id",
        "domain",
        "trait_name",
        "decision",
        "species_identity_correct",
        "value_correct",
        "provenance_complete",
        "cultivar_contamination",
        "reviewer",
        "reviewed_at_utc",
    }
    missing = required.difference(audit.columns)
    if missing:
        raise ValueError(f"audit is missing columns: {sorted(missing)}")
    reviewed = audit.loc[audit["decision"].str.casefold().isin({"accept", "reject"})].copy()
    reviewed["_accepted_correct"] = (
        reviewed["decision"].str.casefold().eq("accept")
        & reviewed["species_identity_correct"].map(_bool)
        & reviewed["value_correct"].map(_bool)
        & reviewed["provenance_complete"].map(_bool)
        & ~reviewed["cultivar_contamination"].map(_bool)
    )
    reviewed["_cultivar"] = reviewed["cultivar_contamination"].map(_bool)

    scopes: list[dict[str, Any]] = []
    for (domain, trait), group in reviewed.groupby(["domain", "trait_name"], sort=True):
        accepted_correct = int(group["_accepted_correct"].sum())
        cultivar_rate = float(group["_cultivar"].mean()) if len(group) else None
        precision = accepted_correct / len(group) if len(group) else None
        approved = bool(
            len(group) >= MIN_REVIEWED_PER_SCOPE
            and precision is not None
            and precision >= MIN_PRECISION
            and cultivar_rate is not None
            and cultivar_rate <= MAX_CULTIVAR_CONTAMINATION
        )
        scopes.append(
            {
                "domain": domain,
                "trait_name": trait,
                "reviewed": len(group),
                "accepted_correct": accepted_correct,
                "precision": precision,
                "cultivar_contamination_rate": cultivar_rate,
                "production_approved": approved,
            }
        )
    scope_table = pd.DataFrame(scopes)
    accepted_correct = int(reviewed["_accepted_correct"].sum())
    precision = accepted_correct / len(reviewed) if len(reviewed) else None
    cultivar_rate = float(reviewed["_cultivar"].mean()) if len(reviewed) else None
    reviewed_domains = int(reviewed["domain"].replace("", pd.NA).nunique())
    approved_scope_count = (
        int(scope_table["production_approved"].sum()) if not scope_table.empty else 0
    )
    summary = {
        "precision_contract": "accepted_correct/all_reviewed",
        "reviewed_pages_or_candidates": len(reviewed),
        "reviewed_domains": reviewed_domains,
        "accepted_correct": accepted_correct,
        "precision": precision,
        "cultivar_contamination_rate": cultivar_rate,
        "minimum_reviewed_pages": MIN_REVIEWED_PAGES,
        "minimum_reviewed_domains": MIN_REVIEWED_DOMAINS,
        "minimum_reviewed_per_domain_trait": MIN_REVIEWED_PER_SCOPE,
        "minimum_precision": MIN_PRECISION,
        "maximum_cultivar_contamination": MAX_CULTIVAR_CONTAMINATION,
        "approved_domain_trait_scopes": approved_scope_count,
        "pilot_gate_passed": bool(
            len(reviewed) >= MIN_REVIEWED_PAGES
            and reviewed_domains >= MIN_REVIEWED_DOMAINS
            and precision is not None
            and precision >= MIN_PRECISION
            and cultivar_rate is not None
            and cultivar_rate <= MAX_CULTIVAR_CONTAMINATION
            and approved_scope_count
        ),
    }
    return scope_table, summary


def accepted_review_ids(
    audit: pd.DataFrame,
    approved_scopes: pd.DataFrame,
) -> set[str]:
    """Return only correct accepted reviews in production-approved scopes."""

    if approved_scopes.empty:
        return set()
    approved = approved_scopes.loc[
        approved_scopes["production_approved"].astype(bool),
        ["domain", "trait_name"],
    ].drop_duplicates()
    reviewed = audit.merge(
        approved,
        on=["domain", "trait_name"],
        how="inner",
        validate="many_to_one",
    )
    good = (
        reviewed["decision"].str.casefold().eq("accept")
        & reviewed["species_identity_correct"].map(_bool)
        & reviewed["value_correct"].map(_bool)
        & reviewed["provenance_complete"].map(_bool)
        & ~reviewed["cultivar_contamination"].map(_bool)
    )
    return set(reviewed.loc[good, "candidate_id"].map(_text))


def coverage_change_counts(
    before_axes: set[tuple[str, str]],
    after_axes: set[tuple[str, str]],
) -> dict[str, Any]:
    """Separate gross fills from invalidations so ``net_change`` is not overstated."""

    gross_gain = after_axes - before_axes
    loss = before_axes - after_axes
    by_axis: dict[str, dict[str, int]] = {}
    for axis in AXES:
        gained = sum(1 for _, value in gross_gain if value == axis)
        lost = sum(1 for _, value in loss if value == axis)
        by_axis[axis] = {
            "gross_gain": gained,
            "loss": lost,
            "net_change": gained - lost,
        }
    return {
        "gross_gain": len(gross_gain),
        "loss": len(loss),
        "net_change": len(after_axes) - len(before_axes),
        "by_axis": by_axis,
    }


def compare_incremental_validated_low(
    before_low: pd.DataFrame,
    after_low: pd.DataFrame,
    incremental_direct: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Classify the incremental pilot's trait-specific effect on Validated Low."""

    def states_by_pair(frame: pd.DataFrame) -> dict[tuple[str, str], set[str]]:
        if frame.empty:
            return {}
        grouped: dict[tuple[str, str], set[str]] = {}
        for row in frame.fillna("").itertuples(index=False):
            key = (_text(row.accepted_species), _text(row.trait_name))
            grouped.setdefault(key, set()).add(_text(row.state_set))
        return grouped

    before = states_by_pair(before_low)
    after = states_by_pair(after_low)
    incremental_direct_pairs = (
        set(
            zip(
                incremental_direct["accepted_species"].map(_text),
                incremental_direct["trait_name"].map(_text),
                strict=True,
            )
        )
        if not incremental_direct.empty
        else set()
    )
    rows: list[dict[str, str]] = []
    for accepted_species, trait_name in sorted(set(before) | set(after)):
        key = (accepted_species, trait_name)
        before_states = before.get(key, set())
        after_states = after.get(key, set())
        if before_states and key in incremental_direct_pairs:
            status = "upgraded_to_incremental_direct"
        elif not before_states and after_states:
            status = "added"
        elif before_states and not after_states:
            status = "invalidated"
        elif before_states != after_states:
            status = "inferred_value_changed"
        else:
            status = "retained"
        rows.append(
            {
                "accepted_species": accepted_species,
                "trait_name": trait_name,
                "axis": trait_axis(trait_name),
                "before_state_sets": json.dumps(
                    sorted(before_states),
                    ensure_ascii=False,
                    separators=(",", ":"),
                ),
                "after_state_sets": json.dumps(
                    sorted(after_states),
                    ensure_ascii=False,
                    separators=(",", ":"),
                ),
                "pilot_change_status": status,
            }
        )
    detail = pd.DataFrame(
        rows,
        columns=[
            "accepted_species",
            "trait_name",
            "axis",
            "before_state_sets",
            "after_state_sets",
            "pilot_change_status",
        ],
    )
    counts = (
        detail["pilot_change_status"].value_counts().to_dict()
        if not detail.empty
        else {}
    )
    summary = {
        "added": int(counts.get("added", 0)),
        "invalidated": int(counts.get("invalidated", 0)),
        "upgraded_to_direct": int(
            counts.get("upgraded_to_incremental_direct", 0)
        ),
        "inferred_value_changed": int(counts.get("inferred_value_changed", 0)),
        "retained": int(counts.get("retained", 0)),
    }
    return detail, summary


def promoted_to_common_evidence(promoted: pd.DataFrame) -> pd.DataFrame:
    """Convert reviewed Web rows to the common species-direct evidence schema."""

    rows: list[dict[str, str]] = []
    for record in promoted.fillna("").to_dict("records"):
        trait = canonical_trait_name(record.get("trait_name"))
        axis = trait_axis(trait)
        if not axis:
            continue
        quality = _text(
            record.get("evidence_quality") or record.get("confidence") or record.get("quality")
        ).casefold()
        if quality not in {"high", "medium"}:
            continue
        rows.append(
            {
                "accepted_species": _text(record.get("accepted_species")),
                "axis": axis,
                "trait_name": trait,
                "normalized_value": _text(record.get("normalized_value")),
                "quality": quality,
                "source_group": "latest_public_web",
                "source_provider": _text(record.get("source_provider") or record.get("domain")),
                "source_url": _text(record.get("source_url")),
                "source_record_id": _text(
                    record.get("candidate_id") or record.get("source_record_id")
                ),
                "source_citation": _text(record.get("source_citation") or record.get("page_title")),
                "source_excerpt": _text(
                    record.get("source_excerpt") or record.get("supporting_excerpt")
                ),
                "evidence_scope": _text(record.get("evidence_scope") or "species_direct"),
                "name_match_method": _text(record.get("name_match_method")),
                "source_lineage": _text(record.get("source_lineage")),
                "lineage_method": _text(
                    record.get("source_lineage_method")
                    or record.get("lineage_method")
                    or "artifact_source_lineage"
                ),
                "source_run_id": _text(record.get("source_run_id")),
                "source_artifact": _text(record.get("source_artifact")),
                "source_file": _text(record.get("source_file")),
                "acceptance_contract": "reviewed_open_web_species_direct_v2",
            }
        )
    return pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)


def _baseline_direct(evidence: pd.DataFrame) -> pd.DataFrame:
    data = evidence.copy().fillna("")
    for column in EVIDENCE_COLUMNS:
        if column not in data:
            data[column] = ""
    data["quality"] = data["quality"].str.casefold()
    data["trait_name"] = data["trait_name"].map(canonical_trait_name)
    data["axis"] = data["trait_name"].map(trait_axis)
    return data.loc[
        data["quality"].isin({"high", "medium"})
        & data["evidence_scope"].isin(DIRECT_SCOPES)
        & data["axis"].ne(""),
        EVIDENCE_COLUMNS,
    ].copy()


def _old_low(
    evidence: pd.DataFrame,
    ontology: dict[str, set[str]],
) -> pd.DataFrame:
    data = evidence.copy().fillna("")
    if "quality" not in data:
        return pd.DataFrame()
    data = data.loc[data["quality"].str.casefold().eq("low")].copy()
    if data.empty:
        return data
    data["trait_name"] = data["trait_name"].map(canonical_trait_name)
    data["axis"] = data["trait_name"].map(trait_axis)
    data["genus"] = data["accepted_species"].str.split().str[0]
    states = [
        normalise_state_set(trait, value, ontology)[0]
        for trait, value in zip(
            data["trait_name"],
            data["normalized_value"],
            strict=True,
        )
    ]
    data["state_set"] = [
        json.dumps(list(value), ensure_ascii=False, separators=(",", ":")) for value in states
    ]
    return data.loc[data["axis"].ne("") & data["state_set"].ne("[]")].drop_duplicates(
        ["accepted_species", "trait_name", "state_set"]
    )


def _master_for_common(
    coverage: pd.DataFrame,
    master: pd.DataFrame | None,
) -> pd.DataFrame:
    universe = set(coverage["accepted_species"].map(_text))
    if master is None:
        data = pd.DataFrame({"accepted_species": sorted(universe)})
    else:
        data = master.copy().fillna("")
        data = data.loc[data["accepted_species"].isin(universe)].copy()
    data = data.drop_duplicates("accepted_species")
    if len(data) != 106_295:
        raise ValueError(f"master has {len(data)} species, expected 106295")
    data["genus"] = data.get(
        "genus",
        data["accepted_species"].str.split().str[0],
    )
    data["genus"] = data["genus"].where(
        data["genus"].ne(""),
        data["accepted_species"].str.split().str[0],
    )
    for column in ("n_islands", "n_records"):
        if column not in data:
            data[column] = "0"
    return data


def rebuild_with_common_all_evidence(
    *,
    coverage: pd.DataFrame,
    evidence: pd.DataFrame,
    promoted: pd.DataFrame,
    prior_promoted: pd.DataFrame | None = None,
    master: pd.DataFrame | None = None,
    ontology_yaml: Path = Path("config/trait_ontology.yml"),
) -> tuple[dict[str, Any], dict[str, pd.DataFrame]]:
    """Run PR #131's shared trait-specific implementation on reviewed evidence."""

    if len(coverage) != 318_885:
        raise ValueError(f"coverage has {len(coverage)} rows, expected 318885")
    common_master = _master_for_common(coverage, master)
    ontology = load_ontology(ontology_yaml)
    base_raw = _baseline_direct(evidence)
    prior_raw = (
        promoted_to_common_evidence(prior_promoted)
        if prior_promoted is not None
        else pd.DataFrame(columns=EVIDENCE_COLUMNS)
    )
    promoted_raw = promoted_to_common_evidence(promoted)

    baseline_raw = pd.concat([base_raw, prior_raw], ignore_index=True).fillna("")
    base_lineages, _ = dedupe_direct_lineages(baseline_raw, ontology)
    base_cells, _ = resolve_direct_cells(base_lineages)
    if not base_cells.empty:
        base_cells["genus"] = base_cells["accepted_species"].str.split().str[0]
    if not base_lineages.empty:
        base_lineages["genus"] = base_lineages["accepted_species"].str.split().str[0]

    raw = pd.concat([baseline_raw, promoted_raw], ignore_index=True).fillna("")
    lineages, lineage_duplicates = dedupe_direct_lineages(raw, ontology)
    direct_cells, cell_audit = resolve_direct_cells(lineages)
    if not direct_cells.empty:
        direct_cells["genus"] = direct_cells["accepted_species"].str.split().str[0]
    if not lineages.empty:
        lineages["genus"] = lineages["accepted_species"].str.split().str[0]

    old_low = _old_low(evidence, ontology)
    base_rules = build_rule_audit(base_cells, base_lineages, old_low)
    base_rebuilt_low = apply_genus_rules(
        common_master,
        base_cells,
        base_rules,
        "current_min3",
    )
    rules = build_rule_audit(direct_cells, lineages, old_low)
    rebuilt_low = apply_genus_rules(
        common_master,
        direct_cells,
        rules,
        "current_min3",
    )
    before = species_axis_coverage(common_master, base_cells, base_rebuilt_low)
    after = species_axis_coverage(common_master, direct_cells, rebuilt_low)
    comparison, comparison_summary = compare_old_low(
        old_low,
        pd.DataFrame(),
        rebuilt_low,
        direct_cells,
        rules,
    )
    queue = acquisition_queue(common_master, after, rules)

    base_direct_pairs = set(zip(base_cells["accepted_species"], base_cells["trait_name"]))
    direct_pairs = set(zip(direct_cells["accepted_species"], direct_cells["trait_name"]))
    incremental_direct = direct_cells.loc[
        [
            (species, trait) not in base_direct_pairs
            for species, trait in zip(
                direct_cells["accepted_species"],
                direct_cells["trait_name"],
                strict=True,
            )
        ]
    ].copy()
    base_axes = set(
        zip(
            before.loc[before["quality"].isin(QUALITY_RANK), "accepted_species"],
            before.loc[before["quality"].isin(QUALITY_RANK), "axis"],
        )
    )
    after_axes = set(
        zip(
            after.loc[after["quality"].isin(QUALITY_RANK), "accepted_species"],
            after.loc[after["quality"].isin(QUALITY_RANK), "axis"],
        )
    )
    base_direct_axes = set(zip(base_cells["accepted_species"], base_cells["axis"]))
    direct_axes = set(zip(direct_cells["accepted_species"], direct_cells["axis"]))
    old_low_keys = set(zip(old_low["accepted_species"], old_low["trait_name"]))
    new_low_keys = set(zip(rebuilt_low["accepted_species"], rebuilt_low["trait_name"]))
    direct_keys = set(zip(direct_cells["accepted_species"], direct_cells["trait_name"]))
    pilot_low_comparison, pilot_low_summary = compare_incremental_validated_low(
        base_rebuilt_low,
        rebuilt_low,
        incremental_direct,
    )
    coverage_change = coverage_change_counts(base_axes, after_axes)
    report = {
        "contract": "open_web_shared_all_evidence_rebuild_v1",
        "implementation": "island_v2.all_evidence_trait_audit",
        "trait_specific_join": "genus_x_trait_name",
        "axis_join_for_inference": False,
        "strict_axis_traits": STRICT_AXIS_TRAITS,
        "separate_non_axis_traits": list(NON_AXIS_TRAITS),
        "denominator": {
            "species": 106_295,
            "axes": len(AXES),
            "species_axis": 318_885,
        },
        "new_direct_species_trait": len(direct_pairs - base_direct_pairs),
        "new_direct_species_axis": len(direct_axes - base_direct_axes),
        "prior_formal_public_web_evidence_rows": len(prior_raw),
        "lineage_duplicates_before_after": {
            "raw_rows": len(raw),
            "deduplicated_lineages": len(lineages),
            "duplicate_rows": len(lineage_duplicates),
        },
        "direct_conflict_classification": {
            str(key): int(value)
            for key, value in cell_audit["classification"].value_counts().items()
        }
        if not cell_audit.empty
        else {},
        "validated_low": {
            **comparison_summary,
            "added": len(new_low_keys - old_low_keys),
            "invalidated": len(old_low_keys - new_low_keys - direct_keys),
            "upgraded_to_direct": len(old_low_keys & direct_keys),
        },
        "pilot_validated_low_change": pilot_low_summary,
        "coverage_before": coverage_snapshot(before),
        "coverage_after": coverage_snapshot(after),
        "coverage_change_species_axis": coverage_change,
        "coverage_increment_species_axis": coverage_change["net_change"],
        "family_inference": False,
        "global_fallback": False,
    }
    return report, {
        "formal_evidence": promoted_raw,
        "lineages": lineages,
        "lineage_duplicates": lineage_duplicates,
        "direct_cells": direct_cells,
        "direct_conflicts": cell_audit,
        "rules": rules,
        "prior_rules": base_rules,
        "rebuilt_low": rebuilt_low,
        "prior_rebuilt_low": base_rebuilt_low,
        "pilot_low_comparison": pilot_low_comparison,
        "old_low_comparison": comparison,
        "strict_coverage": after,
        "next_queue": queue,
    }
