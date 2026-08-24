"""Recover trait-specific morphology measurements from the frozen WFO package.

The original WFO rules emitted some corolla-tube and floral-part dimensions as
``flower_size_class`` candidates.  This recovery pass starts only from the
already acquired, identity-gated WFO evidence, attaches the exact extracted
measurement to its nearest botanical subject, and keeps whole-flower/corolla
size separate from corolla-tube depth.  Other floral-part measurements fail
closed.

The script is deliberately network-free.  It writes a deterministic 200-row
audit queue that excludes the audit used for the first WFO package.  A separate
review decision file is required before it will write a promotion audit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import unicodedata
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS
from island_v2.trait_measurements import derive_class, load_rules

RUN_ID = "wfo-official-bulk-morphology-recovery-20260808"
AUDIT_REVIEWER = "codex_gpt5_wfo_morphology_quote_audit_20260808"
AUDIT_SIZE = 200

_WHOLE_SUBJECT = re.compile(
    r"\b(?:flower(?![- ](?:bearing|stalk|bract|pedicel|peduncle|head|zone))s?|"
    r"corollas?|perianths?|flor(?:es)?|corola)\b",
    re.IGNORECASE,
)
_TUBE_SUBJECT = re.compile(
    r"\b(?:(?:corolla|floral|perianth)\s+tube|"
    r"tube\s+of\s+(?:the\s+)?corolla|tube|"
    r"tubo(?:\s+de\s+la)?\s+corola|tubus\s+corollae)\b",
    re.IGNORECASE,
)
_NON_TARGET_SUBJECT = re.compile(
    r"\b(?:petals?|tepals?|sepals?|calyx|calycine|bracts?|bracteoles?|"
    r"bractlets?|stalks?|"
    r"pedicels?|peduncles?|pistils?|rachis|rachillae?|branches?|branched|"
    r"inflorescences?|heads?|hypanthium|lobes?|lips?|labellum|"
    r"standards?|banners?|keels?|wings?|filaments?|styles?|stigmas?|"
    r"stamens?|anthers?|ovar(?:y|ies)|glumes?|lemmas?|spikelets?|"
    r"throats?|limbs?|flower\s+zone|ligulate\s+extension|c[aá]liz|"
    r"s[eé]palas?|p[eé]talas?|br[aá]cteas?|ped[uú]nculos?|pedicelos?)\b",
    re.IGNORECASE,
)
_AMBIGUOUS_COMPARISON = re.compile(
    r"(?:bigger|larger|longer|more|smaller|shorter|less)\s+than\s*$|"
    r"(?:larger|bigger)\s+or\s+equal\s+to\s*$|"
    r"(?:at\s+least|minimum|above|over|below|under)\s*$|[><≥≤]\s*$",
    re.IGNORECASE,
)
_BAD_AFTER_MEASUREMENT = re.compile(
    r"^\s*(?:long|wide)?\s*(?:pedunculate\s+)?(?:heads?|bracts?|"
    r"pedicels?|peduncles?|rachis)\b|^.{0,45}\b(?:in\s+fruit|"
    r"or\s+(?:perhaps\s+)?longer|longest\s+branch)\b|"
    r"^.{0,120}\b(?:stalk\s+of\s+the\s+inflorescence|"
    r"with\s+\d+(?:-\d+)?\s+flowers)\b",
    re.IGNORECASE,
)
_AMBIGUOUS_WHOLE_SHAPE = re.compile(
    r"^\s*flowers?\s*(?:,|are)?[^.;]{0,35}\b(?:ovate|obovate|"
    r"oblanceolate|spatulate|lanceolate|linear)\b",
    re.IGNORECASE,
)
_WHOLE_COROLLA_OVERRIDES_TUBE = re.compile(
    r"\bcorolla\b[^.;]{0,80}\b(?:above|with)\b[^.;]{0,45}\btube\b",
    re.IGNORECASE,
)


def _text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _normalise(value: Any) -> str:
    return (
        unicodedata.normalize("NFKC", _text(value))
        .replace("–", "-")
        .replace("—", "-")
        .replace("−", "-")
        .casefold()
    )


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def classify_measurement_subject(row: pd.Series) -> tuple[str, str]:
    """Return target trait and a reproducible decision reason.

    The raw measurement must occur verbatim after Unicode/dash normalisation.
    The nearest named structure before that exact measurement owns it.  This
    prevents a flower mention early in a sentence from lending its identity to
    a later calyx, petal, inflorescence, or tube measurement.
    """

    quote = _normalise(row.get("exact_supporting_quote", ""))
    raw = _normalise(row.get("raw_value", ""))
    if not quote or not raw:
        return "", "missing_quote_or_raw_measurement"
    position = quote.find(raw)
    if position < 0:
        return "", "raw_measurement_not_in_exact_quote"

    before = quote[:position]
    after = quote[position + len(raw) :]
    if _AMBIGUOUS_COMPARISON.search(before[-35:]):
        return "", "open_ended_comparison_not_classifiable"
    if _BAD_AFTER_MEASUREMENT.search(after[:120]):
        return "", "measurement_belongs_to_display_or_fruiting_context"
    if re.search(
        r"\b(?:double-flowered|long-exserted|antheriferous)\b",
        quote,
        re.IGNORECASE,
    ):
        return "", "cultivar_or_truncated_reproductive_structure_context"
    if re.search(r"\byel-\s*$", before[-20:], re.IGNORECASE):
        return "", "ocr_truncation_immediately_before_measurement"
    if re.match(r"^\s*corolla\s+and\s+up\s+to\b", quote, re.IGNORECASE):
        return "", "truncated_subject_before_corolla_comparison"
    if re.search(
        r"\bflowers?\s+when\s+short\b[^.;]{0,45}\b(?:elongate|branched)\b",
        before[-120:],
        re.IGNORECASE,
    ):
        return "", "inflorescence_growth_not_flower_size"
    if re.match(
        r"^\s*flower,\s*(?:barely|small|ovate|obovate|oblanceolate|"
        r"spatulate|lanceolate|linear)\b",
        quote,
        re.IGNORECASE,
    ):
        return "", "singular_flower_token_is_truncated_part_context"
    if _AMBIGUOUS_WHOLE_SHAPE.search(quote[:position]):
        return "", "shape_grammar_indicates_truncated_floral_part_subject"
    if re.search(r"\bflowers?\s*\)", before[-100:], re.IGNORECASE):
        return "", "truncated_subject_before_measurement"
    if re.search(
        r"\bflowers?\b[^.;]{0,75}\b(?:number|[0-9]+\s+in\s+the)\b",
        before[-120:],
        re.IGNORECASE,
    ):
        return "", "flower_count_or_display_precedes_measurement"

    subjects: list[tuple[int, str, str]] = []
    for pattern, kind in (
        (_WHOLE_SUBJECT, "flower_size_class"),
        (_TUBE_SUBJECT, "tube_depth_class"),
        (_NON_TARGET_SUBJECT, ""),
    ):
        subjects.extend((match.end(), kind, match.group(0)) for match in pattern.finditer(before))
    if not subjects:
        return "", "no_botanical_subject_before_measurement"
    subject_end, trait, subject = max(
        subjects,
        key=lambda item: (item[0], item[1] == "tube_depth_class"),
    )
    bridge = before[subject_end - len(subject) :]
    if len(bridge) > 95:
        return "", "measurement_too_far_from_botanical_subject"
    if not trait:
        return "", f"measurement_of_non_target_part:{subject}"
    tube_context = before[max(0, subject_end - len(subject) - 20) :]
    if trait == "tube_depth_class" and re.search(
        r"\b(?:calyx|calycine)\s+tube\b", tube_context, re.IGNORECASE
    ):
        return "", "calyx_tube_is_not_corolla_tube_depth"
    if trait == "tube_depth_class" and re.match(
        r"^\s*(?:wide|diam(?:eter)?\b)", after[:20], re.IGNORECASE
    ):
        return "", "tube_mouth_width_is_not_tube_depth"
    if trait == "tube_depth_class" and _WHOLE_COROLLA_OVERRIDES_TUBE.search(before[-140:]):
        return "flower_size_class", "corolla_total_length_after_nested_tube_description"
    if trait == "flower_size_class" and re.match(
        r"^\s*long\s+(?:corolla\s+)?tube\b", after[:30], re.IGNORECASE
    ):
        return "tube_depth_class", "postposed_corolla_tube_subject"
    if trait == "flower_size_class" and re.match(
        r"^.{0,30}\b(?:calyx|tepals?|bractlets?)\b", after[:45], re.IGNORECASE
    ):
        return "", "postposed_non_target_floral_part_subject"
    return trait, (
        "exact_measurement_attached_to_whole_flower_or_corolla"
        if trait == "flower_size_class"
        else "exact_measurement_attached_to_corolla_or_perianth_tube"
    )


def _measurement_parts(raw_value: str) -> tuple[str, str]:
    value = _normalise(raw_value)
    unit_match = re.search(r"\b(mm|cm|m)\b", value)
    if not unit_match:
        return "", ""
    measurement = value[: unit_match.start()].strip()
    measurement = re.sub(r"^[^0-9]*", "", measurement)
    return measurement, unit_match.group(1)


def recover_candidates(
    original_candidates: pd.DataFrame,
    context_gated_evidence: pd.DataFrame,
    rules: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build common-schema evidence and a complete decision ledger."""

    allowed_ids = set(
        context_gated_evidence.loc[
            context_gated_evidence["trait_name"].eq("flower_size_class"),
            "source_record_id",
        ]
    )
    source = original_candidates.loc[
        original_candidates["candidate_id"].isin(allowed_ids)
        & original_candidates["trait_name"].eq("flower_size_class")
    ].copy()
    decisions = source.copy()
    classified = decisions.apply(classify_measurement_subject, axis=1)
    decisions["recovered_trait_name"] = [item[0] for item in classified]
    decisions["recovery_reason"] = [item[1] for item in classified]

    rows: list[dict[str, str]] = []
    normalized_values: list[str] = []
    for row in decisions.itertuples(index=False):
        trait = _text(row.recovered_trait_name)
        measurement, unit = _measurement_parts(_text(row.raw_value))
        normalized = derive_class(trait, measurement, unit, rules) if trait else ""
        normalized_values.append(normalized)
        if not trait or not normalized:
            continue
        record_key = "|".join(
            (
                _text(row.candidate_id),
                trait,
                normalized,
                str(rules.get("classification_rule_version", "")),
            )
        )
        rows.append(
            {
                "accepted_species": _text(row.accepted_species),
                "axis": "floral_structural_complexity",
                "trait_name": trait,
                "normalized_value": normalized,
                "quality": "high",
                "source_group": "latest_public_web",
                "source_provider": _text(row.provider),
                "source_url": _text(row.source_url),
                "source_record_id": f"WFOM-{_sha(record_key)[:20]}",
                "source_citation": _text(row.source_title),
                "source_excerpt": _text(row.exact_supporting_quote),
                "evidence_scope": "species_direct",
                "name_match_method": _text(row.name_match_method),
                "source_lineage": _text(row.source_lineage),
                "lineage_method": "provider_treatment_wfo_id",
                "source_run_id": RUN_ID,
                "source_artifact": "wfo_bulk_morphology_recovery_20260808",
                "source_file": "wfo_morphology_recovered_evidence.csv.gz",
                "acceptance_contract": "wfo_morphology_subject_gate_plus_quote_audit_v1",
                "_original_candidate_id": _text(row.candidate_id),
                "_raw_measurement": _text(row.raw_value),
                "_measurement_unit": unit,
                "_classification_rule_version": str(rules.get("classification_rule_version", "")),
                "_recovery_reason": _text(row.recovery_reason),
                "_provider": _text(row.provider),
            }
        )
    decisions["recovered_normalized_value"] = normalized_values
    evidence = pd.DataFrame(rows)
    if evidence.empty:
        return pd.DataFrame(columns=EVIDENCE_COLUMNS), decisions
    evidence = evidence.sort_values(
        ["trait_name", "source_provider", "accepted_species", "source_record_id"]
    ).drop_duplicates("source_record_id", keep="first")
    return evidence.reset_index(drop=True), decisions


def deterministic_audit_queue(
    evidence: pd.DataFrame,
    excluded_audits: list[pd.DataFrame],
    *,
    size: int = AUDIT_SIZE,
) -> pd.DataFrame:
    """Select a fixed provider x trait holdout without prior/development rows."""

    excluded_ids: set[str] = set()
    excluded_original_ids: set[str] = set()
    for audit in excluded_audits:
        if "candidate_id" in audit:
            excluded_ids.update(audit["candidate_id"].map(_text))
        if "_original_candidate_id" in audit:
            excluded_original_ids.update(audit["_original_candidate_id"].map(_text))
    pool = evidence.loc[
        ~evidence["source_record_id"].isin(excluded_ids)
        & ~evidence["_original_candidate_id"].isin(excluded_ids.union(excluded_original_ids))
    ].copy()
    if len(pool) < size:
        raise ValueError(f"only {len(pool)} independent rows are available for a {size}-row audit")
    pool["_audit_hash"] = pool["source_record_id"].map(
        lambda value: _sha(f"wfo-morphology-audit-20260808|{value}")
    )

    # Guarantee representation of every provider x recovered-trait stratum,
    # then fill globally by the same immutable hash.
    selected_indices: list[int] = []
    for _, group in pool.groupby(["source_provider", "trait_name"], sort=True):
        take = min(3, len(group))
        selected_indices.extend(group.sort_values("_audit_hash").head(take).index.tolist())
    selected = pool.loc[sorted(set(selected_indices))].copy()
    remaining = pool.drop(index=selected.index)
    if len(selected) > size:
        selected = selected.sort_values("_audit_hash").head(size)
    elif len(selected) < size:
        selected = pd.concat(
            [
                selected,
                remaining.sort_values("_audit_hash").head(size - len(selected)),
            ],
            ignore_index=False,
        )

    audit = selected.sort_values(["trait_name", "source_provider", "_audit_hash"]).copy()
    audit = audit.rename(
        columns={
            "source_record_id": "candidate_id",
            "source_provider": "provider",
            "source_excerpt": "exact_supporting_quote",
            "_raw_measurement": "raw_measurement",
            "_recovery_reason": "recovery_reason",
        }
    )
    audit["accepted_correct"] = ""
    audit["cultivar_status"] = "wild_or_flora_treatment"
    audit["reviewer"] = ""
    audit["reviewed_at_utc"] = ""
    audit["audit_reason"] = ""
    columns = [
        "candidate_id",
        "_original_candidate_id",
        "accepted_species",
        "trait_name",
        "normalized_value",
        "raw_measurement",
        "provider",
        "source_url",
        "source_lineage",
        "exact_supporting_quote",
        "recovery_reason",
        "accepted_correct",
        "cultivar_status",
        "reviewer",
        "reviewed_at_utc",
        "audit_reason",
    ]
    return audit[columns].reset_index(drop=True)


def finalize_audit(queue: pd.DataFrame, reviews: pd.DataFrame) -> pd.DataFrame:
    required = {"candidate_id", "accepted_correct", "audit_reason"}
    missing = required.difference(reviews.columns)
    if missing:
        raise ValueError(f"review file is missing columns: {sorted(missing)}")
    if reviews["candidate_id"].duplicated().any():
        raise ValueError("review file has duplicate candidate_id values")
    expected = set(queue["candidate_id"])
    supplied = set(reviews["candidate_id"])
    if supplied != expected:
        raise ValueError(
            f"review IDs differ from audit queue: missing={len(expected - supplied)}, "
            f"extra={len(supplied - expected)}"
        )
    final = queue.drop(
        columns=["accepted_correct", "reviewer", "reviewed_at_utc", "audit_reason"]
    ).merge(reviews[list(required)], on="candidate_id", validate="one_to_one")
    final["reviewer"] = AUDIT_REVIEWER
    final["reviewed_at_utc"] = datetime.now(UTC).isoformat().replace("+00:00", "Z")
    values = final["accepted_correct"].astype(str).str.casefold()
    if not values.isin({"true", "false"}).all():
        raise ValueError("accepted_correct must contain only true/false")
    return final


def _write_manifest(output: Path) -> None:
    rows = []
    for path in sorted(output.iterdir()):
        if path.is_file() and path.name != "file_manifest.sha256":
            rows.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}")
    (output / "file_manifest.sha256").write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-candidates", type=Path, required=True)
    parser.add_argument("--context-gated-evidence", type=Path, required=True)
    parser.add_argument(
        "--exclude-audit",
        type=Path,
        action="append",
        default=[],
        help="Prior or development audit whose rows cannot enter the holdout.",
    )
    parser.add_argument(
        "--approved-provider",
        action="append",
        default=[],
        help="Keep only providers that passed the independent audit.",
    )
    parser.add_argument(
        "--completed-audit",
        type=Path,
        help="Completed line-by-line audit; skips creation of a new audit queue.",
    )
    parser.add_argument("--measurement-config", type=Path, required=True)
    parser.add_argument("--review-decisions", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    original = pd.read_csv(args.original_candidates, dtype=str).fillna("")
    gated = pd.read_csv(args.context_gated_evidence, dtype=str).fillna("")
    excluded_audits = [pd.read_csv(path, dtype=str).fillna("") for path in args.exclude_audit]
    rules = load_rules(args.measurement_config)
    evidence, decisions = recover_candidates(original, gated, rules)
    if args.approved_provider:
        evidence = evidence.loc[
            evidence["source_provider"].isin(set(args.approved_provider))
        ].copy()
    queue = (
        pd.DataFrame()
        if args.completed_audit
        else deterministic_audit_queue(evidence, excluded_audits)
    )

    evidence[EVIDENCE_COLUMNS].to_csv(
        args.output / "wfo_morphology_recovered_evidence.csv.gz",
        index=False,
    )
    decisions.to_csv(args.output / "wfo_morphology_recovery_decisions.csv.gz", index=False)
    if not queue.empty:
        queue.to_csv(args.output / "wfo_morphology_audit_queue_200.csv", index=False)

    summary: dict[str, Any] = {
        "contract": "wfo_morphology_recovery_v1",
        "source_run_id": RUN_ID,
        "network_requests": 0,
        "search_queries": 0,
        "cost_usd": 0.0,
        "input_context_gated_flower_size_rows": int(
            gated["trait_name"].eq("flower_size_class").sum()
        ),
        "recovered_evidence_rows": len(evidence),
        "recovered_species_traits": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "by_trait": evidence.groupby("trait_name").size().to_dict(),
        "approved_providers": sorted(set(args.approved_provider)),
        "rejected_or_unclassifiable_rows": int(decisions["recovered_trait_name"].eq("").sum()),
        "audit_queue_rows": len(queue),
        "excluded_audit_files": [str(path) for path in args.exclude_audit],
        "audit_independent_of_excluded_candidate_ids": bool(args.completed_audit)
        or all(
            not bool(
                set(queue["candidate_id"]).intersection(set(audit.get("candidate_id", [])))
                or set(queue["_original_candidate_id"]).intersection(
                    set(audit.get("candidate_id", []))
                    | set(audit.get("_original_candidate_id", []))
                )
            )
            for audit in excluded_audits
        ),
        "classification_rule_version": str(rules.get("classification_rule_version", "")),
    }
    if args.completed_audit:
        final_audit = pd.read_csv(args.completed_audit, dtype=str).fillna("")
        required = {
            "candidate_id",
            "trait_name",
            "accepted_correct",
            "cultivar_status",
            "reviewer",
            "reviewed_at_utc",
            "audit_reason",
        }
        missing = required.difference(final_audit.columns)
        if missing:
            raise ValueError(f"completed audit is missing columns: {sorted(missing)}")
        if final_audit["candidate_id"].duplicated().any():
            raise ValueError("completed audit has duplicate candidate_id values")
        accepted_mask = final_audit["accepted_correct"].str.casefold().eq("true")
        missing_accepted = set(final_audit.loc[accepted_mask, "candidate_id"]).difference(
            set(evidence["source_record_id"])
        )
        if missing_accepted:
            raise ValueError(
                f"completed audit has {len(missing_accepted)} accepted rows missing from evidence"
            )
        final_audit.to_csv(args.output / "wfo_morphology_full_census_audit.csv", index=False)
        accepted = accepted_mask
    elif args.review_decisions:
        reviews = pd.read_csv(args.review_decisions, dtype=str).fillna("")
        final_audit = finalize_audit(queue, reviews)
        final_audit.to_csv(args.output / "wfo_morphology_manual_audit_200.csv", index=False)
        accepted = final_audit["accepted_correct"].str.casefold().eq("true")
    else:
        final_audit = pd.DataFrame()
        accepted = pd.Series(dtype=bool)
    if not final_audit.empty:
        summary["audit"] = {
            "reviewed": len(final_audit),
            "accepted_correct": int(accepted.sum()),
            "precision": float(accepted.mean()),
            "cultivar_contamination_rate": float(
                final_audit["cultivar_status"]
                .str.casefold()
                .str.contains(r"cultivar|hybrid|horticultural", regex=True)
                .mean()
            ),
            "by_trait": (
                final_audit.assign(_accepted=accepted)
                .groupby("trait_name")
                .agg(reviewed=("candidate_id", "size"), accepted_correct=("_accepted", "sum"))
                .reset_index()
                .assign(precision=lambda frame: frame["accepted_correct"] / frame["reviewed"])
                .to_dict("records")
            ),
        }
    (args.output / "wfo_morphology_recovery_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    _write_manifest(args.output)


if __name__ == "__main__":
    main()
