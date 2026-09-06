"""Promote source-locked legacy LLM candidates after a fail-closed audit.

The repository contains model extractions from frozen, species-direct source
chunks that were intentionally left as unreviewed candidates.  This checkpoint
does not call a model or browse the Web.  It revalidates every quote against the
original packet, applies trait-specific category-error guards, performs a
deterministic domain x trait audit, and promotes only scopes whose frozen manual
review reaches the declared precision and cultivar-contamination thresholds.

All promoted rows are Medium.  The model is an extractor, never a source.  The
source URL, citation, source-text hash, prompt hash, model ID and exact quote
remain attached to the audit.  Genus inference is deliberately absent; the
shared all-evidence implementation must rebuild genus x trait_name rules.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import pandas as pd

from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS

CREATED_AT = "2026-08-28T06:30:00Z"
REVIEWER = "Codex source-locked LLM candidate audit"
SOURCE_GROUP = "wave38_source_locked_llm_promotion"
CONTRACT = "wave38_source_locked_llm_medium_promotion_v1"
MIN_REVIEWED_PER_SCOPE = 10
MIN_PRECISION = 0.90
MAX_CULTIVAR_CONTAMINATION = 0.02

TRAIT_AXIS = {
    "flower_primary_color": "flower_colour",
    "floral_form": "floral_structural_complexity",
    "floral_symmetry": "floral_structural_complexity",
    "tube_depth_class": "floral_structural_complexity",
    "flower_size_class": "floral_structural_complexity",
    "inflorescence_display": "floral_structural_complexity",
    "self_incompatibility": "reproductive_assurance",
    "autonomous_selfing_capacity": "reproductive_assurance",
    "mating_system": "reproductive_assurance",
    "cleistogamy": "reproductive_assurance",
}

COLOUR_PATTERNS = {
    "white": r"\b(?:white|cream(?:y)?|ivory)\b",
    "yellow_orange": r"\b(?:yellow(?:ish)?|orange|gold(?:en)?)\b",
    "red_pink": r"\b(?:red(?:dish)?|pink(?:ish)?|rose|crimson|scarlet|maroon|carmine)\b",
    "blue_purple": r"\b(?:blue|purple|violet|lavender|lilac|mauve)\b",
    "green_brown_inconspicuous": (
        r"\b(?:green(?:ish)?|brown(?:ish)?|inconspicuous|insignificant)\b"
    ),
}
FLORAL_PART = re.compile(
    r"\b(?:flowers?|florets?|corolla|petals?|tepals?|perianth|rays?)\b",
    re.IGNORECASE,
)
CULTIVAR = re.compile(
    r"\b(?:cultivar|cv\.|garden variety|horticultural variety|hybrid cultivar|"
    r"double-flowered cultivar|named variety)\b|['\u2018\u2019][A-Z][^'\u2018\u2019]{1,40}['\u2018\u2019]",
    re.IGNORECASE,
)
NON_FLORAL_COLOUR = re.compile(
    r"\b(?:fruits?|berries|seeds?|leaves?|foliage|stems?|bark|wood|roots?)\b",
    re.IGNORECASE,
)

FORM_PATTERNS = {
    "open_radial": r"\b(?:rotate|star-shaped|star shaped|open radial)\b",
    "bell_campanulate": r"\b(?:bell-shaped|bell shaped|campanulate)\b",
    "tubular": r"\b(?:tubular|tube-shaped|tube shaped)\b",
    "salverform": r"\b(?:salverform|salver-shaped|salver shaped|hypocrateriform)\b",
    "funnel_trumpet": r"\b(?:funnel-shaped|funnel shaped|funnelform|trumpet-shaped|trumpet shaped)\b",
    "urn_urceolate": r"\b(?:urn-shaped|urn shaped|urceolate)\b",
    "brush_puff": r"\b(?:brush-like|brush like|powder-puff|powder puff|puffball)\b",
    "composite_head": r"\b(?:flower heads?|capitulum|capitula)\b",
    "papilionaceous": r"\b(?:papilionaceous|pea-shaped|pea shaped)\b",
    "bilabiate": r"\b(?:bilabiate|two-lipped|two lipped)\b",
    "spurred": r"\b(?:spurred|nectar spur|floral spur)\b",
    "reduced_wind": r"\b(?:petal-less|petalless|flowers? lacking petals?|flowers? without petals?)\b",
}

INFLORESCENCE_PATTERNS = {
    "solitary": r"\b(?:solitary flower|flowers? solitary|single flower)\b",
    "few_flowered": r"\b(?:[1-9](?:\s*(?:-|to|or)\s*[1-9])?-?flowered|cluster(?:s)? of [1-9](?:\s*(?:-|to|or)\s*[1-9])? flowers?)\b",
    "raceme_spike_panicle": r"\b(?:racemes?|spikes?|panicles?|catkins?|spadix|spadices)\b",
    "umbel_corymb": r"\b(?:umbels?|corymbs?|cymes?|cymose)\b",
    "composite_display": r"\b(?:flower heads?|capitulum|capitula|cymose heads?)\b",
    "brush_puff_display": r"\b(?:brush-like display|brush like display|powder-puff|powder puff)\b",
}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _canonical_text_sha256(path: Path) -> str:
    """Hash text after the same universal-newline normalization as packet creation."""
    return _sha256_bytes(path.read_text(encoding="utf-8").encode("utf-8"))


def _sha256_text(value: str) -> str:
    return _sha256_bytes(value.encode("utf-8"))


def _write_gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    payload = frame.to_csv(index=False, lineterminator="\n").encode("utf-8")
    with path.open("wb") as raw, gzip.GzipFile(
        filename="", mode="wb", fileobj=raw, mtime=0
    ) as compressed:
        compressed.write(payload)


def _normalized(value: object) -> str:
    return _text(value).casefold().replace("\u2013", "-").replace("\u2014", "-")


def _contains_quote(source: object, quote: object) -> bool:
    normalized_quote = _text(quote)
    return bool(normalized_quote) and normalized_quote in _text(source)


def _candidate_id(record: dict[str, Any]) -> str:
    material = "|".join(
        _text(record.get(column))
        for column in (
            "accepted_species",
            "trait_name",
            "provisional_candidate_value",
            "source_chunk_id",
            "evidence_quote",
        )
    )
    return "W38-" + _sha256_text(material)[:24]


def _colour_states(quote: str) -> set[str]:
    clauses = re.split(r"[.;]", quote)
    states: set[str] = set()
    for clause in clauses:
        if not FLORAL_PART.search(clause):
            continue
        # A compact clause can mention a fruit after the flower.  Retain only
        # text through the first non-floral organ marker.
        marker = NON_FLORAL_COLOUR.search(clause)
        floral = clause[: marker.start()] if marker else clause
        for state, pattern in COLOUR_PATTERNS.items():
            if re.search(pattern, floral, re.IGNORECASE):
                states.add(state)
    return states


def _colour_supported(value: str, quote: str) -> tuple[bool, str]:
    states = _colour_states(quote)
    if not states:
        return False, "no_flower_bound_colour_state"
    if value == "multicolored_variable":
        variable = bool(
            len(states) >= 2
            or re.search(
                r"\b(?:variable|variously|ranges? from|\bto\b|\bor\b|most sub-species|subspecies)\b",
                quote,
                re.IGNORECASE,
            )
        )
        return (variable, "accepted_multistate_colour" if variable else "missing_variation")
    return (
        value in states,
        "accepted_flower_bound_colour" if value in states else "mapped_colour_absent",
    )


def _form_supported(value: str, quote: str) -> tuple[bool, str]:
    if value == "other_described":
        return False, "other_described_not_formal_axis_evidence"
    pattern = FORM_PATTERNS.get(value)
    if not pattern:
        return False, "unsupported_form_value"
    if value == "composite_head":
        supported = bool(re.search(pattern, quote, re.IGNORECASE))
    else:
        supported = bool(
            FLORAL_PART.search(quote)
            and re.search(pattern, quote, re.IGNORECASE)
            and not re.search(r"\bcalyx\b", quote, re.IGNORECASE)
        )
    return supported, "accepted_explicit_floral_form" if supported else "form_not_explicit"


def _symmetry_supported(value: str, quote: str) -> tuple[bool, str]:
    patterns = {
        "actinomorphic": r"\b(?:actinomorphic|radially symmetrical|radial symmetry|Flower Shape:\s*Radial)\b",
        "zygomorphic": r"\b(?:zygomorphic|bilaterally symmetrical|bilateral symmetry)\b",
        "asymmetric": r"\b(?:asymmetric|asymmetrical)\b",
        "mixed_or_variable": r"\b(?:variable|both actinomorphic and zygomorphic)\b",
    }
    supported = bool(
        FLORAL_PART.search(quote)
        and re.search(patterns.get(value, r"(?!)"), quote, re.IGNORECASE)
        and not re.search(r"\b(?:androecium|stamens?|style|ovary)\b", quote, re.IGNORECASE)
    )
    return supported, "accepted_explicit_whole_flower_symmetry" if supported else "symmetry_inferred_or_part_only"


def _measurement_mm(quote: str, *, tube: bool) -> float | None:
    subject = (
        r"(?:corolla|floral|perianth)\s+tube"
        if tube
        else r"(?:flowers?|florets?|corolla|perianth)"
    )
    match = re.search(
        subject
        + r"[^.;]{0,80}?(\d+(?:\.\d+)?)\s*(?:[-\u2013\u2014]\s*(\d+(?:\.\d+)?))?\s*(mm|cm|inches|inch|in\.)\b",
        quote,
        re.IGNORECASE,
    )
    if not match:
        return None
    low = float(match.group(1))
    high = float(match.group(2) or match.group(1))
    midpoint = (low + high) / 2
    unit = match.group(3).casefold()
    if unit == "cm":
        midpoint *= 10
    elif unit in {"inches", "inch", "in."}:
        midpoint *= 25.4
    return midpoint


def _size_supported(value: str, quote: str) -> tuple[bool, str]:
    midpoint = _measurement_mm(quote, tube=False)
    if midpoint is not None:
        expected = (
            "very_small"
            if midpoint <= 5
            else "small"
            if midpoint <= 15
            else "medium"
            if midpoint <= 30
            else "large"
            if midpoint <= 60
            else "very_large"
        )
        return value == expected, "accepted_numeric_whole_flower_size" if value == expected else "size_class_mismatch"
    qualitative = {
        "very_small": r"\b(?:minute|tiny)\s+flowers?\b|\bflowers?\s+(?:are\s+)?(?:minute|tiny)\b",
        "small": r"\bsmall\s+flowers?\b|\bflowers?\s+(?:are\s+)?small\b",
        "large": r"\blarge\s+flowers?\b|\bflowers?\s+(?:are\s+)?large\b",
        "very_large": r"\b(?:very large|huge|enormous)\s+flowers?\b",
    }
    supported = bool(re.search(qualitative.get(value, r"(?!)"), quote, re.IGNORECASE))
    return supported, "accepted_explicit_qualitative_flower_size" if supported else "whole_flower_size_not_explicit"


def _tube_supported(value: str, quote: str) -> tuple[bool, str]:
    midpoint = _measurement_mm(quote, tube=True)
    if midpoint is not None:
        expected = "shallow" if midpoint <= 5 else "intermediate" if midpoint <= 20 else "deep"
        return value == expected, "accepted_numeric_floral_tube_depth" if value == expected else "tube_class_mismatch"
    if value == "absent_or_open":
        supported = bool(
            FLORAL_PART.search(quote)
            and re.search(
                r"\b(?:petal-less|petalless|flowers? lacking petals?|flowers? without petals?|segments? free to base)\b",
                quote,
                re.IGNORECASE,
            )
        )
        return supported, "accepted_explicit_absent_or_open_tube" if supported else "tube_absence_inferred"
    qualitative = {
        "shallow": r"\b(?:shallow|short)\s+(?:corolla|floral|perianth)\s+tube\b",
        "deep": r"\b(?:deep|long)\s+(?:corolla|floral|perianth)\s+tube\b",
    }
    supported = bool(re.search(qualitative.get(value, r"(?!)"), quote, re.IGNORECASE))
    return supported, "accepted_explicit_qualitative_tube_depth" if supported else "tube_depth_not_explicit"


def _inflorescence_supported(value: str, quote: str) -> tuple[bool, str]:
    pattern = INFLORESCENCE_PATTERNS.get(value)
    if value == "multistate_variable":
        found = [
            state
            for state, candidate in INFLORESCENCE_PATTERNS.items()
            if re.search(candidate, quote, re.IGNORECASE)
        ]
        supported = len(found) >= 2 or bool(re.search(r"\b(?:variable|variously)\b", quote, re.IGNORECASE))
    else:
        supported = bool(pattern and re.search(pattern, quote, re.IGNORECASE))
    return supported, "accepted_explicit_inflorescence" if supported else "inflorescence_state_not_explicit"


def _reproduction_supported(trait: str, value: str, quote: str) -> tuple[bool, str]:
    if trait == "self_incompatibility":
        if value == "SI":
            supported = bool(re.search(r"\bself[- ]incompatib", quote, re.IGNORECASE))
        elif value == "SC":
            supported = bool(
                re.search(r"\bself[- ]compatib", quote, re.IGNORECASE)
                or re.search(r"\bself[- ]fertile\s*[:.]?\s*yes\b", quote, re.IGNORECASE)
            )
        else:
            supported = bool(re.search(r"\b(?:mixed|variable).{0,30}self[- ]compat", quote, re.IGNORECASE))
        return supported, "accepted_explicit_si_sc" if supported else "self_fertility_category_error"
    if trait == "autonomous_selfing_capacity":
        patterns = {
            "autonomous": r"\b(?:autonomous self|self[- ]pollinat|autogam)",
            "delayed": r"\bdelayed self",
            "facilitated": r"\bfacilitated self",
            "absent": r"\b(?:no|without|absent).{0,35}(?:autonomous self|autogam)",
            "mixed_or_variable": r"\b(?:variable|mixed).{0,35}(?:autonomous self|autogam)",
        }
        supported = bool(re.search(patterns.get(value, r"(?!)"), quote, re.IGNORECASE))
        return supported, "accepted_explicit_autonomous_selfing" if supported else "autonomous_selfing_not_explicit"
    if trait == "mating_system":
        patterns = {
            "predominantly_outcrossing": r"\b(?:predominantly|mainly|largely) outcross",
            "mixed_mating": r"\bmixed mating",
            "predominantly_selfing": r"\b(?:predominantly|mainly|largely) selfing",
            "obligate_selfing": r"\b(?:obligate self|no evident outcrossing|exclusively self)",
        }
        supported = bool(re.search(patterns.get(value, r"(?!)"), quote, re.IGNORECASE))
        return supported, "accepted_explicit_mating_system" if supported else "sex_system_or_pollination_not_mating_system"
    if trait == "cleistogamy":
        if value == "absent":
            supported = bool(re.search(r"\b(?:no|not|without).{0,20}cleistogam", quote, re.IGNORECASE))
        elif value == "obligate":
            supported = bool(re.search(r"\b(?:obligate|only|exclusively).{0,25}cleistogam|cleistogam.{0,40}no evident outcrossing", quote, re.IGNORECASE))
        else:
            supported = bool(re.search(r"\b(?:facultative|both).{0,25}cleistogam", quote, re.IGNORECASE))
        return supported, "accepted_explicit_cleistogamy" if supported else "cleistogamy_state_not_explicit"
    return False, "unsupported_reproductive_trait"


def candidate_supported(record: dict[str, Any]) -> tuple[bool, str, str]:
    """Return accepted, reason and normalized value for one validated claim."""
    trait = _text(record.get("trait_name"))
    value = _text(record.get("provisional_candidate_value"))
    quote = _text(record.get("evidence_quote"))
    if trait not in TRAIT_AXIS:
        return False, "not_a_strict_three_axis_trait", ""
    if CULTIVAR.search(" ".join([quote, _text(record.get("source_url"))])):
        return False, "cultivar_or_named_horticultural_form", ""
    if trait == "flower_primary_color":
        accepted, reason = _colour_supported(value, quote)
        if accepted and value == "multicolored_variable":
            normalized = "|".join(sorted(_colour_states(quote)) + ["multicolored_variable"])
        else:
            normalized = value
    elif trait == "floral_form":
        accepted, reason = _form_supported(value, quote)
        normalized = value
    elif trait == "floral_symmetry":
        accepted, reason = _symmetry_supported(value, quote)
        normalized = value
    elif trait == "flower_size_class":
        accepted, reason = _size_supported(value, quote)
        normalized = value
    elif trait == "tube_depth_class":
        accepted, reason = _tube_supported(value, quote)
        normalized = value
    elif trait == "inflorescence_display":
        accepted, reason = _inflorescence_supported(value, quote)
        normalized = value
    else:
        accepted, reason = _reproduction_supported(trait, value, quote)
        normalized = value
    return accepted, reason, normalized if accepted else ""


def _provider_lineage(domain: str, citation: str) -> tuple[str, str]:
    normalized_citation = _normalized(citation)
    if domain == "www.gbif.org" and normalized_citation:
        return (
            "cited-source:" + _sha256_text(normalized_citation)[:24],
            "gbif_description_original_citation_dedup",
        )
    provider = domain.removeprefix("www.")
    return f"provider-treatment:{provider}", "provider_wide_redistribution_lineage"


def _load_packets(root: Path) -> tuple[pd.DataFrame, dict[str, dict[str, Any]], list[dict[str, str]]]:
    frames: list[pd.DataFrame] = []
    sources: dict[str, dict[str, Any]] = {}
    manifests: list[dict[str, str]] = []
    for batch_dir in sorted(path for path in root.iterdir() if path.is_dir()):
        candidate_path = batch_dir / "v2_llm_reported_candidates.csv"
        packet_path = batch_dir / "packet_input.json"
        packet_manifest_path = batch_dir / "packet_manifest.json"
        ingest_path = batch_dir / "llm_ingest_manifest.json"
        if not candidate_path.is_file():
            continue
        for required in (packet_path, packet_manifest_path, ingest_path):
            if not required.is_file():
                raise ValueError(f"candidate batch missing frozen contract file: {required}")
        packet_manifest = json.loads(packet_manifest_path.read_text(encoding="utf-8"))
        if _canonical_text_sha256(packet_path) != _text(
            packet_manifest.get("packet_input_sha256")
        ):
            raise ValueError(f"packet hash mismatch: {batch_dir}")
        packet = json.loads(packet_path.read_text(encoding="utf-8"))
        for source in packet.get("sources", []):
            source_id = _text(source.get("source_chunk_id"))
            if source_id and source_id in sources and sources[source_id] != source:
                raise ValueError(f"source chunk collision: {source_id}")
            if source_id:
                sources[source_id] = source
        frame = pd.read_csv(candidate_path, dtype=str).fillna("")
        frame["batch_dir"] = batch_dir.name
        frame["candidate_file_sha256"] = _sha256(candidate_path)
        frames.append(frame)
        manifests.append(
            {
                "batch_dir": batch_dir.name,
                "candidate_csv_sha256": _sha256(candidate_path),
                "packet_input_sha256": _sha256(packet_path),
                "packet_manifest_sha256": _sha256(packet_manifest_path),
                "ingest_manifest_sha256": _sha256(ingest_path),
            }
        )
    if not frames:
        raise ValueError(f"no frozen LLM candidate batches found under {root}")
    return pd.concat(frames, ignore_index=True).fillna(""), sources, manifests


def _audit_sample(candidates: pd.DataFrame) -> pd.DataFrame:
    """Select ten reproducible records per largest domain x trait scope."""
    frame = candidates.copy()
    frame["scope"] = frame["domain"] + "|" + frame["trait_name"]
    frame["audit_order"] = frame["candidate_id"].map(
        lambda value: _sha256_text("wave38-manual-audit|" + value)
    )
    scope_sizes = frame.groupby("scope").size().sort_values(ascending=False)
    eligible_scopes = list(scope_sizes.loc[scope_sizes.ge(MIN_REVIEWED_PER_SCOPE)].index)
    selected_scopes = eligible_scopes[:20]
    sampled = pd.concat(
        [
            frame.loc[frame["scope"].eq(scope)]
            .sort_values("audit_order", kind="stable")
            .head(MIN_REVIEWED_PER_SCOPE)
            for scope in selected_scopes
        ],
        ignore_index=True,
    )
    return sampled.sort_values(["scope", "audit_order"], kind="stable")


def _bool(value: object) -> bool:
    return _text(value).casefold() in {"1", "true", "yes", "y"}


def validated_candidate_tables(
    *,
    candidate_root: Path,
    master_csv: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict[str, str]]]:
    """Return all candidates, the full audit, deduplicated supported rows and inputs."""
    candidates, sources, input_manifests = _load_packets(candidate_root)
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    accepted_species = set(master["accepted_species"])
    audit_rows: list[dict[str, Any]] = []
    validated_rows: list[dict[str, Any]] = []
    for record in candidates.to_dict("records"):
        record["candidate_id"] = _candidate_id(record)
        source = sources.get(_text(record.get("source_chunk_id")))
        validation_reason = ""
        if source is None:
            validation_reason = "missing_frozen_source_chunk"
        elif _text(source.get("accepted_species")) != _text(record.get("accepted_species")):
            validation_reason = "source_species_mismatch"
        elif _text(source.get("evidence_scope")) != "species_direct":
            validation_reason = "source_not_species_direct"
        elif _text(record.get("accepted_species")) not in accepted_species:
            validation_reason = "species_outside_fixed_master"
        elif not _contains_quote(source.get("chunk_text"), record.get("evidence_quote")):
            validation_reason = "quote_not_in_frozen_source"
        elif _text(source.get("source_text_sha256")) != _text(record.get("source_text_sha256")):
            validation_reason = "source_text_hash_mismatch"
        if validation_reason:
            supported, rule_reason, normalized_value = False, validation_reason, ""
        else:
            supported, rule_reason, normalized_value = candidate_supported(record)
        domain = urlparse(_text(record.get("source_url"))).netloc.casefold()
        citation = _text(source.get("source_citation")) if source else ""
        lineage, lineage_method = _provider_lineage(domain, citation)
        row = {
            **record,
            "candidate_id": record["candidate_id"],
            "domain": domain,
            "source_citation": citation,
            "source_chunk_text": _text(source.get("chunk_text")) if source else "",
            "rule_supported": supported,
            "rule_reason": rule_reason,
            "normalized_value_after_audit": normalized_value,
            "source_lineage_after_audit": lineage,
            "lineage_method_after_audit": lineage_method,
            "normalized_excerpt_sha256": _sha256_text(_normalized(record.get("evidence_quote"))),
            "content_fingerprint": _sha256_text(
                "|".join(
                    [
                        _text(record.get("accepted_species")),
                        _text(record.get("trait_name")),
                        _normalized(record.get("evidence_quote")),
                    ]
                )
            ),
        }
        audit_rows.append(row)
        if supported:
            validated_rows.append(row)
    audit = pd.DataFrame(audit_rows).fillna("")
    validated = pd.DataFrame(validated_rows).fillna("")
    validated = validated.sort_values(
        ["candidate_id", "domain", "source_url"], kind="stable"
    ).drop_duplicates("candidate_id")
    duplicate_excerpt = validated.duplicated(
        ["accepted_species", "trait_name", "content_fingerprint"], keep="first"
    )
    audit.loc[
        audit["candidate_id"].isin(validated.loc[duplicate_excerpt, "candidate_id"]),
        ["rule_supported", "rule_reason"],
    ] = [False, "mirror_or_republication_excerpt_duplicate"]
    validated = validated.loc[~duplicate_excerpt].copy()
    duplicate_lineage = validated.duplicated(
        ["accepted_species", "trait_name", "source_lineage_after_audit"], keep="first"
    )
    audit.loc[
        audit["candidate_id"].isin(validated.loc[duplicate_lineage, "candidate_id"]),
        ["rule_supported", "rule_reason"],
    ] = [False, "same_original_source_lineage_duplicate"]
    validated = validated.loc[~duplicate_lineage].copy()

    return candidates, audit, validated, input_manifests


def prepare_review_sample(
    *, candidate_root: Path, master_csv: Path, output_csv: Path
) -> dict[str, Any]:
    candidates, audit, validated, input_manifests = validated_candidate_tables(
        candidate_root=candidate_root,
        master_csv=master_csv,
    )
    sample = _audit_sample(validated)
    columns = [
        "candidate_id",
        "domain",
        "trait_name",
        "accepted_species",
        "provisional_candidate_value",
        "normalized_value_after_audit",
        "source_url",
        "source_citation",
        "evidence_quote",
        "source_chunk_id",
        "source_text_sha256",
        "model_provider",
        "model_id",
        "prompt_version",
        "prompt_sha256",
        "packet_input_sha256",
        "rule_reason",
        "content_fingerprint",
        "audit_order",
    ]
    sample.reindex(columns=columns).to_csv(output_csv, index=False, lineterminator="\n")
    report = {
        "input_candidate_rows": len(candidates),
        "rule_supported_rows": int(audit["rule_supported"].astype(bool).sum()),
        "deduplicated_supported_rows": len(validated),
        "manual_sample_rows": len(sample),
        "manual_sample_scopes": int(sample[["domain", "trait_name"]].drop_duplicates().shape[0]),
        "candidate_batches": len(input_manifests),
        "output_csv": output_csv.as_posix(),
        "output_sha256": _sha256(output_csv),
    }
    return report


def build_checkpoint(
    *,
    candidate_root: Path,
    decisions_csv: Path,
    master_csv: Path,
    output_dir: Path,
    source_run_id: str,
    source_artifact: str,
) -> dict[str, Any]:
    candidates, audit, validated, input_manifests = validated_candidate_tables(
        candidate_root=candidate_root,
        master_csv=master_csv,
    )

    sample = _audit_sample(validated)
    decisions = pd.read_csv(decisions_csv, dtype=str).fillna("")
    required_decisions = {
        "candidate_id",
        "decision",
        "species_identity_correct",
        "value_correct",
        "provenance_complete",
        "cultivar_contamination",
        "decision_reason",
        "reviewer",
        "reviewed_at_utc",
    }
    if missing := required_decisions.difference(decisions.columns):
        raise ValueError(f"manual decisions missing columns: {sorted(missing)}")
    if decisions["candidate_id"].duplicated().any():
        raise ValueError("manual decisions contain duplicate candidate IDs")
    if set(decisions["candidate_id"]) != set(sample["candidate_id"]):
        missing = sorted(set(sample["candidate_id"]) - set(decisions["candidate_id"]))
        extra = sorted(set(decisions["candidate_id"]) - set(sample["candidate_id"]))
        raise ValueError(f"manual decisions do not match deterministic sample; missing={missing}, extra={extra}")
    manual = sample.merge(decisions, on="candidate_id", validate="one_to_one")
    manual["accepted_correct"] = (
        manual["decision"].str.casefold().eq("accept")
        & manual["species_identity_correct"].map(_bool)
        & manual["value_correct"].map(_bool)
        & manual["provenance_complete"].map(_bool)
        & ~manual["cultivar_contamination"].map(_bool)
    )
    manual["scope"] = manual["domain"] + "|" + manual["trait_name"]
    scope_metrics = (
        manual.groupby(["domain", "trait_name"], sort=True)
        .agg(
            reviewed=("candidate_id", "size"),
            accepted_correct=("accepted_correct", "sum"),
            cultivar_contamination=("cultivar_contamination", lambda values: sum(map(_bool, values))),
        )
        .reset_index()
    )
    scope_metrics["precision"] = scope_metrics["accepted_correct"] / scope_metrics["reviewed"]
    scope_metrics["cultivar_contamination_rate"] = (
        scope_metrics["cultivar_contamination"] / scope_metrics["reviewed"]
    )
    scope_metrics["approved"] = (
        scope_metrics["reviewed"].ge(MIN_REVIEWED_PER_SCOPE)
        & scope_metrics["precision"].ge(MIN_PRECISION)
        & scope_metrics["cultivar_contamination_rate"].le(MAX_CULTIVAR_CONTAMINATION)
    )
    approved_scopes = set(
        zip(
            scope_metrics.loc[scope_metrics["approved"], "domain"],
            scope_metrics.loc[scope_metrics["approved"], "trait_name"],
            strict=True,
        )
    )
    promoted = validated.loc[
        [
            (domain, trait) in approved_scopes
            for domain, trait in zip(validated["domain"], validated["trait_name"], strict=True)
        ]
    ].copy()
    rejected_manual_ids = set(
        manual.loc[~manual["accepted_correct"], "candidate_id"]
    )
    promoted = promoted.loc[~promoted["candidate_id"].isin(rejected_manual_ids)].copy()

    evidence_rows: list[dict[str, str]] = []
    for record in promoted.to_dict("records"):
        evidence_rows.append(
            {
                "accepted_species": _text(record["accepted_species"]),
                "axis": TRAIT_AXIS[_text(record["trait_name"])],
                "trait_name": _text(record["trait_name"]),
                "normalized_value": _text(record["normalized_value_after_audit"]),
                "quality": "medium",
                "source_group": SOURCE_GROUP,
                "source_provider": _text(record["domain"]),
                "source_url": _text(record["source_url"]),
                "source_record_id": _text(record["candidate_id"]),
                "source_citation": _text(record["source_citation"]),
                "source_excerpt": _text(record["evidence_quote"]),
                "evidence_scope": "species_direct",
                "name_match_method": "accepted_name_exact_source_packet",
                "source_lineage": _text(record["source_lineage_after_audit"]),
                "lineage_method": _text(record["lineage_method_after_audit"]),
                "source_run_id": _text(source_run_id),
                "source_artifact": _text(source_artifact),
                "source_file": (
                    "data/v2/staging/traits/llm_evidence_extracted/"
                    + _text(record["batch_dir"])
                    + "/v2_llm_reported_candidates.csv"
                ),
                "acceptance_contract": CONTRACT,
            }
        )
    evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS).fillna("")
    if evidence["source_record_id"].duplicated().any():
        raise ValueError("promoted source record IDs are not unique")

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "wave38_reviewed_direct_evidence.csv.gz"
    full_audit_path = output_dir / "wave38_full_candidate_audit.csv.gz"
    manual_path = output_dir / "wave38_manual_audit_200.csv.gz"
    metrics_path = output_dir / "wave38_scope_precision.csv.gz"
    manifest_path = output_dir / "wave38_input_manifest.json"
    _write_gzip_csv(evidence, evidence_path)
    _write_gzip_csv(audit, full_audit_path)
    _write_gzip_csv(manual, manual_path)
    _write_gzip_csv(scope_metrics, metrics_path)
    manifest_path.write_text(
        json.dumps(
            {
                "contract": CONTRACT,
                "created_at_utc": CREATED_AT,
                "candidate_root": candidate_root.as_posix(),
                "decisions_csv": decisions_csv.as_posix(),
                "decisions_sha256": _sha256(decisions_csv),
                "master_sha256": _sha256(master_csv),
                "source_batches": input_manifests,
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    summary = {
        "contract": CONTRACT,
        "created_at_utc": CREATED_AT,
        "source_run_id": _text(source_run_id),
        "source_artifact": _text(source_artifact),
        "denominator": {"species": 106_295, "species_axis": 318_885},
        "input_candidate_rows": len(candidates),
        "strict_trait_candidate_rows": int(candidates["trait_name"].isin(TRAIT_AXIS).sum()),
        "rule_supported_before_lineage_dedup": int(
            audit["rule_supported"].astype(bool).sum()
        ),
        "rule_supported_after_lineage_dedup": len(validated),
        "manual_reviewed_rows": len(manual),
        "approved_domain_trait_scopes": len(approved_scopes),
        "promoted_direct_rows": len(evidence),
        "promoted_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "promoted_by_axis": evidence["axis"].value_counts().sort_index().to_dict(),
        "promoted_by_trait": evidence["trait_name"].value_counts().sort_index().to_dict(),
        "promoted_by_domain": evidence["source_provider"].value_counts().sort_index().to_dict(),
        "guardrails": {
            "source_free_model_answers": False,
            "exact_quote_reverified": True,
            "model_output_is_source": False,
            "all_promoted_quality": "medium",
            "trait_specific_category_guards": True,
            "self_fertile_no_to_si": False,
            "dioecy_to_mating_system": False,
            "pollen_vector_or_reward_in_strict_axis": False,
            "genus_rule_join": "downstream_genus_x_trait_name",
            "family_inference": False,
            "global_fallback": False,
        },
        "artifact_sha256": {
            path.name: _sha256(path)
            for path in (
                evidence_path,
                full_audit_path,
                manual_path,
                metrics_path,
                manifest_path,
            )
        },
    }
    summary_path = output_dir / "wave38_checkpoint_summary.json"
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", required=True, type=Path)
    parser.add_argument("--decisions-csv", type=Path)
    parser.add_argument("--master-csv", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--source-run-id", default="")
    parser.add_argument("--source-artifact", default="")
    parser.add_argument("--prepare-review-sample", type=Path)
    args = parser.parse_args()
    if args.prepare_review_sample:
        report = prepare_review_sample(
            candidate_root=args.candidate_root,
            master_csv=args.master_csv,
            output_csv=args.prepare_review_sample,
        )
    else:
        if not args.decisions_csv or not args.output_dir:
            parser.error("--decisions-csv and --output-dir are required outside prepare mode")
        if not args.source_run_id or not args.source_artifact:
            parser.error("--source-run-id and --source-artifact are required outside prepare mode")
        report = build_checkpoint(
            candidate_root=args.candidate_root,
            decisions_csv=args.decisions_csv,
            master_csv=args.master_csv,
            output_dir=args.output_dir,
            source_run_id=args.source_run_id,
            source_artifact=args.source_artifact,
        )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
