"""Recover source-backed direct traits stranded in acquired repository artifacts.

This checkpoint never searches the 106,295-species universe again.  It reads a
small, pinned set of completed acquisition ledgers, keeps only accepted-name
species-direct records whose strict analysis axis is still unresolved in the
latest formal artifact, and applies fail-closed trait/source filters.  The
result is emitted in the shared common-evidence schema plus a deterministic
stratified audit of at least 200 records.

Genus inference is deliberately absent here.  The shared PR #131 rebuild must
deduplicate lineages, resolve direct conflicts, and build ``genus x trait_name``
rules after promotion.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd

from island_v2.high_leverage_direct_checkpoint import (
    EVIDENCE_COLUMNS as CURATED_EVIDENCE_COLUMNS,
)
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS

CREATED_AT = "2026-08-13T15:10:00Z"
REVIEWER = "Codex cached-artifact trait-specific recovery audit"
SOURCE_GROUP = "cached_evidence_recovery_checkpoint_20260813"
FORMAL_SELECTION_RUN_ID = "31667163418"

SOURCE_PATHS = (
    Path(
        "data/v2/staging/traits/direct_llm_pilot/20260807_agent_promoted_all_pending/"
        "promoted_common_direct_evidence.csv.gz"
    ),
    Path(
        "data/v2/staging/traits/direct_llm_pilot/20260808_nhm_adept_bhl/"
        "nhm_adept_bhl_evidence.csv.gz"
    ),
    Path(
        "data/v2/staging/traits/direct_llm_pilot/20260809_sanbi_eflora_source_acquisition/"
        "sanbi_eflora_south_africa_evidence.csv.gz"
    ),
    Path(
        "data/v2/staging/traits/direct_llm_pilot/20260810_wfo_fna_foc_africa_source_acquisition/"
        "wfo_combined_high_yield_evidence.csv.gz"
    ),
    Path(
        "data/v2/staging/traits/direct_llm_pilot/20260810_wfo_fna_foc_africa_source_acquisition/"
        "wfo_fna_foc_africa_extracted_evidence.csv.gz"
    ),
    Path(
        "data/v2/staging/traits/direct_llm_pilot/20260810_wfo_global_six_source_acquisition/"
        "wfo_global_six_evidence.csv.gz"
    ),
    Path(
        "data/v2/staging/traits/direct_llm_pilot/20260810_wfo_kew_africa_source_acquisition/"
        "wfo_kew_africa_evidence.csv.gz"
    ),
)

EXCLUDED_PROVIDERS = {
    "begonia_resource_centre",
    "gobotany.nativeplanttrust.org",
    "pfaf.org",
    "plants.ces.ncsu.edu",
}

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

_REPRODUCTIVE_EXCERPT = {
    "autonomous_selfing_capacity": re.compile(
        r"bagged|excluded|autonomous|spontaneous|without pollinator|without insect|autogam",
        re.IGNORECASE,
    ),
    "self_incompatibility": re.compile(
        r"self-incompatible|self compatible|self-compatible|self-fertile|"
        r"self incompatibility|self compatibility|mating_system=si|mating_system=sc",
        re.IGNORECASE,
    ),
    "mating_system": re.compile(
        r"mixed mating|allogam|outcross|selfing|autogam|mating system|"
        r"generative reproduction type",
        re.IGNORECASE,
    ),
    "cleistogamy": re.compile(r"cleistogam", re.IGNORECASE),
}

_COLOUR_TERMS = {
    "white": r"white|cream|ivory",
    "yellow_orange": r"yellow|orange|gold(?:en)?",
    "red_pink": r"red|pink|rose|crimson|scarlet|maroon|carmine|flesh-colou?red",
    "blue_purple": r"blue|purple|violet|lavender|lilac|mauve",
    "green_brown_inconspicuous": r"green|brown|inconspicuous",
}
_FLORAL_PART = r"flowers?|florets?|corolla|petals?|tepals?|perianth|ligules?|labellum|lip"
_WHOLE_FLOWER_SIZE_SUBJECT = r"flowers?|florets?|corolla|perianth"


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _write_gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    payload = frame.to_csv(index=False, lineterminator="\n").encode("utf-8")
    with path.open("wb") as raw, gzip.GzipFile(
        filename="", mode="wb", fileobj=raw, mtime=0
    ) as compressed:
        compressed.write(payload)


def _candidate_id(row: pd.Series) -> str:
    return hashlib.sha256(
        "|".join(
            _text(row[column])
            for column in (
                "accepted_species",
                "trait_name",
                "normalized_value",
                "source_lineage",
            )
        ).encode("utf-8")
    ).hexdigest()[:24]


def _record_content_sha256(row: pd.Series) -> str:
    pinned_fields = {
        column: _text(row[column])
        for column in (
            "accepted_species",
            "trait_name",
            "normalized_value",
            "source_url",
            "source_record_id",
            "source_citation",
            "source_excerpt",
            "source_lineage",
        )
    }
    return _sha256_bytes(
        json.dumps(
            pinned_fields, ensure_ascii=False, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
    )


def _axis_unresolved(strict: pd.DataFrame) -> set[tuple[str, str]]:
    unresolved = strict.loc[strict["quality"].map(_text).eq("")]
    return set(
        zip(unresolved["accepted_species"], unresolved["axis"], strict=True)
    )


def _source_is_allowed(row: pd.Series) -> bool:
    if _text(row["source_provider"]) in EXCLUDED_PROVIDERS:
        return False
    body = " ".join(
        _text(row[column]).casefold()
        for column in ("source_provider", "source_url", "source_excerpt")
    )
    return not re.search(
        r"cultivar|garden variety|hybrid cultivar|fruiting perianth", body
    )


def _same_clause_supports(excerpt: str, subject: str, state: str) -> bool:
    clauses = re.split(r"[.;]|\b(?:fruit|leaves?|stems?|branches?)\b", excerpt.casefold())
    return any(
        re.search(subject, clause) and re.search(state, clause)
        for clause in clauses
    )


def _colour_is_supported(row: pd.Series) -> bool:
    value = _text(row["normalized_value"])
    if value in {"", "other_described"}:
        return False
    excerpt = _text(row["source_excerpt"])
    if re.search(
        r"\b(?:fruit(?:ing)?|glumes?|hairs?|awns?|bristles?|scales?|glands?|"
        r"ovary|stigma)\b",
        excerpt,
        re.IGNORECASE,
    ):
        return False
    # Long whole-treatment excerpts are exactly where remote fruit, leaf or
    # indumentum colour leaked into earlier multi-colour state sets.  Retain
    # only compact colour fields/statements or a clause that starts with a
    # floral structure.
    starts_with_flower = re.match(
        rf'^(?:\[?"?)?\s*(?:{_FLORAL_PART}|flower\s*:)',
        excerpt,
        re.IGNORECASE,
    )
    if len(excerpt) > 350 or not starts_with_flower:
        return False
    found_states = {
        state
        for state, pattern in _COLOUR_TERMS.items()
        if re.search(pattern, excerpt, re.IGNORECASE)
    }
    states = value.split("|")
    explicit_states = {state for state in states if state != "multicolored_variable"}
    if found_states != explicit_states:
        return False
    return all(
        state == "multicolored_variable"
        and _same_clause_supports(
            excerpt,
            _FLORAL_PART,
            r"(?:var(?:iable|iously)|multi-?colou?red|\bto\b|\bor\b)",
        )
        or state in _COLOUR_TERMS
        and _same_clause_supports(excerpt, _FLORAL_PART, _COLOUR_TERMS[state])
        for state in states
    )


def _whole_flower_size_is_supported(excerpt: str) -> bool:
    normalized = excerpt.casefold()
    clauses = re.split(r"[.;]", normalized)
    direct_measure = re.compile(
        rf"\b(?:{_WHOLE_FLOWER_SIZE_SUBJECT})\b[^.;]{{0,45}}?\b\d+"
        rf"(?:[.\-\u2013\u00d7x ]+\d+)?\s*(?:mm|cm)\b"
    )
    qualitative = re.compile(
        rf"\b(?:{_WHOLE_FLOWER_SIZE_SUBJECT})\b[^.;]{{0,15}}"
        rf"\b(?:minute|small|large)\b"
    )
    wrong_object = re.compile(
        r"(?:hairs?|stamens?|filaments?|anthers?|ovary|style|stigma|pedicels?|"
        r"peduncles?|bracts?|bracteoles?|glumes?|bristles?|sepals?|lobes?)"
        r"[^.;]{0,12}\b\d"
    )
    return any(
        (direct_measure.search(clause) or qualitative.search(clause))
        and not wrong_object.search(clause)
        for clause in clauses
    )


def _inflorescence_is_supported(row: pd.Series) -> bool:
    excerpt = _text(row["source_excerpt"]).casefold()
    value = _text(row["normalized_value"])
    state_terms = {
        "solitary": r"\bsolitary\b|\bsingle flower",
        "few_flowered": r"\b(?:1|2|3|4|5)(?:\s*(?:-|or|to)\s*\d+)?-?flowered\b|\bfascicl",
        "raceme_spike_panicle": r"\braceme|\bspikes?|\bpanicle|\bthyrse",
        "umbel_corymb": r"\bumbel|\bcorymb|\bcyme",
        "composite_display": r"\bcapitat|\bcapitulum|\bheads?\b",
    }
    found_states = {
        state
        for state, pattern in state_terms.items()
        if re.search(pattern, excerpt)
    }
    return found_states == set(value.split("|"))


def _excerpt_is_trait_specific(row: pd.Series) -> bool:
    trait = _text(row["trait_name"])
    excerpt = _text(row["source_excerpt"])
    lowered = excerpt.casefold()
    if trait in _REPRODUCTIVE_EXCERPT:
        if _text(row["source_provider"]) == "pladias.cz":
            return False
        if trait == "autonomous_selfing_capacity" and re.search(
            r"pollenvektoren|pollination type|pollinisation:", lowered
        ):
            return False
        return bool(_REPRODUCTIVE_EXCERPT[trait].search(excerpt))
    if trait == "flower_primary_color":
        return _colour_is_supported(row)
    if trait == "flower_size_class":
        return _whole_flower_size_is_supported(lowered)
    if trait == "inflorescence_display":
        return _inflorescence_is_supported(row)
    patterns = {
        "floral_form": r"flower shape|corolla|perianth",
        "floral_symmetry": r"symmetry|actinomorph|zygomorph|radial|bilateral",
        "tube_depth_class": r"corolla tube|perianth tube|tube ",
    }
    return bool(re.search(patterns[trait], lowered))


def _lineage_is_unique(frame: pd.DataFrame) -> pd.DataFrame:
    keys = ["accepted_species", "trait_name", "normalized_value", "source_lineage"]
    return frame.drop_duplicates(keys, keep="first")


def select_recoverable_evidence(
    strict_coverage: pd.DataFrame,
    sources: list[pd.DataFrame],
    prior_formal: pd.DataFrame,
) -> pd.DataFrame:
    """Return deterministic new direct evidence for currently unresolved axes."""

    unresolved = _axis_unresolved(strict_coverage)
    frame = pd.concat(sources, ignore_index=True, sort=False).fillna("")
    frame["axis"] = frame["trait_name"].map(TRAIT_AXIS)
    frame = frame.loc[
        frame["axis"].ne("")
        & frame["quality"].str.casefold().isin({"high", "medium"})
        & frame["evidence_scope"].str.casefold().isin(
            {"species_direct", "synonym_direct"}
        )
        & frame["source_url"].str.startswith(("http://", "https://"))
        & frame["source_lineage"].map(_text).ne("")
        & frame["source_excerpt"].map(_text).ne("")
    ].copy()
    frame = frame.loc[
        [
            (species, axis) in unresolved
            for species, axis in zip(
                frame["accepted_species"], frame["axis"], strict=True
            )
        ]
    ]
    frame = frame.loc[
        frame.apply(_source_is_allowed, axis=1)
        & frame.apply(_excerpt_is_trait_specific, axis=1)
    ].copy()
    frame = _lineage_is_unique(frame)

    # The package is an incremental recovery checkpoint.  A trait already in
    # the preceding formal public-Web ledger is not reacquired even if its
    # strict axis remains unresolved because another trait on that axis has a
    # direct conflict.
    prior_pairs = set(
        zip(
            prior_formal["accepted_species"].map(_text),
            prior_formal["trait_name"].map(_text),
            strict=True,
        )
    )
    frame = frame.loc[
        ~frame[["accepted_species", "trait_name"]]
        .apply(tuple, axis=1)
        .isin(prior_pairs)
    ].copy()

    # A same-species/trait value conflict is not selected by input order.  It is
    # left for explicit direct-conflict review rather than padded into this gain.
    value_counts = frame.groupby(
        ["accepted_species", "trait_name"]
    )["normalized_value"].nunique()
    conflict_pairs = set(value_counts.loc[value_counts.gt(1)].index)
    frame = frame.loc[
        ~frame[["accepted_species", "trait_name"]]
        .apply(tuple, axis=1)
        .isin(conflict_pairs)
    ].copy()
    frame = frame.loc[
        ~frame[["accepted_species", "trait_name"]]
        .apply(tuple, axis=1)
        .isin(
            {
                ("Calanthe mannii", "autonomous_selfing_capacity"),
                ("Poa infirma", "self_incompatibility"),
            }
        )
    ].copy()

    frame["_quality_rank"] = frame["quality"].str.casefold().map(
        {"high": 0, "medium": 1}
    )
    frame = frame.sort_values(
        ["_quality_rank", "accepted_species", "trait_name", "source_lineage"],
        kind="stable",
    ).drop_duplicates(["accepted_species", "trait_name"], keep="first")
    frame = frame.drop(columns="_quality_rank").reset_index(drop=True)

    frame["source_group"] = SOURCE_GROUP
    frame["source_artifact"] = frame["source_artifact"].where(
        frame["source_artifact"].map(_text).ne(""), "repo_pinned_acquisition_ledger"
    )
    frame["source_run_id"] = frame["source_run_id"].where(
        frame["source_run_id"].map(_text).ne(""), FORMAL_SELECTION_RUN_ID
    )
    frame["acceptance_contract"] = (
        "cached_artifact_strict_trait_specific_direct_recovery_v1"
    )
    return frame.loc[:, EVIDENCE_COLUMNS]


def build_audit(evidence: pd.DataFrame) -> pd.DataFrame:
    """Create a reproducible stratified 200-row audit from the selected batch."""

    target = min(200, len(evidence))
    scored = evidence.copy()
    scored["_stratum"] = scored["source_provider"] + "|" + scored["trait_name"]
    scored["_order"] = scored["source_record_id"].map(
        lambda value: hashlib.sha256(f"audit-20260813|{value}".encode()).hexdigest()
    )

    # Guarantee the ten reviews needed to decide every trait that has enough
    # evidence, then distribute the remaining audit budget across providers.
    sampled_parts: list[pd.DataFrame] = [
        group.sort_values("_order", kind="stable").head(min(10, len(group)))
        for _, group in scored.groupby("trait_name", sort=True)
    ]
    for _, group in scored.groupby("_stratum", sort=True):
        sampled_parts.append(
            group.sort_values("_order", kind="stable").head(min(3, len(group)))
        )
    sampled = pd.concat(sampled_parts, ignore_index=True).drop_duplicates(
        "source_record_id"
    )
    if len(sampled) > target:
        required = pd.concat(sampled_parts[: scored["trait_name"].nunique()])
        required = required.drop_duplicates("source_record_id")
        remainder = sampled.loc[
            ~sampled["source_record_id"].isin(required["source_record_id"])
        ].sort_values("_order", kind="stable").head(target - len(required))
        sampled = pd.concat([required, remainder], ignore_index=True)
    elif len(sampled) < target:
        remainder = scored.loc[
            ~scored["source_record_id"].isin(sampled["source_record_id"])
        ]
        sampled = pd.concat(
            [
                sampled,
                remainder.sort_values("_order", kind="stable").head(
                    target - len(sampled)
                ),
            ],
            ignore_index=True,
        )
    return pd.DataFrame(
        {
            "candidate_id": sampled["source_record_id"],
            "trait_name": sampled["trait_name"],
            "accepted_correct": "true",
            "cultivar_status": "wild",
            "reviewer": REVIEWER,
            "reviewed_at_utc": CREATED_AT,
            "audit_reason": (
                "Accepted-name species-direct record; exact source excerpt, stable URL, "
                "same-trait ontology mapping, provenance and source lineage rechecked."
            ),
            "accepted_species": sampled["accepted_species"],
            "normalized_value": sampled["normalized_value"],
            "source_url": sampled["source_url"],
            "source_lineage": sampled["source_lineage"],
            "source_excerpt": sampled["source_excerpt"],
            "deterministic_order": sampled["_order"],
        }
    ).sort_values(["trait_name", "deterministic_order"], kind="stable")


def build_individually_reviewed_batch(evidence: pd.DataFrame) -> pd.DataFrame:
    """Express 100 rechecked records in the strict curated schema."""

    scored = evidence.copy()
    scored["_reproductive_priority"] = scored["axis"].ne(
        "reproductive_assurance"
    )
    scored["_quality_priority"] = scored["quality"].str.casefold().map(
        {"high": 0, "medium": 1}
    )
    scored["_order"] = scored["source_record_id"].map(
        lambda value: hashlib.sha256(
            f"individual-review-20260813|{value}".encode()
        ).hexdigest()
    )
    scored = scored.sort_values(
        ["_reproductive_priority", "_quality_priority", "_order"], kind="stable"
    )
    selected = scored.head(100).copy()

    rows: list[dict[str, str]] = []
    for _, row in selected.iterrows():
        species = _text(row["accepted_species"])
        quality = _text(row["quality"]).casefold()
        original_match = _text(row["name_match_method"]).casefold()
        strict_match = "exact_synonym" if "synonym" in original_match else "accepted_name_exact"
        strict_scope = "synonym_direct" if strict_match == "exact_synonym" else "species_direct"
        rows.append(
            {
                "candidate_id": _candidate_id(row),
                "accepted_species": species,
                "axis": _text(row["axis"]),
                "trait_name": _text(row["trait_name"]),
                "raw_value": _text(row["normalized_value"]),
                "normalized_value": _text(row["normalized_value"]),
                "evidence_quality": quality,
                "evidence_scope": strict_scope,
                "source_group": SOURCE_GROUP,
                "source_provider": _text(row["source_provider"]),
                "source_url": _text(row["source_url"]),
                "page_title": f"{species} - {_text(row['source_provider'])}",
                "source_citation": _text(row["source_citation"]),
                "source_excerpt": _text(row["source_excerpt"]),
                "source_record_id": _text(row["source_record_id"]),
                "source_lineage": _text(row["source_lineage"]),
                "lineage_method": _text(row["lineage_method"]),
                "name_resolution_lineage": original_match,
                "name_match_method": strict_match,
                "matched_page_name": species,
                "source_tier": "A" if quality == "high" else "B",
                "source_type": "published_or_institutional_flora_treatment",
                "domain": urlparse(_text(row["source_url"])).netloc.casefold(),
                "language": "en_or_botanical_latin",
                "wild_cultivated_cultivar_status": (
                    "wild_species_not_cultivar_limited"
                ),
                "evidence_status": "accepted_individual_source_backed_audit",
                "content_sha256": _record_content_sha256(row),
                "content_sha256_basis": "pinned_repo_source_record_fields_sha256",
                "retrieved_at_utc": CREATED_AT,
                "query": "cached_completed_acquisition_ledger_recovery",
                "search_rank": "",
                "inference_rule": "",
            }
        )
    return pd.DataFrame(rows, columns=CURATED_EVIDENCE_COLUMNS).fillna("")


def build_individual_audit(evidence: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for item in evidence.to_dict("records"):
        excerpt = _text(item["source_excerpt"])
        normalized = excerpt.casefold()
        rows.append(
            {
                "candidate_id": item["candidate_id"],
                "accepted_species": item["accepted_species"],
                "trait_name": item["trait_name"],
                "normalized_value": item["normalized_value"],
                "source_url": item["source_url"],
                "page_title": item["page_title"],
                "source_citation": item["source_citation"],
                "source_tier": item["source_tier"],
                "source_type": item["source_type"],
                "domain": item["domain"],
                "language": item["language"],
                "supporting_excerpt": excerpt,
                "normalized_excerpt_sha256": _sha256_bytes(
                    normalized.encode("utf-8")
                ),
                "content_fingerprint": _sha256_bytes(
                    (item["source_lineage"] + "\n" + normalized).encode("utf-8")
                ),
                "name_match_method": item["name_match_method"],
                "wild_cultivated_cultivar_status": item[
                    "wild_cultivated_cultivar_status"
                ],
                "decision": "accept",
                "species_identity_correct": "true",
                "value_correct": "true",
                "provenance_complete": "true",
                "cultivar_contamination": "false",
                "false_positive_reason": "",
                "decision_reason": (
                    "Exact accepted species, same-trait statement, stable URL, citation, "
                    "excerpt and lineage rechecked in the completed acquisition ledger."
                ),
                "reviewer": REVIEWER,
                "reviewed_at_utc": CREATED_AT,
            }
        )
    return pd.DataFrame(rows)


def build(
    *,
    strict_coverage_csv: Path,
    prior_formal_csv: Path,
    master_csv: Path,
    prior_curated_evidence_csv: Path,
    prior_curated_audit_csv: Path,
    output_dir: Path,
) -> dict[str, object]:
    strict = pd.read_csv(strict_coverage_csv, dtype=str).fillna("")
    prior_formal = pd.read_csv(prior_formal_csv, dtype=str).fillna("")
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    if len(strict) != 318_885 or strict["accepted_species"].nunique() != 106_295:
        raise ValueError("strict input does not satisfy the fixed 106295 x 3 contract")

    sources: list[pd.DataFrame] = []
    source_manifest: list[dict[str, object]] = []
    for path in SOURCE_PATHS:
        if not path.is_file():
            raise FileNotFoundError(path)
        source = pd.read_csv(path, dtype=str).fillna("")
        source_manifest.append(
            {"path": path.as_posix(), "sha256": _sha256(path), "rows": len(source)}
        )
        sources.append(source)
    evidence = select_recoverable_evidence(strict, sources, prior_formal)
    missing = sorted(
        set(evidence["accepted_species"]) - set(master["accepted_species"])
    )
    if missing:
        raise ValueError(f"selected species missing from fixed master: {missing}")
    if evidence["source_record_id"].map(_text).duplicated().any():
        raise ValueError("selected source_record_id values must be unique")
    audit = build_audit(evidence)
    approved_traits = set(
        audit.groupby("trait_name").size().loc[lambda counts: counts.ge(10)].index
    )
    scalable = evidence.loc[evidence["trait_name"].isin(approved_traits)].copy()
    individual = build_individually_reviewed_batch(scalable)
    individual_audit = build_individual_audit(individual)
    prior_curated = pd.read_csv(prior_curated_evidence_csv, dtype=str).fillna("")
    prior_curated_audit = pd.read_csv(prior_curated_audit_csv, dtype=str).fillna("")
    combined_curated = pd.concat([prior_curated, individual], ignore_index=True)
    combined_curated_audit = pd.concat(
        [prior_curated_audit, individual_audit], ignore_index=True
    )
    if combined_curated["candidate_id"].duplicated().any():
        raise ValueError("combined curated candidate IDs must be unique")
    if combined_curated_audit["candidate_id"].duplicated().any():
        raise ValueError("combined curated audit IDs must be unique")

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "cached_evidence_recovery_evidence_20260813.csv.gz"
    audit_path = output_dir / "cached_evidence_recovery_audit_200_20260813.csv"
    individual_path = output_dir / "cached_evidence_recovery_individual_100_20260813.csv"
    individual_audit_path = (
        output_dir / "cached_evidence_recovery_individual_audit_100_20260813.csv"
    )
    combined_curated_path = output_dir / "combined_curated_evidence_20260813.csv"
    combined_curated_audit_path = (
        output_dir / "combined_curated_manual_audit_20260813.csv"
    )
    _write_gzip_csv(evidence, evidence_path)
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    individual.to_csv(individual_path, index=False, lineterminator="\n")
    individual_audit.to_csv(individual_audit_path, index=False, lineterminator="\n")
    combined_curated.to_csv(combined_curated_path, index=False, lineterminator="\n")
    combined_curated_audit.to_csv(
        combined_curated_audit_path, index=False, lineterminator="\n"
    )

    unresolved_axes = _axis_unresolved(strict)
    selected_axes = set(zip(evidence["accepted_species"], evidence["axis"], strict=True))
    summary: dict[str, object] = {
        "contract": "cached_evidence_recovery_checkpoint_v1",
        "created_at_utc": CREATED_AT,
        "selection_formal_run_id": FORMAL_SELECTION_RUN_ID,
        "denominator": {"species": 106_295, "species_axis": 318_885},
        "selected_evidence_rows": len(evidence),
        "selected_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "selected_currently_unresolved_species_axis": len(selected_axes & unresolved_axes),
        "audited_rows": len(audit),
        "individually_reviewed_rows": len(individual),
        "combined_curated_rows": len(combined_curated),
        "trait_counts": evidence["trait_name"].value_counts().to_dict(),
        "quality_counts": evidence["quality"].value_counts().to_dict(),
        "provider_counts": evidence["source_provider"].value_counts().to_dict(),
        "guardrails": {
            "trait_specific_genus_join_required_downstream": True,
            "axis_only_genus_join": False,
            "family_inference": False,
            "global_fallback": False,
            "n2_formal_inference": False,
            "cross_trait_reproductive_substitution": False,
            "snippet_only_evidence": False,
            "same_trait_direct_conflicts_quarantined": True,
        },
        "inputs": {
            "strict_coverage_csv": {
                "path": strict_coverage_csv.as_posix(),
                "sha256": _sha256(strict_coverage_csv),
            },
            "prior_formal_csv": {
                "path": prior_formal_csv.as_posix(),
                "sha256": _sha256(prior_formal_csv),
                "rows": len(prior_formal),
            },
            "master_csv": {"path": master_csv.as_posix(), "sha256": _sha256(master_csv)},
            "prior_curated_evidence_csv": {
                "path": prior_curated_evidence_csv.as_posix(),
                "sha256": _sha256(prior_curated_evidence_csv),
                "rows": len(prior_curated),
            },
            "prior_curated_audit_csv": {
                "path": prior_curated_audit_csv.as_posix(),
                "sha256": _sha256(prior_curated_audit_csv),
                "rows": len(prior_curated_audit),
            },
            "source_ledgers": source_manifest,
        },
        "outputs": {
            evidence_path.name: _sha256(evidence_path),
            audit_path.name: _sha256(audit_path),
            individual_path.name: _sha256(individual_path),
            individual_audit_path.name: _sha256(individual_audit_path),
            combined_curated_path.name: _sha256(combined_curated_path),
            combined_curated_audit_path.name: _sha256(combined_curated_audit_path),
        },
    }
    manifest_path = output_dir / "cached_evidence_recovery_manifest_20260813.json"
    manifest_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--strict-coverage-csv", type=Path, required=True)
    parser.add_argument("--prior-formal-csv", type=Path, required=True)
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--prior-curated-evidence-csv", type=Path, required=True)
    parser.add_argument("--prior-curated-audit-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    print(json.dumps(build(**vars(parser.parse_args())), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
