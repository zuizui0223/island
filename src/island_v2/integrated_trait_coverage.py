"""Integrate audited direct and validated-low evidence on the three trait axes.

The coverage ledger is deliberately narrower than the general trait-fill
cascade.  It accepts only species/synonym-direct High or Medium evidence plus
the separately validated genus-consensus Low artifact.  Family inference and
global fallback are rejected.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Iterable
from pathlib import Path
from typing import Annotated, Any
from urllib.parse import parse_qsl, urlencode, urlsplit, urlunsplit

import pandas as pd
import typer

app = typer.Typer(add_completion=False)

AXES = (
    "flower_colour",
    "floral_structural_complexity",
    "reproductive_assurance",
)
AXIS_TRAITS = {
    "flower_colour": {
        "flower_primary_color",
        "flower_color",
    },
    "floral_structural_complexity": {
        "floral_form",
        "floral_symmetry",
        "tube_depth_class",
        "flower_size_class",
        "inflorescence_display",
        "flower_shape",
    },
    "reproductive_assurance": {
        "self_incompatibility",
        "autonomous_selfing_capacity",
        "mating_system",
        "cleistogamy",
    },
}
TRAIT_TO_AXIS = {
    trait: axis for axis, traits in AXIS_TRAITS.items() for trait in traits
}
QUALITY_RANK = {"low": 1, "medium": 2, "high": 3}
DIRECT_SCOPES = {"species_direct", "synonym_direct"}
DIRECT_NAME_MATCHES = {
    "",
    "accepted_name_exact",
    "accepted_name_author_stripped",
    "exact_accepted_name_match",
    "gbif_exact_synonym_to_accepted",
    "exact_synonym_to_accepted",
}
GENERIC_BULK_SOURCES = {
    "austraits_all_target_floral_v7.0.0": "austraits",
    "eol_traitbank_all_zenodo_13305577": "eol_traitbank",
    "gift_v3_2_direct": "gift",
    "meyer_galloway_eckert_2026_dryad": "meyer_2026",
    "ferrer_2024_si_database": "ferrer_2024_si",
    "razanajatovo_etal_2016_nature_communications": "razanajatovo_2016",
    "usda_plants_current": "usda_plants",
}
BULK_SOURCE_KINDS = set(GENERIC_BULK_SOURCES) | {
    "efloras",
    "combined_species_direct_si_sc_v1_genus_inference",
}
BULK_ORDER = (
    "eol_traitbank",
    "austraits",
    "efloras",
    "gift",
    "usda_plants",
    "meyer_2026",
    "ferrer_2024_si",
    "razanajatovo_2016",
    "promoted_public_web",
)
THREE_PASS_SHARD_COUNT = 128
EVIDENCE_COLUMNS = [
    "accepted_species",
    "axis",
    "trait_name",
    "normalized_value",
    "quality",
    "source_group",
    "source_provider",
    "source_url",
    "source_record_id",
    "source_citation",
    "source_excerpt",
    "evidence_scope",
    "name_match_method",
    "source_lineage",
    "lineage_method",
    "source_run_id",
    "source_artifact",
    "source_file",
    "acceptance_contract",
]
AUDIT_COLUMNS = [
    "source_group",
    "source_file",
    "row_number",
    "accepted_species",
    "trait_name",
    "normalized_value",
    "accepted",
    "reason",
]
DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:a-z0-9]+", re.IGNORECASE)


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _truthy(value: object) -> bool:
    return _text(value).casefold() in {"1", "true", "yes", "y"}


def _read_many(paths: Iterable[Path]) -> pd.DataFrame:
    paths = list(paths)
    if not paths:
        raise ValueError("no input files matched")
    return pd.concat(
        [pd.read_csv(path, dtype=str).fillna("") for path in paths],
        ignore_index=True,
        sort=False,
    ).fillna("")


def _canonical_url(value: object) -> str:
    text = _text(value)
    if not text:
        return ""
    try:
        parts = urlsplit(text)
    except ValueError:
        return text.casefold()
    host = parts.netloc.casefold()
    host = host.removeprefix("www.")
    query = urlencode(
        [
            (key, item)
            for key, item in parse_qsl(parts.query, keep_blank_values=True)
            if not key.casefold().startswith("utm_")
        ]
    )
    path = re.sub(r"/+$", "", parts.path) or "/"
    return urlunsplit((parts.scheme.casefold() or "https", host, path, query, ""))


def _citation_fingerprint(value: object) -> str:
    text = _text(value).casefold()
    text = re.sub(r"https?://\S+", " ", text)
    text = re.sub(r"\baccessed\b.*$", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return " ".join(text.split())


def source_lineage(record: dict[str, Any]) -> tuple[str, str]:
    """Return an original-source lineage, not merely the redistributing provider."""
    citation = _text(record.get("source_citation"))
    url = _text(record.get("source_url"))
    excerpt = _text(record.get("source_excerpt"))
    source_group = _text(record.get("source_group"))
    record_id = _text(record.get("source_record_id"))
    combined = f"{citation} {url} {excerpt}"
    if source_group == "austraits" and record_id.startswith("austraits:"):
        dataset = record_id.split(":", 2)[1]
        return f"origin:austraits:{dataset.casefold()}", "provider_component_dataset"
    doi = DOI_RE.search(combined)
    if doi:
        return f"doi:{doi.group(0).rstrip('.,;)').casefold()}", "doi"
    folded = combined.casefold()
    if (
        ("plants database" in folded or "plants.usda.gov" in folded)
        and ("united states department of agriculture" in folded or "usda" in folded)
    ):
        return "origin:usda_plants", "known_origin_alias"
    citation_key = _citation_fingerprint(citation)
    if citation_key and source_group in {"gift", "eol_traitbank"}:
        digest = hashlib.sha256(citation_key.encode("utf-8")).hexdigest()[:20]
        return f"citation:{digest}", "citation_fingerprint"
    canonical_url = _canonical_url(url)
    if canonical_url:
        return f"url:{canonical_url}", "canonical_url"
    fallback = f"{source_group}|{record_id}|{citation_key}|{excerpt}"
    digest = hashlib.sha256(fallback.encode("utf-8")).hexdigest()[:20]
    return f"record:{digest}", "record_fallback"


def _axis(trait_name: object) -> str:
    return TRAIT_TO_AXIS.get(_text(trait_name), "")


def _manifest_entry(
    manifest: dict[str, Any], source_group: str, source_file: Path
) -> tuple[str, str]:
    entries = manifest.get("sources", [])
    matching = [
        entry
        for entry in entries
        if _text(entry.get("source_group")) == source_group
    ]
    if not matching:
        matching = [
            entry
            for entry in entries
            if _text(entry.get("role")) == "bulk_materialization"
        ]
    path_text = str(source_file)
    path_matching = [
        entry
        for entry in matching
        if _text(entry.get("run_id")) in path_text
    ]
    if path_matching:
        matching = path_matching
    if not matching:
        return "", ""
    run_ids = sorted({_text(entry.get("run_id")) for entry in matching if entry.get("run_id")})
    artifacts = sorted(
        {
            _text(entry.get("artifact_name"))
            for entry in matching
            if entry.get("artifact_name")
        }
    )
    return "|".join(run_ids), "|".join(artifacts)


def _make_evidence(
    *,
    species: object,
    trait_name: object,
    value: object,
    quality: str,
    source_group: str,
    provider: object,
    url: object,
    record_id: object,
    citation: object,
    excerpt: object,
    scope: object,
    name_match: object,
    source_file: Path,
    contract: str,
    manifest: dict[str, Any],
    explicit_lineage: object = "",
) -> dict[str, str]:
    run_id, artifact = _manifest_entry(manifest, source_group, source_file)
    row = {
        "accepted_species": _text(species),
        "axis": _axis(trait_name),
        "trait_name": _text(trait_name),
        "normalized_value": _text(value),
        "quality": quality,
        "source_group": source_group,
        "source_provider": _text(provider),
        "source_url": _text(url),
        "source_record_id": _text(record_id),
        "source_citation": _text(citation),
        "source_excerpt": _text(excerpt),
        "evidence_scope": _text(scope),
        "name_match_method": _text(name_match),
        "source_run_id": run_id,
        "source_artifact": artifact,
        "source_file": str(source_file),
        "acceptance_contract": contract,
    }
    lineage = _text(explicit_lineage)
    if lineage:
        method = "candidate_explicit_source_lineage"
    else:
        lineage, method = source_lineage(row)
    row["source_lineage"] = lineage
    row["lineage_method"] = method
    return row


def load_three_pass(
    high_dir: Path,
    medium_dir: Path,
    unresolved_dir: Path,
    manifest: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    high_paths = sorted(high_dir.glob("**/accepted_high_evidence.csv.gz"))
    medium_paths = sorted(medium_dir.glob("**/accepted_medium_evidence.csv.gz"))
    unresolved_paths = sorted(unresolved_dir.glob("**/unresolved_after_pass3.csv.gz"))
    expected = (THREE_PASS_SHARD_COUNT,) * 3
    if (len(high_paths), len(medium_paths), len(unresolved_paths)) != expected:
        raise ValueError(
            f"three-pass artifact counts must be {expected}, got "
            f"{len(high_paths)}/{len(medium_paths)}/{len(unresolved_paths)}"
        )
    evidence: list[dict[str, str]] = []
    for paths, quality, group in (
        (high_paths, "high", "three_pass_high"),
        (medium_paths, "medium", "three_pass_medium"),
    ):
        for path in paths:
            frame = pd.read_csv(path, dtype=str).fillna("")
            for record in frame.to_dict("records"):
                value = record.get("reported_value") or record.get("raw_candidate_value")
                if (
                    _text(record.get("accepted_species"))
                    and _axis(record.get("trait_name"))
                    and _text(value)
                    and _text(record.get("evidence_scope")) in DIRECT_SCOPES
                ):
                    evidence.append(
                        _make_evidence(
                            species=record.get("accepted_species"),
                            trait_name=record.get("trait_name"),
                            value=value,
                            quality=quality,
                            source_group=group,
                            provider=record.get("provider"),
                            url=record.get("source_url"),
                            record_id=record.get("stable_source_id"),
                            citation=record.get("source_title"),
                            excerpt=record.get("excerpt"),
                            scope=record.get("evidence_scope"),
                            name_match="",
                            source_file=path,
                            contract=f"{group}_accepted_artifact_v1",
                            manifest=manifest,
                        )
                    )
    unresolved = _read_many(unresolved_paths)
    unresolved = unresolved.loc[
        unresolved["accepted_species"].map(_text).ne("")
        & unresolved["axis"].isin(AXES),
        ["accepted_species", "axis"],
    ].drop_duplicates()
    return pd.DataFrame(evidence, columns=EVIDENCE_COLUMNS), unresolved


def load_direct_artifact(
    path: Path,
    *,
    source_group: str,
    manifest: dict[str, Any],
) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str).fillna("")
    evidence: list[dict[str, str]] = []
    for record in frame.to_dict("records"):
        quality = _text(record.get("evidence_quality")).casefold()
        value = record.get("reported_value") or record.get("raw_candidate_value")
        if (
            quality not in {"high", "medium"}
            or _text(record.get("evidence_scope")) not in DIRECT_SCOPES
            or not _text(record.get("accepted_species"))
            or not _axis(record.get("trait_name"))
            or not _text(value)
        ):
            continue
        evidence.append(
            _make_evidence(
                species=record.get("accepted_species"),
                trait_name=record.get("trait_name"),
                value=value,
                quality=quality,
                source_group=source_group,
                provider=record.get("provider"),
                url=record.get("source_url"),
                record_id=record.get("stable_source_id"),
                citation=record.get("source_title"),
                excerpt=record.get("excerpt"),
                scope=record.get("evidence_scope"),
                name_match="",
                source_file=path,
                contract="species_direct_high_medium_artifact_v1",
                manifest=manifest,
            )
        )
    return pd.DataFrame(evidence, columns=EVIDENCE_COLUMNS)


def load_validated_low(
    root: Path, manifest: dict[str, Any]
) -> pd.DataFrame:
    paths = sorted(root.glob("**/accepted_validated_low_evidence.csv.gz"))
    if len(paths) != 3:
        raise ValueError(f"validated Low must contain three axis files, got {len(paths)}")
    evidence: list[dict[str, str]] = []
    for path in paths:
        frame = pd.read_csv(path, dtype=str).fillna("")
        for row_number, record in enumerate(frame.to_dict("records"), start=2):
            if _truthy(record.get("family_inference_used")):
                raise ValueError(f"family inference in {path}:{row_number}")
            if _truthy(record.get("global_fallback_used")):
                raise ValueError(f"global fallback in {path}:{row_number}")
            if (
                _text(record.get("evidence_quality")).casefold() != "low"
                or _text(record.get("inference_method")) != "validated_genus_consensus"
                or _text(record.get("axis")) not in AXES
            ):
                continue
            row = _make_evidence(
                species=record.get("accepted_species"),
                trait_name=record.get("trait_name"),
                value=record.get("normalized_value"),
                quality="low",
                source_group="validated_low",
                provider="validated_genus_consensus",
                url="",
                record_id=(
                    f"{record.get('genus')}:{record.get('axis')}:"
                    f"{record.get('trait_name')}:{record.get('normalized_value')}"
                ),
                citation="",
                excerpt=(
                    f"n_direct_species={record.get('n_direct_species')};"
                    f"masked_accuracy={record.get('masked_accuracy')}"
                ),
                scope=record.get("evidence_scope"),
                name_match="",
                source_file=path,
                contract="validated_low_genus_consensus_v1",
                manifest=manifest,
            )
            # A validated-low lineage identifies the rule, not an original source.
            row["source_lineage"] = f"validated-low:{row['source_record_id']}"
            row["lineage_method"] = "validated_rule"
            evidence.append(row)
    return pd.DataFrame(evidence, columns=EVIDENCE_COLUMNS)


def _candidate_rejection_reason(record: dict[str, Any]) -> str:
    if _text(record.get("candidate_kind")) != "source_backed":
        return "candidate_kind_not_source_backed"
    if _text(record.get("evidence_scope")) not in DIRECT_SCOPES:
        return "not_species_or_synonym_direct"
    if _text(record.get("name_match_method")) not in DIRECT_NAME_MATCHES:
        return "name_match_not_exact_or_author_stripped"
    if _text(record.get("inference_rule")):
        return "inference_rule_present"
    if _text(record.get("confidence")).casefold() not in {"high", "medium"}:
        return "quality_not_high_or_medium"
    if not _axis(record.get("trait_name")):
        return "trait_outside_three_axes"
    value = record.get("standardized_value") or record.get("raw_value")
    if not _text(record.get("accepted_species")) or not _text(value):
        return "missing_species_or_value"
    if _text(value).casefold() in {"unknown", "not_reported", "nan"}:
        return "unknown_value"
    required = (
        record.get("source_url"),
        record.get("source_record_id") or record.get("evidence_id"),
        record.get("evidence_excerpt"),
        record.get("source_citation"),
    )
    if not all(_text(value) for value in required):
        return "missing_source_provenance"
    return ""


def load_generic_bulk_candidates(
    path: Path, manifest: dict[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = pd.read_csv(path, dtype=str).fillna("")
    evidence: list[dict[str, str]] = []
    audit: list[dict[str, str]] = []
    for row_number, record in enumerate(frame.to_dict("records"), start=2):
        provider = _text(record.get("source_name"))
        source_group = GENERIC_BULK_SOURCES.get(provider, "unknown_bulk")
        reason = (
            "unapproved_bulk_source_contract"
            if source_group == "unknown_bulk"
            else _candidate_rejection_reason(record)
        )
        value = record.get("standardized_value") or record.get("raw_value")
        audit.append(
            {
                "source_group": source_group,
                "source_file": str(path),
                "row_number": str(row_number),
                "accepted_species": _text(record.get("accepted_species")),
                "trait_name": _text(record.get("trait_name")),
                "normalized_value": _text(value),
                "accepted": str(not reason),
                "reason": reason or "accepted_by_strict_source_contract",
            }
        )
        if reason:
            continue
        quality = _text(record.get("confidence")).casefold()
        evidence.append(
            _make_evidence(
                species=record.get("accepted_species"),
                trait_name=record.get("trait_name"),
                value=value,
                quality=quality,
                source_group=source_group,
                provider=provider,
                url=record.get("source_url"),
                record_id=record.get("source_record_id") or record.get("evidence_id"),
                citation=record.get("source_citation"),
                excerpt=record.get("evidence_excerpt"),
                scope=record.get("evidence_scope"),
                name_match=record.get("name_match_method"),
                source_file=path,
                contract="strict_bulk_species_direct_v1",
                manifest=manifest,
                explicit_lineage=record.get("source_lineage"),
            )
        )
    return (
        pd.DataFrame(evidence, columns=EVIDENCE_COLUMNS),
        pd.DataFrame(audit, columns=AUDIT_COLUMNS),
    )


def load_efloras_candidates(
    path: Path, manifest: dict[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = pd.read_csv(path, dtype=str).fillna("")
    evidence: list[dict[str, str]] = []
    audit: list[dict[str, str]] = []
    for row_number, record in enumerate(frame.to_dict("records"), start=2):
        reason = ""
        if _text(record.get("evidence_scope")) != "species_direct":
            reason = "not_species_direct"
        elif _text(record.get("evidence_status")) != "source_backed_rule_extracted_unreviewed":
            reason = "unexpected_evidence_status"
        elif _text(record.get("name_match_method")) not in DIRECT_NAME_MATCHES:
            reason = "name_match_not_exact_or_author_stripped"
        elif not _axis(record.get("trait_name")):
            reason = "trait_outside_three_axes"
        elif not all(
            _text(record.get(column))
            for column in (
                "accepted_species",
                "candidate_value",
                "source_url",
                "source_record_id",
                "source_excerpt",
            )
        ):
            reason = "missing_species_value_or_provenance"
        audit.append(
            {
                "source_group": "efloras",
                "source_file": str(path),
                "row_number": str(row_number),
                "accepted_species": _text(record.get("accepted_species")),
                "trait_name": _text(record.get("trait_name")),
                "normalized_value": _text(record.get("candidate_value")),
                "accepted": str(not reason),
                "reason": reason or "accepted_by_strict_efloras_contract",
            }
        )
        if reason:
            continue
        evidence.append(
            _make_evidence(
                species=record.get("accepted_species"),
                trait_name=record.get("trait_name"),
                value=record.get("candidate_value"),
                quality="medium",
                source_group="efloras",
                provider=f"efloras:{record.get('flora_id')}",
                url=record.get("source_url"),
                record_id=record.get("source_record_id"),
                citation=record.get("flora_name"),
                excerpt=record.get("source_excerpt"),
                scope=record.get("evidence_scope"),
                name_match=record.get("name_match_method"),
                source_file=path,
                contract="strict_efloras_species_direct_v1",
                manifest=manifest,
            )
        )
    return (
        pd.DataFrame(evidence, columns=EVIDENCE_COLUMNS),
        pd.DataFrame(audit, columns=AUDIT_COLUMNS),
    )


def load_promoted_public_web(
    path: Path, manifest: dict[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = pd.read_csv(path, dtype=str).fillna("")
    evidence: list[dict[str, str]] = []
    audit: list[dict[str, str]] = []
    for row_number, record in enumerate(frame.to_dict("records"), start=2):
        source_kind = _text(record.get("source_kind"))
        source_backed = _truthy(record.get("source_backed"))
        quality = _text(record.get("confidence")).casefold()
        reason = ""
        if source_kind in BULK_SOURCE_KINDS:
            reason = "handled_by_bulk_source_contract"
        elif quality not in {"high", "medium"}:
            reason = "not_high_or_medium"
        elif not source_backed or _text(record.get("evidence_type")) == "inference":
            reason = "not_source_backed_direct_evidence"
        elif not _axis(record.get("field")):
            reason = "field_outside_three_axes"
        elif not all(
            _text(record.get(column))
            for column in (
                "species",
                "value",
                "source_url",
                "source_record_id",
                "source_excerpt",
            )
        ):
            reason = "missing_species_value_or_provenance"
        audit.append(
            {
                "source_group": "promoted_public_web",
                "source_file": str(path),
                "row_number": str(row_number),
                "accepted_species": _text(record.get("species")),
                "trait_name": _text(record.get("field")),
                "normalized_value": _text(record.get("value")),
                "accepted": str(not reason),
                "reason": reason or "accepted_all_master_source_backed_medium",
            }
        )
        if reason:
            continue
        evidence.append(
            _make_evidence(
                species=record.get("species"),
                trait_name=record.get("field"),
                value=record.get("value"),
                quality=quality,
                source_group="promoted_public_web",
                provider=source_kind,
                url=record.get("source_url"),
                record_id=record.get("source_record_id"),
                citation="",
                excerpt=record.get("source_excerpt"),
                scope="species_direct",
                name_match="accepted_species_from_promoted_ledger",
                source_file=path,
                contract="all_master_source_backed_direct_v1",
                manifest=manifest,
            )
        )
    return (
        pd.DataFrame(evidence, columns=EVIDENCE_COLUMNS),
        pd.DataFrame(audit, columns=AUDIT_COLUMNS),
    )


def load_bulk_artifact(
    root: Path, manifest: dict[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    evidence_parts: list[pd.DataFrame] = []
    audit_parts: list[pd.DataFrame] = []
    candidates = sorted(root.glob("**/trait_candidates.csv*"))
    if not candidates:
        raise ValueError(f"no bulk candidate files under {root}")
    for path in candidates:
        if "inference" in {part.casefold() for part in path.parts}:
            continue
        frame = pd.read_csv(path, nrows=1, dtype=str).fillna("")
        if "evidence_status" in frame.columns:
            evidence, audit = load_efloras_candidates(path, manifest)
        else:
            evidence, audit = load_generic_bulk_candidates(path, manifest)
        evidence_parts.append(evidence)
        audit_parts.append(audit)
    ledger_paths = sorted(root.glob("**/all_species_trait_evidence.csv.gz"))
    if len(ledger_paths) != 1:
        raise ValueError(
            f"expected one all_species_trait_evidence.csv.gz under {root}, "
            f"got {len(ledger_paths)}"
        )
    promoted, promoted_audit = load_promoted_public_web(ledger_paths[0], manifest)
    evidence_parts.append(promoted)
    audit_parts.append(promoted_audit)
    return (
        pd.concat(evidence_parts, ignore_index=True).fillna(""),
        pd.concat(audit_parts, ignore_index=True).fillna(""),
    )


def dedupe_source_lineages(
    evidence: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    data = evidence.copy()
    data["quality_rank"] = data["quality"].map(QUALITY_RANK)
    data = data.sort_values(
        [
            "accepted_species",
            "trait_name",
            "source_lineage",
            "quality_rank",
            "source_group",
            "source_record_id",
        ],
        ascending=[True, True, True, False, True, True],
    )
    duplicate_mask = data.duplicated(
        ["accepted_species", "trait_name", "source_lineage"], keep=False
    )
    duplicates = data.loc[duplicate_mask].copy()
    duplicates["lineage_record_count"] = duplicates.groupby(
        ["accepted_species", "trait_name", "source_lineage"]
    )["source_lineage"].transform("size")
    duplicates["lineage_value_count"] = duplicates.groupby(
        ["accepted_species", "trait_name", "source_lineage"]
    )["normalized_value"].transform("nunique")
    deduped = data.drop_duplicates(
        ["accepted_species", "trait_name", "source_lineage"], keep="first"
    ).drop(columns=["quality_rank"])
    return deduped, duplicates.drop(columns=["quality_rank"])


def select_best(evidence: pd.DataFrame) -> pd.DataFrame:
    data = evidence.copy()
    data["quality_rank"] = data["quality"].map(QUALITY_RANK)
    data = data.sort_values(
        [
            "accepted_species",
            "axis",
            "quality_rank",
            "source_group",
            "source_lineage",
        ],
        ascending=[True, True, False, True, True],
    )
    return data.drop_duplicates(["accepted_species", "axis"], keep="first")


def _universe(
    direct: pd.DataFrame,
    unresolved: pd.DataFrame,
    expected_species: int,
) -> pd.DataFrame:
    universe = pd.concat(
        [
            direct[["accepted_species", "axis"]],
            unresolved[["accepted_species", "axis"]],
        ],
        ignore_index=True,
    ).drop_duplicates()
    counts = universe["axis"].value_counts().to_dict()
    if set(counts) != set(AXES) or any(counts[axis] != expected_species for axis in AXES):
        raise ValueError(
            f"universe must be {expected_species} species for every axis, got {counts}"
        )
    if universe["accepted_species"].nunique() != expected_species:
        raise ValueError(
            f"universe has {universe['accepted_species'].nunique()} species, "
            f"expected {expected_species}"
        )
    return universe.sort_values(["accepted_species", "axis"]).reset_index(drop=True)


def _in_universe(evidence: pd.DataFrame, universe: pd.DataFrame) -> pd.DataFrame:
    index = pd.MultiIndex.from_frame(universe[["accepted_species", "axis"]])
    evidence_index = pd.MultiIndex.from_frame(
        evidence[["accepted_species", "axis"]]
    )
    return evidence.loc[evidence_index.isin(index)].copy()


def _coverage_snapshot(
    best: pd.DataFrame, universe: pd.DataFrame
) -> tuple[dict[str, Any], pd.DataFrame]:
    best_keys = pd.MultiIndex.from_frame(best[["accepted_species", "axis"]])
    universe_keys = pd.MultiIndex.from_frame(universe[["accepted_species", "axis"]])
    unresolved = universe.loc[~universe_keys.isin(best_keys)].copy()
    by_axis: dict[str, dict[str, Any]] = {}
    for axis in AXES:
        subset = best.loc[best["axis"].eq(axis)]
        counts = {
            quality: int(subset["quality"].eq(quality).sum())
            for quality in ("high", "medium", "low")
        }
        filled = len(subset)
        by_axis[axis] = {
            "denominator": int((universe["axis"] == axis).sum()),
            "quality_counts": counts,
            "filled_species": filled,
            "fill_rate": filled / int((universe["axis"] == axis).sum()),
            "unresolved_species": int((unresolved["axis"] == axis).sum()),
        }
    per_species = (
        best.groupby("accepted_species")["axis"].nunique().rename("filled_axes")
    )
    species = pd.DataFrame(
        {"accepted_species": sorted(universe["accepted_species"].unique())}
    )
    species = species.join(per_species, on="accepted_species").fillna({"filled_axes": 0})
    species["filled_axes"] = species["filled_axes"].astype(int)
    distribution = {
        str(number): int(species["filled_axes"].eq(number).sum())
        for number in range(4)
    }
    quality_counts = {
        quality: int(best["quality"].eq(quality).sum())
        for quality in ("high", "medium", "low")
    }
    return (
        {
            "quality_counts": quality_counts,
            "filled_species_axis": len(best),
            "unresolved_species_axis": len(unresolved),
            "by_axis": by_axis,
            "species_by_filled_axis_count": distribution,
        },
        unresolved,
    )


def _source_increments(
    before_best: pd.DataFrame,
    bulk: pd.DataFrame,
    universe: pd.DataFrame,
) -> dict[str, dict[str, Any]]:
    before_keys = pd.MultiIndex.from_frame(before_best[["accepted_species", "axis"]])
    current = before_best.copy()
    result: dict[str, dict[str, Any]] = {}
    for source_group in BULK_ORDER:
        source = bulk.loc[bulk["source_group"].eq(source_group)].copy()
        source_best = select_best(source) if len(source) else source
        prior_best = select_best(current)
        if len(source_best):
            source_keys = pd.MultiIndex.from_frame(
                source_best[["accepted_species", "axis"]]
            )
            pure_gain = source_best.loc[~source_keys.isin(before_keys)].copy()
            prior_keys = pd.MultiIndex.from_frame(
                prior_best[["accepted_species", "axis"]]
            )
            sequential_gain = source_best.loc[
                ~source_keys.isin(prior_keys)
            ].copy()
        else:
            pure_gain = source_best.copy()
            sequential_gain = source_best.copy()
        result[source_group] = {
            "accepted_lineages": len(source),
            "unique_species_axis": len(source_best),
            "pure_gain_vs_before": len(pure_gain),
            "pure_gain_vs_before_by_axis": {
                axis: int(pure_gain["axis"].eq(axis).sum()) for axis in AXES
            },
            "sequential_marginal_gain": len(sequential_gain),
            "sequential_marginal_gain_by_axis": {
                axis: int(sequential_gain["axis"].eq(axis).sum()) for axis in AXES
            },
        }
        current = pd.concat([current, source], ignore_index=True)
    return result


def _direct_run_contributions(
    evidence: pd.DataFrame, manifest: dict[str, Any]
) -> dict[str, dict[str, Any]]:
    run_ids = sorted(
        {
            int(entry["run_id"])
            for entry in manifest.get("sources", [])
            if entry.get("role") == "baseline_direct"
        }
    )
    overall = select_best(evidence)
    result: dict[str, dict[str, Any]] = {}
    for run_id in run_ids:
        run_text = str(run_id)
        owned = evidence.loc[evidence["source_run_id"].eq(run_text)].copy()
        without = select_best(
            evidence.loc[~evidence["source_run_id"].eq(run_text)].copy()
        )
        without_quality = without[
            ["accepted_species", "axis", "quality"]
        ].rename(columns={"quality": "without_quality"})
        comparison = overall.merge(
            without_quality, on=["accepted_species", "axis"], how="left"
        )
        comparison["quality_rank"] = comparison["quality"].map(QUALITY_RANK)
        comparison["without_quality_rank"] = comparison["without_quality"].map(
            QUALITY_RANK
        ).fillna(0)
        contributed = comparison.loc[
            comparison["source_run_id"].eq(run_text)
            & (
                comparison["quality_rank"]
                > comparison["without_quality_rank"]
            )
        ].copy()
        pure_gain = contributed.loc[contributed["without_quality"].isna()].copy()
        upgrades = contributed.loc[contributed["without_quality"].notna()].copy()
        result[run_text] = {
            "accepted_lineages": len(owned),
            "unique_species_axis": int(
                owned[["accepted_species", "axis"]].drop_duplicates().shape[0]
            ),
            "pure_gain_species_axis": len(pure_gain),
            "pure_gain_by_axis": {
                axis: int(pure_gain["axis"].eq(axis).sum()) for axis in AXES
            },
            "quality_upgrades": len(upgrades),
            "quality_upgrades_by_transition": {
                f"{old}_to_{new}": int(
                    (
                        upgrades["without_quality"].eq(old)
                        & upgrades["quality"].eq(new)
                    ).sum()
                )
                for old, new in (
                    ("low", "medium"),
                    ("low", "high"),
                    ("medium", "high"),
                )
            },
        }
    return result


def aggregate(
    *,
    three_pass_high_dir: Path,
    three_pass_medium_dir: Path,
    three_pass_unresolved_dir: Path,
    expanded_artifact_dirs: list[Path],
    targeted_dir: Path,
    validated_low_dir: Path,
    all_master_dir: Path,
    source_manifest: dict[str, Any],
    expected_species: int = 106_295,
) -> tuple[dict[str, Any], dict[str, pd.DataFrame]]:
    three_pass, unresolved = load_three_pass(
        three_pass_high_dir,
        three_pass_medium_dir,
        three_pass_unresolved_dir,
        source_manifest,
    )
    direct_parts = [three_pass]
    for root in expanded_artifact_dirs:
        candidates = sorted(
            list(root.glob("**/direct_expanded_trait_evidence.csv.gz"))
            + list(root.glob("**/direct_expanded_trait_pilot_evidence.csv.gz"))
        )
        if len(candidates) != 1:
            raise ValueError(f"expected one expanded evidence file under {root}")
        direct_parts.append(
            load_direct_artifact(
                candidates[0],
                source_group="expanded_wave",
                manifest=source_manifest,
            )
        )
    targeted_paths = sorted(targeted_dir.glob("**/targeted_new_evidence.csv.gz"))
    if len(targeted_paths) != 1:
        raise ValueError(f"expected one targeted evidence file under {targeted_dir}")
    direct_parts.append(
        load_direct_artifact(
            targeted_paths[0],
            source_group="targeted_reacquisition",
            manifest=source_manifest,
        )
    )
    direct_parts.append(load_validated_low(validated_low_dir, source_manifest))
    before_raw = pd.concat(direct_parts, ignore_index=True).fillna("")
    universe = _universe(three_pass, unresolved, expected_species)
    before_in_universe = _in_universe(before_raw, universe)
    before_lineages, before_duplicates = dedupe_source_lineages(before_in_universe)
    before_best = select_best(before_lineages)

    bulk_raw, bulk_audit = load_bulk_artifact(all_master_dir, source_manifest)
    bulk_in_universe = _in_universe(bulk_raw, universe)
    bulk_lineages, bulk_duplicates = dedupe_source_lineages(bulk_in_universe)
    all_lineages, cross_duplicates = dedupe_source_lineages(
        pd.concat([before_lineages, bulk_lineages], ignore_index=True)
    )
    after_best = select_best(all_lineages)

    before_summary, before_unresolved = _coverage_snapshot(before_best, universe)
    after_summary, after_unresolved = _coverage_snapshot(after_best, universe)
    before_quality = before_best[
        ["accepted_species", "axis", "quality"]
    ].rename(columns={"quality": "before_quality"})
    transitions = after_best.merge(
        before_quality, on=["accepted_species", "axis"], how="left"
    )
    low_upgrades = transitions.loc[
        transitions["before_quality"].eq("low")
        & transitions["quality"].isin({"high", "medium"})
    ].copy()
    net_by_axis: dict[str, dict[str, Any]] = {}
    for axis in AXES:
        before_axis = before_summary["by_axis"][axis]
        after_axis = after_summary["by_axis"][axis]
        net_by_axis[axis] = {
            "filled_species": (
                after_axis["filled_species"] - before_axis["filled_species"]
            ),
            "fill_rate": after_axis["fill_rate"] - before_axis["fill_rate"],
            "quality_counts": {
                quality: (
                    after_axis["quality_counts"][quality]
                    - before_axis["quality_counts"][quality]
                )
                for quality in ("high", "medium", "low")
            },
            "validated_low_upgraded_to_direct": int(
                low_upgrades["axis"].eq(axis).sum()
            ),
        }
    coverage = universe.merge(
        before_best[
            [
                "accepted_species",
                "axis",
                "quality",
                "normalized_value",
                "source_group",
                "source_lineage",
            ]
        ].rename(
            columns={
                "quality": "before_quality",
                "normalized_value": "before_value",
                "source_group": "before_source_group",
                "source_lineage": "before_source_lineage",
            }
        ),
        on=["accepted_species", "axis"],
        how="left",
    ).merge(
        after_best[
            [
                "accepted_species",
                "axis",
                "quality",
                "normalized_value",
                "source_group",
                "source_lineage",
            ]
        ].rename(
            columns={
                "quality": "after_quality",
                "normalized_value": "after_value",
                "source_group": "after_source_group",
                "source_lineage": "after_source_lineage",
            }
        ),
        on=["accepted_species", "axis"],
        how="left",
    ).fillna("")
    summary = {
        "contract": "integrated_species_axis_coverage_v1",
        "denominator": {
            "species": expected_species,
            "axes": len(AXES),
            "species_axis": expected_species * len(AXES),
            "axis_denominators": {axis: expected_species for axis in AXES},
        },
        "quality_precedence": ["high", "medium", "validated_low"],
        "acceptance_policy": {
            "bulk_pending_auto_promoted": False,
            "bulk_contract": (
                "approved source, exact accepted/synonym mapping, species-direct, "
                "mapped controlled value, complete provenance, no inference rule"
            ),
            "family_inference_allowed": False,
            "global_fallback_allowed": False,
            "source_lineage_deduplicated": True,
        },
        "source_manifest": source_manifest,
        "before_bulk_integration": before_summary,
        "after_bulk_integration": after_summary,
        "net_change": {
            "filled_species_axis": (
                after_summary["filled_species_axis"]
                - before_summary["filled_species_axis"]
            ),
            "unresolved_species_axis": (
                after_summary["unresolved_species_axis"]
                - before_summary["unresolved_species_axis"]
            ),
            "quality_counts": {
                quality: (
                    after_summary["quality_counts"][quality]
                    - before_summary["quality_counts"][quality]
                )
                for quality in ("high", "medium", "low")
            },
            "by_axis": net_by_axis,
            "validated_low_upgraded_to_direct": len(low_upgrades),
        },
        "source_increments": _source_increments(
            before_best, bulk_lineages, universe
        ),
        "direct_run_contributions": _direct_run_contributions(
            before_lineages, source_manifest
        ),
        "input_audit": {
            "baseline_evidence_rows": len(before_raw),
            "baseline_lineages": len(before_lineages),
            "bulk_candidate_rows_audited": len(bulk_audit),
            "bulk_candidate_rows_contract_accepted": int(
                bulk_audit["accepted"].eq("True").sum()
            ),
            "bulk_candidate_rows_rejected": int(
                bulk_audit["accepted"].ne("True").sum()
            ),
            "bulk_accepted_lineages": len(bulk_lineages),
            "lineage_duplicate_rows": len(cross_duplicates),
            "lineage_duplicate_groups": int(
                cross_duplicates[
                    ["accepted_species", "trait_name", "source_lineage"]
                ].drop_duplicates().shape[0]
            ),
            "out_of_universe_bulk_rows": int(len(bulk_raw) - len(bulk_in_universe)),
        },
        "checks": {
            "denominator_exact": len(universe) == expected_species * 3,
            "quality_counts_sum_before": (
                sum(before_summary["quality_counts"].values())
                == before_summary["filled_species_axis"]
            ),
            "quality_counts_sum_after": (
                sum(after_summary["quality_counts"].values())
                == after_summary["filled_species_axis"]
            ),
            "before_filled_plus_unresolved": (
                before_summary["filled_species_axis"]
                + before_summary["unresolved_species_axis"]
                == expected_species * 3
            ),
            "after_filled_plus_unresolved": (
                after_summary["filled_species_axis"]
                + after_summary["unresolved_species_axis"]
                == expected_species * 3
            ),
            "family_inference_absent": True,
            "global_fallback_absent": True,
        },
    }
    frames = {
        "coverage": coverage,
        "evidence_lineages": all_lineages,
        "lineage_duplicates": cross_duplicates,
        "candidate_audit": bulk_audit,
        "before_unresolved": before_unresolved,
        "after_unresolved": after_unresolved,
        "validated_low_upgrades": low_upgrades,
        "before_duplicates": before_duplicates,
        "bulk_duplicates": bulk_duplicates,
    }
    return summary, frames


def write_outputs(
    summary: dict[str, Any],
    frames: dict[str, pd.DataFrame],
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "integrated_trait_coverage_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (output_dir / "source_run_manifest.json").write_text(
        json.dumps(
            summary["source_manifest"],
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    names = {
        "coverage": "integrated_species_axis_coverage.csv.gz",
        "evidence_lineages": "integrated_evidence_lineage.csv.gz",
        "lineage_duplicates": "source_lineage_duplicates.csv.gz",
        "candidate_audit": "bulk_candidate_acceptance_audit.csv.gz",
        "before_unresolved": "unresolved_before_bulk.csv.gz",
        "after_unresolved": "unresolved_after_bulk.csv.gz",
        "validated_low_upgrades": "validated_low_upgraded_by_direct.csv.gz",
    }
    for key, name in names.items():
        frames[key].to_csv(output_dir / name, index=False)


@app.command()
def build(
    three_pass_high_dir: Annotated[
        Path, typer.Option(exists=True, file_okay=False)
    ],
    three_pass_medium_dir: Annotated[
        Path, typer.Option(exists=True, file_okay=False)
    ],
    three_pass_unresolved_dir: Annotated[
        Path, typer.Option(exists=True, file_okay=False)
    ],
    expanded_artifact_dir: Annotated[
        list[Path], typer.Option(exists=True, file_okay=False)
    ],
    targeted_dir: Annotated[Path, typer.Option(exists=True, file_okay=False)],
    validated_low_dir: Annotated[
        Path, typer.Option(exists=True, file_okay=False)
    ],
    all_master_dir: Annotated[Path, typer.Option(exists=True, file_okay=False)],
    source_manifest_json: Annotated[
        Path, typer.Option(exists=True, dir_okay=False)
    ],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
    expected_species: Annotated[int, typer.Option(min=1)] = 106_295,
) -> None:
    """Build the reproducible before/after species-axis coverage ledger."""
    manifest = json.loads(source_manifest_json.read_text(encoding="utf-8"))
    summary, frames = aggregate(
        three_pass_high_dir=three_pass_high_dir,
        three_pass_medium_dir=three_pass_medium_dir,
        three_pass_unresolved_dir=three_pass_unresolved_dir,
        expanded_artifact_dirs=expanded_artifact_dir,
        targeted_dir=targeted_dir,
        validated_low_dir=validated_low_dir,
        all_master_dir=all_master_dir,
        source_manifest=manifest,
        expected_species=expected_species,
    )
    write_outputs(summary, frames, output_dir)
    typer.echo(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    app()
