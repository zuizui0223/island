"""Build the exact nine-column trait export for every Island master species.

The export is deliberately exhaustive: source-backed database/flora candidates,
public-web rule matches, evidence-locked LLM extractions, and validated LLM-only
inference can contribute values, but no master species is dropped. A companion
long table retains provenance that the requested presentation table cannot carry.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.llm_evidence_extraction import VALIDATED_COLUMNS as LLM_EXTRACTED_COLUMNS
from island_v2.search_enabled_llm_campaign import (
    JOB_HISTORY_COLUMNS,
    TRAIT_FIELDS,
    VALIDATED_RESULT_COLUMNS,
    _job_record,
    _validate_row,
    canonical_result_sha,
)
from island_v2.v1_category_traits import OUTPUT_COLUMNS, validate_result_table

app = typer.Typer(add_completion=False, no_args_is_help=True)

FIELDS = (
    "flower_color",
    "flower_shape",
    "pollination_guild",
    "mating_system",
    "self_incompatibility",
)
REPRODUCTIVE_FIELDS = {"pollination_guild", "mating_system", "self_incompatibility"}
EVIDENCE_ORDER = ("field_study", "review", "flora", "horticulture", "inference")
CONFIDENCE_ORDER = {"low": 0, "medium": 1, "high": 2}
EVIDENCE_COLUMNS = [
    "species",
    "field",
    "value",
    "source_kind",
    "source_url",
    "source_record_id",
    "source_excerpt",
    "evidence_type",
    "confidence",
    "source_rank",
    "source_backed",
]
AUDIT_COLUMNS = [
    "input_file",
    "row_number",
    "species",
    "trait_name",
    "raw_value",
    "reason",
]
LLM_PROVENANCE_COLUMNS = {
    "job_id",
    "job_sha256",
    "target_fields",
    "model_provider",
    "model_id",
    "prompt_version",
    "prompt_sha256",
    "result_sha256",
}
LLM_EXTRACTED_PROVENANCE_COLUMNS = {
    "model_provider",
    "model_id",
    "prompt_version",
    "prompt_sha256",
    "packet_input_sha256",
    "ontology_sha256",
    "source_text_sha256",
    "extraction_run_id",
}

COLOR_MAP = {
    "white": "white",
    "green_brown_inconspicuous": "green/brown/inconspicuous",
    "yellow_orange": "yellow/orange",
    "red_pink": "red/pink",
    "blue_purple": "blue/purple",
    "multicolored_variable": "multicolored/variable",
}
RAW_COLORS = {
    "black",
    "blue",
    "brown",
    "cream",
    "green",
    "orange",
    "pink",
    "purple",
    "red",
    "violet",
    "white",
    "yellow",
}
FORM_MAP = {
    "bilabiate": "bilabiate / two-lipped",
    "campanulate": "bell-shaped",
    "bell_campanulate": "bell-shaped",
    "funnelform": "funnel-shaped",
    "funnel_trumpet": "funnel / trumpet-shaped",
    "rotate": "rotate / open",
    "open_radial": "open / bowl / cup-shaped",
    "salverform": "salverform",
    "spurred": "spurred",
    "tubular": "tubular",
    "urceolate": "urn-shaped",
    "urn_urceolate": "urn-shaped",
    "composite_head": "composite head / capitulum",
    "papilionaceous": "papilionaceous",
}
SYMMETRY_MAP = {
    "radial": "actinomorphic / radially symmetric",
    "actinomorphic": "actinomorphic / radially symmetric",
    "bilateral": "zygomorphic / bilaterally symmetric",
    "zygomorphic": "zygomorphic / bilaterally symmetric",
}
GUILD_MAP = {
    "other_bees": "bees",
    "bee": "bees",
    "bumblebees": "bumblebees",
    "flies": "flies",
    "fly": "flies",
    "butterflies": "butterflies",
    "moths": "moths",
    "butterfly_moth": "mixed",
    "birds": "birds",
    "bird": "birds",
    "mixed_or_generalist": "mixed",
    "self_only": "self",
}
POLLEN_VECTOR_MAP = {
    "abiotic_wind": "wind",
    "wind": "wind",
    "mixed": "mixed",
}
MATING_MAP = {
    "obligate_outcrossing": "obligate_outcrossing",
    "mixed_mating": "mixed_mating",
    "predominantly_selfing": "mainly_selfing",
    "mainly_selfing": "mainly_selfing",
    "obligate_selfing": "obligate_selfing",
}
SI_MAP = {
    "SI": "SI",
    "SC": "SC",
    "likely_SI": "likely_SI",
    "likely_SC": "likely_SC",
    "self_incompatible": "SI",
    "self_compatible": "SC",
}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _is_sha256(value: object) -> bool:
    text = _text(value)
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _load_json_object(path: Path, label: str) -> dict[str, Any]:
    if not path.is_file():
        raise typer.BadParameter(f"{label} does not exist: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise typer.BadParameter(f"{label} is not valid UTF-8 JSON: {path}") from exc
    if not isinstance(value, dict):
        raise typer.BadParameter(f"{label} must contain a JSON object: {path}")
    return value


def _tokens(value: str) -> list[str]:
    return list(dict.fromkeys(token.strip() for token in value.split("|") if token.strip()))


def map_candidate_value(trait_name: str, raw_value: str) -> list[tuple[str, str]]:
    """Map a direct canonical/raw source value to the requested nine-column vocabulary."""
    trait = _text(trait_name)
    value = _text(raw_value)
    tokens = _tokens(value)
    if trait == "flower_primary_color":
        if value in COLOR_MAP:
            return [("flower_color", COLOR_MAP[value])]
        if tokens and set(tokens).issubset(RAW_COLORS):
            return [("flower_color", ", ".join(tokens))]
    if trait == "floral_form":
        mapped = list(dict.fromkeys(FORM_MAP[token] for token in tokens if token in FORM_MAP))
        if mapped and len(mapped) == len(tokens):
            return [("flower_shape", item) for item in mapped]
    if trait == "floral_symmetry":
        mapped = list(
            dict.fromkeys(SYMMETRY_MAP[token] for token in tokens if token in SYMMETRY_MAP)
        )
        if mapped and len(mapped) == len(tokens):
            return [("flower_shape", item) for item in mapped]
    if trait in {"pollination_functional_guild", "pollination_vector_reported"}:
        if value in GUILD_MAP:
            return [("pollination_guild", GUILD_MAP[value])]
    if trait == "pollen_vector_mode" and value in POLLEN_VECTOR_MAP:
        return [("pollination_guild", POLLEN_VECTOR_MAP[value])]
    if trait in {"mating_system", "mating_system_reported"} and value in MATING_MAP:
        return [("mating_system", MATING_MAP[value])]
    if trait in {"self_incompatibility", "self_incompatibility_reported"}:
        if value in SI_MAP:
            return [("self_incompatibility", SI_MAP[value])]
    if trait == "self_compatibility_reported" and value in SI_MAP:
        return [("self_incompatibility", SI_MAP[value])]
    return []


def _source_profile(row: dict[str, str], input_path: Path) -> tuple[str, str, str, int]:
    source_type = _text(row.get("source_type")).casefold()
    citation = _text(row.get("source_citation")).casefold()
    path_text = str(input_path).casefold()
    llm_extracted = _text(row.get("extraction_method")).startswith("provider_neutral_llm_")
    source_kind_prefix = "llm_extracted_" if llm_extracted else ""
    if _text(row.get("flora_id")) or "efloras" in path_text:
        return f"{source_kind_prefix}efloras", "flora", "medium", 0
    if "flora_or_monograph" in source_type or "gbif_species_description" in source_type:
        return f"{source_kind_prefix}{source_type}", "flora", "medium", 0
    if "horticulture" in source_type:
        return f"{source_kind_prefix}{source_type}", "horticulture", "medium", 2
    if "openalex" in source_type:
        if any(term in citation for term in ("review", "meta-analysis", "systematic review")):
            return f"{source_kind_prefix}{source_type}", "review", "medium", 1
        if any(
            term in citation for term in ("field study", "field experiment", "hand pollination")
        ):
            return f"{source_kind_prefix}{source_type}", "field_study", "medium", 1
        return f"{source_kind_prefix}{source_type}", "inference", "low", 3
    if "wikipedia" in source_type or "wikidata" in source_type:
        return f"{source_kind_prefix}{source_type or 'wikimedia'}", "review", "low", 2
    if source_type == "curated_trait_database" or _text(row.get("source_reliability")):
        source_kind = _text(row.get("source_name")) or source_type
        return f"{source_kind_prefix}{source_kind}", "review", "medium", 1
    source_kind = source_type or input_path.stem
    return f"{source_kind_prefix}{source_kind}", "inference", "low", 3


def candidate_evidence(
    frame: pd.DataFrame,
    input_path: Path,
    *,
    allow_validated_llm_extraction: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Convert one materialized source or machine-candidate table to long evidence."""
    required = {"accepted_species", "trait_name"}
    missing = required.difference(frame.columns)
    if missing:
        raise typer.BadParameter(f"{input_path} is missing candidate columns: {sorted(missing)}")
    value_column = next(
        (
            column
            for column in (
                "candidate_value",
                "standardized_value",
                "provisional_candidate_value",
            )
            if column in frame.columns
        ),
        None,
    )
    if value_column is None:
        raise typer.BadParameter(f"{input_path} has no candidate value column")
    llm_extracted = frame.get("extraction_method", pd.Series("", index=frame.index)).map(
        lambda value: _text(value).startswith("provider_neutral_llm_")
    )
    if llm_extracted.any() and not allow_validated_llm_extraction:
        raise typer.BadParameter(
            f"{input_path} contains LLM-extracted rows; pass a validated extraction bundle "
            "with --llm-extracted-bundle instead of generic --candidate-csv"
        )
    evidence: list[dict[str, object]] = []
    audit: list[dict[str, str]] = []
    for row_number, record in enumerate(frame.fillna("").to_dict("records"), start=2):
        species = _text(record.get("accepted_species"))
        trait = _text(record.get("trait_name"))
        raw_value = _text(record.get(value_column))
        mapped = map_candidate_value(trait, raw_value)
        if not mapped:
            audit.append(
                {
                    "input_file": str(input_path),
                    "row_number": str(row_number),
                    "species": species,
                    "trait_name": trait,
                    "raw_value": raw_value,
                    "reason": "not_mappable_to_requested_vocabulary_without_inference",
                }
            )
            continue
        source_kind, evidence_type, confidence, rank = _source_profile(record, input_path)
        supplied_confidence = _text(record.get("confidence"))
        if supplied_confidence in CONFIDENCE_ORDER:
            confidence = min(
                (confidence, supplied_confidence),
                key=lambda item: CONFIDENCE_ORDER[item],
            )
        excerpt = _text(
            record.get("source_excerpt")
            or record.get("evidence_excerpt")
            or record.get("evidence_quote")
            or record.get("raw_description")
        )
        source_url = _text(record.get("source_url"))
        source_record_id = _text(
            record.get("source_record_id")
            or record.get("evidence_id")
            or record.get("source_chunk_id")
            or record.get("extraction_run_id")
            or (
                f"{_text(record.get('first_wave_id'))}:{species}:{trait}:{raw_value}"
                if _text(record.get("first_wave_id"))
                else ""
            )
        )
        source_backed = bool(source_url and excerpt and source_record_id)
        if not source_backed:
            audit.append(
                {
                    "input_file": str(input_path),
                    "row_number": str(row_number),
                    "species": species,
                    "trait_name": trait,
                    "raw_value": raw_value,
                    "reason": "mapped_candidate_missing_url_excerpt_or_record_id",
                }
            )
            evidence_type, confidence, rank = "inference", "low", 4
        for field, value in mapped:
            evidence.append(
                {
                    "species": species,
                    "field": field,
                    "value": value,
                    "source_kind": source_kind,
                    "source_url": source_url,
                    "source_record_id": source_record_id,
                    "source_excerpt": excerpt or f"{trait} = {raw_value}",
                    "evidence_type": evidence_type,
                    "confidence": confidence,
                    "source_rank": rank,
                    "source_backed": source_backed,
                }
            )
    return (
        pd.DataFrame(evidence, columns=EVIDENCE_COLUMNS),
        pd.DataFrame(audit, columns=AUDIT_COLUMNS),
    )


def web_evidence(frame: pd.DataFrame, input_path: Path) -> pd.DataFrame:
    """Normalize the existing public-web long evidence contract."""
    required = {"species", "field", "value"}
    missing = required.difference(frame.columns)
    if missing:
        raise typer.BadParameter(f"{input_path} is missing web evidence columns: {sorted(missing)}")
    rows: list[dict[str, object]] = []
    for record in frame.fillna("").to_dict("records"):
        field = _text(record.get("field"))
        value = _text(record.get("value"))
        if field not in FIELDS or not value or value == "unknown":
            continue
        source_kind, evidence_type, confidence, rank = _source_profile(record, input_path)
        supplied_type = _text(record.get("evidence_type"))
        supplied_confidence = _text(record.get("confidence"))
        if supplied_type in EVIDENCE_ORDER:
            evidence_type = supplied_type
        if supplied_confidence in CONFIDENCE_ORDER:
            confidence = supplied_confidence
        source_url = _text(record.get("source_url"))
        source_record_id = _text(record.get("rule_id"))
        source_excerpt = _text(record.get("source_excerpt"))
        source_backed = bool(source_url and source_record_id and source_excerpt)
        if not source_backed:
            evidence_type, confidence, rank = "inference", "low", 4
        rows.append(
            {
                "species": _text(record.get("species")),
                "field": field,
                "value": value,
                "source_kind": source_kind,
                "source_url": source_url,
                "source_record_id": source_record_id,
                "source_excerpt": source_excerpt,
                "evidence_type": evidence_type,
                "confidence": confidence,
                "source_rank": rank,
                "source_backed": source_backed,
            }
        )
    return pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)


def llm_extracted_evidence(bundle_dirs: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Revalidate evidence-locked LLM bundles against their frozen source packets."""

    evidence_parts: list[pd.DataFrame] = []
    audit_parts: list[pd.DataFrame] = []
    for bundle_dir in bundle_dirs:
        candidates_path = bundle_dir / "v2_llm_reported_candidates.csv"
        ingest_manifest_path = bundle_dir / "llm_ingest_manifest.json"
        packet_input_path = bundle_dir / "packet_input.json"
        packet_manifest_path = bundle_dir / "packet_manifest.json"
        prompt_path = bundle_dir / "prompt.txt"
        result_path = bundle_dir / "result.jsonl"
        for path, label in (
            (candidates_path, "validated LLM candidate table"),
            (prompt_path, "frozen LLM prompt"),
            (result_path, "raw provider result"),
        ):
            if not path.is_file():
                raise typer.BadParameter(f"{label} does not exist: {path}")

        ingest = _load_json_object(ingest_manifest_path, "LLM ingest manifest")
        packet = _load_json_object(packet_input_path, "LLM packet input")
        packet_manifest = _load_json_object(packet_manifest_path, "LLM packet manifest")
        if _text(ingest.get("validation_status")) != "success":
            raise typer.BadParameter(
                f"LLM extraction bundle is not marked validation_status=success: {bundle_dir}"
            )
        expected_hashes = {
            "accepted_csv_sha256": _sha256(candidates_path),
            "packet_input_sha256": _sha256(packet_input_path),
            "packet_manifest_sha256": _sha256(packet_manifest_path),
            "prompt_sha256": _sha256(prompt_path),
            "result_jsonl_sha256": _sha256(result_path),
        }
        for key, observed in expected_hashes.items():
            if not _is_sha256(ingest.get(key)) or _text(ingest.get(key)) != observed:
                raise typer.BadParameter(f"LLM extraction bundle {key} mismatch: {bundle_dir}")
        for key in ("packet_input_sha256", "prompt_sha256", "ontology_sha256"):
            if not _is_sha256(packet_manifest.get(key)):
                raise typer.BadParameter(f"LLM packet manifest has invalid {key}: {bundle_dir}")
            if _text(packet_manifest.get(key)) != _text(ingest.get(key)):
                raise typer.BadParameter(f"LLM ingest/packet {key} mismatch: {bundle_dir}")
        if _text(packet.get("ontology_sha256")) != _text(ingest.get("ontology_sha256")):
            raise typer.BadParameter(f"LLM packet ontology hash mismatch: {bundle_dir}")

        run_material = "|".join(
            [
                _text(ingest.get("model_provider")),
                _text(ingest.get("model_id")),
                _text(ingest.get("prompt_sha256")),
                _text(ingest.get("packet_input_sha256")),
                _text(ingest.get("result_jsonl_sha256")),
            ]
        )
        expected_run_id = f"llm_{hashlib.sha256(run_material.encode('utf-8')).hexdigest()[:20]}"
        if _text(ingest.get("extraction_run_id")) != expected_run_id:
            raise typer.BadParameter(f"LLM extraction_run_id mismatch: {bundle_dir}")
        if int(ingest.get("n_rejected_claims", -1)) != 0:
            raise typer.BadParameter(
                f"LLM extraction bundle contains rejected claims: {bundle_dir}"
            )

        candidates = pd.read_csv(candidates_path, dtype=str).fillna("")
        missing = set(LLM_EXTRACTED_COLUMNS).difference(candidates.columns)
        if missing:
            raise typer.BadParameter(
                f"{candidates_path} is missing validated LLM columns: {sorted(missing)}"
            )
        if _text(ingest.get("accepted_csv_filename")) != candidates_path.name:
            raise typer.BadParameter(f"LLM extraction candidate filename mismatch: {bundle_dir}")
        if int(ingest.get("accepted_csv_row_count", -1)) != len(candidates):
            raise typer.BadParameter(f"LLM extraction accepted row-count mismatch: {bundle_dir}")
        if int(ingest.get("n_validated_claims", -1)) != len(candidates):
            raise typer.BadParameter(f"LLM extraction row-count mismatch: {bundle_dir}")

        tasks = {
            _text(task.get("accepted_species")): task
            for task in packet.get("tasks", [])
            if isinstance(task, dict)
        }
        sources = {
            _text(source.get("source_chunk_id")): source
            for source in packet.get("sources", [])
            if isinstance(source, dict)
        }
        for row_number, row in enumerate(candidates.to_dict("records"), start=2):
            species = _text(row.get("accepted_species"))
            trait = _text(row.get("trait_name"))
            value = _text(row.get("provisional_candidate_value"))
            chunk_id = _text(row.get("source_chunk_id"))
            quote = _text(row.get("evidence_quote"))
            task = tasks.get(species)
            source = sources.get(chunk_id)
            try:
                if task is None or source is None:
                    raise ValueError("species or source chunk is absent from frozen packet")
                if trait not in task.get("target_traits", []):
                    raise ValueError("trait was not requested in frozen packet")
                if value not in (task.get("allowed_values") or {}).get(trait, []):
                    raise ValueError("value is outside the frozen ontology")
                if _text(source.get("accepted_species")) != species:
                    raise ValueError("source chunk species mismatch")
                if _text(source.get("evidence_scope")) != "species_direct":
                    raise ValueError("source chunk is not species-direct")
                if not quote or quote not in _text(source.get("chunk_text")):
                    raise ValueError("evidence quote is absent from frozen source chunk")
                if _text(row.get("source_url")) != _text(source.get("source_url")):
                    raise ValueError("source URL differs from frozen source chunk")
                if _text(row.get("source_text_sha256")) != _text(source.get("source_text_sha256")):
                    raise ValueError("source text hash differs from frozen source chunk")
                for column in LLM_EXTRACTED_PROVENANCE_COLUMNS:
                    if not _text(row.get(column)):
                        raise ValueError(f"blank LLM provenance column: {column}")
                for column in (
                    "prompt_sha256",
                    "packet_input_sha256",
                    "ontology_sha256",
                    "source_text_sha256",
                ):
                    if not _is_sha256(row.get(column)):
                        raise ValueError(f"invalid {column}")
                for column in (
                    "model_provider",
                    "model_id",
                    "prompt_version",
                    "prompt_sha256",
                    "packet_input_sha256",
                    "ontology_sha256",
                    "extraction_run_id",
                ):
                    if _text(row.get(column)) != _text(ingest.get(column)):
                        raise ValueError(f"row/ingest manifest {column} mismatch")
            except (TypeError, ValueError) as exc:
                raise typer.BadParameter(
                    f"invalid evidence-locked LLM row {candidates_path}:{row_number}: {exc}"
                ) from exc

        evidence, audit = candidate_evidence(
            candidates,
            candidates_path,
            allow_validated_llm_extraction=True,
        )
        evidence_parts.append(evidence)
        audit_parts.append(audit)

    evidence = (
        pd.concat(evidence_parts, ignore_index=True)
        if evidence_parts
        else pd.DataFrame(columns=EVIDENCE_COLUMNS)
    )
    audit = (
        pd.concat(audit_parts, ignore_index=True)
        if audit_parts
        else pd.DataFrame(columns=AUDIT_COLUMNS)
    )
    return evidence, audit


def llm_inference_evidence(bundle_dirs: list[Path], master: pd.DataFrame) -> pd.DataFrame:
    """Revalidate model-only rows against job history and the current taxonomy.

    Raw provider responses are intentionally rejected here. The upstream campaign
    must first bind every row to a known job and prompt hash, record the exact model,
    and publish a hash-locked validation bundle.
    """
    rows: list[dict[str, object]] = []
    master_by_species = {
        _text(record.get("accepted_species")): record for record in master.to_dict("records")
    }
    for bundle_dir in bundle_dirs:
        path = bundle_dir / "validated_trait_results.csv"
        history_path = bundle_dir / "job_history.csv"
        manifest_path = bundle_dir / "validated_llm_manifest.json"
        for required_path, label in (
            (path, "validated LLM result table"),
            (history_path, "validated LLM job history"),
        ):
            if not required_path.is_file():
                raise typer.BadParameter(f"{label} does not exist: {required_path}")
        manifest = _load_json_object(manifest_path, "validated LLM manifest")
        if _text(manifest.get("contract_version")) != "validated_llm_traits_v1":
            raise typer.BadParameter(f"unsupported validated LLM bundle: {bundle_dir}")
        if _text(manifest.get("validation_status")) != "success":
            raise typer.BadParameter(f"validated LLM bundle is not successful: {bundle_dir}")
        for key, file_path in (
            ("validated_trait_results_sha256", path),
            ("job_history_sha256", history_path),
        ):
            if not _is_sha256(manifest.get(key)) or _text(manifest.get(key)) != _sha256(file_path):
                raise typer.BadParameter(f"validated LLM bundle {key} mismatch: {bundle_dir}")

        frame = pd.read_csv(path, dtype=str).fillna("")
        required = set(VALIDATED_RESULT_COLUMNS).union(LLM_PROVENANCE_COLUMNS)
        missing = required.difference(frame.columns)
        if missing:
            raise typer.BadParameter(
                f"{path} is not a validated LLM inference table; missing {sorted(missing)}"
            )
        if frame["job_id"].duplicated().any():
            raise typer.BadParameter(f"validated LLM table contains duplicate job_id: {path}")
        history = pd.read_csv(history_path, dtype=str).fillna("")
        history_missing = set(JOB_HISTORY_COLUMNS).difference(history.columns)
        if history_missing:
            raise typer.BadParameter(
                f"{history_path} is missing job-history columns: {sorted(history_missing)}"
            )
        if history["job_id"].duplicated().any():
            raise typer.BadParameter(f"job history contains duplicate job_id: {history_path}")
        history_by_job = {
            _text(record.get("job_id")): record for record in history.to_dict("records")
        }
        if int(manifest.get("n_validated_results", -1)) != len(frame):
            raise typer.BadParameter(f"validated LLM result count mismatch: {bundle_dir}")

        for row_number, payload in enumerate(frame.to_dict("records"), start=2):
            try:
                parsed = {column: _text(payload.get(column)) for column in OUTPUT_COLUMNS}
                _validate_row(parsed, parsed["species"])
                species = parsed["species"]
                source = master_by_species.get(species)
                if source is None:
                    raise ValueError("species is absent from current master")
                for column in LLM_PROVENANCE_COLUMNS:
                    if not _text(payload.get(column)):
                        raise ValueError(f"blank validated provenance field: {column}")
                for column in ("job_sha256", "prompt_sha256", "result_sha256"):
                    value_hash = _text(payload.get(column))
                    if not _is_sha256(value_hash):
                        raise ValueError(f"invalid {column}")
                target_fields = [
                    field for field in _text(payload.get("target_fields")).split("|") if field
                ]
                if not target_fields or not set(target_fields).issubset(TRAIT_FIELDS):
                    raise ValueError("invalid target_fields")
                for field in TRAIT_FIELDS:
                    if field not in target_fields and parsed[field] != "unknown":
                        raise ValueError(f"non-target field {field} must be unknown")
                expected_job = _job_record(
                    {
                        "species": species,
                        "genus": source.get("genus"),
                        "family": source.get("family"),
                        "target_fields": "|".join(target_fields),
                    },
                    0,
                    1,
                )
                for column in (
                    "job_id",
                    "job_sha256",
                    "prompt_version",
                    "prompt_sha256",
                ):
                    if _text(payload.get(column)) != _text(expected_job.get(column)):
                        raise ValueError(f"current job {column} mismatch")
                if _text(payload.get("result_sha256")) != canonical_result_sha(payload):
                    raise ValueError("canonical result_sha256 mismatch")
                historical = history_by_job.get(_text(payload.get("job_id")))
                if historical is None or _text(historical.get("status")) != "completed":
                    raise ValueError("job is absent or not completed in job history")
                for column in (
                    "job_id",
                    "job_sha256",
                    "species",
                    "target_fields",
                    "prompt_version",
                    "prompt_sha256",
                    "result_sha256",
                ):
                    row_value = species if column == "species" else _text(payload.get(column))
                    if row_value != _text(historical.get(column)):
                        raise ValueError(f"result/job-history {column} mismatch")
                note = _text(parsed.get("pollination_notes"))
                for field in FIELDS:
                    value = _text(parsed.get(field))
                    if not value or value == "unknown":
                        continue
                    if field == "self_incompatibility":
                        value = {"SI": "likely_SI", "SC": "likely_SC"}.get(value, value)
                    rows.append(
                        {
                            "species": parsed["species"],
                            "field": field,
                            "value": value,
                            "source_kind": "llm_only",
                            "source_url": "",
                            "source_record_id": _text(payload.get("job_id")),
                            "source_excerpt": (
                                note if note and note != "unknown" else "LLM-only inference"
                            ),
                            "evidence_type": "inference",
                            "confidence": "low",
                            "source_rank": 4,
                            "source_backed": False,
                        }
                    )
            except Exception as exc:
                raise typer.BadParameter(
                    f"invalid validated LLM row {path}:{row_number}: {exc}"
                ) from exc
    return pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)


def _unique(values: pd.Series) -> list[str]:
    return list(dict.fromkeys(_text(value) for value in values if _text(value)))


def _collapse_value(field: str, values: list[str]) -> str:
    if not values:
        return "unknown"
    if field in {"flower_color", "flower_shape"}:
        return ", ".join(values)
    if field == "pollination_guild":
        distinct = set(values)
        if distinct == {"bees", "bumblebees"}:
            return "bumblebees"
        return values[0] if len(values) == 1 else "mixed"
    if field == "mating_system":
        return values[0] if len(values) == 1 else "mixed_mating"
    explicit = [value for value in values if value in {"SI", "SC"}]
    if len(set(explicit)) > 1:
        return "unknown"
    return explicit[0] if explicit else values[0]


def _short_notes(excerpts: list[str], limit: int = 360) -> str:
    text = " | ".join(excerpts[:3])
    if not text:
        return "unknown"
    return text if len(text) <= limit else text[: limit - 1].rstrip() + "…"


def collapse_all_species(master_species: list[str], evidence: pd.DataFrame) -> pd.DataFrame:
    """Choose the strongest available evidence independently for each requested field."""
    evidence = evidence.copy()
    evidence["source_rank"] = pd.to_numeric(evidence["source_rank"], errors="coerce").fillna(9)
    if evidence.empty:
        selected = evidence
    else:
        best_rank = evidence.groupby(["species", "field"], sort=False)["source_rank"].transform(
            "min"
        )
        selected = evidence.loc[evidence["source_rank"].eq(best_rank)].copy()

    value_by_key: dict[tuple[str, str], str] = {}
    for (species, field), group in selected.groupby(["species", "field"], sort=False):
        value_by_key[(_text(species), _text(field))] = _collapse_value(
            _text(field), _unique(group["value"])
        )

    metadata_by_species: dict[str, tuple[str, str, str]] = {}
    for species, group in selected.groupby("species", sort=False):
        evidence_types = set(_unique(group["evidence_type"]))
        ordered_types = [item for item in EVIDENCE_ORDER if item in evidence_types]
        confidences = _unique(group["confidence"])
        confidence = (
            min(confidences, key=lambda item: CONFIDENCE_ORDER.get(item, -1))
            if confidences
            else "low"
        )
        reproductive = group.loc[group["field"].isin(REPRODUCTIVE_FIELDS)]
        notes = _short_notes(_unique(reproductive["source_excerpt"]))
        metadata_by_species[_text(species)] = (
            "|".join(ordered_types) if ordered_types else "inference",
            confidence,
            notes,
        )

    output: list[dict[str, str]] = []
    for species in master_species:
        values = {field: value_by_key.get((species, field), "unknown") for field in FIELDS}
        evidence_type, confidence, notes = metadata_by_species.get(
            species, ("inference", "low", "unknown")
        )
        output.append(
            {
                "species": species,
                "flower_color": values["flower_color"],
                "flower_shape": values["flower_shape"],
                "pollination_guild": values["pollination_guild"],
                "pollination_notes": notes,
                "mating_system": values["mating_system"],
                "self_incompatibility": values["self_incompatibility"],
                "evidence_type": evidence_type,
                "confidence": confidence,
            }
        )
    return validate_result_table(
        pd.DataFrame(output, columns=OUTPUT_COLUMNS), expected_species=master_species
    )


def build_export(
    *,
    master_csv: Path,
    output_dir: Path,
    applicability_csv: Path | None = None,
    candidate_csvs: list[Path] | None = None,
    web_evidence_csvs: list[Path] | None = None,
    llm_inference_bundle_dirs: list[Path] | None = None,
    llm_extracted_bundle_dirs: list[Path] | None = None,
) -> dict[str, Any]:
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    if "accepted_species" not in master.columns:
        raise typer.BadParameter("master CSV must contain accepted_species")
    if "genus" not in master.columns:
        master["genus"] = master["accepted_species"].map(
            lambda value: _text(value).split()[0] if _text(value) else ""
        )
    if "family" not in master.columns:
        master["family"] = ""
    master_species = master["accepted_species"].map(_text).tolist()
    if any(not species for species in master_species):
        raise typer.BadParameter("master CSV contains a blank accepted_species")
    duplicates = pd.Series(master_species)[pd.Series(master_species).duplicated()].unique().tolist()
    if duplicates:
        raise typer.BadParameter(f"master CSV contains duplicate species: {duplicates[:10]}")

    eligible_species = set(master_species)
    if applicability_csv is not None:
        applicability = pd.read_csv(applicability_csv, dtype=str).fillna("")
        required = {"accepted_species", "angiosperm_analysis_eligible"}
        missing = required.difference(applicability.columns)
        if missing:
            raise typer.BadParameter(
                f"{applicability_csv} is missing applicability columns: {sorted(missing)}"
            )
        if applicability["accepted_species"].duplicated().any():
            raise typer.BadParameter("applicability CSV contains duplicate accepted_species")
        observed_species = set(applicability["accepted_species"].map(_text))
        if observed_species != set(master_species):
            raise typer.BadParameter("applicability CSV species do not exactly match the master")
        eligible = applicability["angiosperm_analysis_eligible"].map(
            lambda value: _text(value).casefold() in {"true", "1", "yes"}
        )
        eligible_species = set(applicability.loc[eligible, "accepted_species"].map(_text))

    evidence_parts: list[pd.DataFrame] = []
    audit_parts: list[pd.DataFrame] = []
    for path in candidate_csvs or []:
        frame = pd.read_csv(path, dtype=str).fillna("")
        evidence, audit = candidate_evidence(frame, path)
        evidence_parts.append(evidence)
        audit_parts.append(audit)
    for path in web_evidence_csvs or []:
        evidence_parts.append(web_evidence(pd.read_csv(path, dtype=str).fillna(""), path))
    llm_extracted, llm_extracted_audit = llm_extracted_evidence(llm_extracted_bundle_dirs or [])
    evidence_parts.append(llm_extracted)
    audit_parts.append(llm_extracted_audit)
    evidence_parts.append(llm_inference_evidence(llm_inference_bundle_dirs or [], master))

    evidence = (
        pd.concat(evidence_parts, ignore_index=True)
        if evidence_parts
        else pd.DataFrame(columns=EVIDENCE_COLUMNS)
    )
    evidence = evidence.loc[evidence["species"].isin(eligible_species)].drop_duplicates()
    audit = (
        pd.concat(audit_parts, ignore_index=True)
        if audit_parts
        else pd.DataFrame(columns=AUDIT_COLUMNS)
    )
    result = collapse_all_species(master_species, evidence)

    output_dir.mkdir(parents=True, exist_ok=True)
    output_csv = output_dir / "all_species_traits.csv"
    evidence_csv = output_dir / "all_species_trait_evidence.csv.gz"
    audit_csv = output_dir / "all_species_trait_unmapped.csv"
    result.to_csv(output_csv, index=False)
    evidence.to_csv(evidence_csv, index=False, compression="gzip")
    audit.to_csv(audit_csv, index=False)

    coverage = {field: int(result[field].ne("unknown").sum()) for field in FIELDS}
    any_known = result[list(FIELDS)].ne("unknown").any(axis=1)
    all_known = result[list(FIELDS)].ne("unknown").all(axis=1)
    manifest = {
        "version": "1.0",
        "n_master_species": len(master_species),
        "n_output_rows": len(result),
        "n_trait_applicable_species": len(eligible_species),
        "n_trait_not_applicable_species": len(master_species) - len(eligible_species),
        "output_columns": OUTPUT_COLUMNS,
        "n_evidence_rows": len(evidence),
        "n_source_backed_evidence_rows": int(evidence["source_backed"].astype(bool).sum()),
        "n_llm_only_evidence_rows": int(evidence["source_kind"].eq("llm_only").sum()),
        "n_unmapped_or_invalid_rows": len(audit),
        "coverage": coverage,
        "n_species_any_known": int(any_known.sum()),
        "n_species_all_five_known": int(all_known.sum()),
        "n_species_all_unknown": int((~any_known).sum()),
        "output_sha256": _sha256(output_csv),
        "candidate_csvs": [str(path) for path in candidate_csvs or []],
        "applicability_csv": str(applicability_csv) if applicability_csv is not None else None,
        "web_evidence_csvs": [str(path) for path in web_evidence_csvs or []],
        "llm_inference_bundle_dirs": [str(path) for path in llm_inference_bundle_dirs or []],
        "llm_extracted_bundle_dirs": [str(path) for path in llm_extracted_bundle_dirs or []],
        "policy": (
            "Source-backed candidates, including evidence-locked LLM extraction, outrank "
            "public-web rule matches, which outrank validated LLM-only inference. Raw LLM "
            "responses are never accepted. LLM-only rows are always inference/low. Missing "
            "fields remain unknown; every master species remains exactly once."
        ),
    }
    (output_dir / "all_species_traits_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("run")
def run_command(
    master_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    applicability_csv: Path | None = typer.Option(None, exists=True),
    candidate_csv: list[Path] | None = typer.Option(None, "--candidate-csv", exists=True),
    web_evidence_csv: list[Path] | None = typer.Option(None, "--web-evidence-csv", exists=True),
    llm_inference_bundle: list[Path] | None = typer.Option(
        None,
        "--llm-inference-bundle",
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    llm_extracted_bundle: list[Path] | None = typer.Option(
        None,
        "--llm-extracted-bundle",
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
) -> None:
    report = build_export(
        master_csv=master_csv,
        output_dir=output_dir,
        applicability_csv=applicability_csv,
        candidate_csvs=candidate_csv,
        web_evidence_csvs=web_evidence_csv,
        llm_inference_bundle_dirs=llm_inference_bundle,
        llm_extracted_bundle_dirs=llm_extracted_bundle,
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
