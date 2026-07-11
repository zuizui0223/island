"""Provider-neutral, evidence-locked LLM extraction for cached plant trait text.

This module does not call an LLM. It freezes source text into deterministic
packets, then validates externally produced JSONL claims against the frozen
text and the v2 trait ontology. A claim is accepted only when its evidence
quote is present in the cited source chunk and its value is allowed by the
ontology.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

PROMPT_VERSION = "v2_llm_evidence_v1"
DEFAULT_PROMPT_PATH = Path("config/v2_llm_evidence_prompt.txt")
DEFAULT_ONTOLOGY_PATH = Path("config/trait_ontology.yml")

SOURCE_REQUIRED_COLUMNS = {
    "accepted_species",
    "source_text",
    "source_url",
    "source_citation",
    "source_type",
    "evidence_scope",
}
TARGET_TRAITS = [
    "flower_primary_color",
    "floral_form",
    "floral_symmetry",
    "autonomous_selfing_capacity",
    "mating_system",
    "self_incompatibility",
    "cleistogamy",
    "pollen_vector_mode",
    "pollination_functional_guild",
]
CLAIM_REQUIRED_FIELDS = {
    "accepted_species",
    "trait_name",
    "standardized_value",
    "source_chunk_id",
    "evidence_quote",
    "confidence",
}
CONFIDENCE_VALUES = {"high", "medium", "low"}

VALIDATED_COLUMNS = [
    "accepted_species",
    "trait_layer",
    "trait_name",
    "provisional_candidate_value",
    "matched_term",
    "candidate_class",
    "source",
    "source_type",
    "source_url",
    "raw_description",
    "evidence_scope",
    "wild_or_cultivated",
    "review_status",
    "evidence_quote",
    "source_chunk_id",
    "source_text_sha256",
    "extraction_method",
    "model_provider",
    "model_id",
    "prompt_version",
    "prompt_sha256",
    "packet_input_sha256",
    "ontology_sha256",
    "extraction_run_id",
    "confidence",
]


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_text(text: str) -> str:
    return _sha256_bytes(text.encode("utf-8"))


def _file_sha256(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _json_dump(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def load_ontology(path: Path = DEFAULT_ONTOLOGY_PATH) -> tuple[dict[str, Any], str]:
    if not path.exists():
        raise typer.BadParameter(f"ontology file does not exist: {path}")
    raw = path.read_text(encoding="utf-8")
    ontology = yaml.safe_load(raw) or {}
    traits = ontology.get("traits")
    if not isinstance(traits, dict) or not traits:
        raise typer.BadParameter("trait ontology must contain a non-empty traits mapping")
    return ontology, _sha256_text(raw)


def allowed_values(ontology: dict[str, Any], trait_name: str) -> list[str]:
    trait = (ontology.get("traits") or {}).get(trait_name) or {}
    values = trait.get("allowed_values") or []
    if not isinstance(values, list):
        return []
    return [_text(value) for value in values if _text(value)]


def trait_layer(ontology: dict[str, Any], trait_name: str) -> str:
    domain = _text(((ontology.get("traits") or {}).get(trait_name) or {}).get("domain"))
    if domain in {"floral_signal", "floral_architecture"}:
        return "M0_floral_phenotype"
    if domain in {"reproductive_assurance", "pollen_vector"}:
        return "M1_reproductive_assurance_pollen_vector"
    if domain == "pollination":
        return "M2_pollinator_functional_evidence"
    return "M0_other_reported_trait"


def _candidate_value_column(table: pd.DataFrame) -> str | None:
    for column in (
        "provisional_candidate_value",
        "candidate_value",
        "standardized_value",
        "proxy_value",
    ):
        if column in table.columns:
            return column
    return None


def existing_coverage(paths: Iterable[Path]) -> set[tuple[str, str]]:
    covered: set[tuple[str, str]] = set()
    for path in paths:
        table = pd.read_csv(path, dtype=str).fillna("")
        if "accepted_species" not in table.columns or "trait_name" not in table.columns:
            continue
        value_column = _candidate_value_column(table)
        for row in table.to_dict("records"):
            species = _text(row.get("accepted_species"))
            trait = _text(row.get("trait_name"))
            value = _text(row.get(value_column)) if value_column else "reported"
            if species and trait and value and value != "unresolved":
                covered.add((species, trait))
    return covered


def _species_order(
    sources: pd.DataFrame,
    species_csv: Path | None,
    species_column: str,
) -> list[str]:
    if species_csv is None:
        values = sources["accepted_species"]
    else:
        table = pd.read_csv(species_csv, dtype=str).fillna("")
        if species_column not in table.columns:
            raise typer.BadParameter(
                f"species column {species_column!r} not found in {species_csv}"
            )
        values = table[species_column]
    ordered: list[str] = []
    seen: set[str] = set()
    for value in values:
        species = _text(value)
        if species and species not in seen:
            seen.add(species)
            ordered.append(species)
    return ordered


def _chunk_text(text: str, chunk_chars: int, overlap_chars: int) -> list[str]:
    if chunk_chars < 1000:
        raise typer.BadParameter("chunk_chars must be at least 1000")
    if overlap_chars < 0 or overlap_chars >= chunk_chars:
        raise typer.BadParameter("overlap_chars must be >= 0 and < chunk_chars")
    normalized = _text(text)
    if not normalized:
        return []
    if len(normalized) <= chunk_chars:
        return [normalized]
    chunks: list[str] = []
    start = 0
    while start < len(normalized):
        end = min(len(normalized), start + chunk_chars)
        if end < len(normalized):
            boundary = max(
                normalized.rfind(". ", start + chunk_chars // 2, end),
                normalized.rfind("; ", start + chunk_chars // 2, end),
            )
            if boundary > start:
                end = boundary + 1
        chunks.append(normalized[start:end].strip())
        if end >= len(normalized):
            break
        start = max(start + 1, end - overlap_chars)
    return [chunk for chunk in chunks if chunk]


def _source_chunks(
    sources: pd.DataFrame,
    selected_species: set[str],
    chunk_chars: int,
    overlap_chars: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for record in sources.fillna("").to_dict("records"):
        species = _text(record.get("accepted_species"))
        text = _text(record.get("source_text"))
        scope = _text(record.get("evidence_scope"))
        if species not in selected_species or not text or scope != "species_direct":
            continue
        source_sha = _sha256_text(text)
        chunks = _chunk_text(text, chunk_chars, overlap_chars)
        for index, chunk in enumerate(chunks, start=1):
            identity = "|".join(
                [species, _text(record.get("source_url")), source_sha, str(index), chunk]
            )
            chunk_id = f"src_{_sha256_text(identity)[:20]}"
            if chunk_id in seen:
                continue
            seen.add(chunk_id)
            rows.append(
                {
                    "source_chunk_id": chunk_id,
                    "accepted_species": species,
                    "source_url": _text(record.get("source_url")),
                    "source_citation": _text(record.get("source_citation")),
                    "source_type": _text(record.get("source_type")) or "cached_source_text",
                    "evidence_scope": scope,
                    "source_text_sha256": source_sha,
                    "chunk_index": index,
                    "chunk_text": chunk,
                }
            )
    return rows


def prepare_packets(
    *,
    source_csv: Path,
    output_dir: Path,
    prompt_path: Path = DEFAULT_PROMPT_PATH,
    ontology_path: Path = DEFAULT_ONTOLOGY_PATH,
    species_csv: Path | None = None,
    species_column: str = "accepted_species",
    candidate_csvs: Iterable[Path] = (),
    batch_size: int = 10,
    chunk_chars: int = 12000,
    overlap_chars: int = 500,
    target_traits: Iterable[str] = TARGET_TRAITS,
) -> dict[str, Any]:
    if batch_size < 1 or batch_size > 100:
        raise typer.BadParameter("batch_size must be between 1 and 100")
    if not source_csv.exists():
        raise typer.BadParameter(f"source CSV does not exist: {source_csv}")
    if not prompt_path.exists():
        raise typer.BadParameter(f"prompt file does not exist: {prompt_path}")

    sources = pd.read_csv(source_csv, dtype=str).fillna("")
    missing = SOURCE_REQUIRED_COLUMNS.difference(sources.columns)
    if missing:
        raise typer.BadParameter(f"source CSV is missing columns: {sorted(missing)}")

    ontology, ontology_sha = load_ontology(ontology_path)
    requested_traits = [_text(value) for value in target_traits if _text(value)]
    unknown_traits = [trait for trait in requested_traits if trait not in (ontology.get("traits") or {})]
    if unknown_traits:
        raise typer.BadParameter(f"target traits not present in ontology: {unknown_traits}")

    order = _species_order(sources, species_csv, species_column)
    source_species = set(sources.loc[sources["evidence_scope"].eq("species_direct"), "accepted_species"])
    coverage = existing_coverage(candidate_csvs)

    task_specs: list[dict[str, Any]] = []
    for species in order:
        if species not in source_species:
            continue
        gaps = [trait for trait in requested_traits if (species, trait) not in coverage]
        if not gaps:
            continue
        task_specs.append(
            {
                "accepted_species": species,
                "target_traits": gaps,
                "allowed_values": {
                    trait: [value for value in allowed_values(ontology, trait) if value != "unresolved"]
                    for trait in gaps
                },
            }
        )

    selected_species = {task["accepted_species"] for task in task_specs}
    chunks = _source_chunks(sources, selected_species, chunk_chars, overlap_chars)
    chunks_by_species: dict[str, list[dict[str, Any]]] = {}
    for chunk in chunks:
        chunks_by_species.setdefault(chunk["accepted_species"], []).append(chunk)

    task_specs = [
        {
            **task,
            "source_chunk_ids": [
                chunk["source_chunk_id"]
                for chunk in chunks_by_species.get(task["accepted_species"], [])
            ],
        }
        for task in task_specs
        if chunks_by_species.get(task["accepted_species"])
    ]

    prompt_text = prompt_path.read_text(encoding="utf-8")
    prompt_sha = _sha256_text(prompt_text)
    output_dir.mkdir(parents=True, exist_ok=True)
    packet_rows: list[dict[str, Any]] = []

    for batch_number, start in enumerate(range(0, len(task_specs), batch_size), start=1):
        batch_tasks = task_specs[start : start + batch_size]
        batch_species = {task["accepted_species"] for task in batch_tasks}
        batch_chunks = [chunk for chunk in chunks if chunk["accepted_species"] in batch_species]
        packet_id = f"batch_{batch_number:05d}"
        batch_dir = output_dir / packet_id
        batch_dir.mkdir(parents=True, exist_ok=True)
        packet_input = {
            "contract_version": "1.0",
            "prompt_version": PROMPT_VERSION,
            "ontology_sha256": ontology_sha,
            "tasks": batch_tasks,
            "sources": batch_chunks,
        }
        input_path = batch_dir / "packet_input.json"
        _json_dump(input_path, packet_input)
        (batch_dir / "prompt.txt").write_text(prompt_text, encoding="utf-8")
        (batch_dir / "result_template.jsonl").write_text("", encoding="utf-8")
        input_sha = _file_sha256(input_path)
        manifest = {
            "contract_version": "1.0",
            "packet_id": packet_id,
            "prompt_version": PROMPT_VERSION,
            "prompt_sha256": prompt_sha,
            "packet_input_sha256": input_sha,
            "ontology_sha256": ontology_sha,
            "source_csv_sha256": _file_sha256(source_csv),
            "n_species": len(batch_tasks),
            "n_source_chunks": len(batch_chunks),
            "n_requested_trait_gaps": sum(len(task["target_traits"]) for task in batch_tasks),
        }
        _json_dump(batch_dir / "packet_manifest.json", manifest)
        packet_rows.append(
            {
                "packet_id": packet_id,
                "n_species": len(batch_tasks),
                "n_source_chunks": len(batch_chunks),
                "n_requested_trait_gaps": manifest["n_requested_trait_gaps"],
                "first_species": batch_tasks[0]["accepted_species"],
                "last_species": batch_tasks[-1]["accepted_species"],
                "packet_dir": str(batch_dir),
                "packet_input_sha256": input_sha,
            }
        )

    pd.DataFrame(packet_rows).to_csv(output_dir / "packet_index.csv", index=False)
    report = {
        "contract_version": "1.0",
        "method": "provider-neutral evidence-locked LLM extraction packets",
        "prompt_version": PROMPT_VERSION,
        "prompt_sha256": prompt_sha,
        "ontology_sha256": ontology_sha,
        "source_csv_sha256": _file_sha256(source_csv),
        "n_species_with_source_text": len(source_species),
        "n_species_with_unresolved_target_traits": len(task_specs),
        "n_packets": len(packet_rows),
        "n_source_chunks": len(chunks),
        "n_requested_trait_gaps": sum(len(task["target_traits"]) for task in task_specs),
        "policy": (
            "Packets contain frozen species-direct source text only. The LLM may emit a claim only "
            "when it can cite a source chunk and an evidence quote. Missing claims remain unresolved."
        ),
    }
    _json_dump(output_dir / "llm_packet_manifest.json", report)
    return report


def _read_jsonl(path: Path) -> list[tuple[int, dict[str, Any] | None, str]]:
    rows: list[tuple[int, dict[str, Any] | None, str]] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not raw.strip():
            continue
        try:
            value = json.loads(raw)
        except json.JSONDecodeError:
            rows.append((line_number, None, raw))
            continue
        rows.append((line_number, value if isinstance(value, dict) else None, raw))
    return rows


def _normalised_contains(haystack: str, needle: str) -> bool:
    normalized_haystack = _text(haystack)
    normalized_needle = _text(needle)
    return bool(normalized_needle) and normalized_needle in normalized_haystack


def validate_results(
    *,
    packet_dir: Path,
    result_jsonl: Path,
    output_dir: Path,
    model_id: str,
    model_provider: str = "user_supplied",
    strict: bool = True,
) -> dict[str, Any]:
    if not model_id.strip():
        raise typer.BadParameter("model_id is required")
    input_path = packet_dir / "packet_input.json"
    prompt_path = packet_dir / "prompt.txt"
    manifest_path = packet_dir / "packet_manifest.json"
    for path in (input_path, prompt_path, manifest_path, result_jsonl):
        if not path.exists():
            raise typer.BadParameter(f"required file does not exist: {path}")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if _file_sha256(input_path) != manifest.get("packet_input_sha256"):
        raise typer.BadParameter("packet_input.json hash does not match packet_manifest.json")
    if _file_sha256(prompt_path) != manifest.get("prompt_sha256"):
        raise typer.BadParameter("prompt.txt hash does not match packet_manifest.json")

    packet = json.loads(input_path.read_text(encoding="utf-8"))
    task_by_species = {
        _text(task.get("accepted_species")): task for task in packet.get("tasks", [])
    }
    source_by_id = {
        _text(source.get("source_chunk_id")): source for source in packet.get("sources", [])
    }
    result_sha = _file_sha256(result_jsonl)
    run_material = "|".join(
        [
            _text(model_provider),
            _text(model_id),
            _text(manifest.get("prompt_sha256")),
            _text(manifest.get("packet_input_sha256")),
            result_sha,
        ]
    )
    extraction_run_id = f"llm_{_sha256_text(run_material)[:20]}"

    accepted: list[dict[str, str]] = []
    rejected: list[dict[str, str]] = []
    seen: set[tuple[str, str, str, str, str]] = set()

    for line_number, claim, raw in _read_jsonl(result_jsonl):
        reason = ""
        if claim is None:
            reason = "invalid_json_object"
        else:
            missing = CLAIM_REQUIRED_FIELDS.difference(claim)
            if missing:
                reason = f"missing_fields:{','.join(sorted(missing))}"
        if reason:
            rejected.append({"line_number": str(line_number), "reason": reason, "raw_json": raw})
            continue

        assert claim is not None
        species = _text(claim.get("accepted_species"))
        trait = _text(claim.get("trait_name"))
        value = _text(claim.get("standardized_value"))
        source_chunk_id = _text(claim.get("source_chunk_id"))
        quote = _text(claim.get("evidence_quote"))
        confidence = _text(claim.get("confidence")).lower()

        task = task_by_species.get(species)
        source = source_by_id.get(source_chunk_id)
        if task is None:
            reason = "species_not_in_packet"
        elif trait not in task.get("target_traits", []):
            reason = "trait_not_requested_for_species"
        elif value not in (task.get("allowed_values") or {}).get(trait, []):
            reason = "value_outside_packet_ontology"
        elif confidence not in CONFIDENCE_VALUES:
            reason = "invalid_confidence"
        elif source is None:
            reason = "unknown_source_chunk_id"
        elif _text(source.get("accepted_species")) != species:
            reason = "source_species_mismatch"
        elif _text(source.get("evidence_scope")) != "species_direct":
            reason = "source_not_species_direct"
        elif not _normalised_contains(_text(source.get("chunk_text")), quote):
            reason = "evidence_quote_not_in_source_chunk"

        if reason:
            rejected.append({"line_number": str(line_number), "reason": reason, "raw_json": raw})
            continue

        key = (species, trait, value, source_chunk_id, quote)
        if key in seen:
            continue
        seen.add(key)
        accepted.append(
            {
                "accepted_species": species,
                "trait_layer": trait_layer_from_packet(packet, trait),
                "trait_name": trait,
                "provisional_candidate_value": value,
                "matched_term": quote,
                "candidate_class": "reported",
                "source": "llm_evidence_extraction",
                "source_type": _text(source.get("source_type")),
                "source_url": _text(source.get("source_url")),
                "raw_description": quote,
                "evidence_scope": "species_direct",
                "wild_or_cultivated": "unknown",
                "review_status": "unreviewed_llm_extraction",
                "evidence_quote": quote,
                "source_chunk_id": source_chunk_id,
                "source_text_sha256": _text(source.get("source_text_sha256")),
                "extraction_method": "provider_neutral_llm_from_frozen_source_text",
                "model_provider": _text(model_provider),
                "model_id": _text(model_id),
                "prompt_version": _text(packet.get("prompt_version")),
                "prompt_sha256": _text(manifest.get("prompt_sha256")),
                "packet_input_sha256": _text(manifest.get("packet_input_sha256")),
                "ontology_sha256": _text(packet.get("ontology_sha256")),
                "extraction_run_id": extraction_run_id,
                "confidence": confidence,
            }
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    accepted_frame = pd.DataFrame(accepted, columns=VALIDATED_COLUMNS)
    rejected_frame = pd.DataFrame(rejected, columns=["line_number", "reason", "raw_json"])
    accepted_frame.to_csv(output_dir / "v2_llm_reported_candidates.csv", index=False)
    rejected_frame.to_csv(output_dir / "llm_rejected_claims.csv", index=False)
    report = {
        "contract_version": "1.0",
        "extraction_run_id": extraction_run_id,
        "model_provider": _text(model_provider),
        "model_id": _text(model_id),
        "result_jsonl_sha256": result_sha,
        "prompt_version": _text(packet.get("prompt_version")),
        "prompt_sha256": _text(manifest.get("prompt_sha256")),
        "packet_input_sha256": _text(manifest.get("packet_input_sha256")),
        "ontology_sha256": _text(packet.get("ontology_sha256")),
        "n_validated_claims": len(accepted_frame),
        "n_rejected_claims": len(rejected_frame),
        "n_species_with_validated_claims": (
            int(accepted_frame["accepted_species"].nunique()) if len(accepted_frame) else 0
        ),
        "policy": (
            "Validated rows are source-backed unreviewed candidates, not curated truth. "
            "A missing LLM claim is unresolved and never biological absence."
        ),
    }
    _json_dump(output_dir / "llm_ingest_manifest.json", report)
    if strict and len(rejected_frame):
        raise typer.BadParameter(
            f"{len(rejected_frame)} LLM claim(s) failed validation; see llm_rejected_claims.csv"
        )
    return report


def trait_layer_from_packet(packet: dict[str, Any], trait_name: str) -> str:
    """Infer layer from the packet's ontology hash contract and known v2 trait names."""
    if trait_name in {"flower_primary_color", "floral_form", "floral_symmetry"}:
        return "M0_floral_phenotype"
    if trait_name in {
        "autonomous_selfing_capacity",
        "mating_system",
        "self_incompatibility",
        "cleistogamy",
        "pollen_vector_mode",
    }:
        return "M1_reproductive_assurance_pollen_vector"
    if trait_name == "pollination_functional_guild":
        return "M2_pollinator_functional_evidence"
    return "M0_other_reported_trait"


@app.command("prepare")
def prepare_command(
    source_csv: Path = typer.Option(..., exists=True, help="Cached source_texts.csv."),
    output_dir: Path = typer.Option(...),
    prompt_path: Path = typer.Option(DEFAULT_PROMPT_PATH, exists=True),
    ontology_path: Path = typer.Option(DEFAULT_ONTOLOGY_PATH, exists=True),
    species_csv: Path | None = typer.Option(None, exists=True),
    species_column: str = typer.Option("accepted_species"),
    candidate_csv: list[Path] | None = typer.Option(
        None,
        "--candidate-csv",
        exists=True,
        help="Existing candidate CSV; may be repeated to skip already-covered traits.",
    ),
    batch_size: int = typer.Option(10, min=1, max=100),
    chunk_chars: int = typer.Option(12000, min=1000),
    overlap_chars: int = typer.Option(500, min=0),
) -> None:
    report = prepare_packets(
        source_csv=source_csv,
        output_dir=output_dir,
        prompt_path=prompt_path,
        ontology_path=ontology_path,
        species_csv=species_csv,
        species_column=species_column,
        candidate_csvs=candidate_csv or [],
        batch_size=batch_size,
        chunk_chars=chunk_chars,
        overlap_chars=overlap_chars,
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


@app.command("ingest")
def ingest_command(
    packet_dir: Path = typer.Option(..., exists=True),
    result_jsonl: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    model_id: str = typer.Option(..., help="Exact model identifier recorded for this extraction run."),
    model_provider: str = typer.Option("user_supplied"),
    strict: bool = typer.Option(True, "--strict/--allow-rejected"),
) -> None:
    report = validate_results(
        packet_dir=packet_dir,
        result_jsonl=result_jsonl,
        output_dir=output_dir,
        model_id=model_id,
        model_provider=model_provider,
        strict=strict,
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
