"""Prepare and parse offline OpenAI Batch files without making API calls.

The adapter converts hash-locked ``jobs_shard_*.jsonl`` campaign jobs into
deterministic ``/v1/responses`` Batch JSONL.  It is deliberately an LLM-only
inference lane: no tools are attached, no API key is read, and no network call
is implemented.  Batch output parsing accepts only raw Responses API
``output_text`` content parts and reconciles every result to the frozen request
manifest.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import tempfile
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any

import typer

from island_v2.search_enabled_llm_campaign import (
    JOB_HASH_COLUMNS,
    canonical_job_sha,
)

MODEL_ID = "gpt-5.4-nano-2026-03-17"
RESPONSES_ENDPOINT = "/v1/responses"
CONTRACT_VERSION = "openai-batch-llm-v1"

MODEL_CONTEXT_TOKENS = 400_000
DEFAULT_MAX_OUTPUT_TOKENS = 1_024
MAX_ALLOWED_OUTPUT_TOKENS = 8_192

MAX_BATCH_REQUESTS = 50_000
MAX_BATCH_FILE_BYTES = 200_000_000

BATCH_INPUT_PRICE_PER_MILLION = Decimal("0.20")
BATCH_OUTPUT_PRICE_PER_MILLION = Decimal("1.25")
ONE_MILLION = Decimal(1_000_000)

REQUEST_MANIFEST_FILENAME = "request_manifest.json"
PROVIDER_RESULTS_FILENAME = "provider_results.jsonl"
PARSE_AUDIT_FILENAME = "parse_audit.csv"
PARSE_REPORT_FILENAME = "parse_report.json"

REQUIRED_JOB_FIELDS = (
    "job_id",
    "job_sha256",
    *JOB_HASH_COLUMNS,
    "prompt",
)

AUDIT_COLUMNS = (
    "source_file",
    "line_number",
    "custom_id",
    "job_id",
    "error_code",
    "message",
    "http_status",
    "response_model",
    "expected_model",
)

LLM_ONLY_INSTRUCTIONS = """This is an LLM-only inference lane with no web access or tools.
The source campaign prompt may ask you to search the web; do not do so and do not claim that
you searched, retrieved, cited, or verified a source. Use model knowledge or biological
inference only, keep unsupported fields unknown, set evidence_type to inference, and use
appropriately low confidence. Return the requested single CSV row in the structured `result`
field and no prose."""

RESULT_TEXT_FORMAT: dict[str, Any] = {
    "type": "json_schema",
    "name": "island_trait_inference_result",
    "strict": True,
    "schema": {
        "type": "object",
        "properties": {"result": {"type": "string"}},
        "required": ["result"],
        "additionalProperties": False,
    },
}

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Prepare and parse offline OpenAI Batch files; this command never calls an API.",
)


@dataclass(frozen=True)
class ValidatedJob:
    payload: dict[str, str]
    source_file: str
    source_line: int


@dataclass(frozen=True)
class PreparedRequest:
    job: ValidatedJob
    payload: dict[str, Any]
    line: bytes
    request_sha256: str
    prompt_utf8_bytes: int
    input_tokens_upper_bound: int


@dataclass(frozen=True)
class RawOutputLine:
    source_file: str
    line_number: int
    payload: dict[str, Any]
    custom_id: str


class OutputParseError(ValueError):
    def __init__(self, code: str, message: str) -> None:
        super().__init__(message)
        self.code = code


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_text(value: str) -> str:
    return _sha256_bytes(value.encode("utf-8"))


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def _normalise_text(value: str) -> str:
    return " ".join(value.strip().split())


def _require_string(
    payload: Mapping[str, Any],
    field: str,
    source: str,
    *,
    preserve_whitespace: bool = False,
) -> str:
    value = payload.get(field)
    if not isinstance(value, str):
        raise ValueError(f"{source}: {field} must be a string")
    if not value.strip():
        raise ValueError(f"{source}: {field} must not be blank")
    return value if preserve_whitespace else _normalise_text(value)


def _validate_hash(value: str, field: str, source: str) -> None:
    if re.fullmatch(r"[0-9a-f]{64}", value) is None:
        raise ValueError(f"{source}: {field} must be a lowercase SHA-256")


def _validate_job(payload: object, source_file: str, line_number: int) -> ValidatedJob:
    source = f"{source_file}:{line_number}"
    if not isinstance(payload, dict):
        raise ValueError(f"{source}: each JSONL line must be an object")
    values = {
        field: _require_string(
            payload,
            field,
            source,
            preserve_whitespace=field == "prompt",
        )
        for field in REQUIRED_JOB_FIELDS
    }
    _validate_hash(values["prompt_sha256"], "prompt_sha256", source)
    _validate_hash(values["prompt_policy_sha256"], "prompt_policy_sha256", source)
    _validate_hash(values["job_sha256"], "job_sha256", source)

    expected_prompt_sha = _sha256_text(values["prompt"])
    if values["prompt_sha256"] != expected_prompt_sha:
        raise ValueError(f"{source}: prompt_sha256 mismatch")
    expected_job_sha = canonical_job_sha(values)
    if values["job_sha256"] != expected_job_sha:
        raise ValueError(f"{source}: job_sha256 mismatch")
    expected_job_id = f"trait_{expected_job_sha[:20]}"
    if values["job_id"] != expected_job_id:
        raise ValueError(f"{source}: job_id does not match job_sha256")

    for field in ("shard_index", "shard_count"):
        if field in payload and not isinstance(payload[field], str):
            raise ValueError(f"{source}: {field} must be a string when present")
    if "shard_index" in payload and "shard_count" in payload:
        try:
            shard_index = int(payload["shard_index"])
            shard_count = int(payload["shard_count"])
        except ValueError as exc:
            raise ValueError(f"{source}: shard fields must be integers") from exc
        if shard_count < 1 or not 0 <= shard_index < shard_count:
            raise ValueError(f"{source}: invalid shard_index/shard_count")

    return ValidatedJob(values, source_file, line_number)


def load_campaign_jobs(paths: Sequence[Path]) -> tuple[ValidatedJob, ...]:
    """Load, independently re-hash, deduplicate, and deterministically order jobs."""

    if not paths:
        raise ValueError("at least one jobs_shard_*.jsonl input is required")
    resolved = [Path(path).resolve() for path in paths]
    if len(set(resolved)) != len(resolved):
        raise ValueError("duplicate jobs JSONL path")

    jobs: list[ValidatedJob] = []
    for path in sorted(resolved, key=lambda value: value.as_posix().casefold()):
        if re.fullmatch(r"jobs_shard_[0-9]+\.jsonl", path.name) is None:
            raise ValueError(f"input file must be named jobs_shard_*.jsonl: {path.name}")
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            found = False
            for line_number, raw in enumerate(handle, start=1):
                found = True
                text = raw.rstrip("\r\n")
                if not text.strip():
                    raise ValueError(f"{path.name}:{line_number}: blank JSONL line")
                try:
                    payload = json.loads(text)
                except json.JSONDecodeError as exc:
                    raise ValueError(f"{path.name}:{line_number}: invalid JSON") from exc
                jobs.append(_validate_job(payload, path.name, line_number))
            if not found:
                raise ValueError(f"{path.name}: jobs JSONL is empty")

    duplicate_checks = {
        "job_id": [job.payload["job_id"] for job in jobs],
        "job_sha256": [job.payload["job_sha256"] for job in jobs],
        "species": [job.payload["species"] for job in jobs],
    }
    for field, values in duplicate_checks.items():
        duplicates = sorted(value for value, count in Counter(values).items() if count > 1)
        if duplicates:
            raise ValueError(f"duplicate {field} values: {duplicates[:5]}")
    return tuple(sorted(jobs, key=lambda job: job.payload["job_id"]))


def _responses_body(job: ValidatedJob, max_output_tokens: int) -> dict[str, Any]:
    return {
        "input": job.payload["prompt"],
        "instructions": LLM_ONLY_INSTRUCTIONS,
        "max_output_tokens": max_output_tokens,
        "model": MODEL_ID,
        "reasoning": {"effort": "none"},
        "store": False,
        "text": {"format": RESULT_TEXT_FORMAT},
    }


def _prepare_request(job: ValidatedJob, max_output_tokens: int) -> PreparedRequest:
    body = _responses_body(job, max_output_tokens)
    body_bytes = _json_bytes(body)
    input_tokens_upper_bound = len(body_bytes)
    if input_tokens_upper_bound + max_output_tokens > MODEL_CONTEXT_TOKENS:
        raise ValueError(
            f"{job.payload['job_id']}: conservative input bound plus output cap "
            f"exceeds {MODEL_CONTEXT_TOKENS:,}-token context"
        )
    request = {
        "body": body,
        "custom_id": job.payload["job_id"],
        "method": "POST",
        "url": RESPONSES_ENDPOINT,
    }
    encoded = _json_bytes(request)
    return PreparedRequest(
        job=job,
        payload=request,
        line=encoded + b"\n",
        request_sha256=_sha256_bytes(encoded),
        prompt_utf8_bytes=len(job.payload["prompt"].encode("utf-8")),
        input_tokens_upper_bound=input_tokens_upper_bound,
    )


def _parse_max_usd(value: str | Decimal) -> Decimal:
    try:
        parsed = Decimal(str(value))
    except InvalidOperation as exc:
        raise ValueError("max_usd must be a finite nonnegative decimal") from exc
    if not parsed.is_finite() or parsed < 0:
        raise ValueError("max_usd must be a finite nonnegative decimal")
    return parsed


def _money_text(value: Decimal) -> str:
    return format(value, "f")


def _atomic_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{path.name}.",
            suffix=".tmp",
            dir=path.parent,
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    encoded = (
        json.dumps(
            payload,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        ).encode("utf-8")
        + b"\n"
    )
    _atomic_bytes(path, encoded)


def _atomic_jsonl(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    payload = b"".join(_json_bytes(row) + b"\n" for row in rows)
    _atomic_bytes(path, payload)


def _atomic_csv(
    path: Path,
    columns: Sequence[str],
    rows: Sequence[Mapping[str, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            prefix=f".{path.name}.",
            suffix=".tmp",
            dir=path.parent,
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            writer = csv.DictWriter(
                handle,
                fieldnames=columns,
                extrasaction="raise",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def _partition_requests(
    requests: Sequence[PreparedRequest],
    *,
    max_requests_per_file: int,
    max_file_bytes: int,
) -> tuple[tuple[PreparedRequest, ...], ...]:
    if not 1 <= max_requests_per_file <= MAX_BATCH_REQUESTS:
        raise ValueError(f"max_requests_per_file must be between 1 and {MAX_BATCH_REQUESTS:,}")
    if not 2 <= max_file_bytes <= MAX_BATCH_FILE_BYTES:
        raise ValueError(f"max_file_bytes must be between 2 and {MAX_BATCH_FILE_BYTES:,}")

    partitions: list[tuple[PreparedRequest, ...]] = []
    current: list[PreparedRequest] = []
    current_bytes = 0
    for request in requests:
        line_bytes = len(request.line)
        if line_bytes >= max_file_bytes:
            raise ValueError(
                f"{request.job.payload['job_id']}: one Batch request is not below "
                f"the {max_file_bytes:,}-byte file cap"
            )
        would_overflow = current_bytes + line_bytes >= max_file_bytes
        would_overcount = len(current) >= max_requests_per_file
        if current and (would_overflow or would_overcount):
            partitions.append(tuple(current))
            current = []
            current_bytes = 0
        current.append(request)
        current_bytes += line_bytes
    if current:
        partitions.append(tuple(current))
    return tuple(partitions)


def prepare_batch_requests(
    jobs_jsonl: Sequence[Path],
    output_dir: Path,
    *,
    max_usd: str | Decimal,
    max_output_tokens: int = DEFAULT_MAX_OUTPUT_TOKENS,
    max_requests_per_file: int = MAX_BATCH_REQUESTS,
    max_file_bytes: int = MAX_BATCH_FILE_BYTES,
) -> dict[str, Any]:
    """Validate jobs, enforce a worst-case budget, then atomically freeze Batch files."""

    if not 1 <= max_output_tokens <= MAX_ALLOWED_OUTPUT_TOKENS:
        raise ValueError(f"max_output_tokens must be between 1 and {MAX_ALLOWED_OUTPUT_TOKENS:,}")
    budget = _parse_max_usd(max_usd)
    jobs = load_campaign_jobs(jobs_jsonl)
    requests = tuple(_prepare_request(job, max_output_tokens) for job in jobs)

    prompt_utf8_bytes = sum(request.prompt_utf8_bytes for request in requests)
    input_tokens_upper_bound = sum(request.input_tokens_upper_bound for request in requests)
    output_tokens_upper_bound = len(requests) * max_output_tokens
    input_cost = Decimal(input_tokens_upper_bound) * BATCH_INPUT_PRICE_PER_MILLION / ONE_MILLION
    output_cost = Decimal(output_tokens_upper_bound) * BATCH_OUTPUT_PRICE_PER_MILLION / ONE_MILLION
    worst_case_cost = input_cost + output_cost
    if worst_case_cost > budget:
        raise ValueError(
            "worst-case Batch cost exceeds max_usd: "
            f"{_money_text(worst_case_cost)} > {_money_text(budget)}"
        )

    partitions = _partition_requests(
        requests,
        max_requests_per_file=max_requests_per_file,
        max_file_bytes=max_file_bytes,
    )
    destination = Path(output_dir)
    if destination.exists() and not destination.is_dir():
        raise ValueError("output_dir exists and is not a directory")
    managed_existing = []
    if destination.exists():
        managed_existing = [
            *destination.glob("batch_requests_*.jsonl"),
            *(path for path in [destination / REQUEST_MANIFEST_FILENAME] if path.exists()),
        ]
    if managed_existing:
        raise ValueError("output_dir already contains managed Batch preparation files")

    batch_records: list[dict[str, Any]] = []
    request_records: list[dict[str, Any]] = []
    for batch_index, partition in enumerate(partitions, start=1):
        filename = f"batch_requests_{batch_index:04d}.jsonl"
        payload = b"".join(request.line for request in partition)
        if len(partition) > MAX_BATCH_REQUESTS or len(payload) >= MAX_BATCH_FILE_BYTES:
            raise RuntimeError("internal Batch partition exceeded an official hard limit")
        batch_records.append(
            {
                "filename": filename,
                "request_count": len(partition),
                "bytes": len(payload),
                "sha256": _sha256_bytes(payload),
            }
        )
        for line_number, request in enumerate(partition, start=1):
            job = request.job.payload
            request_records.append(
                {
                    "batch_file": filename,
                    "batch_line_number": line_number,
                    "custom_id": job["job_id"],
                    "job_id": job["job_id"],
                    "job_sha256": job["job_sha256"],
                    "prompt_policy_sha256": job["prompt_policy_sha256"],
                    "prompt_sha256": job["prompt_sha256"],
                    "prompt_utf8_bytes": request.prompt_utf8_bytes,
                    "input_tokens_upper_bound": request.input_tokens_upper_bound,
                    "max_output_tokens": max_output_tokens,
                    "model": MODEL_ID,
                    "request_sha256": request.request_sha256,
                }
            )

    source_files = [
        {
            "filename": path.name,
            "bytes": path.stat().st_size,
            "sha256": _file_sha256(path),
        }
        for path in sorted(
            (Path(path).resolve() for path in jobs_jsonl),
            key=lambda value: value.as_posix().casefold(),
        )
    ]
    manifest: dict[str, Any] = {
        "contract_version": CONTRACT_VERSION,
        "endpoint": RESPONSES_ENDPOINT,
        "model": MODEL_ID,
        "lane": "llm_only_inference",
        "api_calls_made": False,
        "api_key_handled": False,
        "web_search_enabled": False,
        "structured_output": True,
        "reasoning_effort": "none",
        "store": False,
        "n_requests": len(requests),
        "n_batches": len(partitions),
        "source_job_files": source_files,
        "limits": {
            "official_max_requests_per_batch": MAX_BATCH_REQUESTS,
            "official_max_batch_file_bytes": MAX_BATCH_FILE_BYTES,
            "configured_max_requests_per_file": max_requests_per_file,
            "configured_max_file_bytes_exclusive": max_file_bytes,
            "model_context_tokens": MODEL_CONTEXT_TOKENS,
            "max_output_tokens_per_request": max_output_tokens,
        },
        "cost": {
            "currency": "USD",
            "pricing_basis": "OpenAI Batch API, fixed model snapshot",
            "pricing_source": ("https://developers.openai.com/api/docs/models/gpt-5.4-nano"),
            "input_usd_per_million_tokens": _money_text(BATCH_INPUT_PRICE_PER_MILLION),
            "output_usd_per_million_tokens": _money_text(BATCH_OUTPUT_PRICE_PER_MILLION),
            "input_token_bound_method": (
                "UTF-8 bytes of canonical Responses request body, counted as one "
                "input token per byte"
            ),
            "source_prompt_utf8_bytes": prompt_utf8_bytes,
            "input_tokens_upper_bound": input_tokens_upper_bound,
            "output_tokens_upper_bound": output_tokens_upper_bound,
            "input_worst_case_usd": _money_text(input_cost),
            "output_worst_case_usd": _money_text(output_cost),
            "worst_case_usd": _money_text(worst_case_cost),
            "max_usd": _money_text(budget),
        },
        "batches": batch_records,
        "requests": request_records,
    }

    destination.mkdir(parents=True, exist_ok=True)
    for batch, partition in zip(batch_records, partitions, strict=True):
        _atomic_bytes(
            destination / str(batch["filename"]),
            b"".join(request.line for request in partition),
        )
    _atomic_json(destination / REQUEST_MANIFEST_FILENAME, manifest)
    return manifest


def _audit_row(
    *,
    source_file: str,
    line_number: int,
    custom_id: str = "",
    job_id: str = "",
    error_code: str,
    message: str,
    http_status: str = "",
    response_model: str = "",
    expected_model: str = MODEL_ID,
) -> dict[str, str]:
    return {
        "source_file": source_file,
        "line_number": str(line_number),
        "custom_id": custom_id,
        "job_id": job_id,
        "error_code": error_code,
        "message": " ".join(str(message).split())[:1000],
        "http_status": http_status,
        "response_model": response_model,
        "expected_model": expected_model,
    }


def _compact_error(value: object) -> str:
    if isinstance(value, Mapping):
        code = value.get("code")
        message = value.get("message")
        return " ".join(str(part) for part in (code, message) if part)
    return str(value)


def _load_request_manifest(path: Path) -> tuple[dict[str, Any], tuple[dict[str, Any], ...]]:
    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError("request manifest is not valid JSON") from exc
    if not isinstance(payload, dict):
        raise ValueError("request manifest must be a JSON object")
    if payload.get("contract_version") != CONTRACT_VERSION:
        raise ValueError("request manifest contract_version mismatch")
    if payload.get("endpoint") != RESPONSES_ENDPOINT:
        raise ValueError("request manifest endpoint mismatch")
    if payload.get("model") != MODEL_ID:
        raise ValueError("request manifest model mismatch")
    if payload.get("api_calls_made") is not False:
        raise ValueError("request manifest does not declare offline preparation")
    if payload.get("web_search_enabled") is not False:
        raise ValueError("request manifest is not an LLM-only lane")
    requests = payload.get("requests")
    if not isinstance(requests, list) or not requests:
        raise ValueError("request manifest requests must be a nonempty list")

    validated: list[dict[str, Any]] = []
    custom_ids: list[str] = []
    for index, request in enumerate(requests, start=1):
        if not isinstance(request, dict):
            raise ValueError(f"request manifest entry {index} is not an object")
        source = f"request manifest entry {index}"
        custom_id = _require_string(request, "custom_id", source)
        job_id = _require_string(request, "job_id", source)
        prompt_sha = _require_string(request, "prompt_sha256", source)
        prompt_policy_sha = _require_string(request, "prompt_policy_sha256", source)
        job_sha = _require_string(request, "job_sha256", source)
        model = _require_string(request, "model", source)
        _validate_hash(prompt_sha, "prompt_sha256", source)
        _validate_hash(prompt_policy_sha, "prompt_policy_sha256", source)
        _validate_hash(job_sha, "job_sha256", source)
        if custom_id != job_id:
            raise ValueError(f"{source}: custom_id/job_id mismatch")
        if job_id != f"trait_{job_sha[:20]}":
            raise ValueError(f"{source}: job_id/job_sha256 mismatch")
        if model != MODEL_ID:
            raise ValueError(f"{source}: model mismatch")
        custom_ids.append(custom_id)
        validated.append(request)
    duplicates = [value for value, count in Counter(custom_ids).items() if count > 1]
    if duplicates:
        raise ValueError(f"request manifest contains duplicate custom_id: {duplicates[:5]}")
    return payload, tuple(validated)


def _read_output_lines(
    paths: Sequence[Path],
) -> tuple[list[RawOutputLine], list[dict[str, str]]]:
    if not paths:
        raise ValueError("at least one Batch output JSONL is required")
    resolved = [Path(path).resolve() for path in paths]
    if len(set(resolved)) != len(resolved):
        raise ValueError("duplicate Batch output JSONL path")
    rows: list[RawOutputLine] = []
    audit: list[dict[str, str]] = []
    for path in sorted(resolved, key=lambda value: value.as_posix().casefold()):
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            for line_number, raw in enumerate(handle, start=1):
                text = raw.rstrip("\r\n")
                if not text.strip():
                    audit.append(
                        _audit_row(
                            source_file=path.name,
                            line_number=line_number,
                            error_code="blank_line",
                            message="blank Batch output JSONL line",
                        )
                    )
                    continue
                try:
                    payload = json.loads(text)
                except json.JSONDecodeError:
                    audit.append(
                        _audit_row(
                            source_file=path.name,
                            line_number=line_number,
                            error_code="invalid_json",
                            message="invalid Batch output JSON",
                        )
                    )
                    continue
                if not isinstance(payload, dict):
                    audit.append(
                        _audit_row(
                            source_file=path.name,
                            line_number=line_number,
                            error_code="invalid_output_shape",
                            message="Batch output line must be an object",
                        )
                    )
                    continue
                custom_id = payload.get("custom_id")
                if not isinstance(custom_id, str) or not custom_id.strip():
                    audit.append(
                        _audit_row(
                            source_file=path.name,
                            line_number=line_number,
                            error_code="missing_custom_id",
                            message="Batch output line has no custom_id",
                        )
                    )
                    continue
                rows.append(RawOutputLine(path.name, line_number, payload, custom_id.strip()))
    return rows, audit


def _response_output_text(body: Mapping[str, Any]) -> str:
    output = body.get("output")
    if not isinstance(output, list):
        raise OutputParseError("invalid_response_shape", "response output must be a list")
    parts: list[str] = []
    for item in output:
        if not isinstance(item, Mapping):
            raise OutputParseError(
                "invalid_response_shape", "response output item must be an object"
            )
        item_type = item.get("type")
        if item_type == "reasoning":
            continue
        if item_type != "message":
            raise OutputParseError(
                "unexpected_output_item",
                f"response contained prohibited output item {item_type!r}",
            )
        if item.get("status") != "completed":
            raise OutputParseError(
                "incomplete_output_item",
                f"response message status is {item.get('status')!r}",
            )
        content = item.get("content")
        if not isinstance(content, list):
            raise OutputParseError(
                "invalid_response_shape", "response message content must be a list"
            )
        for content_part in content:
            if not isinstance(content_part, Mapping):
                raise OutputParseError(
                    "invalid_response_shape", "response content part must be an object"
                )
            part_type = content_part.get("type")
            if part_type == "refusal":
                raise OutputParseError("refusal", "model returned a refusal")
            if part_type != "output_text":
                raise OutputParseError(
                    "unexpected_content_type",
                    f"response message contained {part_type!r}, not output_text",
                )
            text = content_part.get("text")
            if not isinstance(text, str):
                raise OutputParseError(
                    "invalid_response_shape", "output_text text must be a string"
                )
            parts.append(text)
    if not parts:
        raise OutputParseError("missing_output_text", "response has no output_text")
    return "".join(parts)


def _structured_result(output_text: str) -> str:
    try:
        payload = json.loads(output_text)
    except json.JSONDecodeError as exc:
        raise OutputParseError(
            "invalid_structured_output", "output_text is not valid JSON"
        ) from exc
    if not isinstance(payload, dict) or set(payload) != {"result"}:
        raise OutputParseError(
            "invalid_structured_output",
            "output_text must contain exactly the structured result field",
        )
    result = payload.get("result")
    if not isinstance(result, str) or not result.strip():
        raise OutputParseError(
            "invalid_structured_output", "structured result must be a nonblank string"
        )
    return result


def _parse_one_output(
    row: RawOutputLine,
    expected: Mapping[str, Any],
) -> dict[str, Any]:
    payload = row.payload
    top_error = payload.get("error")
    if top_error is not None:
        raise OutputParseError("batch_error", _compact_error(top_error))
    response = payload.get("response")
    if not isinstance(response, Mapping):
        raise OutputParseError("missing_response", "Batch output has no response object")
    status_code = response.get("status_code")
    if not isinstance(status_code, int) or isinstance(status_code, bool):
        raise OutputParseError("invalid_http_status", "response status_code is not an integer")
    if status_code != 200:
        raise OutputParseError("http_error", f"Responses endpoint returned HTTP {status_code}")
    body = response.get("body")
    if not isinstance(body, Mapping):
        raise OutputParseError("invalid_response_shape", "response body must be an object")
    body_error = body.get("error")
    if body_error is not None:
        raise OutputParseError("response_error", _compact_error(body_error))
    status = body.get("status")
    if status != "completed":
        raise OutputParseError("response_not_completed", f"Responses body status is {status!r}")
    response_model = body.get("model")
    if response_model != expected["model"]:
        raise OutputParseError(
            "model_mismatch",
            f"response model {response_model!r} does not match request manifest",
        )
    output_text = _response_output_text(body)
    result = _structured_result(output_text)
    return {
        "job_id": expected["job_id"],
        "prompt_sha256": expected["prompt_sha256"],
        "result": result,
    }


def parse_batch_outputs(
    request_manifest: Path,
    batch_output_jsonl: Sequence[Path],
    output_dir: Path,
    *,
    strict: bool = True,
) -> dict[str, Any]:
    """Reconcile raw Batch outputs and atomically freeze accepted rows and audit."""

    manifest, requests = _load_request_manifest(Path(request_manifest))
    expected_by_id = {str(request["custom_id"]): request for request in requests}
    raw_rows, audit = _read_output_lines(batch_output_jsonl)
    counts = Counter(row.custom_id for row in raw_rows)
    seen_expected_ids = {row.custom_id for row in raw_rows if row.custom_id in expected_by_id}

    accepted_by_id: dict[str, dict[str, Any]] = {}
    for row in raw_rows:
        expected = expected_by_id.get(row.custom_id)
        if expected is None:
            audit.append(
                _audit_row(
                    source_file=row.source_file,
                    line_number=row.line_number,
                    custom_id=row.custom_id,
                    error_code="unknown_custom_id",
                    message="custom_id is absent from the request manifest",
                )
            )
            continue
        if counts[row.custom_id] > 1:
            audit.append(
                _audit_row(
                    source_file=row.source_file,
                    line_number=row.line_number,
                    custom_id=row.custom_id,
                    job_id=str(expected["job_id"]),
                    error_code="duplicate_output",
                    message="custom_id occurs more than once across Batch outputs",
                )
            )
            continue
        try:
            accepted_by_id[row.custom_id] = _parse_one_output(row, expected)
        except OutputParseError as exc:
            response = row.payload.get("response")
            body = response.get("body") if isinstance(response, Mapping) else None
            status = response.get("status_code") if isinstance(response, Mapping) else ""
            response_model = body.get("model") if isinstance(body, Mapping) else ""
            audit.append(
                _audit_row(
                    source_file=row.source_file,
                    line_number=row.line_number,
                    custom_id=row.custom_id,
                    job_id=str(expected["job_id"]),
                    error_code=exc.code,
                    message=str(exc),
                    http_status=str(status) if status != "" else "",
                    response_model=(str(response_model) if response_model is not None else ""),
                    expected_model=str(expected["model"]),
                )
            )

    for custom_id, expected in expected_by_id.items():
        if custom_id not in seen_expected_ids:
            audit.append(
                _audit_row(
                    source_file="",
                    line_number=0,
                    custom_id=custom_id,
                    job_id=str(expected["job_id"]),
                    error_code="missing_output",
                    message="request manifest entry has no Batch output line",
                    expected_model=str(expected["model"]),
                )
            )

    accepted = [
        accepted_by_id[custom_id] for custom_id in expected_by_id if custom_id in accepted_by_id
    ]
    audit.sort(
        key=lambda row: (
            row["custom_id"],
            row["source_file"],
            int(row["line_number"]),
            row["error_code"],
        )
    )
    strict_rejected = bool(strict and audit)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    result_path = destination / PROVIDER_RESULTS_FILENAME
    audit_path = destination / PARSE_AUDIT_FILENAME
    report_path = destination / PARSE_REPORT_FILENAME

    if strict_rejected:
        result_path.unlink(missing_ok=True)
    else:
        _atomic_jsonl(result_path, accepted)
    _atomic_csv(audit_path, AUDIT_COLUMNS, audit)
    report: dict[str, Any] = {
        "contract_version": CONTRACT_VERSION,
        "request_manifest_sha256": _file_sha256(Path(request_manifest)),
        "model": manifest["model"],
        "strict": strict,
        "strict_rejected": strict_rejected,
        "n_expected_requests": len(requests),
        "n_output_lines_with_custom_id": len(raw_rows),
        "n_accepted": 0 if strict_rejected else len(accepted),
        "n_would_accept_without_strict": len(accepted),
        "n_audit_errors": len(audit),
        "audit_error_counts": dict(sorted(Counter(row["error_code"] for row in audit).items())),
        "files": {
            PARSE_AUDIT_FILENAME: {
                "rows": len(audit),
                "sha256": _file_sha256(audit_path),
            }
        },
    }
    if not strict_rejected:
        report["files"][PROVIDER_RESULTS_FILENAME] = {
            "rows": len(accepted),
            "sha256": _file_sha256(result_path),
        }
    _atomic_json(report_path, report)
    return report


@app.command("prepare")
def prepare_command(
    jobs_jsonl: list[Path] = typer.Option(
        ...,
        "--jobs-jsonl",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Repeat for each jobs_shard_*.jsonl input.",
    ),
    output_dir: Path = typer.Option(..., "--output-dir", file_okay=False),
    max_usd: str = typer.Option(
        ...,
        "--max-usd",
        help="Required fail-closed worst-case Batch spend ceiling in USD.",
    ),
    max_output_tokens: int = typer.Option(
        DEFAULT_MAX_OUTPUT_TOKENS,
        min=1,
        max=MAX_ALLOWED_OUTPUT_TOKENS,
    ),
    max_requests_per_file: int = typer.Option(
        MAX_BATCH_REQUESTS,
        min=1,
        max=MAX_BATCH_REQUESTS,
    ),
    max_file_bytes: int = typer.Option(
        MAX_BATCH_FILE_BYTES,
        min=2,
        max=MAX_BATCH_FILE_BYTES,
    ),
) -> None:
    """Prepare deterministic offline Batch JSONL and a costed request manifest."""

    try:
        manifest = prepare_batch_requests(
            jobs_jsonl,
            output_dir,
            max_usd=max_usd,
            max_output_tokens=max_output_tokens,
            max_requests_per_file=max_requests_per_file,
            max_file_bytes=max_file_bytes,
        )
    except (OSError, UnicodeError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from None
    typer.echo(json.dumps(manifest, ensure_ascii=False, sort_keys=True))


@app.command("parse-output")
def parse_output_command(
    request_manifest: Path = typer.Option(
        ...,
        "--request-manifest",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
    ),
    batch_output_jsonl: list[Path] = typer.Option(
        ...,
        "--batch-output-jsonl",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Repeat for each downloaded OpenAI Batch output or error JSONL.",
    ),
    output_dir: Path = typer.Option(..., "--output-dir", file_okay=False),
    strict: bool = typer.Option(
        True,
        "--strict/--no-strict",
        help="Suppress the accepted-results artifact when any audit error exists.",
    ),
) -> None:
    """Parse downloaded raw Batch output without contacting OpenAI."""

    try:
        report = parse_batch_outputs(
            request_manifest,
            batch_output_jsonl,
            output_dir,
            strict=strict,
        )
    except (OSError, UnicodeError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from None
    typer.echo(json.dumps(report, ensure_ascii=False, sort_keys=True))
    if report["strict_rejected"]:
        raise typer.Exit(code=2)
