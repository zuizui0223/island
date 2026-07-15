from __future__ import annotations

import csv
import hashlib
import json
import tomllib
from collections import Counter
from decimal import Decimal
from pathlib import Path

import pytest
from typer.testing import CliRunner

from island_v2.openai_batch_llm import (
    AUDIT_COLUMNS,
    BATCH_INPUT_PRICE_PER_MILLION,
    BATCH_OUTPUT_PRICE_PER_MILLION,
    MODEL_ID,
    PARSE_AUDIT_FILENAME,
    PROVIDER_RESULTS_FILENAME,
    REQUEST_MANIFEST_FILENAME,
    app,
    parse_batch_outputs,
    prepare_batch_requests,
)
from island_v2.search_enabled_llm_campaign import (
    JOB_HASH_COLUMNS,
    _job_record,
    app as campaign_app,
    canonical_job_sha,
)


def _job(species: str, index: int = 0) -> dict[str, str]:
    genus = species.split()[0]
    return _job_record(
        {
            "species": species,
            "genus": genus,
            "family": f"{genus}aceae",
            "target_fields": "flower_color|pollination_guild",
        },
        index,
        16,
    )


def _write_jobs(path: Path, jobs: list[dict[str, str]]) -> None:
    path.write_text(
        "".join(json.dumps(job, ensure_ascii=False, sort_keys=True) + "\n" for job in jobs),
        encoding="utf-8",
    )


def _batch_output(
    custom_id: str,
    *,
    result: str = "Plantus alba,white,tubular,bees,inferred,mixed_mating,likely_SC,inference,low",
    model: str = MODEL_ID,
    status_code: int = 200,
    top_error: object = None,
    body_error: object = None,
) -> dict[str, object]:
    return {
        "id": f"batch_req_{custom_id}",
        "custom_id": custom_id,
        "response": {
            "status_code": status_code,
            "request_id": f"req_{custom_id}",
            "body": {
                "id": f"resp_{custom_id}",
                "object": "response",
                "status": "completed",
                "error": body_error,
                "model": model,
                "output": [
                    {
                        "type": "message",
                        "role": "assistant",
                        "status": "completed",
                        "content": [
                            {
                                "type": "output_text",
                                "text": json.dumps({"result": result}),
                                "annotations": [],
                            }
                        ],
                    }
                ],
            },
        },
        "error": top_error,
    }


def test_prepare_is_deterministic_offline_responses_batch_with_cost_manifest(
    tmp_path: Path,
) -> None:
    jobs_path = tmp_path / "jobs_shard_0000.jsonl"
    jobs = [_job("Plantus alba"), _job("Plantus rubra", 1)]
    _write_jobs(jobs_path, jobs)
    first = tmp_path / "first"
    second = tmp_path / "second"

    manifest = prepare_batch_requests(
        [jobs_path],
        first,
        max_usd="10",
        max_output_tokens=128,
        max_requests_per_file=1,
    )
    second_manifest = prepare_batch_requests(
        [jobs_path],
        second,
        max_usd=Decimal("10"),
        max_output_tokens=128,
        max_requests_per_file=1,
    )

    assert manifest == second_manifest
    assert manifest["model"] == MODEL_ID
    assert manifest["endpoint"] == "/v1/responses"
    assert manifest["n_requests"] == 2
    assert manifest["n_batches"] == 2
    assert manifest["api_calls_made"] is False
    assert manifest["api_key_handled"] is False
    assert manifest["web_search_enabled"] is False
    assert (first / REQUEST_MANIFEST_FILENAME).read_bytes() == (
        second / REQUEST_MANIFEST_FILENAME
    ).read_bytes()

    request_bodies = []
    for batch in manifest["batches"]:
        first_bytes = (first / batch["filename"]).read_bytes()
        second_bytes = (second / batch["filename"]).read_bytes()
        assert first_bytes == second_bytes
        assert len(first_bytes) == batch["bytes"] < 200_000_000
        assert hashlib.sha256(first_bytes).hexdigest() == batch["sha256"]
        line = json.loads(first_bytes)
        assert line["method"] == "POST"
        assert line["url"] == "/v1/responses"
        assert line["custom_id"].startswith("trait_")
        body = line["body"]
        request_bodies.append(body)
        assert body["model"] == MODEL_ID
        assert body["store"] is False
        assert body["reasoning"] == {"effort": "none"}
        assert body["max_output_tokens"] == 128
        assert body["text"]["format"]["type"] == "json_schema"
        assert body["text"]["format"]["strict"] is True
        assert body["text"]["format"]["schema"]["required"] == ["result"]
        assert "tools" not in body
        assert "web_search" not in json.dumps(body)

    prompt_by_id = {job["job_id"]: job["prompt"] for job in jobs}
    for batch in manifest["batches"]:
        request = json.loads((first / batch["filename"]).read_text(encoding="utf-8"))
        assert request["body"]["input"] == prompt_by_id[request["custom_id"]]

    input_bound = sum(
        len(
            json.dumps(
                body,
                ensure_ascii=False,
                separators=(",", ":"),
                sort_keys=True,
            ).encode("utf-8")
        )
        for body in request_bodies
    )
    assert manifest["cost"]["input_tokens_upper_bound"] == input_bound
    expected_cost = Decimal(input_bound) * BATCH_INPUT_PRICE_PER_MILLION / Decimal(
        1_000_000
    ) + Decimal(2 * 128) * BATCH_OUTPUT_PRICE_PER_MILLION / Decimal(1_000_000)
    assert Decimal(manifest["cost"]["worst_case_usd"]) == expected_cost
    assert not list(first.glob("*.tmp"))


def test_real_campaign_artifact_round_trips_through_adapter_and_ingest(
    tmp_path: Path,
) -> None:
    master_path = tmp_path / "master.csv"
    master_path.write_text(
        "accepted_species,genus,family\nPlantus alba,Plantus,Plantaceae\n",
        encoding="utf-8",
    )
    campaign_dir = tmp_path / "campaign"
    result = CliRunner().invoke(
        campaign_app,
        [
            "prepare",
            "--master-csv",
            str(master_path),
            "--campaign-dir",
            str(campaign_dir),
            "--shard-index",
            "0",
            "--shard-count",
            "1",
            "--batch-size",
            "1",
        ],
    )
    assert result.exit_code == 0, result.output

    jobs_path = campaign_dir / "jobs_shard_0000.jsonl"
    job = json.loads(jobs_path.read_text(encoding="utf-8"))
    assert set(JOB_HASH_COLUMNS).issubset(job)
    assert canonical_job_sha(job) == job["job_sha256"]

    output_dir = tmp_path / "openai"
    manifest = prepare_batch_requests(
        [jobs_path],
        output_dir,
        max_usd="1",
        max_output_tokens=128,
    )
    assert manifest["n_requests"] == 1
    request_record = manifest["requests"][0]
    assert request_record["job_sha256"] == job["job_sha256"]
    assert request_record["prompt_sha256"] == job["prompt_sha256"]
    assert request_record["prompt_policy_sha256"] == job["prompt_policy_sha256"]
    batch_request = json.loads(
        (output_dir / manifest["batches"][0]["filename"]).read_text(encoding="utf-8")
    )
    assert batch_request["body"]["input"] == job["prompt"]

    expected_result = (
        "Plantus alba,white,tubular,bees,inference only,mixed_mating,likely_SC,inference,low"
    )
    raw_output = tmp_path / "batch_output.jsonl"
    raw_output.write_text(
        json.dumps(
            _batch_output(job["job_id"], result=expected_result),
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    provider_dir = tmp_path / "provider"
    parse_report = parse_batch_outputs(
        output_dir / REQUEST_MANIFEST_FILENAME,
        [raw_output],
        provider_dir,
    )
    provider_results = provider_dir / PROVIDER_RESULTS_FILENAME
    assert parse_report["n_accepted"] == 1
    assert provider_results.exists()

    ingest = CliRunner().invoke(
        campaign_app,
        [
            "ingest",
            "--campaign-dir",
            str(campaign_dir),
            "--results-jsonl",
            str(provider_results),
            "--model-provider",
            "OpenAI",
            "--model-id",
            MODEL_ID,
        ],
    )
    assert ingest.exit_code == 0, ingest.output
    with (campaign_dir / "validated_trait_results.csv").open(
        encoding="utf-8", newline=""
    ) as handle:
        validated = list(csv.DictReader(handle))
    assert len(validated) == 1
    assert validated[0]["species"] == "Plantus alba"
    assert validated[0]["job_sha256"] == job["job_sha256"]
    assert validated[0]["prompt_sha256"] == job["prompt_sha256"]
    assert validated[0]["model_provider"] == "OpenAI"
    assert validated[0]["model_id"] == MODEL_ID


def test_prepare_rejects_hash_tampering_duplicates_and_budget_before_writes(
    tmp_path: Path,
) -> None:
    tampered_path = tmp_path / "jobs_shard_0000.jsonl"
    tampered = _job("Plantus alba")
    tampered["prompt_sha256"] = "0" * 64
    _write_jobs(tampered_path, [tampered])
    tampered_output = tmp_path / "tampered-output"
    with pytest.raises(ValueError, match="prompt_sha256 mismatch"):
        prepare_batch_requests([tampered_path], tampered_output, max_usd="1")
    assert not tampered_output.exists()

    duplicate_path = tmp_path / "jobs_shard_0001.jsonl"
    duplicate = _job("Plantus rubra", 1)
    _write_jobs(duplicate_path, [duplicate, duplicate])
    duplicate_output = tmp_path / "duplicate-output"
    with pytest.raises(ValueError, match="duplicate job_id"):
        prepare_batch_requests([duplicate_path], duplicate_output, max_usd="1")
    assert not duplicate_output.exists()

    valid_path = tmp_path / "jobs_shard_0002.jsonl"
    _write_jobs(valid_path, [_job("Plantus caerulea", 2)])
    budget_output = tmp_path / "budget-output"
    with pytest.raises(ValueError, match="exceeds max_usd"):
        prepare_batch_requests([valid_path], budget_output, max_usd="0")
    assert not budget_output.exists()


def test_cli_requires_max_usd_and_never_reads_an_api_key(tmp_path: Path) -> None:
    jobs_path = tmp_path / "jobs_shard_0000.jsonl"
    _write_jobs(jobs_path, [_job("Plantus alba")])
    runner = CliRunner()
    missing = runner.invoke(
        app,
        ["prepare", "--jobs-jsonl", str(jobs_path), "--output-dir", str(tmp_path / "o")],
        env={"OPENAI_API_KEY": "must-not-be-read-or-rendered"},
    )
    assert missing.exit_code == 2
    assert "--max-usd" in missing.output
    help_result = runner.invoke(
        app,
        ["prepare", "--help"],
        env={"OPENAI_API_KEY": "must-not-be-read-or-rendered"},
    )
    assert help_result.exit_code == 0
    assert "must-not-be-read-or-rendered" not in help_result.output
    assert "OPENAI_API_KEY" not in help_result.output

    success = runner.invoke(
        app,
        [
            "prepare",
            "--jobs-jsonl",
            str(jobs_path),
            "--output-dir",
            str(tmp_path / "cli-output"),
            "--max-usd",
            "1",
            "--max-output-tokens",
            "64",
        ],
        env={"OPENAI_API_KEY": "must-not-be-read-or-rendered"},
    )
    assert success.exit_code == 0, success.output
    assert "must-not-be-read-or-rendered" not in success.output
    assert (tmp_path / "cli-output" / REQUEST_MANIFEST_FILENAME).exists()


def test_parse_output_accepts_only_manifest_matched_responses_output_text(
    tmp_path: Path,
) -> None:
    jobs_path = tmp_path / "jobs_shard_0000.jsonl"
    job = _job("Plantus alba")
    _write_jobs(jobs_path, [job])
    prepared = tmp_path / "prepared"
    prepare_batch_requests([jobs_path], prepared, max_usd="1", max_output_tokens=128)

    expected_result = (
        "Plantus alba,white,tubular,bees,inference only,mixed_mating,likely_SC,inference,low"
    )
    batch_output = tmp_path / "batch_output.jsonl"
    batch_output.write_text(
        json.dumps(
            _batch_output(job["job_id"], result=expected_result),
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    parsed_dir = tmp_path / "parsed"
    report = parse_batch_outputs(
        prepared / REQUEST_MANIFEST_FILENAME,
        [batch_output],
        parsed_dir,
        strict=True,
    )

    accepted = [
        json.loads(line)
        for line in (parsed_dir / PROVIDER_RESULTS_FILENAME)
        .read_text(encoding="utf-8")
        .splitlines()
    ]
    assert accepted == [
        {
            "job_id": job["job_id"],
            "prompt_sha256": job["prompt_sha256"],
            "result": expected_result,
        }
    ]
    assert report["strict_rejected"] is False
    assert report["n_accepted"] == 1
    with (parsed_dir / PARSE_AUDIT_FILENAME).open(encoding="utf-8", newline="") as handle:
        assert list(csv.DictReader(handle)) == []
    assert not list(parsed_dir.glob("*.tmp"))


def test_parse_output_audits_every_failure_and_strict_suppresses_valid_rows(
    tmp_path: Path,
) -> None:
    jobs_path = tmp_path / "jobs_shard_0000.jsonl"
    jobs = [
        _job("Plantus alba"),
        _job("Plantus rubra", 1),
        _job("Plantus caerulea", 2),
        _job("Plantus lutea", 3),
        _job("Plantus nigra", 4),
        _job("Plantus viridis", 5),
    ]
    _write_jobs(jobs_path, jobs)
    prepared = tmp_path / "prepared"
    manifest = prepare_batch_requests([jobs_path], prepared, max_usd="10")
    custom_ids = [request["custom_id"] for request in manifest["requests"]]

    rows = [
        _batch_output(custom_ids[0], status_code=500),
        _batch_output(
            custom_ids[1],
            top_error={"code": "batch_expired", "message": "not completed"},
        ),
        _batch_output(custom_ids[2]),
        _batch_output(custom_ids[2]),
        _batch_output(custom_ids[3], model="gpt-wrong"),
        _batch_output(custom_ids[5], result="valid inference row"),
    ]
    batch_output = tmp_path / "batch_output.jsonl"
    batch_output.write_text(
        "".join(json.dumps(row) + "\n" for row in rows),
        encoding="utf-8",
    )

    strict_dir = tmp_path / "strict"
    strict_report = parse_batch_outputs(
        prepared / REQUEST_MANIFEST_FILENAME,
        [batch_output],
        strict_dir,
        strict=True,
    )
    assert strict_report["strict_rejected"] is True
    assert strict_report["n_accepted"] == 0
    assert strict_report["n_would_accept_without_strict"] == 1
    assert not (strict_dir / PROVIDER_RESULTS_FILENAME).exists()
    with (strict_dir / PARSE_AUDIT_FILENAME).open(encoding="utf-8", newline="") as handle:
        audit = list(csv.DictReader(handle))
    assert list(audit[0]) == list(AUDIT_COLUMNS)
    codes = Counter(row["error_code"] for row in audit)
    assert codes == Counter(
        {
            "http_error": 1,
            "batch_error": 1,
            "duplicate_output": 2,
            "model_mismatch": 1,
            "missing_output": 1,
        }
    )

    permissive_dir = tmp_path / "permissive"
    permissive_report = parse_batch_outputs(
        prepared / REQUEST_MANIFEST_FILENAME,
        [batch_output],
        permissive_dir,
        strict=False,
    )
    assert permissive_report["strict_rejected"] is False
    assert permissive_report["n_accepted"] == 1
    accepted = [
        json.loads(line)
        for line in (permissive_dir / PROVIDER_RESULTS_FILENAME)
        .read_text(encoding="utf-8")
        .splitlines()
    ]
    assert accepted[0]["job_id"] == custom_ids[5]
    assert accepted[0]["result"] == "valid inference row"


def test_pyproject_registers_openai_batch_llm_entry_point() -> None:
    project = tomllib.loads(
        (Path(__file__).parents[1] / "pyproject.toml").read_text(encoding="utf-8")
    )
    assert project["project"]["scripts"]["island-v2-openai-batch-llm"] == (
        "island_v2.openai_batch_llm:app"
    )
