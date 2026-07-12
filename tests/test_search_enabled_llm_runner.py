from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from island_v2 import search_enabled_llm_runner as runner


def test_extract_output_text_from_message_content() -> None:
    payload = {
        "output": [
            {"type": "web_search_call", "id": "ws_1"},
            {
                "type": "message",
                "content": [{"type": "output_text", "text": "Species,row,here"}],
            },
        ]
    }
    assert runner._extract_output_text(payload) == "Species,row,here"


def test_runner_resumes_existing_job(monkeypatch, tmp_path: Path) -> None:
    jobs = tmp_path / "jobs.jsonl"
    results = tmp_path / "results.jsonl"
    jobs.write_text(
        "\n".join(
            [
                json.dumps({"job_id": "job_1", "prompt": "one"}),
                json.dumps({"job_id": "job_2", "prompt": "two"}),
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    results.write_text(json.dumps({"job_id": "job_1", "result": "already"}) + "\n", encoding="utf-8")
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")

    calls: list[str] = []

    def fake_call(*args, **kwargs):
        calls.append(kwargs["prompt"])
        return "new,result", "resp_2"

    monkeypatch.setattr(runner, "_call_openai", fake_call)
    result = CliRunner().invoke(
        runner.app,
        [
            "run",
            "--jobs-jsonl",
            str(jobs),
            "--results-jsonl",
            str(results),
            "--max-jobs",
            "10",
            "--pause-seconds",
            "0",
        ],
    )
    assert result.exit_code == 0, result.output
    assert calls == ["two"]
    rows = [json.loads(line) for line in results.read_text(encoding="utf-8").splitlines()]
    assert [row["job_id"] for row in rows] == ["job_1", "job_2"]
