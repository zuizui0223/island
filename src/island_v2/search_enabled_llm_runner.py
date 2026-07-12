"""Execute search-enabled LLM campaign jobs with checkpointed OpenAI Responses API calls."""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import httpx
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

RESPONSES_URL = "https://api.openai.com/v1/responses"


def _load_completed(path: Path) -> set[str]:
    if not path.exists():
        return set()
    completed: set[str] = set()
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw.strip():
            continue
        try:
            payload = json.loads(raw)
            job_id = str(payload.get("job_id", "")).strip()
            if job_id:
                completed.add(job_id)
        except json.JSONDecodeError:
            continue
    return completed


def _extract_output_text(payload: dict[str, object]) -> str:
    direct = payload.get("output_text")
    if isinstance(direct, str) and direct.strip():
        return direct.strip()
    chunks: list[str] = []
    output = payload.get("output", [])
    if isinstance(output, list):
        for item in output:
            if not isinstance(item, dict) or item.get("type") != "message":
                continue
            content = item.get("content", [])
            if not isinstance(content, list):
                continue
            for part in content:
                if not isinstance(part, dict):
                    continue
                text = part.get("text")
                if isinstance(text, str) and text.strip():
                    chunks.append(text.strip())
    if not chunks:
        raise ValueError("Responses API payload contained no output text")
    return "\n".join(chunks)


def _call_openai(
    client: httpx.Client,
    *,
    api_key: str,
    model: str,
    prompt: str,
    timeout_seconds: float,
) -> tuple[str, str]:
    response = client.post(
        RESPONSES_URL,
        headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
        json={
            "model": model,
            "tools": [{"type": "web_search"}],
            "input": prompt,
        },
        timeout=timeout_seconds,
    )
    response.raise_for_status()
    payload = response.json()
    return _extract_output_text(payload), str(payload.get("id", ""))


@app.command("run")
def run(
    jobs_jsonl: Path = typer.Option(..., exists=True),
    results_jsonl: Path = typer.Option(...),
    model: str = typer.Option("gpt-5.5"),
    max_jobs: int = typer.Option(100, min=1, max=10000),
    max_retries: int = typer.Option(4, min=0, max=10),
    timeout_seconds: float = typer.Option(180.0, min=10.0),
    pause_seconds: float = typer.Option(0.5, min=0.0),
) -> None:
    """Run pending jobs and append each successful result immediately for safe resume."""
    api_key = os.environ.get("OPENAI_API_KEY", "").strip()
    if not api_key:
        raise typer.BadParameter("OPENAI_API_KEY is required")

    jobs = [json.loads(raw) for raw in jobs_jsonl.read_text(encoding="utf-8").splitlines() if raw.strip()]
    completed = _load_completed(results_jsonl)
    pending = [job for job in jobs if str(job.get("job_id", "")) not in completed][:max_jobs]
    results_jsonl.parent.mkdir(parents=True, exist_ok=True)

    n_success = 0
    n_failed = 0
    with httpx.Client() as client, results_jsonl.open("a", encoding="utf-8") as handle:
        for position, job in enumerate(pending, start=1):
            job_id = str(job["job_id"])
            prompt = str(job["prompt"])
            last_error = ""
            for attempt in range(max_retries + 1):
                try:
                    result, response_id = _call_openai(
                        client,
                        api_key=api_key,
                        model=model,
                        prompt=prompt,
                        timeout_seconds=timeout_seconds,
                    )
                    record = {
                        "job_id": job_id,
                        "result": result,
                        "response_id": response_id,
                        "model": model,
                    }
                    handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")
                    handle.flush()
                    n_success += 1
                    typer.echo(f"[{position}/{len(pending)}] completed {job_id}")
                    break
                except (httpx.HTTPError, ValueError) as exc:
                    last_error = str(exc)
                    if attempt >= max_retries:
                        n_failed += 1
                        typer.echo(f"[{position}/{len(pending)}] failed {job_id}: {last_error}", err=True)
                        break
                    time.sleep(min(30.0, 2.0**attempt))
            if pause_seconds:
                time.sleep(pause_seconds)

    summary = {
        "n_jobs_in_manifest": len(jobs),
        "n_already_completed": len(completed),
        "n_attempted_this_run": len(pending),
        "n_success": n_success,
        "n_failed": n_failed,
        "results_jsonl": str(results_jsonl),
        "model": model,
    }
    typer.echo(json.dumps(summary, ensure_ascii=False))
    if n_failed:
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
