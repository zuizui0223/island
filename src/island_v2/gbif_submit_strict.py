"""Strict, diagnosable GBIF occurrence-download submission.

This wrapper is used for external-download pilots.  Unlike the legacy submit
command, it records GBIF's HTTP response body in the request manifest and exits
non-zero when any attempted request is rejected.  Thus a failed submission can
never be mistaken for a successful acquisition run.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import httpx
import typer

from island_v2.gbif_flora import GBIF_API, gbif_credentials, read_request_manifest, write_csv_atomic

app = typer.Typer(add_completion=False, no_args_is_help=True)


def response_detail(response: httpx.Response, limit: int = 4_000) -> str:
    """Return a compact GBIF error body suitable for an auditable CSV field."""
    try:
        body = response.text.strip()
    except Exception:  # defensive: an HTTP error should still be recorded
        body = ""
    if len(body) > limit:
        body = body[:limit] + " …[truncated]"
    return f"HTTP {response.status_code} {response.reason_phrase}: {body or '<empty response body>'}"


@app.command("submit")
def submit(
    request_manifest: Path = typer.Option(..., exists=True),
    max_requests: int = typer.Option(3, min=1, max=100),
    sleep_seconds: float = typer.Option(1.0, min=0.0),
) -> None:
    """Submit a bounded batch and fail loudly if GBIF rejects any request."""
    username, password = gbif_credentials()
    table = read_request_manifest(request_manifest)
    todo = table.index[table["request_status"].eq("prepared")].tolist()[:max_requests]
    if not todo:
        typer.echo("No prepared requests remain.")
        return

    submitted = 0
    failed = 0
    with httpx.Client(
        timeout=60.0,
        auth=(username, password),
        headers={"User-Agent": "island-floral-v2/0.2", "Accept": "application/json"},
    ) as client:
        for index in todo:
            payload = json.loads(table.at[index, "payload_json"])
            try:
                response = client.post(f"{GBIF_API}/occurrence/download/request", json=payload)
                if response.is_error:
                    failed += 1
                    table.at[index, "request_status"] = "submit_failed"
                    table.at[index, "error"] = response_detail(response)
                else:
                    download_key = response.text.strip().strip('"')
                    if not download_key:
                        failed += 1
                        table.at[index, "request_status"] = "submit_failed"
                        table.at[index, "error"] = response_detail(response)
                    else:
                        submitted += 1
                        table.at[index, "download_key"] = download_key
                        table.at[index, "request_status"] = "submitted"
                        table.at[index, "error"] = ""
            except Exception as exc:
                failed += 1
                table.at[index, "request_status"] = "submit_failed"
                table.at[index, "error"] = f"Client exception: {exc}"
            write_csv_atomic(table, request_manifest)
            time.sleep(sleep_seconds)

    typer.echo(f"GBIF submissions accepted: {submitted}; rejected: {failed}.")
    if failed:
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
