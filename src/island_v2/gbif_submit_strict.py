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

# GBIF returns 400 for both of these; only the message body distinguishes
# them from an ordinary rejected request. Matching is deliberately loose
# (lower-cased substring) because GBIF does not publish a stable error code
# for either condition.
CAPACITY_ERROR_MARKERS = (
    "simultaneous download",
    "concurrent download",
    "too many active downloads",
    "download limit",
)
GEOMETRY_ERROR_MARKERS = (
    "antimeridian",
    "self-intersect",
    "invalid geometry",
    "unable to parse",
    "geojson",
    "wkt",
)


def response_detail(response: httpx.Response, limit: int = 4_000) -> str:
    """Return a compact GBIF error body suitable for an auditable CSV field."""
    try:
        body = response.text.strip()
    except Exception:  # defensive: an HTTP error should still be recorded
        body = ""
    if len(body) > limit:
        body = body[:limit] + " …[truncated]"
    return f"HTTP {response.status_code} {response.reason_phrase}: {body or '<empty response body>'}"


def classify_rejection(response: httpx.Response) -> str:
    """Distinguish a transient capacity limit from a geometry/permanent error.

    GBIF has no distinct HTTP status or error code for "you already have too
    many active downloads" versus "this request is malformed" -- both come
    back as a 400 with a free-text body. Left unclassified, every rejection
    used to collapse into a generic `submit_failed`, which is why capacity
    limits and antimeridian-invalid blocks previously required a human to
    read the raw response and hand-write a recovery CSV. Recognized classes
    get a status the ledger/reconcile workflows already know how to treat
    (retry later vs. needs a geometry fix); anything else stays
    `submit_failed`.
    """
    try:
        body = response.text.lower()
    except Exception:
        body = ""
    if any(marker in body for marker in CAPACITY_ERROR_MARKERS):
        return "deferred_capacity"
    if any(marker in body for marker in GEOMETRY_ERROR_MARKERS):
        return "geometry_repair_pending"
    return "submit_failed"


@app.command("submit")
def submit(
    request_manifest: Path = typer.Option(..., exists=True),
    max_requests: int = typer.Option(3, min=1, max=100),
    sleep_seconds: float = typer.Option(1.0, min=0.0),
    stop_after_capacity_error: bool = typer.Option(
        True,
        help="Stop submitting further requests in this run once GBIF reports its "
        "concurrent-download limit is reached, instead of burning through the "
        "rest of the batch on requests that will fail the same way.",
    ),
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
    deferred = 0
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
                    status = classify_rejection(response)
                    if status == "deferred_capacity":
                        deferred += 1
                    else:
                        failed += 1
                    table.at[index, "request_status"] = status
                    table.at[index, "error"] = response_detail(response)
                    if status == "deferred_capacity" and stop_after_capacity_error:
                        write_csv_atomic(table, request_manifest)
                        typer.echo(
                            "GBIF reports its concurrent-download limit is reached; "
                            "stopping this batch early rather than rejecting the rest of it."
                        )
                        break
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

    typer.echo(
        f"GBIF submissions accepted: {submitted}; rejected: {failed}; "
        f"deferred (capacity): {deferred}."
    )
    if failed:
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
