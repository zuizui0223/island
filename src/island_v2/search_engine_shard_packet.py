"""Collect raw search-engine evidence packets for one deterministic campaign shard.

This is the high-throughput first pass: exactly the two short exact-name queries used by
``search_engine_trait_recovery.PRIMARY_QUERIES`` are run for every species in a prepared
shard manifest. No page fetching, literature API, trait regex extraction, or genus fallback
is performed here. The goal is to freeze broad search context cheaply before any gap-only
follow-up.
"""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import httpx
import pandas as pd
import typer

from island_v2.search_engine_trait_recovery import (
    PRIMARY_QUERIES,
    USER_AGENT,
    clean,
    search_with_fallback,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _load_jobs(path: Path) -> list[dict[str, str]]:
    jobs: list[dict[str, str]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw.strip():
            continue
        row = json.loads(raw)
        species = clean(row.get("species"))
        if species:
            jobs.append({**row, "species": species})
    return jobs


def _one_species(job: dict[str, str], results_per_query: int) -> dict[str, object]:
    species = str(job["species"])
    headers = {"User-Agent": USER_AGENT, "Accept-Language": "en;q=0.9,*;q=0.5"}
    searches: list[dict[str, object]] = []
    with httpx.Client(headers=headers, follow_redirects=True) as client:
        for query_family, template in PRIMARY_QUERIES.items():
            query = template.format(name=species)
            results: list[dict[str, str]] = []
            for rank, (engine, row) in enumerate(
                search_with_fallback(client, query, results_per_query), start=1
            ):
                results.append(
                    {
                        "rank": rank,
                        "engine": engine,
                        "title": clean(row.get("title")),
                        "snippet": clean(row.get("snippet")),
                        "url": clean(row.get("url")),
                    }
                )
            searches.append(
                {
                    "query_family": query_family,
                    "query": query,
                    "results": results,
                }
            )
    return {
        "job_id": job.get("job_id", ""),
        "species": species,
        "genus": job.get("genus", ""),
        "family": job.get("family", ""),
        "shard_index": job.get("shard_index", ""),
        "shard_count": job.get("shard_count", ""),
        "prompt_version": job.get("prompt_version", ""),
        "prompt_sha256": job.get("prompt_sha256", ""),
        "searches": searches,
    }


@app.command()
def collect(
    jobs_jsonl: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    results_per_query: int = typer.Option(5, min=1, max=20),
    workers: int = typer.Option(8, min=1, max=32),
) -> None:
    jobs = _load_jobs(jobs_jsonl)
    output_dir.mkdir(parents=True, exist_ok=True)
    packets_path = output_dir / "search_packets.jsonl"
    errors_path = output_dir / "search_errors.csv"

    packets: list[dict[str, object]] = []
    errors: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_one_species, job, results_per_query): job for job in jobs}
        for future in as_completed(futures):
            job = futures[future]
            try:
                packets.append(future.result())
            except Exception as exc:
                errors.append(
                    {
                        "job_id": str(job.get("job_id", "")),
                        "species": str(job.get("species", "")),
                        "error": str(exc),
                    }
                )

    order = {str(job["species"]): i for i, job in enumerate(jobs)}
    packets.sort(key=lambda row: order.get(str(row["species"]), 10**9))
    packets_path.write_text(
        "\n".join(json.dumps(row, ensure_ascii=False, sort_keys=True) for row in packets)
        + ("\n" if packets else ""),
        encoding="utf-8",
    )
    pd.DataFrame(errors, columns=["job_id", "species", "error"]).to_csv(errors_path, index=False)

    n_queries = sum(len(row["searches"]) for row in packets)
    n_results = sum(
        len(search["results"])
        for row in packets
        for search in row["searches"]
    )
    n_species_with_results = sum(
        any(search["results"] for search in row["searches"])
        for row in packets
    )
    status = {
        "n_jobs": len(jobs),
        "n_packets": len(packets),
        "n_errors": len(errors),
        "n_queries": n_queries,
        "n_results": n_results,
        "n_species_with_any_search_result": n_species_with_results,
        "query_policy": "two short exact-name searches per species; no page fetch; no literature API; no genus fallback",
    }
    (output_dir / "search_packet_status.json").write_text(
        json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(status, ensure_ascii=False))


if __name__ == "__main__":
    app()
