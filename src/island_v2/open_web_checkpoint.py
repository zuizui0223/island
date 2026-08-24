"""Immutable task/range checkpointing for resumable open-Web production."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _overlap(left: tuple[int, int], right: tuple[int, int]) -> bool:
    return max(left[0], right[0]) < min(left[1], right[1])


def _prior_manifests(root: Path | None) -> list[dict[str, Any]]:
    if root is None or not root.exists():
        return []
    rows: list[dict[str, Any]] = []
    for path in sorted(root.glob("**/checkpoint_manifest.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        payload["_manifest_path"] = str(path)
        rows.append(payload)
    return rows


def prepare_shard(
    *,
    tasks_csv: Path,
    output_dir: Path,
    offset: int,
    limit: int,
    source_run_id: str,
    prior_root: Path | None = None,
    resume_run_id: str = "",
) -> dict[str, Any]:
    """Reject overlapping shards and remove every previously completed task."""

    if offset < 0 or limit <= 0:
        raise ValueError("offset must be nonnegative and limit must be positive")
    tasks = pd.read_csv(tasks_csv, dtype=str).fillna("")
    required = {"task_id", "accepted_species", "trait_name", "query"}
    missing = required.difference(tasks.columns)
    if missing:
        raise ValueError(f"task table missing columns: {sorted(missing)}")
    if tasks["task_id"].eq("").any() or tasks["task_id"].duplicated().any():
        raise ValueError("task_id must be nonempty and globally unique")
    if offset > len(tasks):
        raise ValueError(f"offset {offset} exceeds task count {len(tasks)}")
    end = min(offset + limit, len(tasks))
    requested = (offset, end)
    queue_sha256 = _sha256(tasks_csv)
    priors = _prior_manifests(prior_root)
    completed: set[str] = set()
    resumed_manifest: dict[str, Any] | None = None
    for prior in priors:
        if prior.get("queue_sha256") != queue_sha256:
            raise ValueError("prior checkpoint queue hash differs from the immutable task queue")
        completed.update(str(value) for value in prior.get("completed_task_ids", []))
        prior_range = (
            int(prior.get("range_start", 0)),
            int(prior.get("range_end", 0)),
        )
        if not _overlap(requested, prior_range):
            continue
        exact_resume = bool(
            resume_run_id
            and str(prior.get("source_run_id")) == str(resume_run_id)
            and prior_range == requested
        )
        if not exact_resume:
            raise ValueError(
                "requested range overlaps prior run "
                f"{prior.get('source_run_id')} [{prior_range[0]},{prior_range[1]})"
            )
        resumed_manifest = prior

    selected = tasks.iloc[offset:end].copy()
    full_task_ids = selected["task_id"].tolist()
    pending = selected.loc[~selected["task_id"].isin(completed)].copy()
    if resumed_manifest is not None and not len(pending):
        raise ValueError("resume range is already complete; no task may be re-searched")
    output_dir.mkdir(parents=True, exist_ok=True)
    pending.to_csv(output_dir / "search_tasks_shard.csv", index=False)
    plan = {
        "contract": "open_web_production_checkpoint_v1",
        "source_run_id": str(source_run_id),
        "resume_run_id": str(resume_run_id),
        "queue_sha256": queue_sha256,
        "queue_rows": len(tasks),
        "range_start": offset,
        "range_end": end,
        "range_notation": f"[{offset},{end})",
        "range_task_ids": full_task_ids,
        "inherited_completed_task_ids": sorted(completed),
        "pending_task_ids": pending["task_id"].tolist(),
        "pending_tasks": len(pending),
        "overlap_rejected": False,
        "status": "planned",
    }
    (output_dir / "checkpoint_plan.json").write_text(
        json.dumps(plan, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return plan


def finalize_shard(
    *,
    checkpoint_plan_json: Path,
    discovery_report_json: Path,
    output_json: Path,
) -> dict[str, Any]:
    """Persist attempted task IDs so no-hit and error tasks are also resumable."""

    plan = json.loads(checkpoint_plan_json.read_text(encoding="utf-8"))
    discovery = json.loads(discovery_report_json.read_text(encoding="utf-8"))
    attempted = {str(value) for value in discovery.get("attempted_task_ids", [])}
    pending = {str(value) for value in plan.get("pending_task_ids", [])}
    if not attempted.issubset(pending):
        raise ValueError("discovery attempted a task outside the planned shard")
    completed = {str(value) for value in plan.get("inherited_completed_task_ids", [])} | attempted
    range_tasks = {str(value) for value in plan.get("range_task_ids", [])}
    remaining = sorted(range_tasks - completed)
    result = {
        **plan,
        "completed_task_ids": sorted(completed),
        "attempted_task_ids_this_run": sorted(attempted),
        "remaining_task_ids_in_range": remaining,
        "status": "complete" if not remaining else "partial",
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return result


@app.command("prepare")
def prepare_command(
    tasks_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
    offset: Annotated[int, typer.Option(min=0)],
    limit: Annotated[int, typer.Option(min=1)],
    source_run_id: str,
    prior_root: Annotated[
        Path | None,
        typer.Option(exists=True, file_okay=False),
    ] = None,
    resume_run_id: str = "",
) -> None:
    result = prepare_shard(
        tasks_csv=tasks_csv,
        output_dir=output_dir,
        offset=offset,
        limit=limit,
        source_run_id=source_run_id,
        prior_root=prior_root,
        resume_run_id=resume_run_id,
    )
    typer.echo(json.dumps(result, sort_keys=True))


@app.command("finalize")
def finalize_command(
    checkpoint_plan_json: Annotated[
        Path,
        typer.Option(exists=True, dir_okay=False),
    ],
    discovery_report_json: Annotated[
        Path,
        typer.Option(exists=True, dir_okay=False),
    ],
    output_json: Annotated[Path, typer.Option(dir_okay=False)],
) -> None:
    result = finalize_shard(
        checkpoint_plan_json=checkpoint_plan_json,
        discovery_report_json=discovery_report_json,
        output_json=output_json,
    )
    typer.echo(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    app()
