"""Drive a resumable, outcome-blind global trait-extraction campaign.

The campaign operates on the unique species master rather than selecting islands
or regions from observed floral outcomes. It advances four ordered machine-only
source tasks:

1. reproductive and pollen-vector statements from Wikimedia;
2. reproductive and pollen-vector statements from OpenAlex;
3. floral-access statements for species with a direct biotic-vector candidate;
4. alternative-pollinator guild statements for the same conservative subset.

Every machine hit keeps its source excerpt and remains unreviewed. A species with
no hit is recorded as a completed zero-hit lookup, never as a biological absence.
Transient lookup failures are retried up to a configured limit and then retained
as exhausted audit rows so a small number of inaccessible sources cannot block the
entire global campaign forever.
"""

from __future__ import annotations

import json
import re
import time
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2 import reported_ecology_machine as ecology, trait_web_reported_scout as web_reported
from island_v2.trait_source_discovery import openalex_abstract

app = typer.Typer(add_completion=False, no_args_is_help=True)

MASTER_REQUIRED_COLUMNS = {"accepted_species", "family"}
PROHIBITED_MASTER_COLUMNS = {
    "flower_primary_color",
    "floral_symmetry",
    "floral_form",
    "pollen_vector_mode",
    "pollination_functional_guild",
    "self_incompatibility",
    "autonomous_selfing_capacity",
    "mating_system",
    "bombus_occurrence_evidence",
    "model_results",
    "analysis_included",
}
TERMINAL_STATUSES = {"processed", "exhausted", "not_applicable"}
COMMON_CANDIDATE_COLUMNS = [
    "accepted_species",
    "trait_name",
    "candidate_value",
    "source_url",
    "source_citation",
    "source_excerpt",
    "source_type",
    "evidence_scope",
    "matched_term",
    "confidence",
    "campaign_task",
    "campaign_phase",
    "target_for_task",
]


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise typer.BadParameter("global trait campaign config must be a mapping")
    required = {"task_order", "tasks", "excluded_families", "max_attempts"}
    missing = required.difference(config)
    if missing:
        raise typer.BadParameter(f"global trait campaign config missing: {sorted(missing)}")
    order = [str(value) for value in config["task_order"]]
    tasks = config["tasks"]
    if not isinstance(tasks, dict) or set(order).difference(tasks):
        raise typer.BadParameter("every task_order entry must have a tasks definition")
    return config


def load_species_master(master_csv: Path, config: dict[str, Any]) -> pd.DataFrame:
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    missing = MASTER_REQUIRED_COLUMNS.difference(master.columns)
    if missing:
        raise typer.BadParameter(f"species master missing columns: {sorted(missing)}")
    prohibited = PROHIBITED_MASTER_COLUMNS.intersection(master.columns)
    if prohibited:
        raise typer.BadParameter(
            "species master contains prohibited trait/outcome columns: "
            + ", ".join(sorted(prohibited))
        )

    for column in ("accepted_species", "family"):
        master[column] = master[column].map(_text)
    if "genus" not in master.columns:
        master["genus"] = master["accepted_species"].str.split().str[0].fillna("")
    else:
        master["genus"] = master["genus"].map(_text)
    for column in ("n_islands", "n_records"):
        if column not in master.columns:
            master[column] = 0
        master[column] = pd.to_numeric(master[column], errors="coerce").fillna(0).astype(int)

    excluded = {str(value).strip() for value in config.get("excluded_families") or []}
    valid_name = (
        master["accepted_species"].ne("")
        & master["family"].ne("")
        & ~master["accepted_species"].str.contains("×", regex=False)
        & master["accepted_species"].str.split().str.len().ge(2)
    )
    master = master.loc[valid_name & ~master["family"].isin(excluded)].copy()
    master = master.sort_values(
        ["accepted_species", "n_islands", "n_records"],
        ascending=[True, False, False],
    ).drop_duplicates("accepted_species", keep="first")
    return master[["accepted_species", "genus", "family", "n_islands", "n_records"]].reset_index(
        drop=True
    )


def _task_columns(task: str) -> list[str]:
    return [
        f"{task}_status",
        f"{task}_attempts",
        f"{task}_wave_id",
        f"{task}_candidate_count",
        f"{task}_last_error",
    ]


def reconcile_ledger(
    master: pd.DataFrame,
    existing: pd.DataFrame | None,
    config: dict[str, Any],
) -> pd.DataFrame:
    base = master.copy()
    if existing is not None and not existing.empty:
        old = existing.copy()
        old["accepted_species"] = old["accepted_species"].map(_text)
        preserved = [
            column
            for column in old.columns
            if column not in {"genus", "family", "n_islands", "n_records"}
        ]
        base = base.merge(old[preserved], on="accepted_species", how="left", validate="one_to_one")

    for task in config["task_order"]:
        status, attempts, wave, count, error = _task_columns(str(task))
        if status not in base.columns:
            base[status] = "pending"
        base[status] = base[status].fillna("").replace("", "pending")
        if attempts not in base.columns:
            base[attempts] = 0
        base[attempts] = pd.to_numeric(base[attempts], errors="coerce").fillna(0).astype(int)
        if wave not in base.columns:
            base[wave] = ""
        base[wave] = base[wave].fillna("").astype(str)
        if count not in base.columns:
            base[count] = 0
        base[count] = pd.to_numeric(base[count], errors="coerce").fillna(0).astype(int)
        if error not in base.columns:
            base[error] = ""
        base[error] = base[error].fillna("").astype(str)

    if "machine_biotic_candidate" not in base.columns:
        base["machine_biotic_candidate"] = False
    base["machine_biotic_candidate"] = (
        base["machine_biotic_candidate"]
        .astype("string")
        .fillna("")
        .str.lower()
        .isin({"true", "1", "yes", "y"})
    )
    return base.sort_values("accepted_species").reset_index(drop=True)


def task_is_terminal(ledger: pd.DataFrame, task: str) -> bool:
    return bool(ledger[f"{task}_status"].isin(TERMINAL_STATUSES).all())


def task_eligible_mask(
    ledger: pd.DataFrame,
    task: str,
    config: dict[str, Any],
) -> pd.Series:
    """Return rows eligible to advance now under per-species dependencies."""
    if task not in config["tasks"]:
        raise typer.BadParameter(f"unknown global campaign task: {task}")
    mask = ledger[f"{task}_status"].isin({"pending", "retry"})
    dependencies = [str(value) for value in config["tasks"][task].get("depends_on") or []]
    for dependency in dependencies:
        mask &= ledger[f"{dependency}_status"].isin(TERMINAL_STATUSES)
    rule = str(config["tasks"][task].get("eligibility", "all"))
    if rule == "machine_biotic_candidate":
        mask &= ledger["machine_biotic_candidate"].astype(bool)
    elif rule != "all":
        raise typer.BadParameter(f"unsupported eligibility rule for {task}: {rule}")
    return mask.astype(bool)


def prepare_dependent_statuses(
    ledger: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Close gated rows only when that species' dependencies are terminal."""
    result = ledger.copy()
    for task_value in config["task_order"]:
        task = str(task_value)
        rule = str(config["tasks"][task].get("eligibility", "all"))
        if rule != "machine_biotic_candidate":
            continue
        dependencies_ready = pd.Series(True, index=result.index, dtype=bool)
        dependencies = [str(value) for value in config["tasks"][task].get("depends_on") or []]
        for dependency in dependencies:
            dependencies_ready &= result[f"{dependency}_status"].isin(TERMINAL_STATUSES)
        status_column = f"{task}_status"
        mask = (
            dependencies_ready
            & ~result["machine_biotic_candidate"].astype(bool)
            & result[status_column].isin({"pending", "retry"})
        )
        result.loc[mask, status_column] = "not_applicable"
    return result


def choose_active_task(ledger: pd.DataFrame, config: dict[str, Any]) -> str:
    """Choose the first task with at least one currently eligible species."""
    for task_value in config["task_order"]:
        task = str(task_value)
        if task_eligible_mask(ledger, task, config).any():
            return task
    return "complete"


def family_balanced_batch(
    ledger: pd.DataFrame,
    task: str,
    batch_size: int,
    config: dict[str, Any] | None = None,
) -> pd.DataFrame:
    if batch_size < 1:
        raise typer.BadParameter("batch_size must be at least 1")
    if config is None:
        eligible = ledger[f"{task}_status"].isin({"pending", "retry"})
    else:
        eligible = task_eligible_mask(ledger, task, config)
    pending = ledger.loc[eligible].copy()
    if pending.empty:
        return pending
    pending = pending.sort_values(
        ["family", "n_islands", "n_records", "accepted_species"],
        ascending=[True, False, False, True],
    )
    groups = {
        str(family): group.reset_index(drop=True)
        for family, group in pending.groupby("family", sort=True)
    }
    family_order = sorted(groups, key=lambda family: (-len(groups[family]), family))
    rows: list[pd.Series] = []
    round_index = 0
    while len(rows) < batch_size:
        added = False
        for family in family_order:
            group = groups[family]
            if round_index < len(group):
                rows.append(group.iloc[round_index])
                added = True
                if len(rows) >= batch_size:
                    break
        if not added:
            break
        round_index += 1
    if not rows:
        return pending.head(0)
    batch = pd.DataFrame(rows).reset_index(drop=True)
    return batch[["accepted_species", "genus", "family", "n_islands", "n_records"]]


def _next_wave_id(campaign_dir: Path, task: str) -> str:
    waves = campaign_dir / "waves"
    numbers: list[int] = []
    if waves.exists():
        for path in waves.iterdir():
            match = re.match(r"wave_(\d+)_", path.name)
            if match:
                numbers.append(int(match.group(1)))
    number = max(numbers, default=0) + 1
    return f"wave_{number:05d}_{task}"


def _candidate_lane(trait_name: str, config: dict[str, Any]) -> str:
    for task, rule in config["tasks"].items():
        if trait_name in {str(value) for value in rule.get("target_traits") or []}:
            return str(rule.get("phase") or task)
    return "auxiliary"


def normalize_ecology_candidates(
    frame: pd.DataFrame,
    task: str,
    config: dict[str, Any],
) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=COMMON_CANDIDATE_COLUMNS)
    target = {str(value) for value in config["tasks"][task].get("target_traits") or []}
    rows = []
    for record in frame.fillna("").to_dict("records"):
        trait_name = _text(record.get("trait_name"))
        rows.append(
            {
                "accepted_species": _text(record.get("accepted_species")),
                "trait_name": trait_name,
                "candidate_value": _text(record.get("candidate_value")),
                "source_url": _text(record.get("source_url")),
                "source_citation": _text(record.get("source_citation")),
                "source_excerpt": _text(record.get("source_excerpt")),
                "source_type": _text(record.get("source_type")),
                "evidence_scope": _text(record.get("evidence_scope")) or "species_direct",
                "matched_term": _text(record.get("matched_term")),
                "confidence": _text(record.get("confidence")) or "machine_extracted",
                "campaign_task": task,
                "campaign_phase": _candidate_lane(trait_name, config),
                "target_for_task": trait_name in target,
            }
        )
    return pd.DataFrame(rows, columns=COMMON_CANDIDATE_COLUMNS).drop_duplicates()


def normalize_web_candidates(
    frame: pd.DataFrame,
    task: str,
    config: dict[str, Any],
) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=COMMON_CANDIDATE_COLUMNS)
    target = {str(value) for value in config["tasks"][task].get("target_traits") or []}
    rows = []
    for record in frame.fillna("").to_dict("records"):
        trait_name = _text(record.get("trait_name"))
        excerpt = _text(record.get("raw_description"))
        if trait_name in target and not web_reported.FLORAL_CONTEXT_TERMS.search(excerpt):
            continue
        rows.append(
            {
                "accepted_species": _text(record.get("accepted_species")),
                "trait_name": trait_name,
                "candidate_value": _text(record.get("provisional_candidate_value")),
                "source_url": _text(record.get("source_url")),
                "source_citation": _text(record.get("source")),
                "source_excerpt": excerpt,
                "source_type": _text(record.get("source_type")),
                "evidence_scope": _text(record.get("evidence_scope")) or "species_direct",
                "matched_term": _text(record.get("matched_term")),
                "confidence": "explicit_term_match",
                "campaign_task": task,
                "campaign_phase": _candidate_lane(trait_name, config),
                "target_for_task": trait_name in target,
            }
        )
    return pd.DataFrame(rows, columns=COMMON_CANDIDATE_COLUMNS).drop_duplicates()


def _openalex_terms(config: dict[str, Any], task: str) -> str:
    terms = config["tasks"][task].get("query_terms") or []
    return " ".join(str(value) for value in terms)


def fetch_openalex_sources(
    batch: pd.DataFrame,
    task: str,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, list[str]]:
    import httpx

    terms = _openalex_terms(config, task)
    per_page = int(config["tasks"][task].get("per_page", 5))
    pause_seconds = float(config["tasks"][task].get("pause_seconds", 0.05))
    rows: list[dict[str, str]] = []
    errors: list[str] = []
    with httpx.Client(
        timeout=45.0,
        follow_redirects=True,
        headers={"User-Agent": "island-floral-v2 global-trait-campaign"},
    ) as client:
        for species in batch["accepted_species"].astype(str):
            try:
                response = client.get(
                    ecology.OPENALEX_WORKS_URL,
                    params={"search": f'"{species}" {terms}', "per_page": per_page},
                )
                response.raise_for_status()
                works = response.json().get("results") or []
                for work in works:
                    title = _text(re.sub(r"<[^>]+>", " ", work.get("display_name") or ""))
                    abstract = _text(openalex_abstract(work))
                    source_text = _text(f"{title}. {abstract}")
                    if not source_text:
                        continue
                    scope = (
                        "species_direct"
                        if species.casefold() in source_text.casefold()
                        else "species_indirect"
                    )
                    rows.append(
                        {
                            "accepted_species": species,
                            "source_text": source_text,
                            "source_url": _text(work.get("doi") or work.get("id")),
                            "source_citation": _text(
                                f"{title} {work.get('publication_year') or ''}"
                            ),
                            "source_type": "openalex_title_abstract",
                            "evidence_scope": scope,
                        }
                    )
            except Exception as exc:  # noqa: BLE001
                errors.append(f"openalex:{species}:{exc}")
            if pause_seconds > 0:
                time.sleep(pause_seconds)
    return pd.DataFrame(rows, columns=ecology.TEXT_SOURCE_COLUMNS), errors


def _error_species(errors: list[str]) -> set[str]:
    species: set[str] = set()
    for error in errors:
        parts = str(error).split(":", 2)
        if len(parts) >= 2 and parts[1].strip():
            species.add(parts[1].strip())
    return species


def run_source_task(
    batch: pd.DataFrame,
    task: str,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], set[str]]:
    source_kind = str(config["tasks"][task].get("source_kind"))
    ontology = ecology.load_ontology(Path("config/trait_ontology.yml"))
    source_species: set[str] = set()

    if source_kind == "wikimedia_reported_ecology":
        sources, errors = ecology.wikimedia_text_sources(
            batch,
            web_reported._httpx_getter(  # noqa: SLF001
                pause_seconds=float(config["tasks"][task].get("pause_seconds", 0.05)),
                max_retries=int(config["tasks"][task].get("max_retries", 4)),
            ),
            max_taxa=len(batch),
        )
        source_species = (
            set(sources["accepted_species"].astype(str)) if not sources.empty else set()
        )
        candidates, holdouts = ecology.extract_candidates_from_text_sources(
            sources,
            ontology,
            extraction_model="rule_based_wikimedia_global_campaign_v1",
            prompt_id=f"{task}_rules_v1",
        )
        return (
            normalize_ecology_candidates(candidates, task, config),
            holdouts,
            errors,
            source_species,
        )

    if source_kind == "openalex_reported_ecology":
        sources, errors = fetch_openalex_sources(batch, task, config)
        source_species = (
            set(sources["accepted_species"].astype(str)) if not sources.empty else set()
        )
        candidates, holdouts = ecology.extract_candidates_from_text_sources(
            sources,
            ontology,
            extraction_model="rule_based_openalex_global_campaign_v1",
            prompt_id=f"{task}_rules_v1",
        )
        return (
            normalize_ecology_candidates(candidates, task, config),
            holdouts,
            errors,
            source_species,
        )

    if source_kind == "wikimedia_web_reported":
        paths = {
            "M0": Path(config["tasks"][task]["keyword_config"]),
            "M1": Path("config/m1_reported_keywords.yml"),
        }
        keywords, trait_layers = web_reported.load_keyword_layers(paths)
        raw, report = web_reported.scout_web_reported(
            batch,
            web_reported._httpx_getter(  # noqa: SLF001
                pause_seconds=float(config["tasks"][task].get("pause_seconds", 0.05)),
                max_retries=int(config["tasks"][task].get("max_retries", 4)),
            ),
            keywords,
            trait_layers,
            max_taxa=len(batch),
        )
        errors = [str(value) for value in report.get("lookup_errors") or []]
        failed = _error_species(errors)
        source_species = set(batch["accepted_species"].astype(str)).difference(failed)
        return (
            normalize_web_candidates(raw, task, config),
            pd.DataFrame(columns=["accepted_species", "holdout_reason"]),
            errors,
            source_species,
        )

    raise typer.BadParameter(f"unsupported global campaign source_kind: {source_kind}")


def apply_task_result(
    ledger: pd.DataFrame,
    batch: pd.DataFrame,
    candidates: pd.DataFrame,
    errors: list[str],
    source_species: set[str],
    task: str,
    wave_id: str,
    config: dict[str, Any],
) -> pd.DataFrame:
    result = ledger.copy()
    max_attempts = int(config.get("max_attempts", 3))
    error_by_species: dict[str, list[str]] = {}
    for error in errors:
        parts = str(error).split(":", 2)
        if len(parts) >= 2:
            error_by_species.setdefault(parts[1], []).append(error)

    target = candidates.loc[candidates["target_for_task"].astype(bool)].copy()
    target_counts = target.groupby("accepted_species").size().to_dict() if not target.empty else {}
    direct_biotic = set(
        candidates.loc[
            candidates["trait_name"].eq("pollen_vector_mode")
            & candidates["candidate_value"].eq("biotic")
            & candidates["evidence_scope"].eq("species_direct"),
            "accepted_species",
        ].astype(str)
    )

    for species in batch["accepted_species"].astype(str):
        mask = result["accepted_species"].eq(species)
        attempts_column = f"{task}_attempts"
        attempts = int(result.loc[mask, attempts_column].iloc[0]) + 1
        result.loc[mask, attempts_column] = attempts
        result.loc[mask, f"{task}_wave_id"] = wave_id
        result.loc[mask, f"{task}_candidate_count"] = int(target_counts.get(species, 0))
        if species in source_species or species not in error_by_species:
            result.loc[mask, f"{task}_status"] = "processed"
            result.loc[mask, f"{task}_last_error"] = ""
        else:
            status = "exhausted" if attempts >= max_attempts else "retry"
            result.loc[mask, f"{task}_status"] = status
            result.loc[mask, f"{task}_last_error"] = " | ".join(error_by_species[species])[:1000]
        if species in direct_biotic:
            result.loc[mask, "machine_biotic_candidate"] = True
    return prepare_dependent_statuses(result, config)


def campaign_summary(ledger: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    tasks: dict[str, Any] = {}
    for task_value in config["task_order"]:
        task = str(task_value)
        counts = ledger[f"{task}_status"].value_counts().sort_index()
        tasks[task] = {
            "phase": str(config["tasks"][task].get("phase") or task),
            "status_counts": {str(key): int(value) for key, value in counts.items()},
            "n_candidates": int(ledger[f"{task}_candidate_count"].sum()),
            "n_eligible_pending": int(task_eligible_mask(ledger, task, config).sum()),
            "terminal": task_is_terminal(ledger, task),
        }
    return {
        "version": str(config.get("version", "")),
        "n_global_species": int(len(ledger)),
        "n_families": int(ledger["family"].nunique()),
        "n_machine_biotic_candidates": int(ledger["machine_biotic_candidate"].sum()),
        "active_task": choose_active_task(ledger, config),
        "tasks": tasks,
        "interpretation": (
            "Machine-only source-backed candidates. Zero hits are lookup outcomes, not biological "
            "absence; accepted trait values still require the evidence-review contract."
        ),
    }


def _write_frame_gzip(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, compression="gzip")


@app.command("run")
def run(
    master_csv: Path = typer.Option(..., exists=True),
    campaign_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/global_trait_campaign.yml"), exists=True),
    batch_size: int | None = typer.Option(None, min=1, max=2000),
    task: str | None = typer.Option(
        None,
        "--task",
        help="Advance one named task using per-species dependency eligibility.",
    ),
) -> None:
    """Advance exactly one bounded global campaign wave and persist its audit state."""
    config = load_config(config_path)
    master = load_species_master(master_csv, config)
    ledger_path = campaign_dir / "campaign_ledger.csv.gz"
    existing = pd.read_csv(ledger_path, dtype=str).fillna("") if ledger_path.exists() else None
    ledger = reconcile_ledger(master, existing, config)
    ledger = prepare_dependent_statuses(ledger, config)
    requested_task = _text(task)
    if requested_task and requested_task not in config["tasks"]:
        raise typer.BadParameter(f"unknown global campaign task: {requested_task}")
    task = requested_task or choose_active_task(ledger, config)
    campaign_dir.mkdir(parents=True, exist_ok=True)

    if requested_task and not task_eligible_mask(ledger, task, config).any():
        _write_frame_gzip(ledger, ledger_path)
        summary = campaign_summary(ledger, config)
        status_path = campaign_dir / "campaign_status.json"
        if status_path.exists():
            previous = json.loads(status_path.read_text(encoding="utf-8"))
            if "last_wave" in previous:
                summary["last_wave"] = previous["last_wave"]
        status_path.write_text(
            json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        typer.echo(f"No currently eligible species for requested task {task}.")
        return

    if task == "complete":
        _write_frame_gzip(ledger, ledger_path)
        summary = campaign_summary(ledger, config)
        (campaign_dir / "campaign_status.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        typer.echo("Global trait campaign is complete under the configured machine-source tasks.")
        return

    size = int(batch_size or config.get("default_batch_size", 250))
    batch = family_balanced_batch(ledger, task, size, config)
    if batch.empty:
        raise RuntimeError(f"Task {task} is active but produced no pending batch")
    wave_id = _next_wave_id(campaign_dir, task)
    candidates, holdouts, errors, source_species = run_source_task(batch, task, config)
    ledger = apply_task_result(
        ledger,
        batch,
        candidates,
        errors,
        source_species,
        task,
        wave_id,
        config,
    )

    wave_dir = campaign_dir / "waves" / wave_id
    wave_dir.mkdir(parents=True, exist_ok=False)
    batch.to_csv(wave_dir / "species_batch.csv", index=False)
    _write_frame_gzip(candidates, wave_dir / "machine_candidates.csv.gz")
    _write_frame_gzip(holdouts, wave_dir / "holdouts.csv.gz")
    pd.DataFrame({"error": errors}).to_csv(wave_dir / "lookup_errors.csv", index=False)
    wave_summary = {
        "wave_id": wave_id,
        "task": task,
        "phase": str(config["tasks"][task].get("phase") or task),
        "source_kind": str(config["tasks"][task].get("source_kind")),
        "n_species_attempted": int(len(batch)),
        "n_species_with_source_text": int(len(source_species)),
        "n_candidates_total": int(len(candidates)),
        "n_target_candidates": int(candidates["target_for_task"].astype(bool).sum()),
        "n_lookup_errors": int(len(errors)),
        "policy": "Machine-only; no candidate is curated or used as a biological absence.",
    }
    (wave_dir / "wave_summary.json").write_text(
        json.dumps(wave_summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    _write_frame_gzip(ledger, ledger_path)
    summary = campaign_summary(ledger, config)
    summary["last_wave"] = wave_summary
    (campaign_dir / "campaign_status.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(
        f"Advanced {wave_id}: {len(batch)} species, {len(candidates)} candidate(s), "
        f"{len(errors)} lookup error(s). Next task: {summary['active_task']}."
    )


if __name__ == "__main__":
    app()
