"""Recover multilingual plant traits through English-page Wikidata sitelinks.

Scientific names rarely remain the page title in Japanese, Chinese, or Russian
Wikipedia. This lane resolves an English page to a Wikidata item, retrieves
local-language sitelinks, then downloads the local extract. Direct scientific-
title lookup remains a fallback. All outputs are unreviewed source-backed
machine candidates; missing pages, sitelinks, and matches are not absences.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2 import trait_web_reported_scout as web_reported
from island_v2.global_trait_campaign import _next_wave_id, _write_frame_gzip
from island_v2.multilingual_wikimedia_recovery import (
    DEFAULT_LANGUAGES,
    RECOVERY_PHASE,
    RECOVERY_TASK,
    extract_recovery_candidates,
    fetch_language_extracts,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)
JsonGetter = Callable[[str, dict[str, Any]], dict[str, Any]]

WAVE_TASK = "reproductive_multilingual_sitelink_recovery"
STATUS_COLUMNS = [
    f"{WAVE_TASK}_status",
    f"{WAVE_TASK}_attempts",
    f"{WAVE_TASK}_wave_id",
    f"{WAVE_TASK}_candidate_count",
    f"{WAVE_TASK}_last_error",
]
TERMINAL_PRIMARY = {"processed", "exhausted"}
WIKIDATA_API = "https://www.wikidata.org/w/api.php"


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _normal_title(value: str) -> str:
    return " ".join(str(value or "").replace("_", " ").split()).casefold()


def _follow_alias(title: str, aliases: dict[str, str]) -> str:
    current = _normal_title(title)
    seen: set[str] = set()
    while current in aliases and current not in seen:
        seen.add(current)
        current = aliases[current]
    return current


def _query_aliases(query: dict[str, Any]) -> dict[str, str]:
    aliases: dict[str, str] = {}
    for item in [*(query.get("normalized") or []), *(query.get("redirects") or [])]:
        source = _normal_title(item.get("from") or "")
        target = _normal_title(item.get("to") or "")
        if source and target:
            aliases[source] = target
    return aliases


def fetch_english_wikidata_ids(
    getter: JsonGetter,
    species: list[str],
    *,
    api_batch_size: int = 50,
) -> tuple[dict[str, str], set[str], list[str]]:
    """Resolve scientific-name English pages to Wikidata item IDs."""
    endpoint = web_reported.wikipedia_api_for_language("en")
    qids: dict[str, str] = {}
    completed: set[str] = set()
    errors: list[str] = []
    for start in range(0, len(species), api_batch_size):
        names = species[start : start + api_batch_size]
        try:
            payload = getter(
                endpoint,
                {
                    "action": "query",
                    "prop": "pageprops",
                    "ppprop": "wikibase_item",
                    "redirects": 1,
                    "titles": "|".join(names),
                    "format": "json",
                },
            )
        except Exception as exc:  # noqa: BLE001
            errors.extend(f"enwiki-qid:{name}:transient:{exc}" for name in names)
            continue
        completed.update(names)
        query = payload.get("query") or {}
        aliases = _query_aliases(query)
        pages = {
            _normal_title(page.get("title") or ""): page
            for page in (query.get("pages") or {}).values()
            if "missing" not in page
        }
        for name in names:
            page = pages.get(_follow_alias(name, aliases))
            if page is None:
                continue
            qid = _text((page.get("pageprops") or {}).get("wikibase_item"))
            if qid:
                qids[name] = qid
    return qids, completed, errors


def fetch_wikidata_sitelinks(
    getter: JsonGetter,
    species_to_qid: dict[str, str],
    languages: tuple[str, ...],
    *,
    api_batch_size: int = 50,
) -> tuple[dict[str, dict[str, str]], list[str]]:
    """Return ``species -> language -> local title`` from Wikidata."""
    qid_to_species: dict[str, list[str]] = {}
    for species, qid in species_to_qid.items():
        qid_to_species.setdefault(qid, []).append(species)
    resolved = {species: {} for species in species_to_qid}
    errors: list[str] = []
    qids = sorted(qid_to_species)
    sitefilter = "|".join(f"{language}wiki" for language in languages)
    for start in range(0, len(qids), api_batch_size):
        batch = qids[start : start + api_batch_size]
        try:
            payload = getter(
                WIKIDATA_API,
                {
                    "action": "wbgetentities",
                    "ids": "|".join(batch),
                    "props": "sitelinks",
                    "sitefilter": sitefilter,
                    "format": "json",
                },
            )
        except Exception as exc:  # noqa: BLE001
            for qid in batch:
                for species in qid_to_species.get(qid, []):
                    errors.append(f"wikidata-sitelinks:{species}:transient:{exc}")
            continue
        entities = payload.get("entities") or {}
        for qid in batch:
            sitelinks = (entities.get(qid) or {}).get("sitelinks") or {}
            for species in qid_to_species.get(qid, []):
                for language in languages:
                    title = _text((sitelinks.get(f"{language}wiki") or {}).get("title"))
                    if title:
                        resolved.setdefault(species, {})[language] = title
    return resolved, errors


def fetch_resolved_language_extracts(
    getter: JsonGetter,
    species_to_title: dict[str, str],
    language: str,
    *,
    api_batch_size: int = 50,
) -> tuple[pd.DataFrame, set[str], list[str]]:
    """Fetch local extracts while retaining accepted scientific names."""
    rows: list[dict[str, str]] = []
    completed: set[str] = set()
    errors: list[str] = []
    endpoint = web_reported.wikipedia_api_for_language(language)
    items = list(species_to_title.items())
    for start in range(0, len(items), api_batch_size):
        batch = items[start : start + api_batch_size]
        titles = [title for _, title in batch]
        try:
            payload = getter(
                endpoint,
                {
                    "action": "query",
                    "prop": "extracts",
                    "explaintext": 1,
                    "redirects": 1,
                    "titles": "|".join(titles),
                    "format": "json",
                    "exlimit": "max",
                },
            )
        except Exception as exc:  # noqa: BLE001
            errors.extend(
                f"wikipedia-sitelink-{language}:{species}:transient:{exc}"
                for species, _ in batch
            )
            continue
        completed.update(species for species, _ in batch)
        query = payload.get("query") or {}
        aliases = _query_aliases(query)
        pages = {
            _normal_title(page.get("title") or ""): page
            for page in (query.get("pages") or {}).values()
            if "missing" not in page
        }
        for species, requested_title in batch:
            page = pages.get(_follow_alias(requested_title, aliases))
            if page is None:
                continue
            extract = _text(page.get("extract"))
            if not extract:
                continue
            resolved_title = _text(page.get("title")) or requested_title
            rows.append(
                {
                    "accepted_species": species,
                    "language": language,
                    "source_text": extract,
                    "source_url": f"https://{language}.wikipedia.org/wiki/"
                    + resolved_title.replace(" ", "_"),
                    "source_citation": (
                        f"{resolved_title} Wikipedia full extract ({language}; "
                        "Wikidata sitelink)"
                    ),
                    "source_type": f"wikipedia_sitelink_extract_{language}",
                    "evidence_scope": "species_direct",
                }
            )
    columns = [
        "accepted_species",
        "language",
        "source_text",
        "source_url",
        "source_citation",
        "source_type",
        "evidence_scope",
    ]
    return pd.DataFrame(rows, columns=columns), completed, errors


def _drop_nested_guild_matches(candidates: pd.DataFrame) -> pd.DataFrame:
    """Prefer bumblebee over the nested Japanese broad-bee substring match."""
    if candidates.empty:
        return candidates
    guild = candidates["trait_name"].eq("pollination_functional_guild")
    bumblebee = guild & candidates["candidate_value"].eq("bumblebees")
    other_bee = guild & candidates["candidate_value"].eq("other_bees")
    if not bumblebee.any() or not other_bee.any():
        return candidates
    context_columns = ["accepted_species", "source_url", "source_excerpt"]
    specific_contexts = set(
        candidates.loc[bumblebee, context_columns].itertuples(index=False, name=None)
    )
    broad_context = candidates.loc[other_bee, context_columns].apply(tuple, axis=1)
    nested_indexes = broad_context.index[broad_context.isin(specific_contexts)]
    return candidates.drop(index=nested_indexes).reset_index(drop=True)


def _ensure_status_columns(ledger: pd.DataFrame) -> pd.DataFrame:
    result = ledger.copy()
    defaults: dict[str, object] = {
        STATUS_COLUMNS[0]: "pending",
        STATUS_COLUMNS[1]: 0,
        STATUS_COLUMNS[2]: "",
        STATUS_COLUMNS[3]: 0,
        STATUS_COLUMNS[4]: "",
    }
    for column, default in defaults.items():
        if column not in result.columns:
            result[column] = default
    result[STATUS_COLUMNS[0]] = result[STATUS_COLUMNS[0]].fillna("").replace("", "pending")
    for column in (STATUS_COLUMNS[1], STATUS_COLUMNS[3]):
        result[column] = pd.to_numeric(result[column], errors="coerce").fillna(0).astype(int)
    for column in (STATUS_COLUMNS[2], STATUS_COLUMNS[4]):
        result[column] = result[column].fillna("").astype(str)
    return result


def _family_balanced_batch(ledger: pd.DataFrame, batch_size: int) -> pd.DataFrame:
    eligible = (
        ledger["reproductive_wikimedia_status"].isin(TERMINAL_PRIMARY)
        & ledger[STATUS_COLUMNS[0]].isin({"pending", "retry"})
    )
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
    order = sorted(groups, key=lambda family: (-len(groups[family]), family))
    selected: list[pd.Series] = []
    round_index = 0
    while len(selected) < batch_size:
        added = False
        for family in order:
            group = groups[family]
            if round_index < len(group):
                selected.append(group.iloc[round_index])
                added = True
                if len(selected) >= batch_size:
                    break
        if not added:
            break
        round_index += 1
    return pd.DataFrame(selected).reset_index(drop=True)


def _biotic_species(candidates: pd.DataFrame) -> set[str]:
    if candidates.empty:
        return set()
    selected = candidates.loc[
        candidates["trait_name"].eq("pollen_vector_mode")
        & candidates["candidate_value"].eq("biotic")
        & candidates["evidence_scope"].eq("species_direct")
    ]
    return set(selected["accepted_species"].astype(str))


def _apply_results(
    ledger: pd.DataFrame,
    names: list[str],
    candidates: pd.DataFrame,
    completed_species: set[str],
    errors: list[str],
    wave_id: str,
) -> tuple[pd.DataFrame, set[str]]:
    result = ledger.copy()
    counts = candidates.groupby("accepted_species").size().to_dict() if not candidates.empty else {}
    errors_by_species: dict[str, list[str]] = {}
    for error in errors:
        parts = str(error).split(":", 2)
        if len(parts) >= 2:
            errors_by_species.setdefault(parts[1], []).append(error)
    biotic = _biotic_species(candidates)
    for species in names:
        mask = result["accepted_species"].eq(species)
        attempts = int(result.loc[mask, STATUS_COLUMNS[1]].iloc[0]) + 1
        result.loc[mask, STATUS_COLUMNS[1]] = attempts
        result.loc[mask, STATUS_COLUMNS[2]] = wave_id
        result.loc[mask, STATUS_COLUMNS[3]] = int(counts.get(species, 0))
        result.loc[mask, STATUS_COLUMNS[0]] = (
            "processed" if species in completed_species else "retry"
        )
        result.loc[mask, STATUS_COLUMNS[4]] = (
            "" if species in completed_species else " | ".join(errors_by_species.get(species, []))[:1000]
        )
        if species not in biotic:
            continue
        result.loc[mask, "machine_biotic_candidate"] = "True"
        for dependent in (
            "floral_access_wikimedia_status",
            "alternative_guild_wikimedia_status",
            "alternative_guild_openalex_status",
        ):
            if dependent in result.columns:
                reopen = mask & result[dependent].eq("not_applicable")
                result.loc[reopen, dependent] = "pending"
    return result, biotic


def run_sitelink_recovery_wave(
    *,
    campaign_dir: Path,
    batch_size: int = 1000,
    languages: tuple[str, ...] = DEFAULT_LANGUAGES,
    getter: JsonGetter | None = None,
) -> dict[str, Any]:
    """Advance one Wikidata-sitelink multilingual recovery wave."""
    ledger_path = campaign_dir / "campaign_ledger.csv.gz"
    if not ledger_path.exists():
        raise typer.BadParameter(f"campaign ledger not found: {ledger_path}")
    ledger = _ensure_status_columns(pd.read_csv(ledger_path, dtype=str).fillna(""))
    batch = _family_balanced_batch(ledger, batch_size)
    if batch.empty:
        return {
            "status": "complete_for_current_primary_progress",
            "n_species_attempted": 0,
            "languages": list(languages),
        }

    getter = getter or web_reported._httpx_getter(  # noqa: SLF001
        pause_seconds=0.05,
        max_retries=4,
    )
    names = batch["accepted_species"].astype(str).tolist()
    species_to_qid, qid_completed, errors = fetch_english_wikidata_ids(getter, names)
    sitelinks, sitelink_errors = fetch_wikidata_sitelinks(
        getter,
        species_to_qid,
        languages,
    )
    errors.extend(sitelink_errors)

    source_tables: list[pd.DataFrame] = []
    completed_species: set[str] = set(qid_completed)
    sitelink_source_species: set[str] = set()
    direct_source_species: set[str] = set()
    source_counts_by_language: dict[str, int] = {}
    for language in languages:
        resolved_titles = {
            species: local_titles[language]
            for species, local_titles in sitelinks.items()
            if language in local_titles
        }
        resolved_sources, resolved_completed, resolved_errors = fetch_resolved_language_extracts(
            getter,
            resolved_titles,
            language,
        )
        errors.extend(resolved_errors)
        completed_species.update(resolved_completed)
        if not resolved_sources.empty:
            sitelink_source_species.update(resolved_sources["accepted_species"].astype(str))
        source_tables.append(resolved_sources)

        resolved_with_text = (
            set(resolved_sources["accepted_species"].astype(str))
            if not resolved_sources.empty
            else set()
        )
        fallback_names = [name for name in names if name not in resolved_with_text]
        direct_sources, direct_completed, direct_errors = fetch_language_extracts(
            getter,
            fallback_names,
            language,
        )
        errors.extend(direct_errors)
        completed_species.update(direct_completed)
        if not direct_sources.empty:
            direct_source_species.update(direct_sources["accepted_species"].astype(str))
        source_tables.append(direct_sources)
        combined_language = pd.concat(
            [resolved_sources, direct_sources],
            ignore_index=True,
        )
        source_counts_by_language[language] = int(
            combined_language["accepted_species"].nunique()
        )

    sources = pd.concat(source_tables, ignore_index=True, sort=False).fillna("")
    if not sources.empty:
        sources = sources.drop_duplicates()
    candidates = _drop_nested_guild_matches(extract_recovery_candidates(sources))
    if not candidates.empty:
        candidates["campaign_task"] = RECOVERY_TASK
        candidates["campaign_phase"] = RECOVERY_PHASE

    wave_id = _next_wave_id(campaign_dir, WAVE_TASK)
    ledger, biotic = _apply_results(
        ledger,
        names,
        candidates,
        completed_species,
        errors,
        wave_id,
    )
    wave_dir = campaign_dir / "waves" / wave_id
    wave_dir.mkdir(parents=True, exist_ok=False)
    batch.to_csv(wave_dir / "species_batch.csv", index=False)
    _write_frame_gzip(candidates, wave_dir / "machine_candidates.csv.gz")
    pd.DataFrame(columns=["accepted_species", "holdout_reason"]).to_csv(
        wave_dir / "holdouts.csv.gz",
        index=False,
        compression="gzip",
    )
    pd.DataFrame({"error": errors}).to_csv(wave_dir / "lookup_errors.csv", index=False)

    summary = {
        "wave_id": wave_id,
        "task": WAVE_TASK,
        "candidate_task": RECOVERY_TASK,
        "phase": RECOVERY_PHASE,
        "source_kind": "wikidata_sitelink_multilingual_wikipedia",
        "languages": list(languages),
        "n_species_attempted": len(names),
        "n_species_with_english_qid": len(species_to_qid),
        "n_species_with_any_target_sitelink": sum(bool(values) for values in sitelinks.values()),
        "n_species_with_source_text": int(sources["accepted_species"].nunique())
        if not sources.empty
        else 0,
        "n_species_with_sitelink_source_text": len(sitelink_source_species),
        "n_species_with_direct_fallback_text": len(direct_source_species),
        "source_species_by_language": source_counts_by_language,
        "n_candidates_total": len(candidates),
        "n_species_with_candidates": int(candidates["accepted_species"].nunique())
        if not candidates.empty
        else 0,
        "n_new_biotic_species": len(biotic),
        "n_lookup_errors": len(errors),
        "policy": (
            "Machine-only source-backed recovery; missing QIDs, missing sitelinks, "
            "and zero hits are not biological absences."
        ),
    }
    (wave_dir / "wave_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    _write_frame_gzip(ledger, ledger_path)

    status_path = campaign_dir / "campaign_status.json"
    if status_path.exists():
        status = json.loads(status_path.read_text(encoding="utf-8"))
        status["n_machine_biotic_candidates"] = int(
            ledger["machine_biotic_candidate"]
            .astype("string")
            .fillna("")
            .str.lower()
            .isin({"true", "1", "yes", "y"})
            .sum()
        )
        status["last_multilingual_sitelink_recovery_wave"] = summary
        status_path.write_text(
            json.dumps(status, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
    return summary


@app.command("run")
def run_command(
    campaign_dir: Path = typer.Option(..., exists=True),
    batch_size: int = typer.Option(1000, min=1, max=2000),
    languages: str = typer.Option("ja,zh,ru"),
) -> None:
    selected = tuple(
        value.strip().lower()
        for value in languages.split(",")
        if value.strip()
    )
    if not selected:
        raise typer.BadParameter("at least one Wikipedia language is required")
    report = run_sitelink_recovery_wave(
        campaign_dir=campaign_dir,
        batch_size=batch_size,
        languages=selected,
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
