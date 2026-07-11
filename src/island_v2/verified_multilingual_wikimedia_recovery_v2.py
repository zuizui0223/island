"""P225-verified multilingual recovery with correct TextExtracts pagination.

MediaWiki TextExtracts returns multiple pages only when ``exintro`` is true.
This version therefore retrieves up to 20 lead extracts per request, extracts
conservative candidates, then fetches one full extract at a time only for
verified species that still have no lead candidate.

All rows remain unreviewed source-backed machine candidates. Missing pages,
identity mismatches, uncertainty, and zero hits are not biological absences.
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
)
from island_v2.multilingual_wikimedia_sitelink_recovery import (
    _drop_nested_guild_matches,
    _follow_alias,
    _query_aliases,
    _text,
    fetch_english_wikidata_ids,
    fetch_resolved_language_extracts,
)
from island_v2.verified_multilingual_wikimedia_recovery import (
    fetch_verified_sitelinks,
    filter_uncertain_candidates,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)
JsonGetter = Callable[[str, dict[str, Any]], dict[str, Any]]

WAVE_TASK = "reproductive_multilingual_verified_lead_full_v2"
STATUS_COLUMNS = [
    f"{WAVE_TASK}_status",
    f"{WAVE_TASK}_attempts",
    f"{WAVE_TASK}_wave_id",
    f"{WAVE_TASK}_candidate_count",
    f"{WAVE_TASK}_last_error",
]
TERMINAL_PRIMARY = {"processed", "exhausted"}

SOURCE_COLUMNS = [
    "accepted_species",
    "language",
    "source_text",
    "source_url",
    "source_citation",
    "source_type",
    "evidence_scope",
]


def fetch_lead_language_extracts(
    getter: JsonGetter,
    species_to_title: dict[str, str],
    language: str,
    *,
    api_batch_size: int = 20,
) -> tuple[pd.DataFrame, set[str], list[str]]:
    """Fetch up to 20 lead extracts per request, as required by TextExtracts."""
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
                    "exintro": 1,
                    "exlimit": 20,
                    "redirects": 1,
                    "titles": "|".join(titles),
                    "format": "json",
                },
            )
        except Exception as exc:  # noqa: BLE001 - preserve retry audit
            errors.extend(
                f"wikipedia-lead-{language}:{species}:transient:{exc}"
                for species, _ in batch
            )
            continue

        completed.update(species for species, _ in batch)
        query = payload.get("query") or {}
        aliases = _query_aliases(query)
        pages = {
            " ".join(str(page.get("title") or "").replace("_", " ").split()).casefold(): page
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
                        f"{resolved_title} Wikipedia lead extract "
                        f"({language}; P225-verified Wikidata sitelink)"
                    ),
                    "source_type": f"wikipedia_verified_lead_{language}",
                    "evidence_scope": "species_direct",
                }
            )
    return pd.DataFrame(rows, columns=SOURCE_COLUMNS), completed, errors


def _extract_audited_candidates(
    sources: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = _drop_nested_guild_matches(extract_recovery_candidates(sources))
    retained, rejected = filter_uncertain_candidates(candidates)
    if not retained.empty:
        retained["campaign_task"] = RECOVERY_TASK
        retained["campaign_phase"] = RECOVERY_PHASE
    return retained, rejected


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
    error_by_species: dict[str, list[str]] = {}
    for error in errors:
        parts = str(error).split(":", 2)
        if len(parts) >= 2:
            error_by_species.setdefault(parts[1], []).append(error)
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
            "" if species in completed_species else " | ".join(error_by_species.get(species, []))[:1000]
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


def run_lead_full_recovery_wave(
    *,
    campaign_dir: Path,
    batch_size: int = 1000,
    languages: tuple[str, ...] = DEFAULT_LANGUAGES,
    getter: JsonGetter | None = None,
) -> dict[str, Any]:
    """Advance one verified lead-plus-full-fallback recovery wave."""
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
    sitelinks, identity_audit, identity_errors = fetch_verified_sitelinks(
        getter,
        species_to_qid,
        languages,
    )
    errors.extend(identity_errors)

    lead_tables: list[pd.DataFrame] = []
    lead_by_language: dict[str, int] = {}
    completed_species = set(qid_completed)
    title_maps: dict[str, dict[str, str]] = {}
    for language in languages:
        titles = {
            species: local_titles[language]
            for species, local_titles in sitelinks.items()
            if language in local_titles
        }
        title_maps[language] = titles
        sources, completed, language_errors = fetch_lead_language_extracts(
            getter,
            titles,
            language,
        )
        lead_tables.append(sources)
        completed_species.update(completed)
        errors.extend(language_errors)
        lead_by_language[language] = int(sources["accepted_species"].nunique()) if not sources.empty else 0

    lead_sources = pd.concat(lead_tables, ignore_index=True, sort=False).fillna("")
    if not lead_sources.empty:
        lead_sources = lead_sources.drop_duplicates()
    lead_candidates, lead_rejections = _extract_audited_candidates(lead_sources)
    lead_candidate_species = set(lead_candidates["accepted_species"].astype(str))

    full_tables: list[pd.DataFrame] = []
    full_by_language: dict[str, int] = {}
    full_fallback_species: set[str] = set()
    for language in languages:
        titles = {
            species: title
            for species, title in title_maps[language].items()
            if species not in lead_candidate_species
        }
        full_fallback_species.update(titles)
        sources, completed, language_errors = fetch_resolved_language_extracts(
            getter,
            titles,
            language,
            api_batch_size=1,
        )
        full_tables.append(sources)
        completed_species.update(completed)
        errors.extend(language_errors)
        full_by_language[language] = int(sources["accepted_species"].nunique()) if not sources.empty else 0

    full_sources = pd.concat(full_tables, ignore_index=True, sort=False).fillna("")
    if not full_sources.empty:
        full_sources = full_sources.drop_duplicates()
    full_candidates, full_rejections = _extract_audited_candidates(full_sources)
    candidates = pd.concat(
        [lead_candidates, full_candidates],
        ignore_index=True,
        sort=False,
    ).drop_duplicates()
    rejections = pd.concat(
        [lead_rejections, full_rejections],
        ignore_index=True,
        sort=False,
    ).drop_duplicates()

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
    _write_frame_gzip(identity_audit, wave_dir / "wikidata_identity_audit.csv.gz")
    _write_frame_gzip(rejections, wave_dir / "uncertain_candidate_rejections.csv.gz")
    pd.DataFrame({"error": errors}).to_csv(wave_dir / "lookup_errors.csv", index=False)
    pd.DataFrame(columns=["accepted_species", "holdout_reason"]).to_csv(
        wave_dir / "holdouts.csv.gz",
        index=False,
        compression="gzip",
    )

    identity_counts = (
        identity_audit["verification_status"].value_counts().sort_index().to_dict()
        if not identity_audit.empty
        else {}
    )
    summary = {
        "wave_id": wave_id,
        "task": WAVE_TASK,
        "candidate_task": RECOVERY_TASK,
        "phase": RECOVERY_PHASE,
        "source_kind": "p225_verified_lead_then_full_wikipedia",
        "languages": list(languages),
        "n_species_attempted": len(names),
        "n_species_with_english_qid": len(species_to_qid),
        "wikidata_identity_status_counts": {
            str(key): int(value) for key, value in identity_counts.items()
        },
        "n_species_with_verified_target_sitelink": sum(bool(values) for values in sitelinks.values()),
        "n_species_with_lead_text": int(lead_sources["accepted_species"].nunique()) if not lead_sources.empty else 0,
        "lead_source_species_by_language": lead_by_language,
        "n_species_with_lead_candidates": len(lead_candidate_species),
        "n_species_full_fallback_attempted": len(full_fallback_species),
        "n_species_with_full_text": int(full_sources["accepted_species"].nunique()) if not full_sources.empty else 0,
        "full_source_species_by_language": full_by_language,
        "n_candidates_total": len(candidates),
        "n_species_with_candidates": int(candidates["accepted_species"].nunique()) if not candidates.empty else 0,
        "n_uncertain_candidates_rejected": len(rejections),
        "n_new_biotic_species": len(biotic),
        "n_lookup_errors": len(errors),
        "policy": (
            "P225 identity is required. Batched lead extraction follows the official "
            "TextExtracts multi-page contract; full text is requested individually "
            "only where the lead produced no candidate."
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
        status["last_verified_lead_full_recovery_wave"] = summary
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
    selected = tuple(value.strip().lower() for value in languages.split(",") if value.strip())
    if not selected:
        raise typer.BadParameter("at least one Wikipedia language is required")
    report = run_lead_full_recovery_wave(
        campaign_dir=campaign_dir,
        batch_size=batch_size,
        languages=selected,
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
