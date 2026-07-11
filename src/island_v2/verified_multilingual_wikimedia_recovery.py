"""High-precision multilingual recovery with Wikidata taxon-name verification.

English page redirects are useful for discovery but can land on a genus or other
broader taxon. This lane accepts a Wikidata item only when its P225 scientific
name exactly matches the campaign species, then follows Japanese, Chinese, and
Russian sitelinks. Local extracts use 20-title batches, matching the practical
TextExtracts limit for non-bot clients.

All rows remain unreviewed machine candidates. Mismatches, missing claims, missing
pages, uncertainty statements, and zero hits are acquisition/audit outcomes, not
biological absences.
"""

from __future__ import annotations

import json
import re
import unicodedata
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
    WIKIDATA_API,
    _drop_nested_guild_matches,
    _text,
    fetch_english_wikidata_ids,
    fetch_resolved_language_extracts,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)
JsonGetter = Callable[[str, dict[str, Any]], dict[str, Any]]

WAVE_TASK = "reproductive_multilingual_verified_sitelink_recovery"
STATUS_COLUMNS = [
    f"{WAVE_TASK}_status",
    f"{WAVE_TASK}_attempts",
    f"{WAVE_TASK}_wave_id",
    f"{WAVE_TASK}_candidate_count",
    f"{WAVE_TASK}_last_error",
]
TERMINAL_PRIMARY = {"processed", "exhausted"}

UNCERTAINTY_PATTERNS = {
    "ja": re.compile(r"(?:不明|分かっていない|わかっていない|明らかでない)"),
    "zh": re.compile(r"(?:不明|未知|尚不清楚|尚不明确|尚不明確|未确定|未確定)"),
    "ru": re.compile(r"(?:неизвестн\w*|неясн\w*|не установл\w*|не удалось выяснить)", re.I),
}


def _taxon_key(value: object) -> str:
    text = unicodedata.normalize("NFKC", _text(value)).replace("×", "x")
    return " ".join(text.split()).casefold()


def _p225_values(entity: dict[str, Any]) -> list[str]:
    values: list[str] = []
    for claim in ((entity.get("claims") or {}).get("P225") or []):
        snak = claim.get("mainsnak") or {}
        value = (snak.get("datavalue") or {}).get("value")
        if isinstance(value, str) and _text(value):
            values.append(_text(value))
    return values


def fetch_verified_sitelinks(
    getter: JsonGetter,
    species_to_qid: dict[str, str],
    languages: tuple[str, ...],
    *,
    api_batch_size: int = 50,
) -> tuple[dict[str, dict[str, str]], pd.DataFrame, list[str]]:
    """Resolve local sitelinks only for Wikidata items whose P225 matches."""
    qid_to_species: dict[str, list[str]] = {}
    for species, qid in species_to_qid.items():
        qid_to_species.setdefault(qid, []).append(species)

    resolved: dict[str, dict[str, str]] = {}
    audit_rows: list[dict[str, str]] = []
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
                    "props": "claims|sitelinks",
                    "sitefilter": sitefilter,
                    "format": "json",
                },
            )
        except Exception as exc:  # noqa: BLE001 - preserve retry audit
            for qid in batch:
                for species in qid_to_species.get(qid, []):
                    errors.append(f"wikidata-verified:{species}:transient:{exc}")
            continue

        entities = payload.get("entities") or {}
        for qid in batch:
            entity = entities.get(qid) or {}
            p225_values = _p225_values(entity)
            p225_keys = {_taxon_key(value) for value in p225_values}
            sitelinks = entity.get("sitelinks") or {}
            for species in qid_to_species.get(qid, []):
                species_key = _taxon_key(species)
                if not p225_values:
                    status = "missing_p225"
                elif species_key not in p225_keys:
                    status = "p225_mismatch"
                else:
                    status = "p225_verified"
                    for language in languages:
                        title = _text((sitelinks.get(f"{language}wiki") or {}).get("title"))
                        if title:
                            resolved.setdefault(species, {})[language] = title
                audit_rows.append(
                    {
                        "accepted_species": species,
                        "qid": qid,
                        "p225_values": "|".join(p225_values),
                        "verification_status": status,
                        "n_target_sitelinks": str(
                            sum(bool((sitelinks.get(f"{language}wiki") or {}).get("title")) for language in languages)
                            if status == "p225_verified"
                            else 0
                        ),
                    }
                )

    columns = [
        "accepted_species",
        "qid",
        "p225_values",
        "verification_status",
        "n_target_sitelinks",
    ]
    return resolved, pd.DataFrame(audit_rows, columns=columns), errors


def filter_uncertain_candidates(candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Remove candidates whose local excerpt explicitly marks the state unknown."""
    if candidates.empty:
        return candidates.copy(), pd.DataFrame(columns=[*candidates.columns, "audit_reason"])
    frame = candidates.copy()
    uncertain = pd.Series(False, index=frame.index)
    for language, pattern in UNCERTAINTY_PATTERNS.items():
        source_language = frame["source_type"].astype(str).str.endswith(f"_{language}")
        uncertain |= source_language & frame["source_excerpt"].astype(str).str.contains(
            pattern,
            na=False,
        )
    rejected = frame.loc[uncertain].copy()
    rejected["audit_reason"] = "explicit_uncertainty_context"
    retained = frame.loc[~uncertain].copy()
    return retained.reset_index(drop=True), rejected.reset_index(drop=True)


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
    rows = candidates.loc[
        candidates["trait_name"].eq("pollen_vector_mode")
        & candidates["candidate_value"].eq("biotic")
        & candidates["evidence_scope"].eq("species_direct")
    ]
    return set(rows["accepted_species"].astype(str))


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
        result.loc[mask, STATUS_COLUMNS[1]] = int(result.loc[mask, STATUS_COLUMNS[1]].iloc[0]) + 1
        result.loc[mask, STATUS_COLUMNS[2]] = wave_id
        result.loc[mask, STATUS_COLUMNS[3]] = int(counts.get(species, 0))
        result.loc[mask, STATUS_COLUMNS[0]] = "processed" if species in completed_species else "retry"
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


def run_verified_recovery_wave(
    *,
    campaign_dir: Path,
    batch_size: int = 1000,
    languages: tuple[str, ...] = DEFAULT_LANGUAGES,
    getter: JsonGetter | None = None,
) -> dict[str, Any]:
    """Advance one P225-verified multilingual recovery wave."""
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

    source_tables: list[pd.DataFrame] = []
    source_counts_by_language: dict[str, int] = {}
    completed_species: set[str] = set(qid_completed)
    for language in languages:
        titles = {
            species: local_titles[language]
            for species, local_titles in sitelinks.items()
            if language in local_titles
        }
        sources, completed, language_errors = fetch_resolved_language_extracts(
            getter,
            titles,
            language,
            api_batch_size=20,
        )
        source_tables.append(sources)
        completed_species.update(completed)
        errors.extend(language_errors)
        source_counts_by_language[language] = int(sources["accepted_species"].nunique()) if not sources.empty else 0

    sources = pd.concat(source_tables, ignore_index=True, sort=False).fillna("")
    if not sources.empty:
        sources = sources.drop_duplicates()
    candidates = _drop_nested_guild_matches(extract_recovery_candidates(sources))
    candidates, uncertain_rejections = filter_uncertain_candidates(candidates)
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
    _write_frame_gzip(identity_audit, wave_dir / "wikidata_identity_audit.csv.gz")
    _write_frame_gzip(uncertain_rejections, wave_dir / "uncertain_candidate_rejections.csv.gz")
    pd.DataFrame(columns=["accepted_species", "holdout_reason"]).to_csv(
        wave_dir / "holdouts.csv.gz",
        index=False,
        compression="gzip",
    )
    pd.DataFrame({"error": errors}).to_csv(wave_dir / "lookup_errors.csv", index=False)

    verification_counts = (
        identity_audit["verification_status"].value_counts().sort_index().to_dict()
        if not identity_audit.empty
        else {}
    )
    summary = {
        "wave_id": wave_id,
        "task": WAVE_TASK,
        "candidate_task": RECOVERY_TASK,
        "phase": RECOVERY_PHASE,
        "source_kind": "p225_verified_wikidata_sitelink_wikipedia",
        "languages": list(languages),
        "n_species_attempted": len(names),
        "n_species_with_english_qid": len(species_to_qid),
        "wikidata_identity_status_counts": {
            str(key): int(value) for key, value in verification_counts.items()
        },
        "n_species_with_verified_target_sitelink": sum(bool(values) for values in sitelinks.values()),
        "n_species_with_source_text": int(sources["accepted_species"].nunique()) if not sources.empty else 0,
        "source_species_by_language": source_counts_by_language,
        "n_candidates_total": len(candidates),
        "n_species_with_candidates": int(candidates["accepted_species"].nunique()) if not candidates.empty else 0,
        "n_uncertain_candidates_rejected": len(uncertain_rejections),
        "n_new_biotic_species": len(biotic),
        "n_lookup_errors": len(errors),
        "policy": (
            "Only P225-matched Wikidata taxon items are followed. Missing/mismatched "
            "identity and explicit uncertainty are audit outcomes, not absences."
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
        status["last_verified_multilingual_recovery_wave"] = summary
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
    report = run_verified_recovery_wave(
        campaign_dir=campaign_dir,
        batch_size=batch_size,
        languages=selected,
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
