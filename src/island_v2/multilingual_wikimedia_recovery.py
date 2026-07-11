"""Recover explicit reproductive and pollinator statements from Asian/Russian Wikipedia.

This is a supplemental machine-only lane for species whose primary Wikimedia
screen is already terminal. It batches up to 50 scientific names per Wikipedia
API request, preserves direct source excerpts, and writes a standard campaign
wave so the existing cumulative candidate index can ingest the results.

A zero hit means only that the selected language pages did not contain one of the
explicit patterns. It is never promoted to a biological absence or curated trait.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import pandas as pd
import typer

from island_v2 import reported_ecology_machine as ecology
from island_v2 import trait_web_reported_scout as web_reported
from island_v2.global_trait_campaign import (
    COMMON_CANDIDATE_COLUMNS,
    _next_wave_id,
    _write_frame_gzip,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

JsonGetter = Callable[[str, dict[str, Any]], dict[str, Any]]
RECOVERY_TASK = "reproductive_multilingual_recovery"
RECOVERY_PHASE = "reproductive_and_pollen_vector"
RECOVERY_COLUMNS = [
    f"{RECOVERY_TASK}_status",
    f"{RECOVERY_TASK}_attempts",
    f"{RECOVERY_TASK}_wave_id",
    f"{RECOVERY_TASK}_candidate_count",
    f"{RECOVERY_TASK}_last_error",
]
TERMINAL_PRIMARY = {"processed", "exhausted"}
DEFAULT_LANGUAGES = ("ja", "zh", "ru")


@dataclass(frozen=True)
class RecoveryPattern:
    trait_name: str
    candidate_value: str
    language: str
    expression: str


# CJK pollinator-guild expressions require a nearby pollination word. This is
# deliberately stricter than matching an insect name anywhere in an article.
JA_POLLINATION = r"(?:受粉|送粉|花粉媒介|授粉)"
ZH_POLLINATION = r"(?:授粉|传粉|傳粉|花粉媒介)"
RU_POLLINATION = r"(?:опыля\w*|опылен\w*|опылён\w*)"

PATTERNS: tuple[RecoveryPattern, ...] = (
    # Japanese: reproductive system and vector mode.
    RecoveryPattern("self_incompatibility", "SI", "ja", r"自家不和合(?:性)?"),
    RecoveryPattern("self_incompatibility", "SC", "ja", r"自家和合(?:性)?"),
    RecoveryPattern("autonomous_selfing_capacity", "autonomous", "ja", r"(?:自動)?自家(?:受粉|授粉)"),
    RecoveryPattern("pollen_vector_mode", "abiotic_wind", "ja", r"(?:風媒花|風媒植物|風媒性)"),
    RecoveryPattern("pollen_vector_mode", "biotic", "ja", r"(?:虫媒花|虫媒植物|昆虫媒介)"),
    RecoveryPattern("sex_system", "dioecious", "ja", r"雌雄異株"),
    RecoveryPattern("sex_system", "monoecious", "ja", r"雌雄同株"),
    RecoveryPattern("sex_system", "hermaphroditic", "ja", r"両性花"),
    RecoveryPattern(
        "pollination_functional_guild",
        "bumblebees",
        "ja",
        rf"(?:(?:マルハナバチ).{{0,24}}{JA_POLLINATION}|{JA_POLLINATION}.{{0,24}}(?:マルハナバチ))",
    ),
    RecoveryPattern(
        "pollination_functional_guild",
        "other_bees",
        "ja",
        rf"(?:(?:ハナバチ|ミツバチ).{{0,24}}{JA_POLLINATION}|{JA_POLLINATION}.{{0,24}}(?:ハナバチ|ミツバチ))",
    ),
    RecoveryPattern(
        "pollination_functional_guild",
        "flies",
        "ja",
        rf"(?:(?:ハエ|アブ).{{0,24}}{JA_POLLINATION}|{JA_POLLINATION}.{{0,24}}(?:ハエ|アブ))",
    ),
    RecoveryPattern(
        "pollination_functional_guild",
        "birds",
        "ja",
        rf"(?:(?:鳥|鳥類).{{0,24}}{JA_POLLINATION}|{JA_POLLINATION}.{{0,24}}(?:鳥|鳥類)|鳥媒花)",
    ),
    RecoveryPattern(
        "pollination_functional_guild",
        "bats",
        "ja",
        rf"(?:(?:コウモリ).{{0,24}}{JA_POLLINATION}|{JA_POLLINATION}.{{0,24}}(?:コウモリ)|コウモリ媒)",
    ),
    RecoveryPattern(
        "pollination_functional_guild",
        "butterflies",
        "ja",
        rf"(?:(?:チョウ|蝶).{{0,24}}{JA_POLLINATION}|{JA_POLLINATION}.{{0,24}}(?:チョウ|蝶)|蝶媒花)",
    ),
    RecoveryPattern(
        "pollination_functional_guild",
        "moths",
        "ja",
        rf"(?:(?:ガ|蛾).{{0,24}}{JA_POLLINATION}|{JA_POLLINATION}.{{0,24}}(?:ガ|蛾)|蛾媒花)",
    ),
    # Simplified and traditional Chinese.
    RecoveryPattern("self_incompatibility", "SI", "zh", r"自交不(?:亲|親)和(?:性)?"),
    RecoveryPattern("self_incompatibility", "SC", "zh", r"自交(?:亲|親)和(?:性)?"),
    RecoveryPattern("autonomous_selfing_capacity", "autonomous", "zh", r"(?:自动|自動)?自花(?:授粉|传粉|傳粉)"),
    RecoveryPattern("pollen_vector_mode", "abiotic_wind", "zh", r"(?:风媒花|風媒花|风媒植物|風媒植物)"),
    RecoveryPattern("pollen_vector_mode", "biotic", "zh", r"(?:虫媒花|蟲媒花|昆虫媒介|昆蟲媒介)"),
    RecoveryPattern("sex_system", "dioecious", "zh", r"雌雄(?:异|異)株"),
    RecoveryPattern("sex_system", "monoecious", "zh", r"雌雄同株"),
    RecoveryPattern("sex_system", "hermaphroditic", "zh", r"(?:两|兩)性花"),
    RecoveryPattern(
        "pollination_functional_guild",
        "bumblebees",
        "zh",
        rf"(?:(?:熊蜂).{{0,24}}{ZH_POLLINATION}|{ZH_POLLINATION}.{{0,24}}(?:熊蜂))",
    ),
    RecoveryPattern(
        "pollination_functional_guild",
        "other_bees",
        "zh",
        rf"(?:(?:蜜蜂|蜂类|蜂類).{{0,24}}{ZH_POLLINATION}|{ZH_POLLINATION}.{{0,24}}(?:蜜蜂|蜂类|蜂類))",
    ),
    RecoveryPattern(
        "pollination_functional_guild",
        "flies",
        "zh",
        rf"(?:(?:蝇|蠅|食蚜蝇|食蚜蠅).{{0,24}}{ZH_POLLINATION}|{ZH_POLLINATION}.{{0,24}}(?:蝇|蠅|食蚜蝇|食蚜蠅))",
    ),
    RecoveryPattern(
        "pollination_functional_guild",
        "birds",
        "zh",
        rf"(?:(?:鸟类|鳥類|蜂鸟|蜂鳥).{{0,24}}{ZH_POLLINATION}|{ZH_POLLINATION}.{{0,24}}(?:鸟类|鳥類|蜂鸟|蜂鳥)|(?:鸟媒花|鳥媒花))",
    ),
    RecoveryPattern(
        "pollination_functional_guild",
        "bats",
        "zh",
        rf"(?:(?:蝙蝠).{{0,24}}{ZH_POLLINATION}|{ZH_POLLINATION}.{{0,24}}(?:蝙蝠)|蝙蝠媒)",
    ),
    RecoveryPattern(
        "pollination_functional_guild",
        "butterflies",
        "zh",
        rf"(?:(?:蝴蝶).{{0,24}}{ZH_POLLINATION}|{ZH_POLLINATION}.{{0,24}}(?:蝴蝶)|蝶媒花)",
    ),
    RecoveryPattern(
        "pollination_functional_guild",
        "moths",
        "zh",
        rf"(?:(?:蛾类|蛾類|飞蛾|飛蛾).{{0,24}}{ZH_POLLINATION}|{ZH_POLLINATION}.{{0,24}}(?:蛾类|蛾類|飞蛾|飛蛾)|蛾媒花)",
    ),
    # Russian: morphological stems cover ordinary inflection.
    RecoveryPattern("self_incompatibility", "SI", "ru", r"самонесовместим\w*"),
    RecoveryPattern("self_incompatibility", "SC", "ru", r"самосовместим\w*"),
    RecoveryPattern("autonomous_selfing_capacity", "autonomous", "ru", r"самоопыл\w*"),
    RecoveryPattern("pollen_vector_mode", "abiotic_wind", "ru", r"(?:ветроопыляем\w*|анемофил\w*)"),
    RecoveryPattern("pollen_vector_mode", "biotic", "ru", r"(?:насекомоопыляем\w*|энтомофил\w*)"),
    RecoveryPattern("sex_system", "dioecious", "ru", r"двудомн\w*"),
    RecoveryPattern("sex_system", "monoecious", "ru", r"однодомн\w*"),
    RecoveryPattern("sex_system", "hermaphroditic", "ru", r"обоепол\w*"),
    RecoveryPattern(
        "pollination_functional_guild",
        "bumblebees",
        "ru",
        rf"(?:шмел\w*.{{0,35}}{RU_POLLINATION}|{RU_POLLINATION}.{{0,35}}шмел\w*)",
    ),
    RecoveryPattern(
        "pollination_functional_guild",
        "other_bees",
        "ru",
        rf"(?:пчел\w*.{{0,35}}{RU_POLLINATION}|{RU_POLLINATION}.{{0,35}}пчел\w*)",
    ),
    RecoveryPattern(
        "pollination_functional_guild",
        "flies",
        "ru",
        rf"(?:мух\w*.{{0,35}}{RU_POLLINATION}|{RU_POLLINATION}.{{0,35}}мух\w*)",
    ),
    RecoveryPattern(
        "pollination_functional_guild",
        "birds",
        "ru",
        rf"(?:птиц\w*.{{0,35}}{RU_POLLINATION}|{RU_POLLINATION}.{{0,35}}птиц\w*)",
    ),
    RecoveryPattern(
        "pollination_functional_guild",
        "bats",
        "ru",
        rf"(?:летуч\w* мыш\w*.{{0,35}}{RU_POLLINATION}|{RU_POLLINATION}.{{0,35}}летуч\w* мыш\w*)",
    ),
)


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _normal_title(value: str) -> str:
    return re.sub(r"\s+", " ", str(value or "").replace("_", " ").strip()).casefold()


def _excerpt(text: str, start: int, end: int, pad: int = 100) -> str:
    return re.sub(r"\s+", " ", text[max(0, start - pad): min(len(text), end + pad)]).strip()


def _locally_negated(text: str, start: int, language: str) -> bool:
    local = text[max(0, start - 24):start].casefold()
    if language == "ja":
        return bool(re.search(r"(?:ない|なく|非|無い|認められない)$", local))
    if language == "zh":
        return bool(re.search(r"(?:不|非|无|無|没有|沒有)$", local))
    if language == "ru":
        return bool(re.search(r"(?:\bне\b|\bнет\b|\bбез\b)\s*$", local))
    return False


def _follow_alias(title: str, aliases: dict[str, str]) -> str:
    current = _normal_title(title)
    seen: set[str] = set()
    while current in aliases and current not in seen:
        seen.add(current)
        current = aliases[current]
    return current


def fetch_language_extracts(
    getter: JsonGetter,
    species: list[str],
    language: str,
    *,
    api_batch_size: int = 50,
) -> tuple[pd.DataFrame, set[str], list[str]]:
    """Batch scientific-name page lookups for one Wikipedia language."""
    rows: list[dict[str, str]] = []
    completed: set[str] = set()
    errors: list[str] = []
    endpoint = web_reported.wikipedia_api_for_language(language)

    for start in range(0, len(species), api_batch_size):
        names = species[start:start + api_batch_size]
        try:
            payload = getter(
                endpoint,
                {
                    "action": "query",
                    "prop": "extracts",
                    "explaintext": 1,
                    "redirects": 1,
                    "titles": "|".join(names),
                    "format": "json",
                    "exlimit": "max",
                },
            )
        except Exception as exc:  # noqa: BLE001 - preserve per-name retry state
            errors.extend(f"wikipedia-{language}:{name}:transient:{exc}" for name in names)
            continue

        completed.update(names)
        query = payload.get("query") or {}
        aliases: dict[str, str] = {}
        for item in [*(query.get("normalized") or []), *(query.get("redirects") or [])]:
            source = _normal_title(item.get("from") or "")
            target = _normal_title(item.get("to") or "")
            if source and target:
                aliases[source] = target

        pages_by_title: dict[str, dict[str, Any]] = {}
        for page in (query.get("pages") or {}).values():
            if "missing" in page:
                continue
            pages_by_title[_normal_title(page.get("title") or "")] = page

        for name in names:
            resolved_key = _follow_alias(name, aliases)
            page = pages_by_title.get(resolved_key)
            if page is None:
                continue
            extract = _text(page.get("extract"))
            if not extract:
                continue
            resolved = _text(page.get("title")) or name
            rows.append(
                {
                    "accepted_species": name,
                    "language": language,
                    "source_text": extract,
                    "source_url": f"https://{language}.wikipedia.org/wiki/" + resolved.replace(" ", "_"),
                    "source_citation": f"{resolved} Wikipedia full extract ({language})",
                    "source_type": f"wikipedia_extract_{language}",
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


def extract_recovery_candidates(sources: pd.DataFrame) -> pd.DataFrame:
    """Extract only explicit multilingual terms into standard campaign rows."""
    ontology = ecology.load_ontology(Path("config/trait_ontology.yml"))
    rows: list[dict[str, Any]] = []
    for record in sources.fillna("").to_dict("records"):
        species = _text(record.get("accepted_species"))
        language = _text(record.get("language"))
        text = _text(record.get("source_text"))
        if not species or not text:
            continue
        seen: set[tuple[str, str]] = set()
        for rule in PATTERNS:
            if rule.language != language:
                continue
            for match in re.finditer(rule.expression, text, flags=re.IGNORECASE):
                if _locally_negated(text, match.start(), language):
                    continue
                value = ecology.normalize_value(rule.trait_name, rule.candidate_value, ontology)
                key = (rule.trait_name, value)
                if not value or key in seen:
                    continue
                seen.add(key)
                rows.append(
                    {
                        "accepted_species": species,
                        "trait_name": rule.trait_name,
                        "candidate_value": value,
                        "source_url": _text(record.get("source_url")),
                        "source_citation": _text(record.get("source_citation")),
                        "source_excerpt": _excerpt(text, match.start(), match.end()),
                        "source_type": _text(record.get("source_type")),
                        "evidence_scope": "species_direct",
                        "matched_term": match.group(0),
                        "confidence": "explicit_multilingual_rule_match",
                        "campaign_task": RECOVERY_TASK,
                        "campaign_phase": RECOVERY_PHASE,
                        "target_for_task": True,
                    }
                )

                # A named animal pollinator is also direct evidence for biotic
                # pollen transport, even when the local language lacks a generic
                # insect-pollinated phrase.
                if rule.trait_name == "pollination_functional_guild":
                    biotic_key = ("pollen_vector_mode", "biotic")
                    if biotic_key not in seen:
                        seen.add(biotic_key)
                        rows.append(
                            {
                                "accepted_species": species,
                                "trait_name": "pollen_vector_mode",
                                "candidate_value": "biotic",
                                "source_url": _text(record.get("source_url")),
                                "source_citation": _text(record.get("source_citation")),
                                "source_excerpt": _excerpt(text, match.start(), match.end()),
                                "source_type": _text(record.get("source_type")),
                                "evidence_scope": "species_direct",
                                "matched_term": match.group(0),
                                "confidence": "explicit_multilingual_rule_match",
                                "campaign_task": RECOVERY_TASK,
                                "campaign_phase": RECOVERY_PHASE,
                                "target_for_task": True,
                            }
                        )
                break
    if not rows:
        return pd.DataFrame(columns=COMMON_CANDIDATE_COLUMNS)
    return pd.DataFrame(rows, columns=COMMON_CANDIDATE_COLUMNS).drop_duplicates().reset_index(drop=True)


def _ensure_recovery_columns(ledger: pd.DataFrame) -> pd.DataFrame:
    result = ledger.copy()
    defaults: dict[str, object] = {
        RECOVERY_COLUMNS[0]: "pending",
        RECOVERY_COLUMNS[1]: 0,
        RECOVERY_COLUMNS[2]: "",
        RECOVERY_COLUMNS[3]: 0,
        RECOVERY_COLUMNS[4]: "",
    }
    for column, default in defaults.items():
        if column not in result.columns:
            result[column] = default
    result[RECOVERY_COLUMNS[0]] = result[RECOVERY_COLUMNS[0]].fillna("").replace("", "pending")
    for column in (RECOVERY_COLUMNS[1], RECOVERY_COLUMNS[3]):
        result[column] = pd.to_numeric(result[column], errors="coerce").fillna(0).astype(int)
    for column in (RECOVERY_COLUMNS[2], RECOVERY_COLUMNS[4]):
        result[column] = result[column].fillna("").astype(str)
    return result


def _family_balanced_recovery_batch(ledger: pd.DataFrame, batch_size: int) -> pd.DataFrame:
    eligible = (
        ledger["reproductive_wikimedia_status"].isin(TERMINAL_PRIMARY)
        & ledger[RECOVERY_COLUMNS[0]].isin({"pending", "retry"})
    )
    pending = ledger.loc[eligible].copy()
    if pending.empty:
        return pending
    pending = pending.sort_values(
        ["family", "n_islands", "n_records", "accepted_species"],
        ascending=[True, False, False, True],
    )
    groups = {str(f): g.reset_index(drop=True) for f, g in pending.groupby("family", sort=True)}
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


def run_recovery_wave(
    *,
    campaign_dir: Path,
    batch_size: int = 1000,
    languages: tuple[str, ...] = DEFAULT_LANGUAGES,
    getter: JsonGetter | None = None,
) -> dict[str, Any]:
    """Advance one multilingual recovery wave and return its summary."""
    ledger_path = campaign_dir / "campaign_ledger.csv.gz"
    if not ledger_path.exists():
        raise typer.BadParameter(f"campaign ledger not found: {ledger_path}")
    ledger = _ensure_recovery_columns(pd.read_csv(ledger_path, dtype=str).fillna(""))
    batch = _family_balanced_recovery_batch(ledger, batch_size)
    if batch.empty:
        return {
            "status": "complete_for_current_primary_progress",
            "n_species_attempted": 0,
            "n_candidates": 0,
            "languages": list(languages),
        }

    getter = getter or web_reported._httpx_getter(pause_seconds=0.05, max_retries=4)  # noqa: SLF001
    names = batch["accepted_species"].astype(str).tolist()
    source_tables: list[pd.DataFrame] = []
    completed_by_language: dict[str, set[str]] = {}
    errors: list[str] = []
    for language in languages:
        sources, completed, language_errors = fetch_language_extracts(
            getter,
            names,
            language,
        )
        source_tables.append(sources)
        completed_by_language[language] = completed
        errors.extend(language_errors)

    sources = pd.concat(source_tables, ignore_index=True, sort=False).fillna("")
    if not sources.empty:
        sources = sources.drop_duplicates()
    candidates = extract_recovery_candidates(sources)
    wave_id = _next_wave_id(campaign_dir, RECOVERY_TASK)
    candidate_counts = (
        candidates.groupby("accepted_species").size().to_dict() if not candidates.empty else {}
    )
    error_by_species: dict[str, list[str]] = {}
    for error in errors:
        parts = str(error).split(":", 2)
        if len(parts) >= 2:
            error_by_species.setdefault(parts[1], []).append(error)

    biotic_species = set(
        candidates.loc[
            candidates["trait_name"].eq("pollen_vector_mode")
            & candidates["candidate_value"].eq("biotic")
            & candidates["evidence_scope"].eq("species_direct"),
            "accepted_species",
        ].astype(str)
    )

    for species in names:
        mask = ledger["accepted_species"].eq(species)
        previous_attempts = int(ledger.loc[mask, RECOVERY_COLUMNS[1]].iloc[0])
        completed_languages = sum(species in values for values in completed_by_language.values())
        species_errors = error_by_species.get(species, [])
        ledger.loc[mask, RECOVERY_COLUMNS[1]] = previous_attempts + 1
        ledger.loc[mask, RECOVERY_COLUMNS[2]] = wave_id
        ledger.loc[mask, RECOVERY_COLUMNS[3]] = int(candidate_counts.get(species, 0))
        if completed_languages > 0:
            ledger.loc[mask, RECOVERY_COLUMNS[0]] = "processed"
            ledger.loc[mask, RECOVERY_COLUMNS[4]] = ""
        else:
            ledger.loc[mask, RECOVERY_COLUMNS[0]] = "retry"
            ledger.loc[mask, RECOVERY_COLUMNS[4]] = " | ".join(species_errors)[:1000]

        if species in biotic_species:
            ledger.loc[mask, "machine_biotic_candidate"] = True
            # The primary pass may already have closed these rows as
            # not_applicable. Re-open them now that direct biotic evidence exists.
            for dependent in (
                "floral_access_wikimedia_status",
                "alternative_guild_wikimedia_status",
                "alternative_guild_openalex_status",
            ):
                if dependent in ledger.columns:
                    ledger.loc[mask & ledger[dependent].eq("not_applicable"), dependent] = "pending"

    wave_dir = campaign_dir / "waves" / wave_id
    wave_dir.mkdir(parents=True, exist_ok=False)
    batch.to_csv(wave_dir / "species_batch.csv", index=False)
    _write_frame_gzip(candidates, wave_dir / "machine_candidates.csv.gz")
    pd.DataFrame(columns=["accepted_species", "holdout_reason"]).to_csv(
        wave_dir / "holdouts.csv.gz", index=False, compression="gzip"
    )
    pd.DataFrame({"error": errors}).to_csv(wave_dir / "lookup_errors.csv", index=False)
    summary = {
        "wave_id": wave_id,
        "task": RECOVERY_TASK,
        "phase": RECOVERY_PHASE,
        "source_kind": "batched_direct_title_multilingual_wikipedia",
        "languages": list(languages),
        "n_species_attempted": len(names),
        "n_species_with_source_text": int(sources["accepted_species"].nunique())
        if not sources.empty
        else 0,
        "n_candidates_total": len(candidates),
        "n_species_with_candidates": int(candidates["accepted_species"].nunique())
        if not candidates.empty
        else 0,
        "n_new_biotic_species": len(biotic_species),
        "n_lookup_errors": len(errors),
        "policy": "Machine-only source-backed recovery; zero hits are not biological absences.",
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
        status["last_multilingual_recovery_wave"] = summary
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
    summary = run_recovery_wave(
        campaign_dir=campaign_dir,
        batch_size=batch_size,
        languages=selected,
    )
    typer.echo(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    app()
