"""Measure whether taxonomic aliases unlock additional verified Wikimedia pages.

This is a discovery benchmark, not trait extraction. Accepted names and GBIF
aliases are resolved to English Wikipedia/Wikidata items. A page is retained
only when Wikidata P225 matches the accepted name or the alias canonical name.
The output quantifies genuinely new QIDs reachable only through aliases.
"""

from __future__ import annotations

import json
import unicodedata
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2 import trait_web_reported_scout as web_reported
from island_v2.multilingual_wikimedia_sitelink_recovery import (
    WIKIDATA_API,
    fetch_english_wikidata_ids,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)
JsonGetter = Callable[[str, dict[str, Any]], dict[str, Any]]

OUTPUT_COLUMNS = [
    "accepted_species",
    "discovery_name",
    "name_role",
    "qid",
    "p225_values",
    "verification_status",
    "is_alias_only_qid",
]


def _text(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return " ".join(str(value).strip().split())


def _taxon_key(value: object) -> str:
    text = unicodedata.normalize("NFKC", _text(value)).replace("×", "x")
    return " ".join(text.split()).casefold()


def _p225_values(entity: dict[str, Any]) -> list[str]:
    values: list[str] = []
    for claim in ((entity.get("claims") or {}).get("P225") or []):
        value = (((claim.get("mainsnak") or {}).get("datavalue") or {}).get("value"))
        if isinstance(value, str) and _text(value):
            values.append(_text(value))
    return values


def discover_alias_pages(
    aliases: pd.DataFrame,
    getter: JsonGetter,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    required = {"accepted_species", "search_name", "name_role"}
    missing = required.difference(aliases.columns)
    if missing:
        raise typer.BadParameter(f"alias table missing columns: {sorted(missing)}")

    records: list[dict[str, str]] = []
    for row in aliases.fillna("").to_dict("records"):
        accepted = _text(row.get("accepted_species"))
        search_name = _text(row.get("search_name"))
        canonical = _text(row.get("gbif_canonical_name"))
        role = _text(row.get("name_role"))
        if not accepted or not search_name:
            continue
        discovery_names = [search_name]
        if canonical and canonical.casefold() != search_name.casefold():
            discovery_names.append(canonical)
        for discovery_name in discovery_names:
            records.append(
                {
                    "accepted_species": accepted,
                    "discovery_name": discovery_name,
                    "name_role": role,
                    "verification_names": [accepted, search_name, canonical],
                }
            )

    discovery_names = list(dict.fromkeys(record["discovery_name"] for record in records))
    qids, _, qid_errors = fetch_english_wikidata_ids(getter, discovery_names)
    unique_qids = sorted(set(qids.values()))
    entities: dict[str, Any] = {}
    entity_errors: list[str] = []
    for start in range(0, len(unique_qids), 50):
        batch = unique_qids[start : start + 50]
        try:
            payload = getter(
                WIKIDATA_API,
                {
                    "action": "wbgetentities",
                    "ids": "|".join(batch),
                    "props": "claims",
                    "format": "json",
                },
            )
            entities.update(payload.get("entities") or {})
        except Exception as exc:  # noqa: BLE001
            entity_errors.append(f"wikidata:{'|'.join(batch)}:{exc}")

    rows: list[dict[str, str]] = []
    accepted_qids: dict[str, set[str]] = {}
    verified_intermediate: list[dict[str, str]] = []
    for record in records:
        qid = qids.get(record["discovery_name"], "")
        if not qid:
            continue
        p225 = _p225_values(entities.get(qid) or {})
        p225_keys = {_taxon_key(value) for value in p225}
        valid_keys = {
            _taxon_key(value)
            for value in record["verification_names"]
            if _text(value)
        }
        status = "p225_verified" if p225_keys & valid_keys else "p225_mismatch"
        result = {
            "accepted_species": record["accepted_species"],
            "discovery_name": record["discovery_name"],
            "name_role": record["name_role"],
            "qid": qid,
            "p225_values": "|".join(p225),
            "verification_status": status,
        }
        verified_intermediate.append(result)
        if status == "p225_verified" and record["name_role"] == "accepted_input":
            accepted_qids.setdefault(record["accepted_species"], set()).add(qid)

    for result in verified_intermediate:
        alias_only = (
            result["verification_status"] == "p225_verified"
            and result["name_role"] != "accepted_input"
            and result["qid"] not in accepted_qids.get(result["accepted_species"], set())
        )
        rows.append({**result, "is_alias_only_qid": "true" if alias_only else "false"})

    frame = pd.DataFrame(rows, columns=OUTPUT_COLUMNS).drop_duplicates()
    alias_only_species = set(
        frame.loc[frame["is_alias_only_qid"].eq("true"), "accepted_species"]
    )
    verified_species = set(
        frame.loc[frame["verification_status"].eq("p225_verified"), "accepted_species"]
    )
    report = {
        "n_species_input": int(aliases["accepted_species"].nunique()),
        "n_discovery_names": len(discovery_names),
        "n_names_with_english_qid": len(qids),
        "n_species_with_verified_page": len(verified_species),
        "n_species_with_alias_only_verified_qid": len(alias_only_species),
        "alias_only_species": sorted(alias_only_species),
        "n_qid_lookup_errors": len(qid_errors),
        "n_entity_lookup_errors": len(entity_errors),
        "policy": (
            "Alias-only QIDs are discovery opportunities, not evidence. Trait extraction must still "
            "use species-scoped source text and preserve the discovery alias."
        ),
    }
    return frame, report


def _httpx_getter() -> JsonGetter:
    return web_reported._httpx_getter(pause_seconds=0.02, max_retries=4)  # noqa: SLF001


@app.command()
def run(
    alias_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    aliases = pd.read_csv(alias_csv, dtype=str).fillna("")
    frame, report = discover_alias_pages(aliases, _httpx_getter())
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "alias_wikimedia_pages.csv", index=False)
    (output_dir / "alias_wikimedia_discovery_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
