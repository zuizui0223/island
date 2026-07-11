"""Inventory EOL TraitBank predicates and relevant values without taxon joins.

The complete EOL traits table is large and ordered by broad trait groups. A
head-only smoke test therefore cannot tell whether pollination or floral terms
exist later in the file. This module performs a bounded-memory, vectorised full
scan before the expensive species-level extraction. It records every predicate
count and compact value inventories only for explicitly mapped or biologically
relevant predicates.
"""

from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.eol_traitbank import (
    LITERAL_COLUMNS,
    MEASUREMENT_COLUMNS,
    OBJECT_LABEL_COLUMNS,
    OBJECT_TERM_COLUMNS,
    PREDICATE_COLUMNS,
    PREDICATE_LABEL_COLUMNS,
    UNITS_COLUMNS,
    _column,
    _key,
    _normalise_text,
    _open_archive_member,
    _read_terms,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

DEFAULT_DISCOVERY_KEYWORDS = (
    "pollinat",
    "pollen vector",
    "flower color",
    "flower colour",
    "floral color",
    "floral colour",
    "corolla color",
    "corolla colour",
    "petal color",
    "petal colour",
    "sexual system",
    "sex system",
    "plant sex",
    "mating system",
    "breeding system",
    "reproductive system",
    "self incompat",
    "self compat",
    "self pollin",
    "autogam",
    "dichogam",
    "cleistogam",
    "floral symmetry",
    "flower symmetry",
    "flower form",
    "floral form",
    "flower shape",
    "floral shape",
    "flower size",
    "floral size",
    "flower length",
    "floral tube",
    "flower tube",
    "corolla tube",
    "nectar",
)

PREDICATE_COLUMNS_OUT = [
    "predicate",
    "predicate_label",
    "n_rows",
    "currently_mapped",
    "discovery_keyword_match",
]

VALUE_COLUMNS_OUT = [
    "predicate",
    "predicate_label",
    "value",
    "value_label",
    "n_rows",
]

TERM_COLUMNS_OUT = ["uri", "name", "term_type", "currently_mapped", "discovery_keyword_match"]


def _load_mapping_keys(config_path: Path) -> set[str]:
    payload = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    keys: set[str] = set()
    for rule in payload.get("rules") or []:
        for term in (rule or {}).get("predicate_terms") or []:
            if _normalise_text(term):
                keys.add(_key(term))
    return keys


def _keyword_pattern(keywords: tuple[str, ...]) -> re.Pattern[str]:
    cleaned = sorted({_normalise_text(item).casefold() for item in keywords if _normalise_text(item)})
    if not cleaned:
        raise typer.BadParameter("At least one discovery keyword is required.")
    return re.compile("|".join(re.escape(item) for item in cleaned), flags=re.IGNORECASE)


def _series(table: pd.DataFrame, column: str | None) -> pd.Series:
    if not column or column not in table.columns:
        return pd.Series([""] * len(table), index=table.index, dtype=str)
    return table[column].fillna("").astype(str).str.strip()


def _cap_counter(counter: Counter[tuple[str, ...]], *, keep: int = 100_000) -> Counter[tuple[str, ...]]:
    if len(counter) <= keep * 2:
        return counter
    return Counter(dict(counter.most_common(keep)))


def scan_eol_traitbank_inventory(
    *,
    eol_archive: Path,
    output_dir: Path,
    config_path: Path = Path("config/eol_traitbank_terms.yml"),
    chunksize: int = 500_000,
    max_trait_rows: int | None = None,
    max_values_per_predicate: int = 200,
    discovery_keywords: tuple[str, ...] = DEFAULT_DISCOVERY_KEYWORDS,
) -> dict[str, Any]:
    """Scan all EOL predicates and selected values with bounded memory."""
    if not eol_archive.exists():
        raise typer.BadParameter(f"EOL archive does not exist: {eol_archive}")
    if not config_path.exists():
        raise typer.BadParameter(f"EOL mapping config does not exist: {config_path}")
    if chunksize < 1:
        raise typer.BadParameter("chunksize must be positive.")
    if max_trait_rows is not None and max_trait_rows < 1:
        raise typer.BadParameter("max_trait_rows must be positive when supplied.")
    if max_values_per_predicate < 1:
        raise typer.BadParameter("max_values_per_predicate must be positive.")

    output_dir.mkdir(parents=True, exist_ok=True)
    terms = _read_terms(eol_archive)
    term_names = {uri: term.name for uri, term in terms.items()}
    mapped_keys = _load_mapping_keys(config_path)
    keyword_pattern = _keyword_pattern(discovery_keywords)

    predicate_counts: Counter[tuple[str, str]] = Counter()
    relevant_value_counts: Counter[tuple[str, str, str, str]] = Counter()
    n_scanned = 0

    with _open_archive_member(eol_archive, "traits.csv") as handle:
        reader = pd.read_csv(handle, dtype=str, chunksize=chunksize)
        for chunk in reader:
            if max_trait_rows is not None:
                remaining = max_trait_rows - n_scanned
                if remaining <= 0:
                    break
                chunk = chunk.head(remaining)
            if chunk.empty:
                continue

            columns = list(chunk.columns)
            predicate_col = _column(columns, PREDICATE_COLUMNS)
            predicate_label_col = _column(columns, PREDICATE_LABEL_COLUMNS)
            object_term_col = _column(columns, OBJECT_TERM_COLUMNS)
            object_label_col = _column(columns, OBJECT_LABEL_COLUMNS)
            literal_col = _column(columns, LITERAL_COLUMNS)
            measurement_col = _column(columns, MEASUREMENT_COLUMNS)
            units_col = _column(columns, UNITS_COLUMNS)
            if not predicate_col:
                raise typer.BadParameter("EOL traits.csv must contain a predicate column.")

            predicate = _series(chunk, predicate_col)
            explicit_predicate_label = _series(chunk, predicate_label_col)
            term_predicate_label = predicate.map(term_names).fillna("")
            predicate_label = explicit_predicate_label.where(
                explicit_predicate_label.ne(""), term_predicate_label
            )
            predicate_label = predicate_label.where(predicate_label.ne(""), predicate)

            grouped_predicates = pd.DataFrame(
                {"predicate": predicate, "predicate_label": predicate_label}
            ).groupby(["predicate", "predicate_label"], dropna=False).size()
            predicate_counts.update(
                {
                    (str(pred), str(label)): int(count)
                    for (pred, label), count in grouped_predicates.items()
                }
            )

            predicate_keys = predicate.str.casefold()
            label_keys = predicate_label.str.casefold()
            mapped_mask = predicate_keys.isin(mapped_keys) | label_keys.isin(mapped_keys)
            search_text = predicate_keys + " " + label_keys
            keyword_mask = search_text.str.contains(keyword_pattern, na=False)
            relevant_mask = mapped_mask | keyword_mask

            if relevant_mask.any():
                selected = chunk.loc[relevant_mask].copy()
                selected_predicate = predicate.loc[relevant_mask]
                selected_predicate_label = predicate_label.loc[relevant_mask]
                object_term = _series(selected, object_term_col)
                explicit_object_label = _series(selected, object_label_col)
                object_term_label = object_term.map(term_names).fillna("")
                object_label = explicit_object_label.where(
                    explicit_object_label.ne(""), object_term_label
                )
                literal = _series(selected, literal_col)
                measurement = _series(selected, measurement_col)
                units = _series(selected, units_col)
                measurement_with_units = (measurement + " " + units).str.strip()

                raw_value = object_term.where(object_term.ne(""), literal)
                raw_value = raw_value.where(raw_value.ne(""), measurement_with_units)
                value_label = object_label.where(object_label.ne(""), literal)
                value_label = value_label.where(value_label.ne(""), measurement_with_units)
                value_label = value_label.where(value_label.ne(""), raw_value)

                grouped_values = pd.DataFrame(
                    {
                        "predicate": selected_predicate.to_numpy(),
                        "predicate_label": selected_predicate_label.to_numpy(),
                        "value": raw_value.to_numpy(),
                        "value_label": value_label.to_numpy(),
                    }
                ).groupby(
                    ["predicate", "predicate_label", "value", "value_label"],
                    dropna=False,
                ).size()
                relevant_value_counts.update(
                    {
                        (str(pred), str(label), str(value), str(value_label)): int(count)
                        for (pred, label, value, value_label), count in grouped_values.items()
                    }
                )
                relevant_value_counts = _cap_counter(relevant_value_counts)

            n_scanned += len(chunk)
            if max_trait_rows is not None and n_scanned >= max_trait_rows:
                break

    predicate_rows: list[dict[str, Any]] = []
    observed_mapped = 0
    observed_relevant = 0
    for (predicate, label), count in predicate_counts.most_common():
        mapped = _key(predicate) in mapped_keys or _key(label) in mapped_keys
        keyword = bool(keyword_pattern.search(f"{predicate} {label}"))
        observed_mapped += int(mapped)
        observed_relevant += int(mapped or keyword)
        predicate_rows.append(
            {
                "predicate": predicate,
                "predicate_label": label,
                "n_rows": count,
                "currently_mapped": mapped,
                "discovery_keyword_match": keyword,
            }
        )

    by_predicate: dict[tuple[str, str], list[tuple[str, str, int]]] = defaultdict(list)
    for (predicate, label, value, value_label), count in relevant_value_counts.items():
        by_predicate[(predicate, label)].append((value, value_label, count))
    value_rows: list[dict[str, Any]] = []
    for (predicate, label), values in sorted(by_predicate.items()):
        for value, value_label, count in sorted(values, key=lambda item: item[2], reverse=True)[
            :max_values_per_predicate
        ]:
            value_rows.append(
                {
                    "predicate": predicate,
                    "predicate_label": label,
                    "value": value,
                    "value_label": value_label,
                    "n_rows": count,
                }
            )

    term_rows: list[dict[str, Any]] = []
    for term in terms.values():
        mapped = _key(term.uri) in mapped_keys or _key(term.name) in mapped_keys
        keyword = bool(keyword_pattern.search(f"{term.uri} {term.name}"))
        if mapped or keyword:
            term_rows.append(
                {
                    "uri": term.uri,
                    "name": term.name,
                    "term_type": term.term_type,
                    "currently_mapped": mapped,
                    "discovery_keyword_match": keyword,
                }
            )

    predicate_path = output_dir / "eol_predicate_inventory.csv"
    value_path = output_dir / "eol_relevant_value_inventory.csv"
    term_path = output_dir / "eol_relevant_terms.csv"
    manifest_path = output_dir / "eol_inventory_manifest.json"
    pd.DataFrame(predicate_rows, columns=PREDICATE_COLUMNS_OUT).to_csv(predicate_path, index=False)
    pd.DataFrame(value_rows, columns=VALUE_COLUMNS_OUT).to_csv(value_path, index=False)
    pd.DataFrame(term_rows, columns=TERM_COLUMNS_OUT).to_csv(term_path, index=False)

    report = {
        "version": "1.0",
        "eol_archive": str(eol_archive),
        "config_path": str(config_path),
        "n_trait_rows_scanned": int(n_scanned),
        "n_terms": len(terms),
        "n_observed_predicates": len(predicate_counts),
        "n_currently_mapped_predicates_observed": int(observed_mapped),
        "n_relevant_predicates_observed": int(observed_relevant),
        "n_relevant_value_rows_retained": len(value_rows),
        "max_values_per_predicate": int(max_values_per_predicate),
        "max_trait_rows": max_trait_rows,
        "predicate_inventory_csv": str(predicate_path),
        "relevant_value_inventory_csv": str(value_path),
        "relevant_terms_csv": str(term_path),
        "interpretation": (
            "Predicate and value frequencies are acquisition diagnostics. They are not "
            "species-level trait evidence or biological absence records."
        ),
    }
    manifest_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    report["manifest_json"] = str(manifest_path)
    return report


@app.command("scan")
def scan_command(
    eol_archive: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/eol_traitbank_terms.yml"), exists=True),
    chunksize: int = typer.Option(500_000, min=1),
    max_trait_rows: int | None = typer.Option(None, min=1),
    max_values_per_predicate: int = typer.Option(200, min=1),
) -> None:
    report = scan_eol_traitbank_inventory(
        eol_archive=eol_archive,
        output_dir=output_dir,
        config_path=config_path,
        chunksize=chunksize,
        max_trait_rows=max_trait_rows,
        max_values_per_predicate=max_values_per_predicate,
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
