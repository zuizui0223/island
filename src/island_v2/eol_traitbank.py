"""Prepare EOL TraitBank bulk exports for conservative v2 bulk ingest.

The EOL "All trait data" archive is a graph export, not the canonical v2 long
CSV used by :mod:`island_v2.bulk_trait_ingest`. This module streams the archive,
keeps only source terms that have an explicit reviewed mapping, and writes the
canonical long CSV plus an aggregate unmapped-term audit.
"""

from __future__ import annotations

import hashlib
import json
import re
import unicodedata
import zipfile
from collections import Counter
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO

import pandas as pd
import typer
import yaml

EOL_DATASET_URL = "https://zenodo.org/records/13305577"
EOL_PAGE_URL_PREFIX = "https://eol.org/pages"

LONG_COLUMNS = [
    "scientific_name",
    "trait_name",
    "trait_value",
    "source_record_id",
    "source_citation",
    "source_url",
    "evidence_excerpt",
]

UNMAPPED_COLUMNS = [
    "reason",
    "predicate",
    "predicate_label",
    "value",
    "value_label",
    "n_rows",
]

PAGE_ID_COLUMNS = ("page_id", "pageid", "eol_page_id", "eolid", "id")
PAGE_NAME_COLUMNS = ("canonical", "scientific_name", "scientificname", "name")
TERM_URI_COLUMNS = ("uri", "term_uri", "termuri", "id")
TERM_NAME_COLUMNS = ("name", "label", "term_name", "termname")
TERM_TYPE_COLUMNS = ("type", "term_type", "termtype")
TRAIT_ID_COLUMNS = ("eol_pk", "trait_id", "traitid", "id")
RESOURCE_PK_COLUMNS = ("resource_pk", "resourcepk", "measurementorfactid", "associationid")
TRAIT_PAGE_ID_COLUMNS = ("page_id", "pageid", "eol_page_id", "eolid")
TRAIT_NAME_COLUMNS = ("scientific_name", "scientificname", "canonical")
PREDICATE_COLUMNS = ("predicate", "predicate_uri", "predicateuri", "measurementtype")
PREDICATE_LABEL_COLUMNS = ("predicate_label", "predicate_name", "measurementtypelabel")
OBJECT_TERM_COLUMNS = ("object_term", "object_term_uri", "objectterm", "objecttermuri", "measurementvalueid")
OBJECT_LABEL_COLUMNS = ("object_term_label", "object_term_name", "object_name", "measurementvaluelabel")
LITERAL_COLUMNS = ("literal", "measurementvalue", "value")
MEASUREMENT_COLUMNS = ("normal_measurement", "normalmeasurement", "measurement")
UNITS_COLUMNS = ("normal_units", "normalunits", "units", "measurementunit")
CITATION_COLUMNS = ("citation", "bibliographiccitation", "source_citation")
SOURCE_COLUMNS = ("source", "source_url")


@dataclass(frozen=True)
class EolTerm:
    uri: str
    name: str
    term_type: str = ""


@dataclass(frozen=True)
class EolMappingRule:
    rule_id: str
    target_trait: str
    predicate_terms: tuple[str, ...]
    value_map: dict[str, str]


def _normalise_text(value: Any) -> str:
    text = unicodedata.normalize("NFKC", str(value or "")).strip()
    return re.sub(r"\s+", " ", text)


def _key(value: Any) -> str:
    return _normalise_text(value).casefold()


def _compact_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", _key(value))


def _column(columns: list[str], candidates: tuple[str, ...]) -> str | None:
    by_header = {_compact_header(column): column for column in columns}
    for candidate in candidates:
        match = by_header.get(_compact_header(candidate))
        if match:
            return match
    return None


def _sha256_file(path: Path) -> str:
    if path.is_dir():
        return ""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_yaml(path: Path) -> dict[str, Any]:
    parsed = yaml.safe_load(path.read_text(encoding="utf-8"))
    return parsed if isinstance(parsed, dict) else {}


def _load_ontology(path: Path) -> dict[str, Any]:
    ontology = _load_yaml(path)
    if not isinstance(ontology.get("traits"), dict):
        raise typer.BadParameter(f"Trait ontology must contain a traits mapping: {path}")
    return ontology


def _validate_rules(config: dict[str, Any], ontology: dict[str, Any]) -> list[EolMappingRule]:
    raw_rules = config.get("rules") or []
    if not isinstance(raw_rules, list) or not raw_rules:
        raise typer.BadParameter("EOL TraitBank mapping config must contain a non-empty rules list.")
    traits = ontology.get("traits") or {}
    rules: list[EolMappingRule] = []
    for index, raw in enumerate(raw_rules, start=1):
        if not isinstance(raw, dict):
            raise typer.BadParameter(f"EOL rule #{index} must be a mapping.")
        target_trait = _normalise_text(raw.get("target_trait") or raw.get("trait_name"))
        if target_trait not in traits:
            raise typer.BadParameter(f"EOL rule #{index} targets unknown v2 trait: {target_trait!r}")
        predicates = tuple(_normalise_text(item) for item in (raw.get("predicate_terms") or []) if _normalise_text(item))
        if not predicates:
            raise typer.BadParameter(f"EOL rule #{index} lacks predicate_terms.")
        value_map = {
            _key(source): _normalise_text(target)
            for source, target in (raw.get("value_map") or {}).items()
            if _normalise_text(source)
        }
        if not value_map:
            raise typer.BadParameter(f"EOL rule #{index} lacks value_map.")
        allowed = set((traits.get(target_trait) or {}).get("allowed_values") or [])
        invalid = sorted(set(value_map.values()).difference(allowed))
        if invalid:
            raise typer.BadParameter(
                f"EOL rule #{index} maps to invalid values for {target_trait}: {invalid}"
            )
        rules.append(
            EolMappingRule(
                rule_id=_normalise_text(raw.get("id")) or f"rule_{index}",
                target_trait=target_trait,
                predicate_terms=predicates,
                value_map=value_map,
            )
        )
    return rules


def _predicate_index(rules: list[EolMappingRule]) -> dict[str, list[EolMappingRule]]:
    index: dict[str, list[EolMappingRule]] = {}
    for rule in rules:
        for predicate in rule.predicate_terms:
            index.setdefault(_key(predicate), []).append(rule)
    return index


def _matching_rules(
    predicate: str,
    predicate_label: str,
    predicate_index: dict[str, list[EolMappingRule]],
) -> list[EolMappingRule]:
    selected: dict[str, EolMappingRule] = {}
    for term in (predicate, predicate_label):
        for rule in predicate_index.get(_key(term), []):
            selected[rule.rule_id] = rule
    return list(selected.values())


def _find_archive_member(path: Path, basename: str) -> str | Path:
    wanted = basename.casefold()
    if path.is_dir():
        matches = [
            item
            for item in path.rglob("*")
            if item.is_file() and item.name.casefold() == wanted
        ]
        if not matches:
            raise typer.BadParameter(f"EOL archive directory lacks {basename}: {path}")
        return sorted(matches)[0]
    with zipfile.ZipFile(path) as archive:
        matches = [
            info.filename
            for info in archive.infolist()
            if Path(info.filename).name.casefold() == wanted
        ]
    if not matches:
        raise typer.BadParameter(f"EOL zip archive lacks {basename}: {path}")
    return sorted(matches)[0]


@contextmanager
def _open_archive_member(path: Path, basename: str) -> Iterator[BinaryIO]:
    member = _find_archive_member(path, basename)
    if isinstance(member, Path):
        with member.open("rb") as handle:
            yield handle
        return
    with zipfile.ZipFile(path) as archive:
        with archive.open(member) as handle:
            yield handle


def _read_pages(path: Path) -> dict[str, str]:
    with _open_archive_member(path, "pages.csv") as handle:
        pages = pd.read_csv(handle, dtype=str).fillna("")
    columns = list(pages.columns)
    page_id_col = _column(columns, PAGE_ID_COLUMNS)
    name_col = _column(columns, PAGE_NAME_COLUMNS)
    if not page_id_col or not name_col:
        raise typer.BadParameter(
            "EOL pages.csv must contain page_id and canonical/scientific-name columns."
        )
    pages[page_id_col] = pages[page_id_col].map(_normalise_text)
    pages[name_col] = pages[name_col].map(_normalise_text)
    pages = pages.loc[pages[page_id_col].ne("") & pages[name_col].ne("")]
    return dict(zip(pages[page_id_col], pages[name_col], strict=False))


def _read_terms(path: Path) -> dict[str, EolTerm]:
    with _open_archive_member(path, "terms.csv") as handle:
        terms = pd.read_csv(handle, dtype=str).fillna("")
    columns = list(terms.columns)
    uri_col = _column(columns, TERM_URI_COLUMNS)
    name_col = _column(columns, TERM_NAME_COLUMNS)
    type_col = _column(columns, TERM_TYPE_COLUMNS)
    if not uri_col or not name_col:
        raise typer.BadParameter("EOL terms.csv must contain uri and name/label columns.")
    output: dict[str, EolTerm] = {}
    for row in terms[[uri_col, name_col] + ([type_col] if type_col else [])].itertuples(index=False):
        uri = _normalise_text(row[0])
        name = _normalise_text(row[1])
        term_type = _normalise_text(row[2]) if type_col and len(row) > 2 else ""
        if uri:
            output[uri] = EolTerm(uri=uri, name=name, term_type=term_type)
    return output


def _series_or_empty(table: pd.DataFrame, column: str | None) -> pd.Series:
    if column and column in table.columns:
        return table[column].fillna("").map(_normalise_text)
    return pd.Series([""] * len(table), index=table.index, dtype=str)


def _term_label(raw_term: str, explicit_label: str, terms: dict[str, EolTerm]) -> str:
    if explicit_label:
        return explicit_label
    term = terms.get(raw_term)
    if term and term.name:
        return term.name
    return raw_term


def _first_mapped_value(rule: EolMappingRule, candidates: tuple[str, ...]) -> str | None:
    for candidate in candidates:
        mapped = rule.value_map.get(_key(candidate))
        if mapped:
            return mapped
    return None


def _source_url(page_id: str, fallback: str) -> str:
    if page_id:
        return f"{EOL_PAGE_URL_PREFIX}/{page_id}/data"
    return fallback


def _record_id(eol_pk: str, resource_pk: str, row_number: int) -> str:
    if eol_pk:
        return eol_pk
    if resource_pk:
        return resource_pk
    return str(row_number)


def _write_counter(counter: Counter[tuple[str, str, str, str, str]], path: Path) -> None:
    rows = [
        {
            "reason": reason,
            "predicate": predicate,
            "predicate_label": predicate_label,
            "value": value,
            "value_label": value_label,
            "n_rows": count,
        }
        for (reason, predicate, predicate_label, value, value_label), count in counter.most_common()
    ]
    pd.DataFrame(rows, columns=UNMAPPED_COLUMNS).to_csv(path, index=False)


def prepare_eol_traitbank_long(
    *,
    eol_archive: Path,
    output_dir: Path,
    config_path: Path,
    ontology_path: Path = Path("config/trait_ontology.yml"),
    source_name: str = "eol_traitbank_all",
    chunksize: int = 200_000,
    max_trait_rows: int | None = None,
) -> dict[str, Any]:
    """Stream EOL TraitBank all-traits export into reviewed v2 long rows."""
    if not eol_archive.exists():
        raise typer.BadParameter(f"EOL archive does not exist: {eol_archive}")
    if chunksize < 1:
        raise typer.BadParameter("chunksize must be positive.")
    if max_trait_rows is not None and max_trait_rows < 1:
        raise typer.BadParameter("max_trait_rows must be positive when supplied.")

    output_dir.mkdir(parents=True, exist_ok=True)
    config = _load_yaml(config_path)
    ontology = _load_ontology(ontology_path)
    rules = _validate_rules(config, ontology)
    predicate_index = _predicate_index(rules)
    pages = _read_pages(eol_archive)
    terms = _read_terms(eol_archive)

    long_path = output_dir / "eol_traitbank_v2_long.csv"
    unmapped_path = output_dir / "eol_traitbank_unmapped_terms.csv"
    manifest_path = output_dir / "eol_traitbank_extract_manifest.json"
    source_citation = _normalise_text((config.get("dataset") or {}).get("source_citation")) or source_name
    dataset_url = _normalise_text((config.get("dataset") or {}).get("landing_page_url")) or EOL_DATASET_URL

    unmapped: Counter[tuple[str, str, str, str, str]] = Counter()
    n_scanned = 0
    n_candidates = 0
    n_predicate_matched = 0
    n_unmapped_predicate = 0
    n_unmapped_value = 0
    header_written = False

    with _open_archive_member(eol_archive, "traits.csv") as handle:
        reader = pd.read_csv(handle, dtype=str, chunksize=chunksize)
        for chunk in reader:
            chunk = chunk.fillna("")
            if max_trait_rows is not None:
                remaining = max_trait_rows - n_scanned
                if remaining <= 0:
                    break
                chunk = chunk.head(remaining)
            if chunk.empty:
                continue
            columns = list(chunk.columns)
            trait_id_col = _column(columns, TRAIT_ID_COLUMNS)
            resource_pk_col = _column(columns, RESOURCE_PK_COLUMNS)
            page_id_col = _column(columns, TRAIT_PAGE_ID_COLUMNS)
            scientific_name_col = _column(columns, TRAIT_NAME_COLUMNS)
            predicate_col = _column(columns, PREDICATE_COLUMNS)
            predicate_label_col = _column(columns, PREDICATE_LABEL_COLUMNS)
            object_term_col = _column(columns, OBJECT_TERM_COLUMNS)
            object_label_col = _column(columns, OBJECT_LABEL_COLUMNS)
            literal_col = _column(columns, LITERAL_COLUMNS)
            measurement_col = _column(columns, MEASUREMENT_COLUMNS)
            units_col = _column(columns, UNITS_COLUMNS)
            citation_col = _column(columns, CITATION_COLUMNS)
            source_col = _column(columns, SOURCE_COLUMNS)
            if not predicate_col:
                raise typer.BadParameter("EOL traits.csv must contain a predicate column.")
            if not page_id_col and not scientific_name_col:
                raise typer.BadParameter(
                    "EOL traits.csv must contain page_id or scientific_name for taxon matching."
                )

            work = pd.DataFrame(
                {
                    "eol_pk": _series_or_empty(chunk, trait_id_col),
                    "resource_pk": _series_or_empty(chunk, resource_pk_col),
                    "page_id": _series_or_empty(chunk, page_id_col),
                    "scientific_name": _series_or_empty(chunk, scientific_name_col),
                    "predicate": _series_or_empty(chunk, predicate_col),
                    "predicate_label": _series_or_empty(chunk, predicate_label_col),
                    "object_term": _series_or_empty(chunk, object_term_col),
                    "object_label": _series_or_empty(chunk, object_label_col),
                    "literal": _series_or_empty(chunk, literal_col),
                    "measurement": _series_or_empty(chunk, measurement_col),
                    "units": _series_or_empty(chunk, units_col),
                    "citation": _series_or_empty(chunk, citation_col),
                    "source": _series_or_empty(chunk, source_col),
                }
            )
            rows: list[dict[str, str]] = []
            for offset, row in enumerate(work.itertuples(index=False), start=1):
                row_number = n_scanned + offset
                page_name = pages.get(row.page_id, "")
                scientific_name = row.scientific_name or page_name
                predicate_label = _term_label(row.predicate, row.predicate_label, terms)
                object_label = _term_label(row.object_term, row.object_label, terms)
                value_candidates = (
                    row.object_term,
                    object_label,
                    row.literal,
                    row.measurement,
                    f"{row.measurement} {row.units}".strip(),
                )
                matched_rules = _matching_rules(row.predicate, predicate_label, predicate_index)
                if not matched_rules:
                    n_unmapped_predicate += 1
                    unmapped[
                        (
                            "unmapped_predicate",
                            row.predicate,
                            predicate_label,
                            row.object_term,
                            object_label or row.literal or row.measurement,
                        )
                    ] += 1
                    continue
                n_predicate_matched += 1
                for rule in matched_rules:
                    mapped_value = _first_mapped_value(rule, value_candidates)
                    if mapped_value is None:
                        n_unmapped_value += 1
                        unmapped[
                            (
                                "unmapped_value",
                                row.predicate,
                                predicate_label,
                                row.object_term,
                                object_label or row.literal or row.measurement,
                            )
                        ] += 1
                        continue
                    record_id = _record_id(row.eol_pk, row.resource_pk, row_number)
                    citation = row.citation or source_citation
                    source_url = _source_url(row.page_id, dataset_url)
                    raw_value = object_label or row.literal or row.measurement or row.object_term
                    excerpt = (
                        f"EOL TraitBank record {record_id}: {predicate_label or row.predicate} = "
                        f"{raw_value}; EOL page {row.page_id or 'unknown'}."
                    )
                    rows.append(
                        {
                            "scientific_name": scientific_name,
                            "trait_name": rule.target_trait,
                            "trait_value": mapped_value,
                            "source_record_id": record_id,
                            "source_citation": citation,
                            "source_url": source_url,
                            "evidence_excerpt": excerpt,
                        }
                    )
            if rows:
                frame = pd.DataFrame(rows, columns=LONG_COLUMNS)
                frame.to_csv(long_path, index=False, mode="a", header=not header_written)
                header_written = True
                n_candidates += len(frame)
            n_scanned += len(chunk)
            if max_trait_rows is not None and n_scanned >= max_trait_rows:
                break

    if not header_written:
        pd.DataFrame(columns=LONG_COLUMNS).to_csv(long_path, index=False)
    _write_counter(unmapped, unmapped_path)
    report = {
        "source_name": source_name,
        "dataset_url": dataset_url,
        "eol_archive": str(eol_archive),
        "eol_archive_sha256": _sha256_file(eol_archive),
        "config_path": str(config_path),
        "ontology_path": str(ontology_path),
        "n_pages": len(pages),
        "n_terms": len(terms),
        "n_mapping_rules": len(rules),
        "n_trait_rows_scanned": int(n_scanned),
        "n_rows_with_mapped_predicate": int(n_predicate_matched),
        "n_unmapped_predicate_rows": int(n_unmapped_predicate),
        "n_unmapped_value_rows": int(n_unmapped_value),
        "n_long_rows": int(n_candidates),
        "long_csv": str(long_path),
        "unmapped_terms_csv": str(unmapped_path),
        "direct_trait_only": True,
        "inferred_csv_policy": "EOL inferred.csv is not used for v2 direct evidence candidates.",
        "curation_rule": "Only explicitly mapped EOL predicate/value terms are exported; all others remain in the aggregate audit.",
    }
    manifest_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    report["manifest_json"] = str(manifest_path)
    return report
