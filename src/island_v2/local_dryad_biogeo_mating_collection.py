"""Prepare review-only SI/SC candidates from a pinned Dryad mating dataset.

The provider is a static, one-shot bulk scan.  It never queries per species and
never turns the quantitative ``mean.tm`` field into an Island mating-system
category.  Only the dataset's explicit categorical contract is used:
``si=0`` means self-incompatible and ``si=1`` means self-compatible.
"""

from __future__ import annotations

import csv
import json
import re
from datetime import UTC, date, datetime
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.local_reproductive_assurance_collection import (
    CHECKPOINT_COLUMNS,
    CONTRACT_VERSION,
    TERMINAL_STATUSES,
    _atomic_write_csv,
    _atomic_write_text,
    _bool,
    _file_sha256,
    _nonnegative_int,
    _sha256_text,
    _text,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

PROVIDER = "dryad_biogeo_mating_2018"
PROVIDER_VERSION = "dryad_577q1_version_1_sha256_pinned"
COLLECTOR_CONTRACT = "local_dryad_biogeo_mating_review_only_v1"
TRAIT = "self_incompatibility"
SOURCE_TYPE = "dryad_biogeo_mating_si_category"

DATASET_DOI = "10.5061/dryad.577q1"
DATASET_URL = f"https://doi.org/{DATASET_DOI}"
DATASET_TITLE = "Global biogeography of mating system variation in seed plants"
DATASET_CITATION = (
    "Moeller et al. (2018), Data from: Global biogeography of mating system "
    "variation in seed plants, Dryad, version 1."
)
DATASET_YEAR = "2018"
ARTICLE_DOI = "10.1111/ele.12738"
ARTICLE_URL = f"https://doi.org/{ARTICLE_DOI}"
ARTICLE_PUBLICATION = (
    "Moeller et al. (2017), Global biogeography of mating system variation in "
    "seed plants, Ecology Letters 20:375-384."
)
ARTICLE_YEAR = "2017"

EXPECTED_CSV_SHA256 = "94c427b84dc35c12f5f35cea0ee7aded7c96b9de16231f2efd771bc0035a0902"
EXPECTED_CSV_ROWS = 624
SOURCE_ENCODING = "latin-1"
SI_VALUE_MAP = {"0": "SI", "1": "SC"}
BINOMIAL_PATTERN = re.compile(r"[A-Z][A-Za-z-]+ [a-z][A-Za-z-]+")

REQUIRED_SOURCE_COLUMNS = (
    "year",
    "reference",
    "family",
    "genus",
    "species",
    "m.d.g",
    "latitude",
    "hemisphere",
    "biome",
    "growth",
    "life.history",
    "si",
    "mean.tm",
    "References",
)

SOURCE_COLUMNS = [
    "source_record_id",
    "source_row_number",
    "accepted_species",
    "queue_family",
    "source_family",
    "source_si_raw",
    "source_si_value",
    "source_mean_tm_raw_unmapped",
    "source_reference",
    "source_references",
    "original_literature_reference",
    "original_literature_year",
    "source_text",
    "source_text_hash",
    "source_url",
    "source_citation",
    "dataset_title",
    "dataset_doi",
    "dataset_url",
    "dataset_citation",
    "dataset_year",
    "article_doi",
    "article_url",
    "article_publication",
    "article_year",
    "source_type",
    "evidence_scope",
    "local_file_path",
    "local_file_hash",
    "retrieval_date",
    "source_run_id",
    "mapping_status",
]

LEAD_COLUMNS = [
    "candidate_id",
    "source_row_index",
    "source_record_id",
    "accepted_species",
    "family",
    "trait_name",
    "provisional_candidate_value",
    "provider",
    "provider_version",
    "source_type",
    "source_url",
    "source_citation",
    "raw_description",
    "source_text",
    "evidence_scope",
    "local_file_path",
    "local_file_hash",
    "source_text_hash",
    "retrieval_date",
    "source_run_id",
    "dataset_doi",
    "article_doi",
    "publication",
    "year",
    "original_literature_reference",
    "source_mean_tm_raw_unmapped",
    "mapping_basis",
    "review_status",
    "review_note",
]

ERROR_COLUMNS = [
    "species",
    "family",
    "provider",
    "error_class",
    "message",
    "source_row_numbers",
    "source_run_id",
    "updated_at",
]


@app.callback()
def main() -> None:
    """Scan the hash-pinned Dryad dataset without materializing traits."""


def _utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _validate_retrieval_date(value: str) -> str:
    normalized = _text(value)
    try:
        parsed = date.fromisoformat(normalized)
    except ValueError as exc:
        raise ValueError("retrieval_date must use YYYY-MM-DD") from exc
    if parsed.isoformat() != normalized:
        raise ValueError("retrieval_date must use YYYY-MM-DD")
    return normalized


def _is_exact_binomial(value: object) -> bool:
    return bool(BINOMIAL_PATTERN.fullmatch(_text(value)))


def _validate_queue(frame: pd.DataFrame, *, label: str) -> pd.DataFrame:
    required = {"accepted_species", "family"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{label} queue missing columns: {sorted(missing)}")
    queue = frame[["accepted_species", "family"]].copy().fillna("")
    queue["accepted_species"] = queue["accepted_species"].map(_text)
    queue["family"] = queue["family"].map(_text)
    if (queue["accepted_species"] == "").any():
        raise ValueError(f"{label} queue contains blank accepted_species")
    if (queue["family"] == "").any():
        raise ValueError(f"{label} queue contains blank family")
    if queue["accepted_species"].duplicated().any():
        raise ValueError(f"{label} queue contains duplicate accepted_species")
    return queue.reset_index(drop=True)


def _validate_queue_relationship(original: pd.DataFrame, current: pd.DataFrame) -> None:
    original_family = original.set_index("accepted_species")["family"].to_dict()
    outside = sorted(set(current["accepted_species"]).difference(original_family))
    if outside:
        raise ValueError(f"current queue contains species outside original queue: {outside[:5]}")
    conflicts = [
        row.accepted_species
        for row in current.itertuples(index=False)
        if original_family[row.accepted_species] != row.family
    ]
    if conflicts:
        raise ValueError(f"current/original queue family conflicts: {conflicts[:5]}")


def _load_source_rows(source_csv: Path, *, expected_hash: str) -> tuple[pd.DataFrame, str]:
    observed_hash = _file_sha256(source_csv)
    if observed_hash != expected_hash:
        raise ValueError(
            f"Dryad CSV SHA-256 mismatch: expected {expected_hash}, observed {observed_hash}"
        )
    records: list[dict[str, str]] = []
    with source_csv.open("r", encoding=SOURCE_ENCODING, newline="") as handle:
        reader = csv.DictReader(handle)
        fields = tuple(reader.fieldnames or ())
        missing = set(REQUIRED_SOURCE_COLUMNS).difference(fields)
        if missing:
            raise ValueError(f"Dryad CSV missing columns: {sorted(missing)}")
        for source_row_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"Dryad CSV row {source_row_number} has extra fields")
            if any(value is None for value in row.values()):
                raise ValueError(f"Dryad CSV row {source_row_number} has missing fields")
            preserved = {field: str(row[field]) for field in fields}
            genus = _text(preserved["genus"])
            epithet = _text(preserved["species"])
            source_name = _text(f"{genus} {epithet}")
            raw_text = json.dumps(
                preserved,
                ensure_ascii=False,
                separators=(",", ":"),
            )
            records.append(
                {
                    **preserved,
                    "source_row_number": str(source_row_number),
                    "source_name": source_name,
                    "source_is_exact_binomial": str(_is_exact_binomial(source_name)).lower(),
                    "source_text": raw_text,
                    "source_text_hash": _sha256_text(raw_text),
                }
            )
    if len(records) != EXPECTED_CSV_ROWS and expected_hash == EXPECTED_CSV_SHA256:
        raise ValueError(
            f"Dryad CSV row-count mismatch: expected {EXPECTED_CSV_ROWS}, observed {len(records)}"
        )
    return pd.DataFrame(records).fillna(""), observed_hash


def _empty_checkpoint() -> pd.DataFrame:
    return pd.DataFrame(columns=CHECKPOINT_COLUMNS)


def ensure_dryad_checkpoint_rows(
    checkpoint: pd.DataFrame,
    original_queue: pd.DataFrame,
    current_queue: pd.DataFrame,
) -> pd.DataFrame:
    """Ensure one resumable provider key for every original-priority species."""
    missing = set(CHECKPOINT_COLUMNS).difference(checkpoint.columns)
    if missing:
        raise ValueError(f"provider checkpoint missing columns: {sorted(missing)}")
    frame = checkpoint[CHECKPOINT_COLUMNS].copy().fillna("")
    if frame.duplicated(["species", "trait", "provider"]).any():
        raise ValueError("provider checkpoint contains duplicate keys")

    original = _validate_queue(original_queue, label="original")
    current = _validate_queue(current_queue, label="current")
    _validate_queue_relationship(original, current)
    current_species = set(current["accepted_species"])
    existing = set(frame[["species", "trait", "provider"]].itertuples(index=False, name=None))
    additions: list[dict[str, object]] = []
    for species in original["accepted_species"]:
        key = (species, TRAIT, PROVIDER)
        if key in existing:
            continue
        if not _is_exact_binomial(species):
            status = "invalid_species_name"
            basis = "original_priority_name_is_not_exact_binomial"
            error = "permanent:invalid_species_name"
        elif species not in current_species:
            status = "provider_seed_covered"
            basis = "reproductive_assurance_completed_before_static_provider_scan"
            error = ""
        else:
            status = "pending"
            basis = "hash_pinned_static_provider_unscanned"
            error = ""
        additions.append(
            {
                "species": species,
                "trait": TRAIT,
                "provider": PROVIDER,
                "provider_version": PROVIDER_VERSION,
                "status": status,
                "terminal": status in TERMINAL_STATUSES,
                "attempts": 0,
                "candidate_count": 0,
                "accepted_record_count": 0,
                "legacy_status": "",
                "legacy_enabled": "",
                "last_error": error,
                "updated_at": "",
                "last_wave_id": "",
                "source_run_id": "local_dryad_biogeo_mating_2018",
                "migration_basis": basis,
                "contract_version": CONTRACT_VERSION,
            }
        )
    if additions:
        frame = pd.concat([frame, pd.DataFrame(additions)], ignore_index=True, sort=False)

    frame["terminal"] = frame["terminal"].map(_bool).astype(bool)
    for column in ("attempts", "candidate_count", "accepted_record_count"):
        frame[column] = frame[column].map(_nonnegative_int).astype("int64")

    provider_rows = frame["provider"].eq(PROVIDER)
    wrong_trait = provider_rows & ~frame["trait"].eq(TRAIT)
    if wrong_trait.any():
        raise ValueError("Dryad provider checkpoint contains unsupported trait rows")
    wrong_version = provider_rows & ~frame["provider_version"].eq(PROVIDER_VERSION)
    if wrong_version.any():
        raise ValueError("Dryad provider checkpoint contains a different provider version")

    interrupted = provider_rows & frame["status"].eq("running")
    frame.loc[interrupted, "status"] = "retry"
    frame.loc[interrupted, "terminal"] = False
    frame.loc[interrupted, "last_error"] = "transient:interrupted_static_scan"

    no_longer_current = (
        provider_rows
        & frame["species"].isin(set(original["accepted_species"]).difference(current_species))
        & frame["status"].isin({"pending", "retry"})
    )
    frame.loc[no_longer_current, "status"] = "provider_seed_covered"
    frame.loc[no_longer_current, "terminal"] = True
    frame.loc[no_longer_current, "last_error"] = ""
    frame.loc[no_longer_current, "migration_basis"] = (
        "reproductive_assurance_completed_before_static_provider_scan"
    )

    if frame.duplicated(["species", "trait", "provider"]).any():
        raise AssertionError("extended provider checkpoint contains duplicate keys")
    return frame.sort_values(["species", "trait", "provider"]).reset_index(drop=True)


def _append_unique(
    existing_path: Path,
    new_rows: pd.DataFrame,
    *,
    columns: list[str],
    keys: list[str],
) -> pd.DataFrame:
    if existing_path.exists():
        existing = pd.read_csv(existing_path, dtype=str).fillna("")
        missing = set(columns).difference(existing.columns)
        if missing:
            raise ValueError(f"existing {existing_path.name} missing columns: {sorted(missing)}")
        existing = existing[columns]
    else:
        existing = pd.DataFrame(columns=columns)
    combined = pd.concat([existing, new_rows[columns]], ignore_index=True, sort=False).fillna("")
    combined = combined.drop_duplicates(keys, keep="first").reset_index(drop=True)
    _atomic_write_csv(combined, existing_path)
    return combined


def _source_output_row(
    row: pd.Series,
    *,
    queue_family: str,
    source_csv: Path,
    source_hash: str,
    retrieval_date: str,
    source_run_id: str,
    mapping_status: str,
) -> dict[str, object]:
    reference = str(row["reference"])
    references = str(row["References"])
    original_reference = _text(references) or _text(reference)
    source_row_number = _text(row["source_row_number"])
    source_citation = (
        f"{DATASET_CITATION} DOI:{DATASET_DOI}. Primary article: "
        f"{ARTICLE_PUBLICATION} DOI:{ARTICLE_DOI}."
    )
    if original_reference:
        source_citation += f" Original literature: {original_reference}"
    return {
        "source_record_id": f"dryad_577q1:csv_row_{source_row_number}",
        "source_row_number": source_row_number,
        "accepted_species": _text(row["source_name"]),
        "queue_family": queue_family,
        "source_family": _text(row["family"]),
        "source_si_raw": _text(row["si"]),
        "source_si_value": SI_VALUE_MAP.get(_text(row["si"]), ""),
        "source_mean_tm_raw_unmapped": _text(row["mean.tm"]),
        "source_reference": reference,
        "source_references": references,
        "original_literature_reference": original_reference,
        "original_literature_year": _text(row["year"]),
        "source_text": str(row["source_text"]),
        "source_text_hash": _text(row["source_text_hash"]),
        "source_url": DATASET_URL,
        "source_citation": source_citation,
        "dataset_title": DATASET_TITLE,
        "dataset_doi": DATASET_DOI,
        "dataset_url": DATASET_URL,
        "dataset_citation": DATASET_CITATION,
        "dataset_year": DATASET_YEAR,
        "article_doi": ARTICLE_DOI,
        "article_url": ARTICLE_URL,
        "article_publication": ARTICLE_PUBLICATION,
        "article_year": ARTICLE_YEAR,
        "source_type": SOURCE_TYPE,
        "evidence_scope": "species_direct",
        "local_file_path": str(source_csv.resolve()),
        "local_file_hash": source_hash,
        "retrieval_date": retrieval_date,
        "source_run_id": source_run_id,
        "mapping_status": mapping_status,
    }


def _candidate_row(source: dict[str, object]) -> dict[str, object]:
    species = _text(source["accepted_species"])
    value = _text(source["source_si_value"])
    record_id = _text(source["source_record_id"])
    candidate_key = f"{PROVIDER}|{record_id}|{species}|{TRAIT}|{value}"
    original_reference = _text(source["original_literature_reference"])
    return {
        "candidate_id": f"dryad-biogeo-{_sha256_text(candidate_key)[:24]}",
        "source_row_index": "",
        "source_record_id": record_id,
        "accepted_species": species,
        "family": _text(source["queue_family"]),
        "trait_name": TRAIT,
        "provisional_candidate_value": value,
        "provider": PROVIDER,
        "provider_version": PROVIDER_VERSION,
        "source_type": SOURCE_TYPE,
        "source_url": _text(source["source_url"]),
        "source_citation": _text(source["source_citation"]),
        "raw_description": str(source["source_text"]),
        "source_text": str(source["source_text"]),
        "evidence_scope": "species_direct",
        "local_file_path": _text(source["local_file_path"]),
        "local_file_hash": _text(source["local_file_hash"]),
        "source_text_hash": _text(source["source_text_hash"]),
        "retrieval_date": _text(source["retrieval_date"]),
        "source_run_id": _text(source["source_run_id"]),
        "dataset_doi": DATASET_DOI,
        "article_doi": ARTICLE_DOI,
        "publication": ARTICLE_PUBLICATION,
        "year": _text(source["original_literature_year"]) or ARTICLE_YEAR,
        "original_literature_reference": original_reference,
        "source_mean_tm_raw_unmapped": _text(source["source_mean_tm_raw_unmapped"]),
        "mapping_basis": "Dryad field contract: si=0 means SI; si=1 means SC.",
        "review_status": "unreviewed",
        "review_note": (
            "Review the original species-level literature context before materialization; "
            "mean.tm was preserved but not mapped."
        ),
    }


def _classify_species_rows(
    species: str,
    family: str,
    source_rows: pd.DataFrame,
) -> tuple[str, list[int], str]:
    """Return terminal status, candidate row indexes, and a fail-closed error."""
    exact = source_rows.loc[
        source_rows["source_name"].eq(species) & source_rows["source_is_exact_binomial"].eq("true")
    ]
    with_si = exact.loc[exact["si"].map(_text).ne("")]
    if with_si.empty:
        return "provider_no_result", [], ""

    unsupported = with_si.loc[~with_si["si"].map(_text).isin(SI_VALUE_MAP)]
    if not unsupported.empty:
        rows = ",".join(unsupported["source_row_number"].map(_text))
        return "permanent_failure", [], f"unsupported_si_value:rows={rows}"

    family_mismatch = with_si.loc[~with_si["family"].map(_text).eq(family)]
    if not family_mismatch.empty:
        rows = ",".join(family_mismatch["source_row_number"].map(_text))
        return "permanent_failure", [], f"source_family_mismatch:rows={rows}"

    mapped = with_si["si"].map(_text).map(SI_VALUE_MAP)
    if mapped.nunique() != 1:
        rows = ",".join(with_si["source_row_number"].map(_text))
        values = ",".join(sorted(set(mapped)))
        return "permanent_failure", [], f"conflicting_si_values:{values}:rows={rows}"
    return "provider_pass_complete_with_candidates", list(with_si.index), ""


def _set_running(
    checkpoint: pd.DataFrame,
    species: list[str],
    *,
    batch_id: str,
    updated_at: str,
    source_run_id: str,
) -> pd.DataFrame:
    frame = checkpoint.copy()
    scheduled = set(species)
    mask = (
        frame["species"].isin(scheduled) & frame["trait"].eq(TRAIT) & frame["provider"].eq(PROVIDER)
    )
    observed = frame.loc[mask, "species"].map(_text)
    if len(observed) != len(scheduled) or set(observed) != scheduled:
        missing = sorted(scheduled.difference(observed))
        raise ValueError(f"missing Dryad checkpoint keys: {missing[:5]}")
    frame.loc[mask, "status"] = "running"
    frame.loc[mask, "terminal"] = False
    frame.loc[mask, "attempts"] = frame.loc[mask, "attempts"].map(_nonnegative_int) + 1
    frame.loc[mask, "last_error"] = ""
    frame.loc[mask, "updated_at"] = updated_at
    frame.loc[mask, "last_wave_id"] = batch_id
    frame.loc[mask, "source_run_id"] = source_run_id
    return frame


def collect_dryad_biogeo_candidates(
    *,
    original_queue: pd.DataFrame,
    current_queue: pd.DataFrame,
    provider_checkpoint: pd.DataFrame,
    provider_checkpoint_path: Path,
    source_csv: Path,
    output_dir: Path,
    retrieval_date: str,
    expected_csv_sha256: str = EXPECTED_CSV_SHA256,
) -> dict[str, Any]:
    """Scan the pinned static CSV and persist review-only outputs."""
    retrieval_date = _validate_retrieval_date(retrieval_date)
    original = _validate_queue(original_queue, label="original")
    current = _validate_queue(current_queue, label="current")
    _validate_queue_relationship(original, current)
    source_rows, source_hash = _load_source_rows(
        source_csv,
        expected_hash=expected_csv_sha256,
    )
    source_run_id = f"{PROVIDER}_{source_hash[:16]}"
    checkpoint = ensure_dryad_checkpoint_rows(provider_checkpoint, original, current)

    current_order = current["accepted_species"].tolist()
    eligible = checkpoint.loc[
        checkpoint["provider"].eq(PROVIDER)
        & checkpoint["trait"].eq(TRAIT)
        & checkpoint["status"].isin({"pending", "retry"})
        & ~checkpoint["terminal"].map(_bool)
    ]
    eligible_species = set(eligible["species"].map(_text))
    scheduled = [species for species in current_order if species in eligible_species]
    started_at = _utc_now()
    batch_id = f"dryad_biogeo_{started_at.replace('-', '').replace(':', '').replace('.', '')}"
    checkpoint = _set_running(
        checkpoint,
        scheduled,
        batch_id=batch_id,
        updated_at=started_at,
        source_run_id=source_run_id,
    )
    _atomic_write_csv(checkpoint, provider_checkpoint_path)

    family_by_species = current.set_index("accepted_species")["family"].to_dict()
    new_sources: list[dict[str, object]] = []
    new_candidates: list[dict[str, object]] = []
    errors: list[dict[str, object]] = []
    outcomes: dict[str, tuple[str, int, str]] = {}

    for species in scheduled:
        family = family_by_species[species]
        exact_rows = source_rows.loc[
            source_rows["source_name"].eq(species)
            & source_rows["source_is_exact_binomial"].eq("true")
        ]
        status, candidate_indexes, error = _classify_species_rows(
            species,
            family,
            source_rows,
        )
        candidate_index_set = set(candidate_indexes)
        for source_index, row in exact_rows.iterrows():
            if error:
                mapping_status = "rejected_fail_closed"
            elif source_index in candidate_index_set:
                mapping_status = "review_only_candidate"
            else:
                mapping_status = "no_explicit_si_category_mean_tm_not_mapped"
            source_output = _source_output_row(
                row,
                queue_family=family,
                source_csv=source_csv,
                source_hash=source_hash,
                retrieval_date=retrieval_date,
                source_run_id=source_run_id,
                mapping_status=mapping_status,
            )
            new_sources.append(source_output)
            if source_index in candidate_index_set and not error:
                new_candidates.append(_candidate_row(source_output))
        outcomes[species] = (status, len(candidate_indexes) if not error else 0, error)
        if error:
            row_numbers = ",".join(exact_rows["source_row_number"].map(_text))
            errors.append(
                {
                    "species": species,
                    "family": family,
                    "provider": PROVIDER,
                    "error_class": "permanent",
                    "message": error,
                    "source_row_numbers": row_numbers,
                    "source_run_id": source_run_id,
                    "updated_at": _utc_now(),
                }
            )

    output_dir.mkdir(parents=True, exist_ok=True)
    new_source_frame = pd.DataFrame(new_sources, columns=SOURCE_COLUMNS).fillna("")
    source_path = output_dir / "source_evidence.csv"
    sources = _append_unique(
        source_path,
        new_source_frame,
        columns=SOURCE_COLUMNS,
        keys=["source_record_id"],
    )
    source_indexes = {
        record_id: index for index, record_id in enumerate(sources["source_record_id"].map(_text))
    }
    for candidate in new_candidates:
        candidate["source_row_index"] = source_indexes[_text(candidate["source_record_id"])]

    new_candidate_frame = pd.DataFrame(new_candidates, columns=LEAD_COLUMNS).fillna("")
    candidate_path = output_dir / "candidate_leads.csv"
    candidates = _append_unique(
        candidate_path,
        new_candidate_frame,
        columns=LEAD_COLUMNS,
        keys=["candidate_id"],
    )
    new_error_frame = pd.DataFrame(errors, columns=ERROR_COLUMNS).fillna("")
    error_path = output_dir / "lookup_errors.csv"
    cumulative_errors = _append_unique(
        error_path,
        new_error_frame,
        columns=ERROR_COLUMNS,
        keys=["species", "message", "source_run_id"],
    )

    finished_at = _utc_now()
    outcome_status = {species: status for species, (status, _count, _error) in outcomes.items()}
    outcome_count = {species: count for species, (_status, count, _error) in outcomes.items()}
    outcome_error = {species: error for species, (_status, _count, error) in outcomes.items()}
    outcome_mask = (
        checkpoint["species"].isin(outcomes)
        & checkpoint["trait"].eq(TRAIT)
        & checkpoint["provider"].eq(PROVIDER)
    )
    outcome_species = checkpoint.loc[outcome_mask, "species"].map(_text)
    if (
        len(outcome_species) != len(outcomes)
        or set(outcome_species) != set(outcomes)
        or not checkpoint.loc[outcome_mask, "status"].eq("running").all()
    ):
        raise AssertionError("Dryad checkpoint did not retain every running provider key")
    checkpoint.loc[outcome_mask, "status"] = outcome_species.map(outcome_status)
    checkpoint.loc[outcome_mask, "terminal"] = True
    checkpoint.loc[outcome_mask, "candidate_count"] = outcome_species.map(outcome_count)
    checkpoint.loc[outcome_mask, "accepted_record_count"] = 0
    checkpoint.loc[outcome_mask, "last_error"] = outcome_species.map(
        lambda species: f"permanent:{outcome_error[species]}" if outcome_error[species] else ""
    )
    checkpoint.loc[outcome_mask, "updated_at"] = finished_at
    checkpoint.loc[outcome_mask, "migration_basis"] = outcome_species.map(
        lambda species: (
            "hash_pinned_static_scan_failed_closed"
            if outcome_error[species]
            else "hash_pinned_static_scan_review_only"
        )
    )
    _atomic_write_csv(checkpoint, provider_checkpoint_path)

    original_species = set(original["accepted_species"])
    final_provider = checkpoint.loc[
        checkpoint["provider"].eq(PROVIDER)
        & checkpoint["trait"].eq(TRAIT)
        & checkpoint["species"].isin(original_species)
    ]
    if len(final_provider) != len(original) or not final_provider["terminal"].map(_bool).all():
        raise AssertionError("not every original-priority Dryad provider key is terminal")
    if (final_provider["accepted_record_count"].map(_nonnegative_int) != 0).any():
        raise AssertionError("review-only provider must not mark accepted records")

    report: dict[str, Any] = {
        "collector_contract": COLLECTOR_CONTRACT,
        "provider": PROVIDER,
        "provider_version": PROVIDER_VERSION,
        "source_run_id": source_run_id,
        "batch_id": batch_id,
        "started_at_utc": started_at,
        "finished_at_utc": finished_at,
        "dataset_doi": DATASET_DOI,
        "dataset_url": DATASET_URL,
        "article_doi": ARTICLE_DOI,
        "source_csv_path": str(source_csv.resolve()),
        "source_csv_sha256": source_hash,
        "source_csv_rows": int(len(source_rows)),
        "original_priority_species": int(len(original)),
        "current_queue_species": int(len(current)),
        "scheduled_species": int(len(scheduled)),
        "new_source_rows": int(len(new_source_frame)),
        "new_review_only_candidates": int(len(new_candidate_frame)),
        "new_candidate_species": int(new_candidate_frame["accepted_species"].nunique())
        if len(new_candidate_frame)
        else 0,
        "new_permanent_failures": int(len(new_error_frame)),
        "cumulative_source_rows": int(len(sources)),
        "cumulative_review_only_candidates": int(len(candidates)),
        "cumulative_errors": int(len(cumulative_errors)),
        "mean_tm_mapped": False,
        "biological_values_materialized": False,
        "shared_materialized_outputs_mutated": False,
        "global_fallback_used": False,
        "genus_or_family_inference_used": False,
        "all_original_provider_keys_terminal": True,
        "outputs": {},
    }
    for path in (source_path, candidate_path, error_path, provider_checkpoint_path):
        report["outputs"][path.name] = {
            "path": str(path.resolve()),
            "bytes": path.stat().st_size,
            "sha256": _file_sha256(path),
        }
    report_path = output_dir / "campaign_status.json"
    _atomic_write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        report_path,
    )
    return report


@app.command("scan")
def scan_command(
    original_queue_csv: Path = typer.Option(..., exists=True, readable=True),
    current_queue_csv: Path = typer.Option(..., exists=True, readable=True),
    source_csv: Path = typer.Option(..., exists=True, readable=True),
    input_provider_checkpoint_csv: Path = typer.Option(..., exists=True, readable=True),
    output_provider_checkpoint_csv: Path = typer.Option(...),
    output_dir: Path = typer.Option(...),
    retrieval_date: str = typer.Option(..., help="Pinned source retrieval date, YYYY-MM-DD."),
) -> None:
    """Scan once, writing only to the explicit local output checkpoint."""
    if input_provider_checkpoint_csv.resolve() == output_provider_checkpoint_csv.resolve():
        raise typer.BadParameter(
            "output-provider-checkpoint-csv must differ from the input checkpoint"
        )
    checkpoint_source = (
        output_provider_checkpoint_csv
        if output_provider_checkpoint_csv.exists()
        else input_provider_checkpoint_csv
    )
    report = collect_dryad_biogeo_candidates(
        original_queue=pd.read_csv(original_queue_csv, dtype=str).fillna(""),
        current_queue=pd.read_csv(current_queue_csv, dtype=str).fillna(""),
        provider_checkpoint=pd.read_csv(checkpoint_source, dtype=str).fillna(""),
        provider_checkpoint_path=output_provider_checkpoint_csv,
        source_csv=source_csv,
        output_dir=output_dir,
        retrieval_date=retrieval_date,
    )
    typer.echo(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    app()
