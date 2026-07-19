"""Strict local import for explicitly approved targeted-web review lanes.

This module performs no acquisition.  It binds an explicit approval manifest to
the exact discovery tables, current priority queue bytes, and locally cached raw
source bytes.  Discovery rows are never accepted merely because they appear in
``candidate_records.csv``.  Only approval-selected rows can become standard
source evidence, reviewed candidates, and direct records.

The original ``import-reviewed`` command remains pinned to its exact 20-species
lane.  ``import-reviewed-round`` supports later, explicitly approved rounds by
deriving the inspected set from the approval and requiring the discovery tables
to cover that exact same set.  Later rounds may only append complete, disjoint
species-by-trait blocks to an already-terminal targeted-web checkpoint.

The provider checkpoint is written last.  Consequently, an interrupted import
may leave replaceable output files, but it cannot advertise terminal provider
state before all materialization outputs exist.
"""

from __future__ import annotations

import html
import json
import re
from datetime import UTC, datetime
from io import BytesIO
from pathlib import Path
from urllib.parse import unquote, urlparse

import pandas as pd
import typer
from pypdf import PdfReader

from island_v2.local_reproductive_assurance_collection import (
    ALLOWED_DIRECT_VALUES,
    CHECKPOINT_COLUMNS,
    CONTRACT_VERSION,
    RA_TRAITS,
    REVIEW_BASIS,
    TERMINAL_STATUSES,
    _atomic_write_csv,
    _atomic_write_text,
    _bool,
    _explicit_mapping_supported,
    _file_sha256,
    _nonnegative_int,
    _sha256_text,
    _source_provider,
    _text,
    _validate_queue,
    materialize_reviewed_cache,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

PROVIDER = "targeted_web"
PROVIDER_VERSION = "targeted_web_exact_20_species_explicit_approval_v1"
APPROVAL_CONTRACT = "targeted_web_explicit_approval_v1"
IMPORT_CONTRACT = "targeted_web_local_import_v1"
DISCOVERY_QUEUE_SHA256 = "81477c220df8173f0a805fa588b2343ecc1939c062ed6cbd234b2fb4ba56eda1"

INSPECTED_SPECIES = (
    "Begonia hirtella",
    "Campanula persicifolia",
    "Cuscuta gronovii",
    "Cuscuta pentagona",
    "Gentiana andrewsii",
    "Ipomoea hederifolia",
    "Ipomoea imperati",
    "Lonicera periclymenum",
    "Mimulus gracilis",
    "Penstemon digitalis",
    "Primula denticulata",
    "Ribes aureum",
    "Salvia azurea",
    "Silene antirrhina",
    "Silene dioica",
    "Silene gallica",
    "Solanum melongena",
    "Vaccinium ovatum",
    "Viola cornuta",
    "Viola tricolor",
)

CANDIDATE_REQUIRED_COLUMNS = {
    "status",
    "species",
    "trait",
    "category",
    "raw_text",
    "url",
    "retrieval_url",
    "doi",
    "publication",
    "venue",
    "year",
    "source_type",
    "evidence_basis",
    "local_file",
    "local_file_sha256",
    "retrieval_date",
    "review_note",
}
REJECTED_REQUIRED_COLUMNS = {"status", "species", "reason"}

SOURCE_EVIDENCE_COLUMNS = [
    "source_evidence_id",
    "accepted_species",
    "source_text",
    "evidence_excerpt",
    "source_url",
    "retrieval_url",
    "source_citation",
    "source_type",
    "provider",
    "evidence_scope",
    "DOI",
    "publication",
    "venue",
    "year",
    "local_file_path",
    "local_file_hash",
    "retrieval_date",
    "source_run_id",
    "source_text_hash",
    "discovery_row_index",
    "discovery_file_hash",
]

REVIEWED_CANDIDATE_COLUMNS = [
    "candidate_id",
    "accepted_species",
    "trait_name",
    "provisional_candidate_value",
    "matched_term",
    "candidate_class",
    "provider",
    "source_type",
    "source_url",
    "source_citation",
    "raw_description",
    "source_text",
    "evidence_scope",
    "wild_or_cultivated",
    "review_status",
    "provider_status",
    "provider_attempts",
    "provider_updated_at",
    "provider_policy_version",
    "local_file_path",
    "local_file_hash",
    "retrieval_date",
    "retrieval_date_basis",
    "source_run_id",
    "audit_status",
    "audit_reason",
]


@app.callback()
def main() -> None:
    """Import explicitly approved targeted-web evidence without acquisition."""


def _literal(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _valid_sha256(value: object) -> bool:
    return bool(re.fullmatch(r"[0-9a-f]{64}", _text(value).casefold()))


def _valid_url(value: object) -> bool:
    parsed = urlparse(_text(value))
    return parsed.scheme in {"http", "https"} and bool(parsed.netloc)


def _clean_doi(value: object) -> str:
    doi = unquote(_text(value)).casefold()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if doi.startswith(prefix):
            doi = doi[len(prefix) :]
    return doi.rstrip(".,;)")


def _valid_doi(value: str) -> bool:
    return bool(re.fullmatch(r"10\.\d{4,9}/\S+", value, flags=re.IGNORECASE))


def _valid_date(value: object) -> bool:
    text = _text(value)
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", text):
        return False
    try:
        datetime.strptime(text, "%Y-%m-%d")
    except ValueError:
        return False
    return True


def _valid_timestamp(value: object) -> bool:
    text = _text(value)
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return False
    return parsed.tzinfo is not None


def _exact_binomial(value: object) -> bool:
    return bool(re.fullmatch(r"[A-Z][A-Za-z-]+ [a-z][A-Za-z-]+", _text(value)))


def _exact_species_named(species: str, excerpt: str) -> bool:
    expression = r"\s+".join(re.escape(token) for token in species.split())
    return bool(re.search(rf"(?<!\w){expression}(?!\w)", excerpt, flags=re.IGNORECASE))


def _require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{label} missing columns: {sorted(missing)}")


def _validate_discovery_tables(
    candidates: pd.DataFrame,
    rejected: pd.DataFrame,
    *,
    inspected_species: tuple[str, ...],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    _require_columns(candidates, CANDIDATE_REQUIRED_COLUMNS, "targeted-web candidate table")
    _require_columns(rejected, REJECTED_REQUIRED_COLUMNS, "targeted-web rejected table")
    candidate_frame = candidates.copy().fillna("")
    rejected_frame = rejected.copy().fillna("")
    candidate_frame["species"] = candidate_frame["species"].map(_text)
    rejected_frame["species"] = rejected_frame["species"].map(_text)

    if not candidate_frame["status"].map(_text).eq("candidate").all():
        raise ValueError("targeted-web candidate table contains a non-candidate row")
    if not rejected_frame["status"].map(_text).str.match(r"^(?:rejected|held)_").all():
        raise ValueError("targeted-web rejected table contains an unsupported status")
    if not candidate_frame["trait"].map(_text).eq("reproductive_assurance").all():
        raise ValueError("targeted-web discovery row is outside reproductive assurance")

    candidate_species = set(candidate_frame["species"])
    rejected_species = set(rejected_frame["species"])
    overlap = sorted(candidate_species.intersection(rejected_species))
    if overlap:
        raise ValueError(f"held/rejected species also appear as candidates: {overlap}")
    if rejected_frame["species"].duplicated().any():
        raise ValueError("targeted-web rejected table contains duplicate species")
    observed = candidate_species.union(rejected_species)
    expected = set(inspected_species)
    if observed != expected:
        raise ValueError(
            "targeted-web discovery tables do not cover the exact inspected species set"
        )
    invalid = sorted(species for species in observed if not _exact_binomial(species))
    if invalid:
        raise ValueError(f"targeted-web inspected names are not exact binomials: {invalid}")
    return candidate_frame.reset_index(drop=True), rejected_frame.reset_index(drop=True)


def _validate_approval(
    payload: object,
    *,
    discovery_queue_hash: str,
    current_queue_hash: str,
    candidate_hash: str,
    rejected_hash: str,
    inspected_species: tuple[str, ...],
) -> dict[str, object]:
    if not isinstance(payload, dict):
        raise ValueError("targeted-web approval must contain a JSON object")
    required = {
        "contract_version",
        "approval_id",
        "discovery_queue_sha256",
        "current_queue_sha256",
        "candidate_file_sha256",
        "rejected_file_sha256",
        "inspected_species",
        "completed_at_utc",
        "accepted",
    }
    missing = required.difference(payload)
    if missing:
        raise ValueError(f"targeted-web approval missing fields: {sorted(missing)}")
    if _text(payload["contract_version"]) != APPROVAL_CONTRACT:
        raise ValueError("targeted-web approval contract_version mismatch")
    if not _text(payload["approval_id"]):
        raise ValueError("targeted-web approval has a blank approval_id")
    expected_hashes = {
        "discovery_queue_sha256": discovery_queue_hash,
        "current_queue_sha256": current_queue_hash,
        "candidate_file_sha256": candidate_hash,
        "rejected_file_sha256": rejected_hash,
    }
    for field, expected in expected_hashes.items():
        observed = _text(payload[field]).casefold()
        if not _valid_sha256(observed) or observed != expected:
            raise ValueError(f"targeted-web approval {field} mismatch")
    approved_species = payload["inspected_species"]
    if not isinstance(approved_species, list) or not all(
        isinstance(species, str) for species in approved_species
    ):
        raise ValueError("targeted-web approval inspected_species must be a string list")
    if len(approved_species) != len(set(approved_species)):
        raise ValueError("targeted-web approval inspected_species contains duplicates")
    if set(approved_species) != set(inspected_species):
        raise ValueError("targeted-web approval inspected_species mismatch")
    if not _valid_timestamp(payload["completed_at_utc"]):
        raise ValueError("targeted-web approval completed_at_utc is invalid")
    accepted = payload["accepted"]
    if not isinstance(accepted, list) or not all(isinstance(row, dict) for row in accepted):
        raise ValueError("targeted-web approval accepted must be a list of objects")
    return payload


def _dynamic_inspected_species(
    candidate_records: pd.DataFrame,
    rejected_hits: pd.DataFrame,
    approval_payload: object,
) -> tuple[str, ...]:
    """Return the exact approved/discovered species set for one later round."""
    if not isinstance(approval_payload, dict):
        raise ValueError("targeted-web approval must contain a JSON object")
    approved = approval_payload.get("inspected_species")
    if not isinstance(approved, list) or not all(isinstance(species, str) for species in approved):
        raise ValueError("targeted-web approval inspected_species must be a string list")
    inspected = tuple(_text(species) for species in approved)
    if not inspected or any(not species for species in inspected):
        raise ValueError("targeted-web approval inspected_species must be nonempty")
    if len(inspected) != len(set(inspected)):
        raise ValueError("targeted-web approval inspected_species contains duplicates")
    # This validates exact equality with the candidate+rejected discovery union,
    # exact binomials, statuses, traits, and duplicate/overlap constraints before
    # any checkpoint state is considered.
    _validate_discovery_tables(
        candidate_records,
        rejected_hits,
        inspected_species=inspected,
    )
    return inspected


def _candidate_raw_path(row: pd.Series, candidate_csv: Path) -> Path:
    relative = Path(_text(row["local_file"]))
    if not relative.parts or relative.is_absolute():
        raise ValueError("targeted-web local_file must be relative to the discovery directory")
    raw_root = (candidate_csv.parent / "raw").resolve()
    path = (candidate_csv.parent / relative).resolve()
    try:
        path.relative_to(raw_root)
    except ValueError as exc:
        raise ValueError("targeted-web local_file escapes the discovery raw directory") from exc
    if not path.is_file():
        raise ValueError(f"targeted-web local raw file is unavailable: {path}")
    return path


def _excerpt_is_exact(excerpt: str, discovery_text: str, raw_bytes: bytes) -> bool:
    if excerpt in discovery_text:
        return True
    decoded: list[str] = []
    for encoding in ("utf-8", "utf-8-sig", "latin-1"):
        try:
            decoded.append(raw_bytes.decode(encoding))
        except UnicodeDecodeError:
            continue
    if any(excerpt in text or excerpt in html.unescape(text) for text in decoded):
        return True
    if not raw_bytes.startswith(b"%PDF-"):
        return False
    try:
        extracted = "\n".join(
            page.extract_text() or "" for page in PdfReader(BytesIO(raw_bytes)).pages
        )
    except Exception:
        return False
    return excerpt in extracted


def _validate_source_provenance(
    row: pd.Series,
    *,
    approval_row: dict[str, object],
    candidate_csv: Path,
) -> tuple[Path, str, str]:
    source_url = _text(row["url"])
    if not _valid_url(source_url):
        raise ValueError("targeted-web selected source has an invalid URL")
    if _text(approval_row.get("source_url")) != source_url:
        raise ValueError("targeted-web approval source_url does not match discovery row")
    retrieval_url = _text(row["retrieval_url"])
    if retrieval_url and not _valid_url(retrieval_url):
        raise ValueError("targeted-web selected source has an invalid retrieval URL")
    publication = _text(row["publication"])
    if not publication:
        raise ValueError("targeted-web selected source lacks a publication")
    year = _text(row["year"])
    if not re.fullmatch(r"\d{4}", year) or not 1500 <= int(year) <= datetime.now(UTC).year:
        raise ValueError("targeted-web selected source has an invalid publication year")
    retrieval_date = _text(row["retrieval_date"])
    if not _valid_date(retrieval_date):
        raise ValueError("targeted-web selected source has an invalid retrieval date")
    source_type = _text(row["source_type"])
    if _source_provider(source_type, source_url) != PROVIDER:
        raise ValueError(f"targeted-web selected source_type is not registered: {source_type}")

    doi = _clean_doi(row["doi"])
    if doi and not _valid_doi(doi):
        raise ValueError("targeted-web selected source has an invalid DOI")
    parsed = urlparse(source_url)
    if parsed.netloc.casefold() == "doi.org":
        url_doi = _clean_doi(parsed.path.lstrip("/"))
        if not doi or url_doi != doi:
            raise ValueError("targeted-web DOI URL does not match the DOI field")

    path = _candidate_raw_path(row, candidate_csv)
    claimed_hash = _text(row["local_file_sha256"]).casefold()
    if not _valid_sha256(claimed_hash):
        raise ValueError("targeted-web selected source has an invalid local file hash")
    if _text(approval_row.get("local_file_sha256")).casefold() != claimed_hash:
        raise ValueError("targeted-web approval local file hash does not match discovery row")
    if _file_sha256(path) != claimed_hash:
        raise ValueError(f"targeted-web local raw hash mismatch: {path}")
    return path, claimed_hash, doi


def _source_citation(row: pd.Series, doi: str) -> str:
    citation = f"{_text(row['publication'])} ({_text(row['year'])})."
    venue = _text(row["venue"])
    if venue:
        citation = f"{citation} {venue}."
    if doi:
        citation = f"{citation} DOI:{doi}"
    return citation


def build_targeted_web_review_batch(
    *,
    queue: pd.DataFrame,
    candidate_records: pd.DataFrame,
    rejected_hits: pd.DataFrame,
    approval_payload: object,
    current_queue_hash: str,
    candidate_hash: str,
    rejected_hash: str,
    candidate_csv: Path,
    expected_current_queue_hash: str,
    expected_discovery_queue_hash: str = DISCOVERY_QUEUE_SHA256,
    inspected_species: tuple[str, ...] = INSPECTED_SPECIES,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object], dict[str, object]]:
    """Build approved evidence and candidates; never infer approval from discovery status."""
    if not _valid_sha256(expected_discovery_queue_hash):
        raise ValueError("targeted-web expected discovery queue hash is invalid")
    if not _valid_sha256(expected_current_queue_hash):
        raise ValueError("targeted-web expected current queue hash is invalid")
    if current_queue_hash != expected_current_queue_hash:
        raise ValueError("targeted-web current queue hash does not match the expected hash")
    queue_species = set(_validate_queue(queue))
    candidates, rejected = _validate_discovery_tables(
        candidate_records,
        rejected_hits,
        inspected_species=inspected_species,
    )
    approval = _validate_approval(
        approval_payload,
        discovery_queue_hash=expected_discovery_queue_hash,
        current_queue_hash=current_queue_hash,
        candidate_hash=candidate_hash,
        rejected_hash=rejected_hash,
        inspected_species=inspected_species,
    )
    rejected_species = set(rejected["species"])
    approval_id = _text(approval["approval_id"])
    completed_at = _text(approval["completed_at_utc"])

    source_rows: list[dict[str, object]] = []
    candidate_rows: list[dict[str, object]] = []
    decisions: list[dict[str, str]] = []
    seen_keys: set[tuple[int, str]] = set()
    accepted_rows = sorted(
        approval["accepted"],
        key=lambda item: (
            int(str(item.get("discovery_row_index", "-1"))),
            _text(item.get("accepted_trait")),
            _text(item.get("accepted_value")),
        ),
    )
    for approval_row in accepted_rows:
        try:
            row_index = int(str(approval_row.get("discovery_row_index", "")).strip())
        except ValueError as exc:
            raise ValueError("targeted-web approval has an invalid discovery_row_index") from exc
        if row_index < 0 or row_index not in candidates.index:
            raise ValueError(f"targeted-web approval row is unavailable: {row_index}")
        row = candidates.loc[row_index]
        species = _text(row["species"])
        if _text(approval_row.get("species")) != species:
            raise ValueError("targeted-web approval species does not match discovery row")
        if species in rejected_species:
            raise ValueError(f"held/rejected species cannot be approved: {species}")
        if species not in queue_species:
            raise ValueError(
                f"targeted-web approved species is outside the current queue: {species}"
            )
        trait = _text(approval_row.get("accepted_trait"))
        value = _text(approval_row.get("accepted_value"))
        if trait not in RA_TRAITS or value not in ALLOWED_DIRECT_VALUES.get(trait, frozenset()):
            raise ValueError(f"targeted-web approval has an invalid mapping: {trait}={value}")
        key = (row_index, trait)
        if key in seen_keys:
            raise ValueError("targeted-web approval contains a duplicate source-row trait key")
        seen_keys.add(key)

        excerpt = _literal(approval_row.get("evidence_excerpt"))
        if not excerpt:
            raise ValueError("targeted-web approval has a blank evidence_excerpt")
        reason = _text(approval_row.get("review_reason"))
        if not reason:
            raise ValueError("targeted-web approval has a blank review_reason")
        path, local_hash, doi = _validate_source_provenance(
            row,
            approval_row=approval_row,
            candidate_csv=candidate_csv,
        )
        if not _excerpt_is_exact(excerpt, _literal(row["raw_text"]), path.read_bytes()):
            raise ValueError("targeted-web approval excerpt is not exact source text")
        if not _exact_species_named(species, excerpt):
            raise ValueError(
                "targeted-web approval excerpt must name the exact full accepted species"
            )
        if not _explicit_mapping_supported(trait, value, excerpt):
            raise ValueError(
                f"targeted-web approval excerpt does not explicitly support {trait}={value}"
            )

        source_url = _text(row["url"])
        source_type = _text(row["source_type"])
        citation = _source_citation(row, doi)
        source_text_hash = _sha256_text(excerpt)
        identity = "\x1f".join(
            (
                IMPORT_CONTRACT,
                expected_discovery_queue_hash,
                current_queue_hash,
                candidate_hash,
                str(row_index),
                species,
                trait,
                value,
                source_url,
                local_hash,
                source_text_hash,
            )
        )
        candidate_id = _sha256_text(f"candidate\x1f{identity}")
        source_evidence_id = _sha256_text(f"source\x1f{identity}")
        source_rows.append(
            {
                "source_evidence_id": source_evidence_id,
                "accepted_species": species,
                "source_text": excerpt,
                "evidence_excerpt": excerpt,
                "source_url": source_url,
                "retrieval_url": _text(row["retrieval_url"]),
                "source_citation": citation,
                "source_type": source_type,
                "provider": PROVIDER,
                "evidence_scope": "species_direct",
                "DOI": doi,
                "publication": _text(row["publication"]),
                "venue": _text(row["venue"]),
                "year": _text(row["year"]),
                "local_file_path": str(path),
                "local_file_hash": local_hash,
                "retrieval_date": _text(row["retrieval_date"]),
                "source_run_id": approval_id,
                "source_text_hash": source_text_hash,
                "discovery_row_index": row_index,
                "discovery_file_hash": candidate_hash,
            }
        )
        candidate_rows.append(
            {
                "candidate_id": candidate_id,
                "accepted_species": species,
                "trait_name": trait,
                "provisional_candidate_value": value,
                "matched_term": "explicit_targeted_web_approval",
                "candidate_class": "reported",
                "provider": PROVIDER,
                "source_type": source_type,
                "source_url": source_url,
                "source_citation": citation,
                "raw_description": excerpt,
                "source_text": excerpt,
                "evidence_scope": "species_direct",
                "wild_or_cultivated": "unknown",
                "review_status": "strictly_reviewed",
                "provider_status": "provider_pass_complete_with_candidates",
                "provider_attempts": "1",
                "provider_updated_at": completed_at,
                "provider_policy_version": PROVIDER_VERSION,
                "local_file_path": str(path),
                "local_file_hash": local_hash,
                "retrieval_date": _text(row["retrieval_date"]),
                "retrieval_date_basis": "targeted_web_discovery_record",
                "source_run_id": approval_id,
                "audit_status": "strict_review_candidate",
                "audit_reason": "explicit_approval_bound_to_exact_local_raw_hash",
            }
        )
        decisions.append(
            {
                "candidate_id": candidate_id,
                "accepted_trait": trait,
                "accepted_value": value,
                "review_reason": reason,
                "evidence_excerpt": excerpt,
            }
        )

    source_frame = pd.DataFrame(source_rows, columns=SOURCE_EVIDENCE_COLUMNS).fillna("")
    candidate_frame = pd.DataFrame(candidate_rows, columns=REVIEWED_CANDIDATE_COLUMNS).fillna("")
    if source_frame["source_evidence_id"].duplicated().any():
        raise ValueError("targeted-web source evidence contains duplicate IDs")
    if candidate_frame["candidate_id"].duplicated().any():
        raise ValueError("targeted-web reviewed candidates contain duplicate IDs")
    source_frame = source_frame.sort_values(
        ["accepted_species", "discovery_row_index", "source_evidence_id"]
    ).reset_index(drop=True)
    candidate_frame = candidate_frame.sort_values(
        ["accepted_species", "trait_name", "candidate_id"]
    ).reset_index(drop=True)
    review_payload: dict[str, object] = {
        "contract_version": APPROVAL_CONTRACT,
        "approval_id": approval_id,
        "completed_at_utc": completed_at,
        "discovery_queue_sha256": expected_discovery_queue_hash,
        "current_queue_sha256": current_queue_hash,
        "candidate_file_sha256": candidate_hash,
        "rejected_file_sha256": rejected_hash,
        "accepted": decisions,
    }
    lane_metadata: dict[str, object] = {
        "approval_id": approval_id,
        "completed_at_utc": completed_at,
        "inspected_species": list(inspected_species),
        "held_or_rejected_species": sorted(rejected_species),
        "discovery_queue_sha256": expected_discovery_queue_hash,
        "current_queue_sha256": current_queue_hash,
        "candidate_file_sha256": candidate_hash,
        "rejected_file_sha256": rejected_hash,
    }
    return source_frame, candidate_frame, review_payload, lane_metadata


def _targeted_checkpoint_rows(
    candidates: pd.DataFrame,
    review_payload: dict[str, object],
    *,
    inspected_species: tuple[str, ...],
    already_covered_species: frozenset[str],
) -> pd.DataFrame:
    counts = (
        candidates.groupby(["accepted_species", "trait_name"]).size().astype(int).to_dict()
        if not candidates.empty
        else {}
    )
    completed_at = _text(review_payload["completed_at_utc"])
    approval_id = _text(review_payload["approval_id"])
    rows: list[dict[str, object]] = []
    for species in inspected_species:
        for trait in RA_TRAITS:
            count = int(counts.get((species, trait), 0))
            if count:
                status = "accepted_evidence"
            elif species in already_covered_species:
                status = "provider_seed_covered"
            else:
                status = "provider_pass_complete"
            rows.append(
                {
                    "species": species,
                    "trait": trait,
                    "provider": PROVIDER,
                    "provider_version": PROVIDER_VERSION,
                    "status": status,
                    "terminal": True,
                    "attempts": 1,
                    "candidate_count": count,
                    "accepted_record_count": count,
                    "legacy_status": "",
                    "legacy_enabled": "false",
                    "last_error": "",
                    "updated_at": completed_at,
                    "last_wave_id": "",
                    "source_run_id": approval_id,
                    "migration_basis": "|".join(
                        basis
                        for basis in (
                            "targeted_web_exact_20_species_explicit_approval",
                            REVIEW_BASIS if count else "",
                        )
                        if basis
                    ),
                    "contract_version": CONTRACT_VERSION,
                }
            )
    return pd.DataFrame(rows, columns=CHECKPOINT_COLUMNS)


def _already_covered_inspected_species(
    *,
    queue: pd.DataFrame,
    baseline_checkpoint: pd.DataFrame,
    candidates: pd.DataFrame,
    inspected_species: tuple[str, ...],
) -> frozenset[str]:
    """Allow queue departures only when the latest baseline proves direct coverage."""
    _require_columns(
        baseline_checkpoint,
        {"accepted_species", "reproductive_assurance_direct"},
        "three-axis checkpoint",
    )
    queue_species = set(_validate_queue(queue))
    baseline = baseline_checkpoint.copy().fillna("")
    baseline["accepted_species"] = baseline["accepted_species"].map(_text)
    if baseline["accepted_species"].eq("").any():
        raise ValueError("three-axis checkpoint contains a blank species")
    if baseline["accepted_species"].duplicated().any():
        raise ValueError("three-axis checkpoint contains duplicate species")
    baseline_direct = set(
        baseline.loc[baseline["reproductive_assurance_direct"].map(_bool), "accepted_species"]
    )
    inspected = set(inspected_species)
    missing_baseline = sorted(inspected.difference(set(baseline["accepted_species"])))
    if missing_baseline:
        raise ValueError(
            "targeted-web inspected species are outside the latest three-axis checkpoint: "
            f"{missing_baseline}"
        )
    stale_queue = sorted(inspected.intersection(queue_species, baseline_direct))
    if stale_queue:
        raise ValueError(
            f"targeted-web current queue still contains already-covered species: {stale_queue}"
        )
    outside_queue = inspected.difference(queue_species)
    unsupported = sorted(outside_queue.difference(baseline_direct))
    if unsupported:
        raise ValueError(
            "targeted-web inspected species are outside the current queue without "
            f"direct baseline coverage: {unsupported}"
        )
    approved_species = set(candidates["accepted_species"].map(_text))
    approved_outside = sorted(approved_species.intersection(outside_queue))
    if approved_outside:
        raise ValueError(
            f"targeted-web approved species are outside the current queue: {approved_outside}"
        )
    return frozenset(outside_queue)


def _canonical_checkpoint(frame: pd.DataFrame) -> pd.DataFrame:
    _require_columns(frame, set(CHECKPOINT_COLUMNS), "provider checkpoint")
    result = frame[CHECKPOINT_COLUMNS].copy().fillna("")
    if result.duplicated(["species", "trait", "provider"]).any():
        raise ValueError("provider checkpoint contains duplicate keys")
    for column in CHECKPOINT_COLUMNS:
        if column == "terminal":
            result[column] = result[column].map(lambda value: str(_bool(value)).lower())
        elif column in {"attempts", "candidate_count", "accepted_record_count"}:
            result[column] = result[column].map(lambda value: str(_nonnegative_int(value)))
        else:
            result[column] = result[column].map(_text)
    return result.sort_values(["species", "trait", "provider"]).reset_index(drop=True)


def _require_complete_targeted_blocks(frame: pd.DataFrame, *, label: str) -> None:
    """Require exactly one row for every RA trait in each targeted-web species block."""
    if frame.empty:
        return
    expected_traits = set(RA_TRAITS)
    for species, block in frame.groupby("species", sort=False):
        observed_traits = set(block["trait"])
        if len(block) != len(RA_TRAITS) or observed_traits != expected_traits:
            raise ValueError(f"{label} has an incomplete targeted-web species block: {species}")


def _seed_targeted_checkpoint(
    provider_checkpoint: pd.DataFrame,
    planned: pd.DataFrame,
    *,
    allow_append: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    base = _canonical_checkpoint(provider_checkpoint)
    planned_canonical = _canonical_checkpoint(planned)
    existing = base.loc[base["provider"].eq(PROVIDER)].copy()
    _require_complete_targeted_blocks(planned_canonical, label="planned checkpoint")
    if existing.empty:
        seeded = pd.concat([base, planned_canonical], ignore_index=True)
        return seeded[CHECKPOINT_COLUMNS], existing, int(len(planned_canonical))
    _require_complete_targeted_blocks(existing, label="existing checkpoint")
    if (
        not existing["terminal"].eq("true").all()
        or not existing["status"].isin(TERMINAL_STATUSES).all()
    ):
        raise ValueError("existing targeted-web checkpoint contains nonterminal state")
    existing_keys = set(existing[["species", "trait"]].itertuples(index=False, name=None))
    planned_keys = set(planned_canonical[["species", "trait"]].itertuples(index=False, name=None))
    if not allow_append and existing_keys != planned_keys:
        raise ValueError("existing targeted-web checkpoint does not have the exact 20x4 key set")
    overlap_keys = existing_keys.intersection(planned_keys)
    if overlap_keys:
        existing_overlap = existing.loc[
            existing.set_index(["species", "trait"]).index.isin(overlap_keys)
        ]
        planned_overlap = planned_canonical.loc[
            planned_canonical.set_index(["species", "trait"]).index.isin(overlap_keys)
        ]
        if not _canonical_checkpoint(existing_overlap).equals(
            _canonical_checkpoint(planned_overlap)
        ):
            message = (
                "overlapping targeted-web checkpoint rows are not identical"
                if allow_append
                else "existing terminal targeted-web checkpoint rows are immutable"
            )
            raise ValueError(message)
    new_rows = planned_canonical.loc[
        ~planned_canonical.set_index(["species", "trait"]).index.isin(existing_keys)
    ].copy()
    seeded = pd.concat([base, new_rows], ignore_index=True)
    return seeded[CHECKPOINT_COLUMNS], existing, int(len(new_rows))


def _assert_terminal_immutable(existing: pd.DataFrame, final: pd.DataFrame) -> None:
    if existing.empty:
        return
    final_targeted = final.loc[final["provider"].map(_text).eq(PROVIDER)]
    left = _canonical_checkpoint(existing)
    existing_keys = set(left[["species", "trait"]].itertuples(index=False, name=None))
    right = _canonical_checkpoint(
        final_targeted.loc[final_targeted.set_index(["species", "trait"]).index.isin(existing_keys)]
    )
    if not left.equals(right):
        raise ValueError("existing terminal targeted-web checkpoint rows are immutable")


def _revalidate_source_files(source_evidence: pd.DataFrame) -> None:
    for row in source_evidence.to_dict("records"):
        path = Path(_text(row["local_file_path"]))
        expected = _text(row["local_file_hash"]).casefold()
        if not path.is_file() or _file_sha256(path) != expected:
            raise ValueError(f"targeted-web final raw-response provenance mismatch: {path}")


def materialize_targeted_web_lane(
    *,
    queue: pd.DataFrame,
    baseline_checkpoint: pd.DataFrame,
    provider_checkpoint: pd.DataFrame,
    source_evidence: pd.DataFrame,
    candidates: pd.DataFrame,
    review_payload: dict[str, object],
    prior_records: pd.DataFrame | None = None,
    inspected_species: tuple[str, ...] = INSPECTED_SPECIES,
    allow_checkpoint_append: bool = False,
) -> dict[str, object]:
    """Materialize approved rows and one terminal 4-trait provider block per inspected species."""
    already_covered = _already_covered_inspected_species(
        queue=queue,
        baseline_checkpoint=baseline_checkpoint,
        candidates=candidates,
        inspected_species=inspected_species,
    )
    planned = _targeted_checkpoint_rows(
        candidates,
        review_payload,
        inspected_species=inspected_species,
        already_covered_species=already_covered,
    )
    seeded, existing_targeted, appended_rows = _seed_targeted_checkpoint(
        provider_checkpoint,
        planned,
        allow_append=allow_checkpoint_append,
    )
    _revalidate_source_files(source_evidence)
    result = materialize_reviewed_cache(
        queue=queue,
        baseline_checkpoint=baseline_checkpoint,
        provider_checkpoint=seeded,
        candidates=candidates,
        review_payload=review_payload,
        prior_records=prior_records,
    )
    _revalidate_source_files(source_evidence)
    final_targeted = result["provider_checkpoint"].loc[
        result["provider_checkpoint"]["provider"].map(_text).eq(PROVIDER)
    ]
    round_rows = len(inspected_species) * len(RA_TRAITS)
    expected_total_rows = len(existing_targeted) + appended_rows
    if len(final_targeted) != expected_total_rows:
        raise AssertionError(
            f"targeted-web checkpoint rows={len(final_targeted)}, expected={expected_total_rows}"
        )
    if not final_targeted["terminal"].map(_bool).all():
        raise AssertionError("targeted-web checkpoint contains nonterminal rows")
    existing_species = set(existing_targeted["species"].map(_text))
    expected_species = existing_species.union(inspected_species)
    if set(final_targeted["species"].map(_text)) != expected_species:
        raise AssertionError("targeted-web checkpoint escaped the cumulative inspected set")
    _assert_terminal_immutable(existing_targeted, result["provider_checkpoint"])
    result["source_evidence"] = source_evidence
    result["reviewed_candidates"] = candidates
    result["review_payload"] = review_payload
    result["summary"].update(
        {
            "external_acquisition_started": True,
            # The legacy field remains the total targeted-web rows in the output.
            # Round/appended counts are explicit so later rounds cannot be
            # mistaken for replacement of prior terminal state.
            "targeted_web_provider_rows": expected_total_rows,
            "targeted_web_round_provider_rows": round_rows,
            "targeted_web_existing_provider_rows": int(len(existing_targeted)),
            "targeted_web_appended_provider_rows": appended_rows,
            "targeted_web_total_provider_rows": expected_total_rows,
            "targeted_web_inspected_species": len(inspected_species),
            "targeted_web_round_inspected_species": len(inspected_species),
            "targeted_web_appended_species": appended_rows // len(RA_TRAITS),
            "targeted_web_total_inspected_species": len(expected_species),
            "targeted_web_approved_sources": int(len(source_evidence)),
            "targeted_web_reviewed_candidates": int(len(candidates)),
            "targeted_web_already_covered_species": sorted(already_covered),
            "targeted_web_all_terminal": True,
            "targeted_web_global_fallback_used": False,
            "targeted_web_discovery_auto_accepted": False,
            "targeted_web_import_performed_external_acquisition": False,
        }
    )
    return result


def run_targeted_web_import(
    *,
    queue_csv: Path,
    candidate_csv: Path,
    rejected_csv: Path,
    approval_json: Path,
    baseline_checkpoint_csv: Path,
    provider_checkpoint_csv: Path,
    output_dir: Path,
    expected_current_queue_hash: str,
    prior_records_csv: Path | None = None,
    expected_discovery_queue_hash: str = DISCOVERY_QUEUE_SHA256,
    inspected_species: tuple[str, ...] | None = INSPECTED_SPECIES,
    allow_checkpoint_append: bool = False,
) -> dict[str, object]:
    """Validate, materialize, and write checkpoint state after every other output."""
    input_paths = (
        queue_csv,
        candidate_csv,
        rejected_csv,
        approval_json,
        baseline_checkpoint_csv,
        provider_checkpoint_csv,
        *(() if prior_records_csv is None else (prior_records_csv,)),
    )
    input_metadata = {
        str(path.resolve()): {
            "bytes": path.stat().st_size,
            "sha256": _file_sha256(path),
        }
        for path in input_paths
    }
    queue_hash = _file_sha256(queue_csv)
    candidate_hash = _file_sha256(candidate_csv)
    rejected_hash = _file_sha256(rejected_csv)
    queue = pd.read_csv(queue_csv, dtype=str).fillna("")
    candidate_records = pd.read_csv(candidate_csv, dtype=str).fillna("")
    rejected_hits = pd.read_csv(rejected_csv, dtype=str).fillna("")
    approval_payload = json.loads(approval_json.read_text(encoding="utf-8"))
    if inspected_species is None:
        inspected_species = _dynamic_inspected_species(
            candidate_records,
            rejected_hits,
            approval_payload,
        )
    source_evidence, candidates, review_payload, lane_metadata = build_targeted_web_review_batch(
        queue=queue,
        candidate_records=candidate_records,
        rejected_hits=rejected_hits,
        approval_payload=approval_payload,
        current_queue_hash=queue_hash,
        candidate_hash=candidate_hash,
        rejected_hash=rejected_hash,
        candidate_csv=candidate_csv,
        expected_current_queue_hash=expected_current_queue_hash,
        expected_discovery_queue_hash=expected_discovery_queue_hash,
        inspected_species=inspected_species,
    )
    baseline = pd.read_csv(baseline_checkpoint_csv, dtype=str).fillna("")
    providers = pd.read_csv(provider_checkpoint_csv, dtype=str).fillna("")
    prior = (
        pd.read_csv(prior_records_csv, dtype=str).fillna("")
        if prior_records_csv is not None
        else None
    )
    result = materialize_targeted_web_lane(
        queue=queue,
        baseline_checkpoint=baseline,
        provider_checkpoint=providers,
        source_evidence=source_evidence,
        candidates=candidates,
        review_payload=review_payload,
        prior_records=prior,
        inspected_species=inspected_species,
        allow_checkpoint_append=allow_checkpoint_append,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "source_evidence": output_dir / "source_evidence.csv",
        "reviewed_candidates": output_dir / "reviewed_candidates.csv",
        "review": output_dir / "targeted_web_review.json",
        "records": output_dir / "new_reproductive_assurance_records.csv",
        "axes": output_dir / "three_axis_checkpoint.csv.gz",
        "audit": output_dir / "targeted_web_review_audit.csv",
        "next_queue": output_dir / "reproductive_assurance_priority_next.csv",
        "summary": output_dir / "coverage_summary.json",
        "provider_checkpoint": output_dir / "provider_checkpoint.csv",
    }
    _atomic_write_csv(result["source_evidence"], paths["source_evidence"])
    _atomic_write_csv(result["reviewed_candidates"], paths["reviewed_candidates"])
    _atomic_write_text(
        json.dumps(result["review_payload"], ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        paths["review"],
    )
    _atomic_write_csv(result["records"], paths["records"])
    _atomic_write_csv(result["three_axis_checkpoint"], paths["axes"], gzip=True)
    _atomic_write_csv(result["review_audit"], paths["audit"])
    _atomic_write_csv(result["next_queue"], paths["next_queue"])

    before_checkpoint = [
        paths["source_evidence"],
        paths["reviewed_candidates"],
        paths["review"],
        paths["records"],
        paths["axes"],
        paths["audit"],
        paths["next_queue"],
    ]
    summary = dict(result["summary"])
    summary.update(
        {
            "contract_version": IMPORT_CONTRACT,
            "materialized_at_utc": lane_metadata["completed_at_utc"],
            "lane": lane_metadata,
            "checkpoint_written_last": True,
            "outputs_written_before_terminal_checkpoint": [path.name for path in before_checkpoint],
            "inputs": input_metadata,
            "outputs": {
                path.name: {
                    "path": str(path.resolve()),
                    "bytes": path.stat().st_size,
                    "sha256": _file_sha256(path),
                }
                for path in before_checkpoint
            },
        }
    )
    _atomic_write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        paths["summary"],
    )

    # This is deliberately the final persistent state transition.
    _revalidate_source_files(result["source_evidence"])
    _atomic_write_csv(result["provider_checkpoint"], paths["provider_checkpoint"])
    summary["provider_checkpoint_runtime"] = {
        "path": str(paths["provider_checkpoint"].resolve()),
        "bytes": paths["provider_checkpoint"].stat().st_size,
        "sha256": _file_sha256(paths["provider_checkpoint"]),
    }
    return summary


def run_targeted_web_round_import(
    *,
    queue_csv: Path,
    candidate_csv: Path,
    rejected_csv: Path,
    approval_json: Path,
    baseline_checkpoint_csv: Path,
    provider_checkpoint_csv: Path,
    output_dir: Path,
    expected_current_queue_hash: str,
    expected_discovery_queue_hash: str,
    prior_records_csv: Path,
) -> dict[str, object]:
    """Append one dynamic explicit-approval round without performing acquisition."""
    return run_targeted_web_import(
        queue_csv=queue_csv,
        candidate_csv=candidate_csv,
        rejected_csv=rejected_csv,
        approval_json=approval_json,
        baseline_checkpoint_csv=baseline_checkpoint_csv,
        provider_checkpoint_csv=provider_checkpoint_csv,
        output_dir=output_dir,
        expected_current_queue_hash=expected_current_queue_hash,
        prior_records_csv=prior_records_csv,
        expected_discovery_queue_hash=expected_discovery_queue_hash,
        inspected_species=None,
        allow_checkpoint_append=True,
    )


@app.command("import-reviewed")
def import_reviewed_command(
    queue_csv: Path = typer.Option(..., exists=True),
    candidate_csv: Path = typer.Option(..., exists=True),
    rejected_csv: Path = typer.Option(..., exists=True),
    approval_json: Path = typer.Option(..., exists=True),
    baseline_checkpoint_csv: Path = typer.Option(..., exists=True),
    provider_checkpoint_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    expected_current_queue_hash: str = typer.Option(...),
    prior_records_csv: Path | None = typer.Option(None),
    expected_discovery_queue_hash: str = typer.Option(DISCOVERY_QUEUE_SHA256),
) -> None:
    """Import only explicit approvals; perform no web requests or global fallback."""
    report = run_targeted_web_import(
        queue_csv=queue_csv,
        candidate_csv=candidate_csv,
        rejected_csv=rejected_csv,
        approval_json=approval_json,
        baseline_checkpoint_csv=baseline_checkpoint_csv,
        provider_checkpoint_csv=provider_checkpoint_csv,
        output_dir=output_dir,
        expected_current_queue_hash=expected_current_queue_hash,
        prior_records_csv=prior_records_csv,
        expected_discovery_queue_hash=expected_discovery_queue_hash,
    )
    typer.echo(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))


@app.command("import-reviewed-round")
def import_reviewed_round_command(
    queue_csv: Path = typer.Option(..., exists=True),
    candidate_csv: Path = typer.Option(..., exists=True),
    rejected_csv: Path = typer.Option(..., exists=True),
    approval_json: Path = typer.Option(..., exists=True),
    baseline_checkpoint_csv: Path = typer.Option(..., exists=True),
    provider_checkpoint_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    expected_current_queue_hash: str = typer.Option(...),
    expected_discovery_queue_hash: str = typer.Option(...),
    prior_records_csv: Path = typer.Option(..., exists=True),
) -> None:
    """Append one exact approved/discovered species round; perform no web requests."""
    report = run_targeted_web_round_import(
        queue_csv=queue_csv,
        candidate_csv=candidate_csv,
        rejected_csv=rejected_csv,
        approval_json=approval_json,
        baseline_checkpoint_csv=baseline_checkpoint_csv,
        provider_checkpoint_csv=provider_checkpoint_csv,
        output_dir=output_dir,
        expected_current_queue_hash=expected_current_queue_hash,
        expected_discovery_queue_hash=expected_discovery_queue_hash,
        prior_records_csv=prior_records_csv,
    )
    typer.echo(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
