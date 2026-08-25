"""Bootstrap a local, trait-scoped reproductive-assurance provider ledger.

The public-web v6 campaign persisted provider work at ``species x provider``.
Local continuation needs the narrower ``species x trait x provider`` key so a
completed lookup is never repeated and, equally importantly, completion is not
invented for traits a provider did not search.  This module performs only that
auditable migration.  It does not collect evidence or assign biological values.
"""

from __future__ import annotations

import hashlib
import csv
import json
import os
import re
from datetime import UTC, datetime
from pathlib import Path
from tempfile import NamedTemporaryFile

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT_VERSION = "local_reproductive_assurance_provider_checkpoint_v2"
RA_TRAITS = (
    "self_incompatibility",
    "autonomous_selfing_capacity",
    "mating_system",
    "cleistogamy",
)
PROVIDERS = (
    "gbif",
    "world_flora",
    "web_descriptions",
    "wikimedia",
    "openalex",
    "europe_pmc",
)
RECORD_PROVIDERS = frozenset(
    (
        *PROVIDERS,
        "crossref",
        "semantic_scholar",
        "doaj",
        "openaire",
        "pladias",
        "dryad_biogeo_mating_2018",
        "targeted_web",
    )
)
EXTERNAL_ACQUISITION_PROVIDERS = frozenset(
    {
        "crossref",
        "semantic_scholar",
        "doaj",
        "openaire",
        "pladias",
        "dryad_biogeo_mating_2018",
        "targeted_web",
    }
)

DRYAD_BIOGEO_PROVIDER = "dryad_biogeo_mating_2018"
DRYAD_BIOGEO_SOURCE_TYPE = "dryad_biogeo_mating_si_category"
DRYAD_BIOGEO_CSV_SHA256 = "94c427b84dc35c12f5f35cea0ee7aded7c96b9de16231f2efd771bc0035a0902"
DRYAD_BIOGEO_SOURCE_COLUMNS = (
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

TARGETED_WEB_SOURCE_TYPES = frozenset(
    {
        "journal_article_primary_pollination_experiment",
        "journal_article_comparative_breeding_system",
        "flora_treatment_explicit_mating_system",
        "journal_article_primary_breeding_system_experiment",
        "journal_article_primary_population_genetics",
        "dissertation_primary_pollination_experiment",
    }
)

# These are provider-query scopes, not biological claims. Europe PMC v2 used
# one exact-species query containing selfing, outcrossing, autogamy, SI/SC, and
# cleistogamy terms, so an actual v2 request accounts for all four RA traits.
PROVIDER_QUERY_TRAITS = {
    "gbif": frozenset(),
    "world_flora": frozenset({"self_incompatibility", "mating_system"}),
    "web_descriptions": frozenset({"self_incompatibility", "mating_system"}),
    "wikimedia": frozenset({"self_incompatibility", "mating_system"}),
    "openalex": frozenset(RA_TRAITS),
    "europe_pmc": frozenset(RA_TRAITS),
}
OLD_STATIC_TARGET_TRAITS = {
    "europe_pmc": frozenset({"self_incompatibility", "mating_system"}),
}
TERMINAL_STATUSES = frozenset(
    {
        "accepted_evidence",
        "provider_pass_complete",
        "provider_pass_complete_with_candidates",
        "provider_no_result",
        "provider_seed_covered",
        "permanent_failure",
        "invalid_species_name",
        "not_applicable",
    }
)
RETRYABLE_STATUSES = frozenset({"pending", "retry"})
ALL_STATUSES = TERMINAL_STATUSES | RETRYABLE_STATUSES | {"running"}

CHECKPOINT_COLUMNS = [
    "species",
    "trait",
    "provider",
    "provider_version",
    "status",
    "terminal",
    "attempts",
    "candidate_count",
    "accepted_record_count",
    "legacy_status",
    "legacy_enabled",
    "last_error",
    "updated_at",
    "last_wave_id",
    "source_run_id",
    "migration_basis",
    "contract_version",
]

ALLOWED_DIRECT_VALUES = {
    "self_incompatibility": frozenset({"SI", "SC", "mixed_or_variable"}),
    "autonomous_selfing_capacity": frozenset(
        {"absent", "delayed", "facilitated", "autonomous", "mixed_or_variable"}
    ),
    "mating_system": frozenset(
        {
            "predominantly_outcrossing",
            "mixed_mating",
            "predominantly_selfing",
            "obligate_selfing",
        }
    ),
    "cleistogamy": frozenset({"absent", "facultative", "obligate"}),
}

RECORD_COLUMNS = [
    "record_id",
    "species",
    "trait",
    "value",
    "raw_text",
    "evidence_excerpt",
    "URL",
    "DOI",
    "publication",
    "year",
    "source_type",
    "provider",
    "evidence_scope",
    "local_file_path",
    "local_file_hash",
    "source_text_hash",
    "retrieval_date",
    "candidate_id",
    "source_run_id",
    "direct_evidence",
    "review_basis",
]

REQUIRED_CANDIDATE_COLUMNS = {
    "candidate_id",
    "accepted_species",
    "trait_name",
    "provisional_candidate_value",
    "provider",
    "source_type",
    "source_url",
    "source_citation",
    "raw_description",
    "source_text",
    "evidence_scope",
    "local_file_path",
    "local_file_hash",
    "retrieval_date",
    "source_run_id",
}

REVIEW_BASIS = "strict_species_direct_source_text_review_v1"

SOURCE_PROVIDER_MAP = {
    "europe_pmc_title_abstract": "europe_pmc",
    "crossref_title_abstract": "crossref",
    "semantic_scholar_title_abstract": "semantic_scholar",
    "doaj_title_abstract": "doaj",
    "openaire_title_description": "openaire",
    "pladias_generative_reproduction_type": "pladias",
    DRYAD_BIOGEO_SOURCE_TYPE: DRYAD_BIOGEO_PROVIDER,
    "openalex_title_abstract": "openalex",
    "gbif_species_description": "gbif",
    "wikipedia_extract": "wikimedia",
    "wikipedia_sitelink_extract_ja": "wikimedia",
    "horticulture_site": "web_descriptions",
    **{source_type: "targeted_web" for source_type in TARGETED_WEB_SOURCE_TYPES},
}


@app.callback()
def main() -> None:
    """Build and advance the local reproductive-assurance campaign."""


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _atomic_write_csv(frame: pd.DataFrame, path: Path, *, gzip: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
    try:
        compression: str | dict[str, object] | None = (
            {"method": "gzip", "mtime": 0} if gzip else None
        )
        frame.to_csv(temporary, index=False, compression=compression)
        with temporary.open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_write_text(value: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        handle.write(value)
        handle.flush()
        os.fsync(handle.fileno())
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _nonnegative_int(value: object) -> int:
    parsed = pd.to_numeric(value, errors="coerce")
    if pd.isna(parsed):
        return 0
    integer = int(parsed)
    if integer < 0:
        raise ValueError(f"expected a non-negative integer, got {value!r}")
    return integer


def _bool(value: object) -> bool:
    return _text(value).casefold() in {"1", "true", "t", "yes", "y"}


def _citation_fields(citation: object) -> tuple[str, str, str]:
    publication = _text(citation).split(" | ", maxsplit=1)[0].strip()
    year_match = re.search(r"\((\d{4})\)", publication)
    doi_match = re.search(r"\bDOI:([^;|\s]+)", publication, flags=re.IGNORECASE)
    return (
        publication,
        year_match.group(1) if year_match else "",
        doi_match.group(1).rstrip(".,").casefold() if doi_match else "",
    )


def _source_provider(source_type: object, source_url: object) -> str:
    kind = _text(source_type).casefold()
    url = _text(source_url).casefold()
    if "worldfloraonline.org" in url:
        return "world_flora"
    if kind == "flora_or_monograph":
        return "web_descriptions"
    return SOURCE_PROVIDER_MAP.get(kind, "")


def _match_excerpt(text: str, match: re.Match[str], pad: int = 500) -> str:
    start = max(0, match.start() - pad)
    end = min(len(text), match.end() + pad)
    return _text(text[start:end])


def augment_cleistogamy_cache_candidates(
    *,
    queue: pd.DataFrame,
    search_evidence: pd.DataFrame,
    legacy_provider: pd.DataFrame,
    existing_candidates: pd.DataFrame,
    local_file_path: str,
    local_file_hash: str,
    source_run_id: int = 29637034874,
) -> pd.DataFrame:
    """Add review-only cleistogamy candidates missed by the legacy extractor."""
    queue_species = set(_validate_queue(queue))
    required_search = {
        "accepted_species",
        "source_text",
        "source_url",
        "source_citation",
        "source_type",
        "evidence_scope",
    }
    missing = required_search.difference(search_evidence.columns)
    if missing:
        raise ValueError(f"search evidence missing columns: {sorted(missing)}")
    if not REQUIRED_CANDIDATE_COLUMNS.issubset(existing_candidates.columns):
        missing_candidates = REQUIRED_CANDIDATE_COLUMNS.difference(existing_candidates.columns)
        raise ValueError(f"existing candidate table missing columns: {sorted(missing_candidates)}")

    legacy = legacy_provider.copy().fillna("")
    required_legacy = {"species", "provider", "status", "attempts", "updated_at", "policy_version"}
    missing_legacy = required_legacy.difference(legacy.columns)
    if missing_legacy:
        raise ValueError(f"legacy provider checkpoint missing columns: {sorted(missing_legacy)}")
    legacy["species"] = legacy["species"].map(_text)
    legacy["provider"] = legacy["provider"].map(_text)
    if legacy.duplicated(["species", "provider"]).any():
        raise ValueError("legacy provider checkpoint has duplicate keys")
    legacy = legacy.set_index(["species", "provider"], drop=False)

    rows: list[dict[str, str]] = []
    evidence = search_evidence.copy().fillna("")
    evidence = evidence.loc[evidence["accepted_species"].map(_text).isin(queue_species)]
    for source in evidence.to_dict("records"):
        species = _text(source["accepted_species"])
        source_text = _text(source["source_text"])
        match = re.search(r"\bcleistogam\w*\b", source_text, flags=re.IGNORECASE)
        if match is None:
            continue
        provider = _source_provider(source["source_type"], source["source_url"])
        if not provider or (species, provider) not in legacy.index:
            continue
        folded = source_text.casefold()
        if re.search(r"\b(?:no|not|without|lacks?)\s+cleistogam\w*\b", folded):
            value = "absent"
        elif re.search(r"\bobligate\s+cleistogam\w*\b", folded):
            value = "obligate"
        elif (
            re.search(r"\bfacultative\s+cleistogam\w*\b", folded)
            or "chasmogam" in folded
            or "open flower" in folded
            or re.search(r"\b(?:often|usually|sometimes)\s+cleistogam\w*\b", folded)
        ):
            value = "facultative"
        else:
            value = "unresolved"
        provider_row = legacy.loc[(species, provider)]
        retrieval = _text(provider_row.get("updated_at"))
        citation = _text(source["source_citation"])
        excerpt = _match_excerpt(source_text, match)
        candidate_id = _sha256_text(
            "\x1f".join(
                (
                    "cache_cleistogamy_v1",
                    species,
                    provider,
                    _text(source["source_url"]),
                    _sha256_text(source_text),
                )
            )
        )
        rows.append(
            {
                "candidate_id": candidate_id,
                "accepted_species": species,
                "trait_name": "cleistogamy",
                "provisional_candidate_value": value,
                "matched_term": match.group(0),
                "candidate_class": "reported",
                "provider": provider,
                "source_type": _text(source["source_type"]),
                "source_url": _text(source["source_url"]),
                "source_citation": citation,
                "raw_description": excerpt,
                "source_text": source_text,
                "evidence_scope": _text(source["evidence_scope"]),
                "wild_or_cultivated": "unknown",
                "review_status": "unreviewed",
                "provider_status": _text(provider_row.get("status")),
                "provider_attempts": str(_nonnegative_int(provider_row.get("attempts"))),
                "provider_updated_at": retrieval,
                "provider_policy_version": _text(provider_row.get("policy_version")),
                "local_file_path": local_file_path,
                "local_file_hash": local_file_hash,
                "retrieval_date": retrieval[:10],
                "retrieval_date_basis": "provider_checkpoint.updated_at",
                "source_run_id": str(source_run_id),
                "audit_status": "review_candidate",
                "audit_reason": "explicit_cleistogamy_term_in_species_direct_cached_source",
            }
        )
    augmented = pd.concat(
        [existing_candidates.copy().fillna(""), pd.DataFrame(rows)],
        ignore_index=True,
        sort=False,
    ).fillna("")
    if augmented["candidate_id"].duplicated().any():
        augmented = augmented.drop_duplicates("candidate_id", keep="first")
    return augmented.sort_values(
        ["accepted_species", "trait_name", "provider", "candidate_id"]
    ).reset_index(drop=True)


def _search_selection_rows(payload: object) -> list[dict[str, str | int]]:
    if not isinstance(payload, dict):
        raise ValueError("search selection file must contain a JSON object")
    accepted = payload.get("accepted", [])
    if not isinstance(accepted, list) or not all(isinstance(row, dict) for row in accepted):
        raise ValueError("search selection accepted must be a list of objects")
    rows: list[dict[str, str | int]] = []
    for row in accepted:
        try:
            source_row_index = int(str(row.get("source_row_index", "")).strip())
        except ValueError as exc:
            raise ValueError("search selection has an invalid source_row_index") from exc
        if source_row_index < 0:
            raise ValueError("search selection source_row_index must be non-negative")
        rows.append(
            {
                "source_row_index": source_row_index,
                "accepted_species": _text(row.get("accepted_species")),
                "source_url": _text(row.get("source_url")),
                "accepted_trait": _text(row.get("accepted_trait")),
                "accepted_value": _text(row.get("accepted_value")),
                "evidence_excerpt": _text(row.get("evidence_excerpt")),
                "review_reason": _text(row.get("review_reason")),
            }
        )
    return rows


def build_reviewed_search_batch(
    *,
    queue: pd.DataFrame,
    search_evidence: pd.DataFrame,
    provider_checkpoint: pd.DataFrame,
    selection_payload: object,
    local_file_path: str,
    local_file_hash: str,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Reconstruct reviewed cache candidates from stable source-row selections."""
    queue_species = set(_validate_queue(queue))
    required_search = {
        "accepted_species",
        "source_text",
        "source_url",
        "source_citation",
        "source_type",
        "evidence_scope",
    }
    missing_search = required_search.difference(search_evidence.columns)
    if missing_search:
        raise ValueError(f"search evidence missing columns: {sorted(missing_search)}")
    if not re.fullmatch(r"[0-9a-f]{64}", _text(local_file_hash).casefold()):
        raise ValueError("search evidence has an invalid local file hash")

    required_provider = set(CHECKPOINT_COLUMNS)
    missing_provider = required_provider.difference(provider_checkpoint.columns)
    if missing_provider:
        raise ValueError(f"provider checkpoint missing columns: {sorted(missing_provider)}")
    providers = provider_checkpoint[CHECKPOINT_COLUMNS].copy().fillna("")
    if providers.duplicated(["species", "trait", "provider"]).any():
        raise ValueError("provider checkpoint contains duplicate keys")
    providers = providers.set_index(["species", "trait", "provider"], drop=False)

    selections = pd.DataFrame(
        _search_selection_rows(selection_payload),
        columns=[
            "source_row_index",
            "accepted_species",
            "source_url",
            "accepted_trait",
            "accepted_value",
            "evidence_excerpt",
            "review_reason",
        ],
    ).fillna("")
    if selections.empty:
        return pd.DataFrame(columns=sorted(REQUIRED_CANDIDATE_COLUMNS)), {"accepted": []}
    if selections.duplicated(["source_row_index", "accepted_trait"]).any():
        raise ValueError("search selections contain duplicate source-row trait keys")
    conflicts = (
        selections.groupby(["accepted_species", "accepted_trait"])["accepted_value"]
        .nunique()
        .loc[lambda values: values > 1]
    )
    if not conflicts.empty:
        raise ValueError(f"search selections contain trait conflicts: {conflicts.index.tolist()}")

    candidates: list[dict[str, object]] = []
    decisions: list[dict[str, str]] = []
    for source_row_index, group in selections.groupby("source_row_index", sort=True):
        if source_row_index not in search_evidence.index:
            raise ValueError(f"search selection row is unavailable: {source_row_index}")
        source = search_evidence.loc[source_row_index].fillna("")
        if isinstance(source, pd.DataFrame):
            raise ValueError("search evidence index must be unique")
        species = _text(source["accepted_species"])
        selected_species = set(group["accepted_species"].map(_text))
        if selected_species != {species}:
            raise ValueError(
                f"search selection species does not match source row {source_row_index}"
            )
        if species not in queue_species:
            raise ValueError(f"search selection species is outside the current queue: {species}")
        source_url = _text(source["source_url"])
        selected_urls = set(group["source_url"].map(_text))
        if selected_urls != {source_url}:
            raise ValueError(f"search selection URL does not match source row {source_row_index}")
        source_type = _text(source["source_type"])
        source_text = (
            str(source["source_text"]).strip()
            if source_type == DRYAD_BIOGEO_SOURCE_TYPE
            else _text(source["source_text"])
        )
        citation = _text(source["source_citation"])
        evidence_scope = _text(source["evidence_scope"])
        provider = _source_provider(source_type, source_url)
        if not provider:
            raise ValueError(f"search selection has an unsupported provider: {source_row_index}")
        if evidence_scope not in {"species_direct", "synonym_direct"}:
            raise ValueError(f"search selection is not species-direct: {source_row_index}")
        publication, year, doi = _citation_fields(citation)
        if not source_text or not source_url or not publication or not year:
            raise ValueError(f"search selection lacks dated source provenance: {source_row_index}")
        if "doi:" in citation.casefold() and not doi:
            raise ValueError(f"search selection DOI could not be parsed: {source_row_index}")

        row_local_path = _text(source.get("local_file_path"))
        row_local_hash = _text(source.get("local_file_hash")).casefold()
        if bool(row_local_path) != bool(row_local_hash):
            raise ValueError(
                f"search selection has incomplete row-level provenance: {source_row_index}"
            )
        effective_local_path = _text(local_file_path)
        effective_local_hash = _text(local_file_hash).casefold()
        if row_local_path:
            if not re.fullmatch(r"[0-9a-f]{64}", row_local_hash):
                raise ValueError(
                    f"search selection has invalid row-level source hash: {source_row_index}"
                )
            row_path = Path(row_local_path)
            if not row_path.is_file() or _file_sha256(row_path) != row_local_hash:
                raise ValueError(
                    f"search selection row-level source hash mismatch: {source_row_index}"
                )
            effective_local_path = str(row_path.resolve())
            effective_local_hash = row_local_hash

        source_text_hash = _sha256_text(source_text)
        candidate_id = _sha256_text(
            "\x1f".join(
                (
                    "reviewed_search_cache_v1",
                    effective_local_hash,
                    str(source_row_index),
                    species,
                    provider,
                    source_url,
                    source_text_hash,
                )
            )
        )
        first = group.iloc[0]
        first_trait = _text(first["accepted_trait"])
        first_value = _text(first["accepted_value"])
        provider_key = (species, first_trait, provider)
        if provider_key not in providers.index:
            raise ValueError(f"search selection lacks provider state: {provider_key}")
        provider_row = providers.loc[provider_key]
        retrieval_timestamp = _text(provider_row["updated_at"])
        source_retrieval_date = _text(source.get("retrieval_date"))
        retrieval_date = source_retrieval_date or retrieval_timestamp[:10]
        if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", retrieval_date):
            raise ValueError(f"search selection lacks provider retrieval date: {source_row_index}")
        retrieval_date_basis = (
            "search_evidence.retrieval_date"
            if source_retrieval_date
            else "provider_checkpoint.updated_at"
        )
        source_run_id = _text(source.get("source_run_id")) or _text(provider_row["source_run_id"])

        for selection in group.to_dict("records"):
            trait = _text(selection["accepted_trait"])
            value = _text(selection["accepted_value"])
            excerpt = _text(selection["evidence_excerpt"])
            reason = _text(selection["review_reason"])
            if trait not in RA_TRAITS or value not in ALLOWED_DIRECT_VALUES.get(trait, frozenset()):
                raise ValueError(f"search selection has an invalid mapping: {trait}={value}")
            if provider == "pladias" and not _pladias_field_mapping_supported(
                trait, value, source_text
            ):
                raise ValueError(
                    "Pladias review selection is not an allowed exact-field mapping: "
                    f"{source_row_index} -> {trait}={value}"
                )
            if not excerpt or (
                excerpt != source_text
                if provider == DRYAD_BIOGEO_PROVIDER
                else excerpt not in source_text
            ):
                raise ValueError(
                    f"search selection excerpt is not exact source text: {source_row_index}"
                )
            if not reason:
                raise ValueError(f"search selection lacks a review reason: {source_row_index}")
            if provider == DRYAD_BIOGEO_PROVIDER and not _dryad_biogeo_mapping_supported(
                species=species,
                family=_queue_family(queue, species),
                trait=trait,
                value=value,
                source_type=source_type,
                source_text=source_text,
                local_file_path=effective_local_path,
                local_file_hash=effective_local_hash,
                candidate_queue_family=_text(source.get("queue_family")),
                candidate_source_family=_text(source.get("source_family")),
            ):
                raise ValueError(
                    "Dryad BioGeo review selection is not an exact pinned SI/SC row: "
                    f"{source_row_index} -> {trait}={value}"
                )
            if provider != DRYAD_BIOGEO_PROVIDER and not _explicit_mapping_supported(
                trait, value, source_text
            ):
                raise ValueError(
                    "search selection lacks an explicit source statement: "
                    f"{source_row_index} -> {trait}={value}"
                )
            if (species, trait, provider) not in providers.index:
                raise ValueError(
                    f"search selection lacks provider state: {(species, trait, provider)}"
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

        candidates.append(
            {
                "candidate_id": candidate_id,
                "accepted_species": species,
                "trait_name": first_trait,
                "provisional_candidate_value": first_value,
                "matched_term": "strict_reviewed_search_selection",
                "candidate_class": "reported",
                "provider": provider,
                "source_type": source_type,
                "source_url": source_url,
                "source_citation": citation,
                "raw_description": source_text,
                "source_text": source_text,
                "evidence_scope": evidence_scope,
                "wild_or_cultivated": "unknown",
                "review_status": "strictly_reviewed",
                "provider_status": _text(provider_row["status"]),
                "provider_attempts": str(_nonnegative_int(provider_row["attempts"])),
                "provider_updated_at": retrieval_timestamp,
                "provider_policy_version": _text(provider_row["provider_version"]),
                "local_file_path": effective_local_path,
                "local_file_hash": effective_local_hash,
                "retrieval_date": retrieval_date,
                "retrieval_date_basis": retrieval_date_basis,
                "source_run_id": source_run_id,
                "queue_family": _text(source.get("queue_family")),
                "source_family": _text(source.get("source_family")),
                "source_record_id": _text(source.get("source_record_id")),
                "source_row_number": _text(source.get("source_row_number")),
                "dataset_doi": _text(source.get("dataset_doi")),
                "article_doi": _text(source.get("article_doi")),
                "article_publication": _text(source.get("article_publication")),
                "article_year": _text(source.get("article_year")),
                "original_literature_reference": _text(source.get("original_literature_reference")),
                "source_mean_tm_raw_unmapped": _text(source.get("source_mean_tm_raw_unmapped")),
                "audit_status": "strict_review_candidate",
                "audit_reason": "manual_species_direct_selection_from_local_search_cache",
            }
        )

    candidate_frame = pd.DataFrame(candidates).fillna("")
    if candidate_frame["candidate_id"].duplicated().any():
        raise ValueError("reviewed search batch contains duplicate candidate IDs")
    candidate_frame = candidate_frame.sort_values(
        ["accepted_species", "provider", "candidate_id"]
    ).reset_index(drop=True)
    review_payload: dict[str, object] = {
        "contract_version": "reviewed_search_cache_selection_v1",
        "accepted": decisions,
    }
    return candidate_frame, review_payload


def _review_decision_rows(payload: object) -> list[dict[str, str]]:
    if not isinstance(payload, dict):
        raise ValueError("review decision file must contain a JSON object")
    accepted = payload.get("accepted", [])
    if not isinstance(accepted, list) or not all(isinstance(row, dict) for row in accepted):
        raise ValueError("review decision accepted must be a list of objects")
    rows: list[dict[str, str]] = []
    for row in accepted:
        rows.append(
            {
                "candidate_id": _text(row.get("candidate_id")),
                "accepted_trait": _text(row.get("accepted_trait")),
                "accepted_value": _text(row.get("accepted_value")),
                "review_reason": _text(row.get("review_reason")),
                "evidence_excerpt": _text(row.get("evidence_excerpt")),
            }
        )
    return rows


def _pladias_field_mapping_supported(trait: str, value: str, text: str) -> bool:
    """Restrict Pladias review overrides to its exact field/category contract."""
    folded = re.sub(r"[\u2010-\u2015\u2212]", "-", text.casefold())
    match = re.search(
        r"\bgenerative reproduction type\s*[:=]\s*([^\.\r\n]+)",
        folded,
    )
    if match is None:
        return False
    category = _text(match.group(1)).strip(" .:;,|/").casefold()
    has_si = bool(re.search(r"\bself[- ]incompatib\w*\b", category))
    has_sc = bool(re.search(r"\bself[- ]compatib\w*\b", category))
    base = re.sub(r"\bself[- ](?:incompatib|compatib)\w*\b", " ", category)
    segments = [
        _text(segment).casefold() for segment in re.split(r"[,;|/()]+", base) if _text(segment)
    ]
    mating_values = {
        "allogamy": "predominantly_outcrossing",
        "facultative allogamy": "predominantly_outcrossing",
        "autogamy": "obligate_selfing",
        "facultative autogamy": "predominantly_selfing",
        "mixed mating": "mixed_mating",
    }
    if trait == "mating_system":
        observed = {mating_values[segment] for segment in segments if segment in mating_values}
        if "mixed_mating" in observed:
            return value == "mixed_mating"
        return len(observed) == 1 and value in observed
    if trait != "self_incompatibility":
        return False
    if has_si and has_sc:
        return value == "mixed_or_variable"
    if has_si:
        return value == "SI"
    if has_sc:
        return value == "SC"
    return False


def _dryad_biogeo_mapping_supported(
    *,
    species: str,
    family: str,
    trait: str,
    value: str,
    source_type: str,
    source_text: str,
    local_file_path: str,
    local_file_hash: str,
    candidate_queue_family: str,
    candidate_source_family: str,
) -> bool:
    """Validate one exact categorical row against the pinned Dryad CSV.

    This deliberately does not inspect or interpret ``mean.tm``.  A valid
    candidate must reproduce a complete CSV row as JSON and use only the
    dataset's explicit ``si`` contract.
    """
    species = _text(species)
    family = _text(family)
    if (
        source_type != DRYAD_BIOGEO_SOURCE_TYPE
        or trait != "self_incompatibility"
        or value not in {"SI", "SC"}
        or not re.fullmatch(r"[A-Z][A-Za-z-]+ [a-z][A-Za-z-]+", species)
        or not family
        or _text(candidate_queue_family) != family
        or _text(candidate_source_family) != family
        or _text(local_file_hash).casefold() != DRYAD_BIOGEO_CSV_SHA256
    ):
        return False
    source_path = Path(_text(local_file_path))
    if not source_path.is_file() or _file_sha256(source_path) != DRYAD_BIOGEO_CSV_SHA256:
        return False
    try:
        raw = json.loads(str(source_text).strip())
    except (TypeError, ValueError, json.JSONDecodeError):
        return False
    if (
        not isinstance(raw, dict)
        or set(raw) != set(DRYAD_BIOGEO_SOURCE_COLUMNS)
        or not all(isinstance(raw.get(column), str) for column in DRYAD_BIOGEO_SOURCE_COLUMNS)
    ):
        return False
    source_species = _text(f"{raw['genus']} {raw['species']}")
    si_map = {"0": "SI", "1": "SC"}
    if source_species != species or _text(raw["family"]) != family:
        return False
    if raw["si"] not in si_map or si_map[raw["si"]] != value:
        return False

    exact_row_found = False
    relevant_si_rows: list[dict[str, str]] = []
    try:
        with source_path.open("r", encoding="latin-1", newline="") as handle:
            reader = csv.DictReader(handle)
            if tuple(reader.fieldnames or ()) != DRYAD_BIOGEO_SOURCE_COLUMNS:
                return False
            for source_row in reader:
                if None in source_row or any(item is None for item in source_row.values()):
                    return False
                preserved = {column: str(source_row[column]) for column in reader.fieldnames or ()}
                row_species = _text(f"{preserved['genus']} {preserved['species']}")
                if preserved == raw:
                    exact_row_found = True
                if row_species == species and _text(preserved["si"]):
                    relevant_si_rows.append(preserved)
    except (OSError, UnicodeError, csv.Error):
        return False
    if not exact_row_found or not relevant_si_rows:
        return False
    if any(_text(row["family"]) != family for row in relevant_si_rows):
        return False
    observed_si = {_text(row["si"]) for row in relevant_si_rows}
    return observed_si == {raw["si"]} and observed_si.issubset(si_map)


def _queue_family(queue: pd.DataFrame, species: str) -> str:
    if "family" not in queue.columns:
        return ""
    matches = queue.loc[queue[_species_column(queue)].map(_text).eq(species), "family"].map(_text)
    return matches.iloc[0] if len(matches) == 1 else ""


def _explicit_mapping_supported(trait: str, value: str, text: str) -> bool:
    folded = re.sub(r"[\u2010-\u2015\u2212]", "-", text.casefold())
    patterns = {
        ("self_incompatibility", "SI"): (
            r"\bself[- ]incompatib\w*\b|"
            r"\bpollen tubes?\b.{0,180}\b(?:rarely|failed to|did not)\b.{0,100}"
            r"\b(?:selfing|self[- ]pollination)\b"
        ),
        ("self_incompatibility", "SC"): (
            r"\bself[- ]compatib\w*\b|\b(?:partially )?self[- ]fertile\b|"
            r"\bproduced seeds? by self[- ]fertilization\b|"
            r"\b(?:seeds?|fruits?) (?:were )?(?:obtained|produced|set) "
            r"(?:through|by|after) self[- ]pollination\b|"
            r"\b(?:spontaneous|autonomous) self[- ]pollination\b.{0,180}"
            r"\b\d+(?:\.\d+)?\s*%.{0,80}\b(?:set|formed) fruits?\b|"
            r"\bself[- ]fecundation rate\b.{0,80}\b\d+(?:\.\d+)?\s*%|"
            r"\bfruits? matured through autogamy\b|"
            r"\bafter (?:experimental )?selfing\b.{0,180}"
            r"\b(?:germination|viable|fruit|seed)\w*\b|"
            r"\bafter self[- ]pollination\b.{0,220}"
            r"\b(?:full seeds?|fruit weight|seed viability)\b|"
            r"\b(?:inbred )?lines?\b.{0,120}"
            r"\b(?:generated|obtained|purified)\b.{0,80}"
            r"\b(?:across|for|over|with)\s+"
            r"(?:\d+|one|two|three|four|five|six|seven|eight|nine|ten)"
            r"(?:\s+to\s+(?:\d+|one|two|three|four|five|six|seven|eight|nine|ten))?"
            r"\s+(?:successive|consecutive)?\s*generations?\s+of\s+"
            r"self[- ]pollination\b|"
            r"\b(?:inbred )?lines?\b.{0,120}\b(?:generated|obtained)\b.{0,80}"
            r"\b(?:by|through)\s+self[- ]pollination\b.{0,80}"
            r"\b(?:across|for|over)\s+"
            r"(?:\d+|one|two|three|four|five|six|seven|eight|nine|ten)\s+"
            r"(?:successive|consecutive)?\s*generations?\b|"
            r"\bpopulations?\b.{0,80}\bobtained from lines? with\s+"
            r"(?:\d+|one|two|three|four|five|six|seven|eight|nine|ten)"
            r"(?:\s+to\s+(?:\d+|one|two|three|four|five|six|seven|eight|nine|ten))?"
            r"\s+(?:successive|consecutive)?\s*generations?\s+of\s+"
            r"self[- ]pollination\b|"
            r"\bbagging of \w+ yielded\b.{0,80}\bfruit set\b|"
            r"\bobligate self[- ]pollination system\b"
        ),
        ("self_incompatibility", "mixed_or_variable"): (
            r"\b(?:partial|variable|quantitative|facultative)\w*\b.{0,100}"
            r"\bself[- ]incompatib\w*\b|\bfacultative success of selfing\b|"
            r"\bpredominantly self[- ]incompatible\b.{0,180}"
            r"\b(?:spontaneous selfing|viable offspring)\b|"
            r"\bself[- ]compatibility, self[- ]incompatibility\b"
        ),
        ("autonomous_selfing_capacity", "absent"): (
            r"(?:did not occur|incapable of autonomous self[- ]pollination|"
            r"did not self[- ]pollinate autonomously|"
            r"neither autogamous|no (?:autonomous )?self[- ]pollination|"
            r"not autonomously selfing|not autogamous|non-autogamous|failed to set seed|"
            r"reliance .{0,100} mobile pollen vector for reproduction|"
            r"pollinators? .{0,80} necessary for successful pollination and fruit set|"
            r"no (?:pod|fruit|seed) set .{0,100} absence of (?:insect )?pollinators)"
        ),
        ("autonomous_selfing_capacity", "delayed"): r"\bdelayed selfing\b",
        ("autonomous_selfing_capacity", "facilitated"): r"\bfacilitated selfing\b",
        ("autonomous_selfing_capacity", "autonomous"): (
            r"\b(?:autonomous|spontaneous) self(?:ing|[- ]pollination)\b|"
            r"\bautomatic self[- ]pollination\b|\bself[- ]pollinating\b|"
            r"\bcould be self[- ]?pollinated before anthesis\b|"
            r"\bautogam\w*\b|"
            r"\b(?:closed|cleistogamous|cl) flowers?\b.{0,80}"
            r"\b(?:obligately )?self[- ]pollinated\b|"
            r"\bself[- ]pollinated (?:cl|cleistogamous) flowers?\b|"
            r"\bobligate self[- ]pollination system\b"
        ),
        ("autonomous_selfing_capacity", "mixed_or_variable"): (
            r"\bvariable autonomous selfing\b|\bmixed autonomous selfing\b|"
            r"\blow level of spontaneous selfing\b|"
            r"\blabile pollination strategy\b.{0,100}\bautogamy sometimes occurs\b"
        ),
        ("mating_system", "predominantly_outcrossing"): (
            r"\b(?:predominant(?:ly)?|primarily|mainly|mostly|highly|dominantly|largely) "
            r"(?:outcrossing|outcrossed|outcross[- ]pollinated|xenogamous|allogamous)\b|"
            r"\bobligate (?:outcrossing|outbreeding)\b|"
            r"\bthe outcrossing [a-z]+ [a-z-]+\b|"
            r"\bmating system as outcrossing\b|"
            r"\boutcrossing(?:,? [a-z-]+){0,3} (?:species|crop|plant|aster)\b|"
            r"\boutcrossed.{0,60}\bmating strategy\b|"
            r"\b(?:outcrossing|outbreeding) (?:nature|reproductive system|"
            r"breeding system|mating system)\b|\boutbreeders?\b|\ballogamous crops?\b|"
            r"\bfacultative xenogamous species\b|"
            r"\bconserved by outcrossing\b|\([^)]*\boutcrossing\b[,)]|"
            r"\b(?:primarily|predominantly) an? outcrosser\b|"
            r"\bgenerative reproduction type\s*[:=]\s*"
            r"(?:facultative\s+)?allogamy(?:\s+self[- ]incompatib\w*)?\b"
        ),
        ("mating_system", "mixed_mating"): (
            r"\bmixed[- ]mating(?: system)?\b|"
            r"\bmixed (?:breeding|reproductive) (?:system|strateg(?:y|ies))\b|"
            r"\bboth cross- and self-fertilizations had occurred\b|"
            r"\bboth self[- ]fertilization\b.{0,100}\b(?:as well as |and )?"
            r"cross[- ]pollination\b|"
            r"\buses both compatibility strategies\b.{0,100}\b(?:xenogamy|autogamy)\b"
        ),
        ("mating_system", "predominantly_selfing"): (
            r"\b(?:predominantly|mainly) selfing\b|\bpredominant(?:ly)? autogam\w*\b|"
            r"\b(?:a |considered to be a )selfing (?:annual|species)\b|"
            r"\bgenerative reproduction type\s*[:=]\s*facultative\s+autogamy\b"
        ),
        ("mating_system", "obligate_selfing"): (
            r"\bobligate(?:ly)? selfing\b|\bobligate self[- ]pollination system\b|"
            r"\bgenerative reproduction type\s*[:=]\s*autogamy\b"
        ),
        ("cleistogamy", "absent"): r"\b(?:no|not|without) cleistogam\w*\b",
        ("cleistogamy", "facultative"): (
            r"\bfacultative cleistogam\w*\b|\bcleistogam\w*\b.*\bchasmogam\w*\b|"
            r"\bcleistogam\w*\b.*\bopen flowers?\b|"
            r"\b(?:often|usually|sometimes) cleistogam\w*\b"
        ),
        ("cleistogamy", "obligate"): r"\bobligate cleistogam\w*\b",
    }
    pattern = patterns.get((trait, value))
    return bool(pattern and re.search(pattern, folded, flags=re.IGNORECASE | re.DOTALL))


def _validate_reviewed_candidates(
    queue: pd.DataFrame,
    candidates: pd.DataFrame,
    decision_rows: list[dict[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    missing = REQUIRED_CANDIDATE_COLUMNS.difference(candidates.columns)
    if missing:
        raise ValueError(f"candidate table missing columns: {sorted(missing)}")
    frame = candidates.copy().fillna("")
    if frame["candidate_id"].map(_text).eq("").any():
        raise ValueError("candidate table contains a blank candidate_id")
    if frame["candidate_id"].duplicated().any():
        raise ValueError("candidate table contains duplicate candidate_id values")

    queue_species = set(_validate_queue(queue))
    review = pd.DataFrame(
        decision_rows,
        columns=[
            "candidate_id",
            "accepted_trait",
            "accepted_value",
            "review_reason",
            "evidence_excerpt",
        ],
    ).fillna("")
    if review["candidate_id"].eq("").any():
        raise ValueError("accepted review decision contains a blank candidate_id")
    unknown = sorted(set(review["candidate_id"]).difference(frame["candidate_id"]))
    if unknown:
        raise ValueError(f"review decisions refer to unknown candidates: {unknown[:5]}")
    accepted = review.merge(frame, on="candidate_id", how="left", validate="many_to_one")
    accepted["review_evidence_excerpt"] = accepted["evidence_excerpt"].where(
        accepted["evidence_excerpt"].map(_text).ne(""),
        accepted["raw_description"],
    )
    accepted["accepted_trait"] = accepted["accepted_trait"].where(
        accepted["accepted_trait"].map(_text).ne(""),
        accepted["trait_name"],
    )
    if accepted.duplicated(["candidate_id", "accepted_trait"]).any():
        raise ValueError("accepted review decision contains duplicate candidate-trait keys")

    accepted_ids = set(accepted["candidate_id"])
    mappings_by_id = {
        candidate_id: "|".join(
            f"{_text(row.accepted_trait)}={_text(row.accepted_value)}"
            for row in group.itertuples(index=False)
        )
        for candidate_id, group in accepted.groupby("candidate_id", sort=False)
    }
    reasons_by_id = {
        candidate_id: " | ".join(dict.fromkeys(group["review_reason"].map(_text)))
        for candidate_id, group in accepted.groupby("candidate_id", sort=False)
    }
    audit = frame.copy()
    audit["review_decision"] = audit["candidate_id"].map(
        lambda candidate_id: "accepted" if candidate_id in accepted_ids else "rejected"
    )
    audit["accepted_mappings"] = audit["candidate_id"].map(mappings_by_id).fillna("")
    audit["review_reason"] = audit["candidate_id"].map(
        lambda candidate_id: (
            reasons_by_id.get(candidate_id, "")
            if candidate_id in accepted_ids
            else "fail_closed_not_selected_as_unambiguous_species_direct_evidence"
        )
    )
    audit["review_basis"] = REVIEW_BASIS

    for row in accepted.to_dict("records"):
        species = _text(row["accepted_species"])
        trait = _text(row["accepted_trait"])
        value = _text(row["accepted_value"])
        provisional_trait = _text(row["trait_name"])
        provisional_value = _text(row["provisional_candidate_value"])
        provider = _text(row["provider"])
        if species not in queue_species:
            raise ValueError(f"accepted candidate is outside the exact priority queue: {species}")
        if trait not in RA_TRAITS:
            raise ValueError(f"accepted candidate has unsupported trait: {trait}")
        if value not in ALLOWED_DIRECT_VALUES[trait]:
            raise ValueError(f"accepted candidate has unsupported value: {trait}={value}")
        if _text(row["evidence_scope"]) not in {"species_direct", "synonym_direct"}:
            raise ValueError(f"accepted candidate is not species-direct: {row['candidate_id']}")
        if provider not in RECORD_PROVIDERS:
            raise ValueError(f"accepted candidate has unsupported provider: {row['provider']}")
        if provider == "pladias" and not _pladias_field_mapping_supported(
            trait, value, _text(row["source_text"])
        ):
            raise ValueError(
                "accepted Pladias candidate is not an allowed exact-field mapping: "
                f"{row['candidate_id']} -> {trait}={value}"
            )
        if provider == DRYAD_BIOGEO_PROVIDER and not _dryad_biogeo_mapping_supported(
            species=species,
            family=_queue_family(queue, species),
            trait=trait,
            value=value,
            source_type=_text(row["source_type"]),
            source_text=str(row["source_text"]).strip(),
            local_file_path=_text(row["local_file_path"]),
            local_file_hash=_text(row["local_file_hash"]),
            candidate_queue_family=_text(row.get("queue_family")),
            candidate_source_family=_text(row.get("source_family")),
        ):
            raise ValueError(
                "accepted Dryad BioGeo candidate is not an exact pinned SI/SC row: "
                f"{row['candidate_id']} -> {trait}={value}"
            )
        if not _text(row["review_reason"]):
            raise ValueError(f"accepted candidate lacks a review reason: {row['candidate_id']}")
        required_provenance = (
            "source_url",
            "source_citation",
            "raw_description",
            "source_text",
            "source_type",
            "local_file_hash",
            "retrieval_date",
        )
        blank = [column for column in required_provenance if not _text(row[column])]
        if blank:
            raise ValueError(
                f"accepted candidate lacks required provenance {blank}: {row['candidate_id']}"
            )
        publication, year, doi = _citation_fields(row["source_citation"])
        if not publication or not year:
            raise ValueError(f"accepted candidate lacks a dated publication: {row['candidate_id']}")
        if "doi:" in _text(row["source_citation"]).casefold() and not doi:
            raise ValueError(f"accepted candidate DOI could not be parsed: {row['candidate_id']}")
        if not re.fullmatch(r"[0-9a-f]{64}", _text(row["local_file_hash"]).casefold()):
            raise ValueError(
                f"accepted candidate has an invalid local file hash: {row['candidate_id']}"
            )
        if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", _text(row["retrieval_date"])):
            raise ValueError(
                f"accepted candidate has an invalid retrieval date: {row['candidate_id']}"
            )
        text = f"{_text(row['raw_description'])} {_text(row['source_text'])}".casefold()
        if provider in EXTERNAL_ACQUISITION_PROVIDERS and provider != DRYAD_BIOGEO_PROVIDER:
            literature_provider = provider
            excerpt = _text(row["review_evidence_excerpt"])
            local_path = Path(_text(row["local_file_path"]))
            if (
                not local_path.is_file()
                or _file_sha256(local_path) != _text(row["local_file_hash"]).casefold()
            ):
                raise ValueError(
                    f"{literature_provider} raw-response provenance mismatch: {row['candidate_id']}"
                )
            species_expression = r"\s+".join(re.escape(token) for token in species.split())
            species_named = bool(
                re.search(
                    rf"(?<!\w){species_expression}(?!\w)",
                    excerpt,
                    flags=re.IGNORECASE,
                )
            )
            if not species_named or not _explicit_mapping_supported(trait, value, excerpt):
                raise ValueError(
                    f"{literature_provider} acceptance excerpt must name the exact species "
                    "and state "
                    f"the mapped condition: {row['candidate_id']} -> {trait}={value}"
                )
        if provider != DRYAD_BIOGEO_PROVIDER and not _explicit_mapping_supported(
            trait, value, text
        ):
            if (trait, value) != (provisional_trait, provisional_value):
                raise ValueError(
                    "review may override a machine mapping only for an explicit source "
                    f"statement: {row['candidate_id']} -> {trait}={value}"
                )
            raise ValueError(
                "accepted review mapping lacks an explicit source statement: "
                f"{row['candidate_id']} -> {trait}={value}"
            )

    conflicts = (
        accepted.groupby(["accepted_species", "accepted_trait"])["accepted_value"]
        .nunique()
        .loc[lambda values: values > 1]
    )
    if not conflicts.empty:
        raise ValueError(
            f"accepted review decisions contain trait conflicts: {conflicts.index.tolist()}"
        )
    return accepted, audit


def _records_from_accepted(accepted: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for row in accepted.to_dict("records"):
        species = _text(row["accepted_species"])
        trait = _text(row["accepted_trait"])
        value = _text(row["accepted_value"])
        provider = _text(row["provider"])
        source_url = _text(row["source_url"])
        source_text = (
            str(row["source_text"]).strip()
            if provider == DRYAD_BIOGEO_PROVIDER
            else _text(row["source_text"])
        )
        publication, year, doi = _citation_fields(row["source_citation"])
        source_text_hash = _sha256_text(source_text)
        identity = "\x1f".join((species, trait, value, provider, source_url, source_text_hash))
        rows.append(
            {
                "record_id": _sha256_text(identity),
                "species": species,
                "trait": trait,
                "value": value,
                "raw_text": source_text,
                "evidence_excerpt": _text(row.get("review_evidence_excerpt"))
                or _text(row["raw_description"]),
                "URL": source_url,
                "DOI": doi,
                "publication": publication,
                "year": year,
                "source_type": _text(row["source_type"]),
                "provider": provider,
                "evidence_scope": _text(row["evidence_scope"]),
                "local_file_path": _text(row["local_file_path"]),
                "local_file_hash": _text(row["local_file_hash"]),
                "source_text_hash": source_text_hash,
                "retrieval_date": _text(row["retrieval_date"]),
                "candidate_id": _text(row["candidate_id"]),
                "source_run_id": _text(row["source_run_id"]),
                "direct_evidence": True,
                "review_basis": REVIEW_BASIS,
            }
        )
    return pd.DataFrame(rows, columns=RECORD_COLUMNS)


def _merge_records(current: pd.DataFrame, prior: pd.DataFrame | None) -> pd.DataFrame:
    frames = [current]
    if prior is not None and not prior.empty:
        missing = set(RECORD_COLUMNS).difference(prior.columns)
        if missing:
            raise ValueError(f"existing record table missing columns: {sorted(missing)}")
        frames.insert(0, prior[RECORD_COLUMNS].copy().fillna(""))
    records = (
        pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=RECORD_COLUMNS)
    )
    if records.empty:
        return pd.DataFrame(columns=RECORD_COLUMNS)
    records = records.drop_duplicates("record_id", keep="first")
    conflicts = (
        records.groupby(["species", "trait"])["value"].nunique().loc[lambda values: values > 1]
    )
    if not conflicts.empty:
        raise ValueError(f"accepted records contain trait conflicts: {conflicts.index.tolist()}")
    return records.sort_values(["species", "trait", "provider", "record_id"]).reset_index(drop=True)


def _update_provider_checkpoint(
    provider_checkpoint: pd.DataFrame,
    records: pd.DataFrame,
) -> pd.DataFrame:
    required = set(CHECKPOINT_COLUMNS)
    missing = required.difference(provider_checkpoint.columns)
    if missing:
        raise ValueError(f"provider checkpoint missing columns: {sorted(missing)}")
    checkpoint = provider_checkpoint[CHECKPOINT_COLUMNS].copy().fillna("")
    if checkpoint.duplicated(["species", "trait", "provider"]).any():
        raise ValueError("provider checkpoint contains duplicate keys")
    checkpoint["attempts"] = checkpoint["attempts"].map(_nonnegative_int).astype("int64")
    checkpoint["candidate_count"] = (
        checkpoint["candidate_count"].map(_nonnegative_int).astype("int64")
    )
    checkpoint["accepted_record_count"] = (
        checkpoint["accepted_record_count"].map(_nonnegative_int).astype("int64")
    )
    checkpoint["terminal"] = checkpoint["terminal"].map(_bool).astype(bool)
    counts = (
        records.groupby(["species", "trait", "provider"]).size().to_dict()
        if not records.empty
        else {}
    )
    observed = set(checkpoint[["species", "trait", "provider"]].itertuples(index=False, name=None))
    missing_keys = sorted(set(counts).difference(observed))
    if missing_keys:
        raise ValueError(f"accepted records lack provider checkpoint keys: {missing_keys[:5]}")

    def append_review_basis(value: object) -> str:
        bases = [item for item in _text(value).split("|") if item]
        if REVIEW_BASIS not in bases:
            bases.append(REVIEW_BASIS)
        return "|".join(bases)

    def append_seed_covered_basis(value: object) -> str:
        basis = "reproductive_assurance_completed_by_another_provider"
        bases = [item for item in _text(value).split("|") if item]
        if basis not in bases:
            bases.append(basis)
        return "|".join(bases)

    for key, count in counts.items():
        mask = (
            checkpoint["species"].eq(key[0])
            & checkpoint["trait"].eq(key[1])
            & checkpoint["provider"].eq(key[2])
        )
        checkpoint.loc[mask, "accepted_record_count"] = int(count)
        checkpoint.loc[mask, "status"] = "accepted_evidence"
        checkpoint.loc[mask, "terminal"] = True
        checkpoint.loc[mask, "migration_basis"] = checkpoint.loc[mask, "migration_basis"].map(
            append_review_basis
        )
    covered_species = {species for species in records["species"].map(_text) if species}
    seed_covered = checkpoint["species"].isin(covered_species) & ~checkpoint["terminal"]
    checkpoint.loc[seed_covered, "status"] = "provider_seed_covered"
    checkpoint.loc[seed_covered, "terminal"] = True
    checkpoint.loc[seed_covered, "migration_basis"] = checkpoint.loc[
        seed_covered, "migration_basis"
    ].map(append_seed_covered_basis)
    return checkpoint.sort_values(["species", "trait", "provider"]).reset_index(drop=True)


def _updated_axis_checkpoint(
    baseline: pd.DataFrame,
    queue: pd.DataFrame,
    records: pd.DataFrame,
) -> pd.DataFrame:
    required = {
        "accepted_species",
        "colour_vividness_direct",
        "floral_complexity_direct",
        "reproductive_assurance_direct",
        "reproductive_assurance_evidence_traits",
        "all_three_axes_direct",
        "missing_direct_axes",
        "reproductive_assurance_priority",
        "provider_state_available",
        "checkpoint_scope",
    }
    missing = required.difference(baseline.columns)
    if missing:
        raise ValueError(f"three-axis checkpoint missing columns: {sorted(missing)}")
    checkpoint = baseline.copy().fillna("")
    checkpoint["accepted_species"] = checkpoint["accepted_species"].map(_text)
    if checkpoint["accepted_species"].duplicated().any():
        raise ValueError("three-axis checkpoint contains duplicate species")
    boolean_columns = (
        "colour_vividness_direct",
        "floral_complexity_direct",
        "reproductive_assurance_direct",
        "all_three_axes_direct",
        "reproductive_assurance_priority",
        "provider_state_available",
    )
    for column in boolean_columns:
        checkpoint[column] = checkpoint[column].map(_bool).astype(bool)
    queue_species = set(_validate_queue(queue))
    if not queue_species.issubset(set(checkpoint["accepted_species"])):
        raise ValueError("priority queue contains species outside the three-axis checkpoint")
    checkpoint.loc[
        checkpoint["accepted_species"].isin(queue_species), "provider_state_available"
    ] = True
    checkpoint.loc[checkpoint["accepted_species"].isin(queue_species), "checkpoint_scope"] = (
        "materialized_axis_with_local_ra_continuation"
    )
    traits_by_species = (
        records.groupby("species")["trait"].agg(lambda values: sorted(set(values))).to_dict()
        if not records.empty
        else {}
    )
    for species, traits in traits_by_species.items():
        mask = checkpoint["accepted_species"].eq(species)
        if not mask.any():
            raise ValueError(
                f"accepted record species is outside the three-axis checkpoint: {species}"
            )
        existing = {
            token
            for token in _text(
                checkpoint.loc[mask, "reproductive_assurance_evidence_traits"].iloc[0]
            ).split("|")
            if token
        }
        checkpoint.loc[mask, "reproductive_assurance_evidence_traits"] = "|".join(
            sorted(existing.union(traits))
        )
        checkpoint.loc[mask, "reproductive_assurance_direct"] = True

    colour = checkpoint["colour_vividness_direct"].map(_bool)
    complexity = checkpoint["floral_complexity_direct"].map(_bool)
    reproductive = checkpoint["reproductive_assurance_direct"].map(_bool)
    checkpoint["all_three_axes_direct"] = colour & complexity & reproductive
    checkpoint["missing_direct_axes"] = [
        "|".join(
            axis
            for axis, complete in (
                ("colour_vividness", colour_value),
                ("floral_complexity", complexity_value),
                ("reproductive_assurance", reproductive_value),
            )
            if not complete
        )
        for colour_value, complexity_value, reproductive_value in zip(
            colour, complexity, reproductive, strict=True
        )
    ]
    checkpoint["reproductive_assurance_priority"] = colour & complexity & ~reproductive
    return checkpoint.sort_values("accepted_species").reset_index(drop=True)


def materialize_reviewed_cache(
    *,
    queue: pd.DataFrame,
    baseline_checkpoint: pd.DataFrame,
    provider_checkpoint: pd.DataFrame,
    candidates: pd.DataFrame,
    review_payload: object,
    prior_records: pd.DataFrame | None = None,
) -> dict[str, object]:
    """Materialize only explicitly reviewed species-direct cache records."""
    decisions = _review_decision_rows(review_payload)
    accepted, audit = _validate_reviewed_candidates(queue, candidates, decisions)
    new_records = _records_from_accepted(accepted)
    prior_record_ids = (
        set(prior_records["record_id"].map(_text))
        if prior_records is not None and not prior_records.empty
        else set()
    )
    added_record_count = int(
        new_records.loc[~new_records["record_id"].isin(prior_record_ids), "record_id"].nunique()
    )
    records = _merge_records(new_records, prior_records)
    providers = _update_provider_checkpoint(provider_checkpoint, records)
    axes = _updated_axis_checkpoint(baseline_checkpoint, queue, records)
    next_queue = axes.loc[axes["reproductive_assurance_priority"].map(_bool)].copy()

    baseline_reproductive = baseline_checkpoint["reproductive_assurance_direct"].map(_bool)
    baseline_all_three = baseline_checkpoint["all_three_axes_direct"].map(_bool)
    latest_batch_completed = sorted(
        set(
            axes.loc[axes["reproductive_assurance_direct"].map(_bool), "accepted_species"]
        ).difference(baseline_checkpoint.loc[baseline_reproductive, "accepted_species"].map(_text))
    )
    cumulative_completed = sorted({species for species in records["species"].map(_text) if species})
    summary = {
        "contract_version": "local_reproductive_assurance_materialization_v1",
        "n_species": int(len(axes)),
        "direct_coverage": {
            "colour_vividness": int(axes["colour_vividness_direct"].map(_bool).sum()),
            "floral_complexity": int(axes["floral_complexity_direct"].map(_bool).sum()),
            "reproductive_assurance": int(axes["reproductive_assurance_direct"].map(_bool).sum()),
        },
        "remaining_species": {
            "priority_queue": int(len(next_queue)),
            "reproductive_assurance": int(
                (~axes["reproductive_assurance_direct"].map(_bool)).sum()
            ),
        },
        "newly_completed_species": int(len(cumulative_completed)),
        "newly_completed_species_names": cumulative_completed,
        "cumulative_newly_completed_species": int(len(cumulative_completed)),
        "cumulative_newly_completed_species_names": cumulative_completed,
        "latest_batch_newly_completed_species": int(len(latest_batch_completed)),
        "latest_batch_newly_completed_species_names": latest_batch_completed,
        "total_three_axis_completion": int(axes["all_three_axes_direct"].map(_bool).sum()),
        "baseline_total_three_axis_completion": int(baseline_all_three.sum()),
        "new_record_count": added_record_count,
        "total_record_count": int(len(records)),
        "cache_candidate_count": int(len(candidates)),
        "cache_accepted_candidate_count": int(accepted["candidate_id"].nunique()),
        "cache_accepted_trait_record_count": int(len(accepted)),
        "cache_rejected_candidate_count": int(len(audit) - accepted["candidate_id"].nunique()),
        "provider_terminal_count": int(providers["terminal"].map(_bool).sum()),
        "provider_remaining_count": int((~providers["terminal"].map(_bool)).sum()),
        "external_acquisition_started": bool(
            set(candidates["provider"].map(_text)).intersection(EXTERNAL_ACQUISITION_PROVIDERS)
        ),
        "global_fallback_used": False,
        "genus_or_family_inference_used": False,
        # A reviewed materialization batch is not, by itself, an entire provider
        # pass.  Only the explicit coverage-refresh command may set this stop flag.
        "no_new_direct_evidence_in_entire_pass": False,
    }
    return {
        "records": records,
        "review_audit": audit,
        "provider_checkpoint": providers,
        "three_axis_checkpoint": axes,
        "next_queue": next_queue,
        "summary": summary,
    }


def _species_column(frame: pd.DataFrame) -> str:
    for column in ("accepted_species", "species"):
        if column in frame.columns:
            return column
    raise ValueError("queue must contain accepted_species or species")


def _validate_queue(queue: pd.DataFrame) -> list[str]:
    species_column = _species_column(queue)
    species = queue[species_column].map(_text)
    if species.eq("").any():
        raise ValueError("queue contains a blank species")
    if species.duplicated().any():
        raise ValueError("queue contains duplicate species")
    return species.tolist()


def _validate_legacy_provider(legacy: pd.DataFrame, species: set[str]) -> pd.DataFrame:
    required = {
        "species",
        "provider",
        "policy_version",
        "enabled",
        "status",
        "attempts",
        "last_error",
        "updated_at",
    }
    missing = required.difference(legacy.columns)
    if missing:
        raise ValueError(f"legacy provider checkpoint missing columns: {sorted(missing)}")
    frame = legacy.copy().fillna("")
    frame["species"] = frame["species"].map(_text)
    frame["provider"] = frame["provider"].map(_text)
    frame = frame.loc[frame["species"].isin(species)].copy()
    if frame.duplicated(["species", "provider"]).any():
        raise ValueError("legacy provider checkpoint has duplicate species-provider keys")
    expected = {(name, provider) for name in species for provider in PROVIDERS}
    observed = set(frame[["species", "provider"]].itertuples(index=False, name=None))
    if observed != expected:
        raise ValueError(
            "legacy provider checkpoint does not cover the exact queue and provider set"
        )
    return frame.set_index(["species", "provider"], drop=False)


def _candidate_counts(candidates: pd.DataFrame | None) -> dict[tuple[str, str, str], int]:
    if candidates is None or candidates.empty:
        return {}
    required = {"accepted_species", "trait_name", "provider"}
    missing = required.difference(candidates.columns)
    if missing:
        raise ValueError(f"candidate table missing columns: {sorted(missing)}")
    frame = candidates.copy().fillna("")
    frame = frame.loc[frame["trait_name"].isin(RA_TRAITS)].copy()
    return {
        (_text(species), _text(trait), _text(provider)): int(len(group))
        for (species, trait, provider), group in frame.groupby(
            ["accepted_species", "trait_name", "provider"], sort=False
        )
    }


def _openalex_state(
    species: str,
    ledger: pd.DataFrame | None,
) -> tuple[str, int, str, str, str]:
    if ledger is None or species not in ledger.index:
        return "pending", 0, "", "", "global_openalex_unattempted"
    row = ledger.loc[species]
    status = _text(row.get("reproductive_openalex_status")) or "pending"
    attempts = _nonnegative_int(row.get("reproductive_openalex_attempts"))
    error = _text(row.get("reproductive_openalex_last_error"))
    wave_id = _text(row.get("reproductive_openalex_wave_id"))
    if status == "processed":
        return (
            "provider_pass_complete",
            attempts,
            error,
            wave_id,
            "global_reproductive_openalex_processed",
        )
    if status == "retry":
        return "retry", attempts, error, wave_id, "global_reproductive_openalex_retry"
    if status in {"exhausted", "permanent_failure"}:
        return (
            "permanent_failure",
            attempts,
            error,
            wave_id,
            "global_reproductive_openalex_exhausted",
        )
    return "pending", attempts, error, wave_id, "global_reproductive_openalex_unattempted"


def _legacy_state(
    *,
    provider: str,
    trait: str,
    row: pd.Series,
    candidate_count: int,
) -> tuple[str, int, str]:
    legacy_status = _text(row.get("status"))
    attempts = _nonnegative_int(row.get("attempts"))

    if provider == "europe_pmc" and trait not in OLD_STATIC_TARGET_TRAITS["europe_pmc"]:
        # Only an actual v2 request accounts for autogamy/cleistogamy. A static
        # seed-covered row did not execute the reproductive query.
        if legacy_status == "skipped_covered" or attempts == 0:
            return "pending", 0, "europe_pmc_v2_query_not_executed"

    if legacy_status == "hit":
        status = (
            "provider_pass_complete_with_candidates"
            if candidate_count
            else "provider_pass_complete"
        )
        return status, attempts, "provider_source_hit_not_trait_acceptance"
    if legacy_status == "no_hit":
        return "provider_no_result", attempts, "provider_pass_completed_without_source_result"
    if legacy_status == "skipped_covered":
        return "provider_seed_covered", attempts, "provider_skipped_due_to_static_seed_coverage"
    if legacy_status == "legacy_completed":
        return "provider_pass_complete", attempts, "legacy_provider_completion_preserved"
    if legacy_status == "exhausted":
        return "permanent_failure", attempts, "legacy_provider_exhausted"
    if legacy_status == "running":
        return "retry", attempts, "interrupted_legacy_provider_run"
    if legacy_status in {"pending", "retry"}:
        return legacy_status, attempts, "legacy_provider_unfinished"
    raise ValueError(f"unsupported legacy provider status: {legacy_status!r}")


def bootstrap_provider_checkpoint(
    queue: pd.DataFrame,
    legacy_provider: pd.DataFrame,
    *,
    candidates: pd.DataFrame | None = None,
    global_ledger: pd.DataFrame | None = None,
    source_run_id: int = 29637034874,
) -> pd.DataFrame:
    """Return one conservative row per queue species, RA trait, and provider."""
    species = _validate_queue(queue)
    species_set = set(species)
    legacy = _validate_legacy_provider(legacy_provider, species_set)
    counts = _candidate_counts(candidates)
    ledger: pd.DataFrame | None = None
    if global_ledger is not None and not global_ledger.empty:
        ledger = global_ledger.copy().fillna("")
        column = _species_column(ledger)
        ledger[column] = ledger[column].map(_text)
        if ledger[column].duplicated().any():
            raise ValueError("global campaign ledger has duplicate species")
        ledger = ledger.set_index(column, drop=False)

    rows: list[dict[str, object]] = []
    for name in species:
        valid_species = len(name.split()) >= 2
        for trait in RA_TRAITS:
            for provider in PROVIDERS:
                legacy_row = legacy.loc[(name, provider)]
                candidate_count = counts.get((name, trait, provider), 0)
                query_traits = PROVIDER_QUERY_TRAITS[provider]
                status = "not_applicable"
                attempts = 0
                error = ""
                updated_at = ""
                last_wave_id = ""
                basis = "provider_does_not_query_trait"
                provider_version = _text(legacy_row.get("policy_version"))

                if trait in query_traits:
                    if not valid_species:
                        status = "invalid_species_name"
                        basis = "species_name_has_fewer_than_two_tokens"
                        attempts = _nonnegative_int(legacy_row.get("attempts"))
                        legacy_error = _text(legacy_row.get("last_error"))
                        error = "permanent:invalid_species_name"
                        if legacy_error:
                            error = f"{error} | legacy:{legacy_error}"
                        updated_at = _text(legacy_row.get("updated_at"))
                    elif provider == "openalex":
                        status, attempts, error, last_wave_id, basis = _openalex_state(name, ledger)
                        provider_version = "global_reproductive_openalex_v1"
                    elif provider == "world_flora":
                        status = "pending"
                        basis = "legacy_provider_disabled_and_unattempted"
                    else:
                        status, attempts, basis = _legacy_state(
                            provider=provider,
                            trait=trait,
                            row=legacy_row,
                            candidate_count=candidate_count,
                        )
                        error = _text(legacy_row.get("last_error"))
                        updated_at = _text(legacy_row.get("updated_at"))

                if status not in ALL_STATUSES:
                    raise AssertionError(f"invalid migrated status: {status}")
                rows.append(
                    {
                        "species": name,
                        "trait": trait,
                        "provider": provider,
                        "provider_version": provider_version,
                        "status": status,
                        "terminal": status in TERMINAL_STATUSES,
                        "attempts": attempts,
                        "candidate_count": candidate_count,
                        "accepted_record_count": 0,
                        "legacy_status": _text(legacy_row.get("status")),
                        "legacy_enabled": _text(legacy_row.get("enabled")),
                        "last_error": error,
                        "updated_at": updated_at,
                        "last_wave_id": last_wave_id,
                        "source_run_id": source_run_id,
                        "migration_basis": basis,
                        "contract_version": CONTRACT_VERSION,
                    }
                )

    checkpoint = pd.DataFrame(rows, columns=CHECKPOINT_COLUMNS)
    if checkpoint.duplicated(["species", "trait", "provider"]).any():
        raise AssertionError("migrated provider checkpoint contains duplicate keys")
    expected = len(species) * len(RA_TRAITS) * len(PROVIDERS)
    if len(checkpoint) != expected:
        raise AssertionError(
            f"migrated provider checkpoint rows={len(checkpoint)}, expected={expected}"
        )
    return checkpoint.sort_values(["species", "trait", "provider"]).reset_index(drop=True)


def write_bootstrap_outputs(
    checkpoint: pd.DataFrame,
    output_dir: Path,
    *,
    source_paths: list[Path],
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "provider_checkpoint.csv"
    with NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        dir=output_dir,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary_checkpoint = Path(handle.name)
        checkpoint.to_csv(handle, index=False)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary_checkpoint, path)
    status_counts = {
        status: int(count)
        for status, count in checkpoint["status"].value_counts().sort_index().items()
    }
    report = {
        "contract_version": CONTRACT_VERSION,
        "n_rows": int(len(checkpoint)),
        "n_species": int(checkpoint["species"].nunique()),
        "n_traits": int(checkpoint["trait"].nunique()),
        "n_providers": int(checkpoint["provider"].nunique()),
        "n_duplicate_keys": int(checkpoint.duplicated(["species", "trait", "provider"]).sum()),
        "status_counts": status_counts,
        "n_terminal": int(checkpoint["terminal"].astype(bool).sum()),
        "n_remaining": int((~checkpoint["terminal"].astype(bool)).sum()),
        "checkpoint": {
            "path": str(path.resolve()),
            "bytes": path.stat().st_size,
            "sha256": _file_sha256(path),
        },
        "sources": {
            str(source.resolve()): {
                "bytes": source.stat().st_size,
                "sha256": _file_sha256(source),
            }
            for source in source_paths
        },
        "biological_values_assigned": False,
        "global_fallback_used": False,
    }
    report_path = output_dir / "provider_checkpoint_bootstrap.json"
    report_payload = json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    with NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        dir=output_dir,
        prefix=f".{report_path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary_report = Path(handle.name)
        handle.write(report_payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary_report, report_path)
    return report


@app.command()
def bootstrap(
    queue_csv: Path = typer.Option(..., exists=True),
    legacy_provider_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    candidate_csv: Path | None = typer.Option(None),
    global_ledger_csv: Path | None = typer.Option(None),
    source_run_id: int = typer.Option(29637034874, min=1),
) -> None:
    """Create the local trait-scoped ledger without performing provider requests."""
    queue = pd.read_csv(queue_csv, dtype=str).fillna("")
    legacy = pd.read_csv(legacy_provider_csv, dtype=str).fillna("")
    candidates = (
        pd.read_csv(candidate_csv, dtype=str).fillna("") if candidate_csv is not None else None
    )
    global_ledger = (
        pd.read_csv(global_ledger_csv, dtype=str).fillna("")
        if global_ledger_csv is not None
        else None
    )
    checkpoint = bootstrap_provider_checkpoint(
        queue,
        legacy,
        candidates=candidates,
        global_ledger=global_ledger,
        source_run_id=source_run_id,
    )
    sources = [queue_csv, legacy_provider_csv]
    if candidate_csv is not None:
        sources.append(candidate_csv)
    if global_ledger_csv is not None:
        sources.append(global_ledger_csv)
    report = write_bootstrap_outputs(checkpoint, output_dir, source_paths=sources)
    typer.echo(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))


@app.command("materialize-reviewed-cache")
def materialize_reviewed_cache_command(
    queue_csv: Path = typer.Option(..., exists=True),
    baseline_checkpoint_csv: Path = typer.Option(..., exists=True),
    provider_checkpoint_csv: Path = typer.Option(..., exists=True),
    candidate_csv: Path = typer.Option(..., exists=True),
    review_json: list[Path] = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    prior_records_csv: Path | None = typer.Option(None),
) -> None:
    """Accept a fail-closed review manifest and update local continuation outputs."""
    queue = pd.read_csv(queue_csv, dtype=str).fillna("")
    baseline = pd.read_csv(baseline_checkpoint_csv, dtype=str).fillna("")
    providers = pd.read_csv(provider_checkpoint_csv, dtype=str).fillna("")
    candidates = pd.read_csv(candidate_csv, dtype=str).fillna("")
    accepted_decisions: list[dict[str, str]] = []
    for review_path in review_json:
        payload = json.loads(review_path.read_text(encoding="utf-8"))
        accepted_decisions.extend(_review_decision_rows(payload))
    review_payload = {"accepted": accepted_decisions}
    prior = (
        pd.read_csv(prior_records_csv, dtype=str).fillna("")
        if prior_records_csv is not None
        else None
    )
    source_paths = [
        queue_csv,
        baseline_checkpoint_csv,
        provider_checkpoint_csv,
        candidate_csv,
        *review_json,
    ]
    if prior_records_csv is not None:
        source_paths.append(prior_records_csv)
    source_metadata = {
        str(path.resolve()): {
            "bytes": path.stat().st_size,
            "sha256": _file_sha256(path),
        }
        for path in source_paths
    }
    result = materialize_reviewed_cache(
        queue=queue,
        baseline_checkpoint=baseline,
        provider_checkpoint=providers,
        candidates=candidates,
        review_payload=review_payload,
        prior_records=prior,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    records_path = output_dir / "new_reproductive_assurance_records.csv"
    provider_path = output_dir / "provider_checkpoint.csv"
    axes_path = output_dir / "three_axis_checkpoint.csv.gz"
    audit_path = output_dir / "cached_reproductive_assurance_review_audit.csv"
    next_path = output_dir / "reproductive_assurance_priority_next.csv"
    summary_path = output_dir / "coverage_summary.json"
    _atomic_write_csv(result["records"], records_path)
    _atomic_write_csv(result["provider_checkpoint"], provider_path)
    _atomic_write_csv(result["three_axis_checkpoint"], axes_path, gzip=True)
    _atomic_write_csv(result["review_audit"], audit_path)
    _atomic_write_csv(result["next_queue"], next_path)
    remaining = len(result["next_queue"])
    counted_next_path = output_dir / f"reproductive_assurance_priority_{remaining}.csv"
    _atomic_write_csv(result["next_queue"], counted_next_path)

    output_paths = [
        records_path,
        provider_path,
        axes_path,
        audit_path,
        next_path,
        counted_next_path,
    ]
    summary = dict(result["summary"])
    summary["materialized_at_utc"] = datetime.now(UTC).isoformat().replace("+00:00", "Z")
    summary["sources"] = source_metadata
    summary["outputs"] = {
        path.name: {
            "path": str(path.resolve()),
            "bytes": path.stat().st_size,
            "sha256": _file_sha256(path),
        }
        for path in output_paths
    }
    _atomic_write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        summary_path,
    )
    typer.echo(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


@app.command("build-reviewed-search-batch")
def build_reviewed_search_batch_command(
    queue_csv: Path = typer.Option(..., exists=True),
    search_evidence_csv: Path = typer.Option(..., exists=True),
    provider_checkpoint_csv: Path = typer.Option(..., exists=True),
    selection_json: list[Path] = typer.Option(..., exists=True),
    output_candidate_csv: Path = typer.Option(...),
    output_review_json: Path = typer.Option(...),
) -> None:
    """Build deterministic candidates from reviewed rows in a local cache."""
    queue = pd.read_csv(queue_csv, dtype=str).fillna("")
    search = pd.read_csv(search_evidence_csv, dtype=str).fillna("")
    providers = pd.read_csv(provider_checkpoint_csv, dtype=str).fillna("")
    source_hash = _file_sha256(search_evidence_csv)
    accepted_selections: list[dict[str, object]] = []
    for selection_path in selection_json:
        payload = json.loads(selection_path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError("selection manifest must contain a JSON object")
        claimed_hash = _text(payload.get("source_file_sha256"))
        if claimed_hash and claimed_hash.casefold() != source_hash:
            raise ValueError("selection manifest source hash does not match search cache")
        selected = payload.get("accepted", [])
        if not isinstance(selected, list):
            raise ValueError("selection manifest accepted must be a list")
        accepted_selections.extend(selected)
    selection_payload = {"accepted": accepted_selections}
    candidates, review_payload = build_reviewed_search_batch(
        queue=queue,
        search_evidence=search,
        provider_checkpoint=providers,
        selection_payload=selection_payload,
        local_file_path=str(search_evidence_csv.resolve()),
        local_file_hash=source_hash,
    )
    _validate_reviewed_candidates(
        queue,
        candidates,
        _review_decision_rows(review_payload),
    )
    review_payload["source"] = {
        "path": str(search_evidence_csv.resolve()),
        "bytes": search_evidence_csv.stat().st_size,
        "sha256": source_hash,
    }
    review_payload["selection_manifests"] = {
        str(path.resolve()): {
            "bytes": path.stat().st_size,
            "sha256": _file_sha256(path),
        }
        for path in selection_json
    }
    _atomic_write_csv(candidates, output_candidate_csv)
    _atomic_write_text(
        json.dumps(review_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        output_review_json,
    )
    report = {
        "contract_version": "reviewed_search_cache_batch_v1",
        "n_candidates": int(len(candidates)),
        "n_trait_records": int(len(review_payload["accepted"])),
        "n_species": int(candidates["accepted_species"].nunique()),
        "candidate_output": {
            "path": str(output_candidate_csv.resolve()),
            "bytes": output_candidate_csv.stat().st_size,
            "sha256": _file_sha256(output_candidate_csv),
        },
        "review_output": {
            "path": str(output_review_json.resolve()),
            "bytes": output_review_json.stat().st_size,
            "sha256": _file_sha256(output_review_json),
        },
        "external_acquisition_started": bool(
            set(candidates["provider"].map(_text)).intersection(EXTERNAL_ACQUISITION_PROVIDERS)
        ),
        "biological_values_materialized": False,
    }
    typer.echo(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))


@app.command("augment-cleistogamy-cache")
def augment_cleistogamy_cache_command(
    queue_csv: Path = typer.Option(..., exists=True),
    search_evidence_csv: Path = typer.Option(..., exists=True),
    legacy_provider_csv: Path = typer.Option(..., exists=True),
    existing_candidate_csv: Path = typer.Option(..., exists=True),
    output_csv: Path = typer.Option(...),
    source_run_id: int = typer.Option(29637034874, min=1),
) -> None:
    """Surface cached cleistogamy text for strict review without accepting it."""
    queue = pd.read_csv(queue_csv, dtype=str).fillna("")
    search = pd.read_csv(search_evidence_csv, dtype=str).fillna("")
    legacy = pd.read_csv(legacy_provider_csv, dtype=str).fillna("")
    existing = pd.read_csv(existing_candidate_csv, dtype=str).fillna("")
    augmented = augment_cleistogamy_cache_candidates(
        queue=queue,
        search_evidence=search,
        legacy_provider=legacy,
        existing_candidates=existing,
        local_file_path=str(search_evidence_csv.resolve()),
        local_file_hash=_file_sha256(search_evidence_csv),
        source_run_id=source_run_id,
    )
    _atomic_write_csv(augmented, output_csv)
    report = {
        "contract_version": "local_cleistogamy_cache_augmentation_v1",
        "n_existing_candidates": int(len(existing)),
        "n_augmented_candidates": int(len(augmented)),
        "n_new_candidates": int(len(augmented) - len(existing)),
        "n_new_species": int(
            augmented.loc[
                ~augmented["candidate_id"].isin(existing["candidate_id"]), "accepted_species"
            ].nunique()
        ),
        "output": {
            "path": str(output_csv.resolve()),
            "bytes": output_csv.stat().st_size,
            "sha256": _file_sha256(output_csv),
        },
        "biological_values_accepted": False,
    }
    typer.echo(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))


def refresh_coverage_summary(
    *,
    prior_summary: dict[str, object],
    axes: pd.DataFrame,
    queue: pd.DataFrame,
    providers: pd.DataFrame,
    records: pd.DataFrame,
    acquisition_status_paths: list[Path],
    output_paths: list[Path],
    no_new_direct_evidence_in_entire_pass: bool = False,
) -> dict[str, object]:
    """Refresh coverage/checkpoint metadata without changing biological values."""
    required_axes = {
        "accepted_species",
        "colour_vividness_direct",
        "floral_complexity_direct",
        "reproductive_assurance_direct",
        "all_three_axes_direct",
        "reproductive_assurance_priority",
    }
    missing_axes = required_axes.difference(axes.columns)
    if missing_axes:
        raise ValueError(f"three-axis checkpoint missing columns: {sorted(missing_axes)}")
    if axes["accepted_species"].map(_text).duplicated().any():
        raise ValueError("three-axis checkpoint contains duplicate species")
    queue_species = set(_validate_queue(queue))
    expected_queue = set(
        axes.loc[axes["reproductive_assurance_priority"].map(_bool), "accepted_species"].map(_text)
    )
    if queue_species != expected_queue:
        raise ValueError("priority queue does not match the three-axis checkpoint")
    required_providers = set(CHECKPOINT_COLUMNS)
    missing_providers = required_providers.difference(providers.columns)
    if missing_providers:
        raise ValueError(f"provider checkpoint missing columns: {sorted(missing_providers)}")
    if providers.duplicated(["species", "trait", "provider"]).any():
        raise ValueError("provider checkpoint contains duplicate keys")
    missing_records = set(RECORD_COLUMNS).difference(records.columns)
    if missing_records:
        raise ValueError(f"record table missing columns: {sorted(missing_records)}")
    if records["record_id"].map(_text).duplicated().any():
        raise ValueError("record table contains duplicate record IDs")

    acquisition_batches: list[dict[str, object]] = []
    for path in acquisition_status_paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict) or not _text(payload.get("provider")):
            raise ValueError(f"invalid acquisition status file: {path}")
        acquisition_batches.append(
            {
                "provider": _text(payload.get("provider")),
                "provider_version": _text(payload.get("provider_version")),
                "batch_id": _text(payload.get("batch_id")),
                "requested_species": _nonnegative_int(payload.get("requested_species")),
                "source_records": _nonnegative_int(payload.get("source_records")),
                "candidate_leads": _nonnegative_int(payload.get("candidate_leads")),
                "transient_errors": _nonnegative_int(payload.get("transient_errors")),
                "permanent_errors": _nonnegative_int(payload.get("permanent_errors")),
                "path": str(path.resolve()),
                "bytes": path.stat().st_size,
                "sha256": _file_sha256(path),
            }
        )

    result = dict(prior_summary)
    colour = axes["colour_vividness_direct"].map(_bool)
    complexity = axes["floral_complexity_direct"].map(_bool)
    reproductive = axes["reproductive_assurance_direct"].map(_bool)
    all_three = axes["all_three_axes_direct"].map(_bool)
    terminal = providers["terminal"].map(_bool)
    cumulative_completed = sorted({species for species in records["species"].map(_text) if species})
    latest_completed_names = prior_summary.get(
        "latest_batch_newly_completed_species_names",
        prior_summary.get("newly_completed_species_names", []),
    )
    if not isinstance(latest_completed_names, list):
        latest_completed_names = []
    latest_completed_names = sorted(
        {_text(species) for species in latest_completed_names if _text(species)}
    )
    latest_completed_count = _nonnegative_int(
        prior_summary.get(
            "latest_batch_newly_completed_species",
            prior_summary.get("newly_completed_species", len(latest_completed_names)),
        )
    )
    result.update(
        {
            "contract_version": "local_reproductive_assurance_coverage_summary_v2",
            "n_species": int(len(axes)),
            "direct_coverage": {
                "colour_vividness": int(colour.sum()),
                "floral_complexity": int(complexity.sum()),
                "reproductive_assurance": int(reproductive.sum()),
            },
            "remaining_species": {
                "priority_queue": int(len(queue)),
                "reproductive_assurance": int((~reproductive).sum()),
            },
            "total_three_axis_completion": int(all_three.sum()),
            "total_record_count": int(len(records)),
            "newly_completed_species": int(len(cumulative_completed)),
            "newly_completed_species_names": cumulative_completed,
            "cumulative_newly_completed_species": int(len(cumulative_completed)),
            "cumulative_newly_completed_species_names": cumulative_completed,
            "latest_batch_newly_completed_species": latest_completed_count,
            "latest_batch_newly_completed_species_names": latest_completed_names,
            "provider_checkpoint_rows": int(len(providers)),
            "provider_terminal_count": int(terminal.sum()),
            "provider_remaining_count": int((~terminal).sum()),
            "external_acquisition_started": bool(
                prior_summary.get("external_acquisition_started")
                or acquisition_batches
                or (
                    providers["provider"].isin(EXTERNAL_ACQUISITION_PROVIDERS)
                    & providers["attempts"].map(_nonnegative_int).gt(0)
                ).any()
            ),
            "no_new_direct_evidence_in_entire_pass": bool(no_new_direct_evidence_in_entire_pass),
            "global_fallback_used": False,
            "genus_or_family_inference_used": False,
            "coverage_refreshed_at_utc": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
            "acquisition_batches": acquisition_batches,
        }
    )
    result["provider_status_counts"] = [
        {
            "provider": _text(provider),
            "status": _text(status),
            "terminal": bool(is_terminal),
            "count": int(len(group)),
        }
        for (provider, status, is_terminal), group in providers.assign(
            terminal_bool=terminal
        ).groupby(["provider", "status", "terminal_bool"], sort=True)
    ]
    outputs = dict(prior_summary.get("outputs") or {})
    for path in output_paths:
        outputs[path.name] = {
            "path": str(path.resolve()),
            "bytes": path.stat().st_size,
            "sha256": _file_sha256(path),
        }
    result["outputs"] = outputs
    return result


@app.command("refresh-coverage-summary")
def refresh_coverage_summary_command(
    summary_json: Path = typer.Option(..., exists=True),
    axes_csv: Path = typer.Option(..., exists=True),
    queue_csv: Path = typer.Option(..., exists=True),
    provider_checkpoint_csv: Path = typer.Option(..., exists=True),
    records_csv: Path = typer.Option(..., exists=True),
    acquisition_status_json: list[Path] | None = typer.Option(
        None, "--acquisition-status-json", exists=True
    ),
    output_path: Path | None = typer.Option(None),
    no_new_direct_evidence_in_entire_pass: bool = typer.Option(False),
) -> None:
    """Recompute coverage metadata after a provider-only local batch."""
    summary = json.loads(summary_json.read_text(encoding="utf-8"))
    if not isinstance(summary, dict):
        raise typer.BadParameter("coverage summary must contain a JSON object")
    axes = pd.read_csv(axes_csv, dtype=str).fillna("")
    queue = pd.read_csv(queue_csv, dtype=str).fillna("")
    providers = pd.read_csv(provider_checkpoint_csv, dtype=str).fillna("")
    records = pd.read_csv(records_csv, dtype=str).fillna("")
    paths = [axes_csv, queue_csv, provider_checkpoint_csv, records_csv]
    refreshed = refresh_coverage_summary(
        prior_summary=summary,
        axes=axes,
        queue=queue,
        providers=providers,
        records=records,
        acquisition_status_paths=list(acquisition_status_json or []),
        output_paths=paths,
        no_new_direct_evidence_in_entire_pass=no_new_direct_evidence_in_entire_pass,
    )
    destination = output_path or summary_json
    _atomic_write_text(
        json.dumps(refreshed, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        destination,
    )
    typer.echo(json.dumps(refreshed, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    app()
