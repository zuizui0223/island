"""Acquire conservative species-direct SI evidence from Ferrer et al.

The published workbook mixes directly supported classifications with 725 SI
assignments transferred from congeners.  This adapter therefore reconstructs
the direct subset from the structured source fields instead of importing the
headline ``SI analyses`` column wholesale:

* SC requires an explicit species row plus at least one cited reference.
* SI and SC-SI require a phenotype row with an SI reaction and phenotype
  reference.
* SS and SC-SS are never mapped to self-incompatibility.
* conflicting rows that resolve to one accepted species are quarantined.

No genus inference is emitted here.  The common all-evidence, trait-specific
Validated Low implementation is the only downstream inference path.
"""

from __future__ import annotations

import hashlib
import html
import json
import re
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import requests
import typer
import yaml
from bs4 import BeautifulSoup

from island_v2.angiosperm_scope import classify_scope
from island_v2.angiosperm_scope import load_config as load_scope_config
from island_v2.meyer_2026_reproductive_traits import (
    CANDIDATE_COLUMNS as BASE_CANDIDATE_COLUMNS,
)
from island_v2.meyer_2026_reproductive_traits import (
    fetch_gbif_records,
    strict_gbif_resolution,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

CANDIDATE_COLUMNS = [*BASE_CANDIDATE_COLUMNS, "source_lineage"]
REFERENCE_RE = re.compile(r"\s*;\s*|\s*\|\s*")


@app.callback()
def main() -> None:
    """Download and prepare the Ferrer species-level SI evidence package."""


def _text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _key(value: Any) -> str:
    return _text(value).casefold()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise typer.BadParameter(f"Expected a mapping in {path}")
    return value


def _write_csv_gzip(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})


def _normalise_reference(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", " ", _text(value).casefold()).strip()


def _unique_references(values: Iterable[Any]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        for item in REFERENCE_RE.split(_text(value)):
            item = _text(item)
            key = _normalise_reference(item)
            if item and key and key not in seen:
                seen.add(key)
                result.append(item)
    return result


def _citation_set_lineage(citations: Iterable[str]) -> str:
    normalized = sorted({_normalise_reference(value) for value in citations if _text(value)})
    digest = hashlib.sha256("|".join(normalized).encode("utf-8")).hexdigest()[:24]
    return f"citation-set:{digest}"


def download_source_workbook(
    destination: Path,
    *,
    config_path: Path = Path("config/ferrer_2024_si_traits.yml"),
    timeout_seconds: float = 90.0,
) -> dict[str, str]:
    """Resolve the expiring OUP attachment URL from a stable article page."""
    config = _load_yaml(config_path)
    source = config["source"]
    landing_url = _text(source["attachment_resolver_url"])
    filename = _text(source["attachment_filename"])
    headers = {"User-Agent": "island-floral-v2/0.1 (academic trait acquisition)"}
    page = requests.get(landing_url, headers=headers, timeout=timeout_seconds)
    page.raise_for_status()
    soup = BeautifulSoup(page.text, "html.parser")
    urls = [
        html.unescape(str(link.get("href", "")))
        for link in soup.find_all("a")
        if filename in str(link.get("href", ""))
    ]
    if len(set(urls)) != 1:
        raise typer.BadParameter(
            f"Expected one {filename} attachment link at {landing_url}, found {len(set(urls))}"
        )
    response = requests.get(urls[0], headers=headers, timeout=timeout_seconds)
    response.raise_for_status()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(response.content)
    observed = _sha256(destination)
    expected = _text(source["expected_workbook_sha256"]).casefold()
    if observed != expected:
        destination.unlink(missing_ok=True)
        raise typer.BadParameter(
            f"Source SHA-256 mismatch: expected {expected}, observed {observed}"
        )
    return {
        "destination": str(destination),
        "sha256": observed,
        "resolver_url": landing_url,
        "resolved_attachment_url": urls[0],
    }


def _safe_source_rows(workbook: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    references = pd.read_excel(workbook, sheet_name="References", dtype=str).fillna("")
    phenotype = pd.read_excel(
        workbook, sheet_name="SI Phenotype data", dtype=str
    ).fillna("")
    required_references = {"SI_orden", "Family", "Genera", "work_species", "SI analyses"}
    required_phenotype = {
        "SI_orden",
        "SI status",
        "SI reaction",
        "Notes BS",
        "Ref BS",
        "Ref Gen",
    }
    if missing := required_references.difference(references.columns):
        raise typer.BadParameter(f"References sheet lacks columns: {sorted(missing)}")
    if missing := required_phenotype.difference(phenotype.columns):
        raise typer.BadParameter(f"SI Phenotype data sheet lacks columns: {sorted(missing)}")

    citation_columns = [
        column
        for column in references.columns
        if str(column).startswith("Review") or str(column).startswith("Reference")
    ]
    phenotype_by_id = {
        _text(record["SI_orden"]): record for record in phenotype.to_dict("records")
    }
    accepted: list[dict[str, Any]] = []
    exclusions: list[dict[str, Any]] = []
    for row_number, record in enumerate(references.to_dict("records"), start=2):
        source_id = _text(record["SI_orden"])
        source_name = _text(record["work_species"])
        source_value = _text(record["SI analyses"])
        phenotype_record = phenotype_by_id.get(source_id, {})
        cited = _unique_references(record.get(column) for column in citation_columns)
        phenotype_cited = _unique_references(
            (phenotype_record.get("Ref BS"), phenotype_record.get("Ref Gen"))
        )
        phenotype_status = _text(phenotype_record.get("SI status"))
        reaction = _text(phenotype_record.get("SI reaction"))
        phenotype_support = phenotype_cited or cited

        candidates: list[tuple[str, str, list[str]]] = []
        if source_value == "SC" and cited:
            candidates.append(("SC", "reference_supported_SC", cited))
        if phenotype_status == "SI" and reaction and phenotype_support:
            candidates.append(("SI", "phenotype_confirmed_SI", phenotype_support))
        if phenotype_status == "SC-SI" and reaction and phenotype_support:
            candidates.append(
                ("mixed_or_variable", "phenotype_confirmed_partial_SI", phenotype_support)
            )

        if not candidates:
            if source_value in {"SS", "SC-SS"}:
                reason = "self_sterility_not_self_incompatibility"
            elif source_value == "SI":
                reason = "dataset_level_SI_without_confirming_phenotype"
            elif source_value in {"SC-SI", "SI-SC"}:
                reason = "mixed_status_without_confirming_phenotype"
            elif source_value == "SC":
                reason = "SC_without_cited_reference"
            else:
                reason = "unsupported_source_classification"
            exclusions.append(
                {
                    "row_number": row_number,
                    "source_record_id": source_id,
                    "source_scientific_name": source_name,
                    "source_value": source_value,
                    "phenotype_status": phenotype_status,
                    "phenotype_reaction": reaction,
                    "reason": reason,
                }
            )
            continue

        for value, evidence_class, direct_citations in candidates:
            citations = _unique_references([*cited, *direct_citations])
            excerpt_fields = [
                f"work_species={source_name}",
                f"SI analyses={source_value}",
            ]
            if phenotype_status:
                excerpt_fields.extend(
                    [
                        f"SI status={phenotype_status}",
                        f"SI reaction={reaction}",
                    ]
                )
                notes = _text(phenotype_record.get("Notes BS"))
                if notes:
                    excerpt_fields.append(f"Notes BS={notes}")
            excerpt_fields.append("references=" + " | ".join(citations))
            accepted.append(
                {
                    "row_number": row_number,
                    "source_record_id": source_id,
                    "source_scientific_name": source_name,
                    "source_family": _text(record["Family"]),
                    "source_genus": _text(record["Genera"]),
                    "source_value": source_value,
                    "standardized_value": value,
                    "evidence_class": evidence_class,
                    "source_citation": " | ".join(citations),
                    "evidence_excerpt": "; ".join(excerpt_fields),
                    "source_lineage": _citation_set_lineage(citations),
                }
            )
    return pd.DataFrame(accepted), pd.DataFrame(exclusions)


def _meyer_values(path: Path | None) -> dict[str, str]:
    if path is None or not path.is_file():
        return {}
    frame = pd.read_csv(path, dtype=str).fillna("")
    values: dict[str, str] = {}
    for record in frame.to_dict("records"):
        species = _text(record.get("accepted_species"))
        value = _text(record.get("standardized_value"))
        if species and value in {"SC", "SI"}:
            values[_key(species)] = value
    return values


def _candidate(
    record: Mapping[str, Any],
    *,
    accepted_species: str,
    genus: str,
    family: str,
    match_method: str,
    config: Mapping[str, Any],
    meyer_values: Mapping[str, str],
) -> dict[str, Any]:
    source = config["source"]
    value = _text(record["standardized_value"])
    lineage = _text(record["source_lineage"])
    if meyer_values.get(_key(accepted_species)) == value:
        lineage = "doi:10.5061/dryad.cc2fqz6hr"
    source_id = _text(record["source_record_id"])
    return {
        "evidence_id": f"ferrer2024:{source_id}:{accepted_species}:{value}",
        "source_name": source["name"],
        "source_record_id": f"supplementary_data.xlsx:{source_id}",
        "accepted_species": accepted_species,
        "genus": genus,
        "family": family,
        "name_match_method": match_method,
        "trait_name": "self_incompatibility",
        "raw_value": _text(record["source_value"]),
        "standardized_value": value,
        "candidate_kind": "source_backed",
        "evidence_scope": (
            "species_direct" if match_method == "accepted_name_exact" else "synonym_direct"
        ),
        "source_type": "curated_literature_database",
        "source_reliability": "B_curated_database_or_institution",
        "source_citation": (
            f"{source['citation']} Underlying source(s): "
            f"{_text(record['source_citation'])}"
        ),
        "source_url": source["dataset_url"],
        "evidence_excerpt": _text(record["evidence_excerpt"]),
        "supporting_taxa": "[]",
        "inference_rule": "",
        "confidence": "medium",
        "inference_note": (
            f"{record['evidence_class']}; source-level genus-transferred SI, SS and SC-SS "
            "were excluded before name reconciliation."
        ),
        "needs_human_review": True,
        "review_status": "pending",
        "analysis_role": "confirmatory_species_direct_after_conflict_resolution",
        "source_lineage": lineage,
    }


def prepare_ferrer_2024(
    *,
    source_xlsx: Path,
    master_csv: Path,
    output_dir: Path,
    config_path: Path = Path("config/ferrer_2024_si_traits.yml"),
    scope_config_path: Path = Path("config/angiosperm_scope.yml"),
    meyer_candidates_csv: Path | None = Path(
        "data/v2/staging/traits/bulk/meyer_2026_reproductive_traits/trait_candidates.csv.gz"
    ),
    resolve_gbif: bool = True,
    gbif_workers: int = 8,
) -> dict[str, Any]:
    """Validate, reconcile and materialize the conservative direct subset."""
    config = _load_yaml(config_path)
    observed_hash = _sha256(source_xlsx)
    expected_hash = _text(config["source"]["expected_workbook_sha256"]).casefold()
    if observed_hash != expected_hash:
        raise typer.BadParameter(
            f"Source SHA-256 mismatch: expected {expected_hash}, observed {observed_hash}"
        )

    safe, exclusions = _safe_source_rows(source_xlsx)
    expected_rows = int(config["source"]["expected_source_rows"])
    source_rows = pd.read_excel(
        source_xlsx, sheet_name="References", usecols=["SI_orden"], dtype=str
    )
    if len(source_rows) != expected_rows:
        raise typer.BadParameter(f"Expected {expected_rows} source rows, found {len(source_rows)}")

    master = pd.read_csv(master_csv, dtype=str).fillna("")
    required_master = {"accepted_species", "genus", "family"}
    if missing := required_master.difference(master.columns):
        raise typer.BadParameter(f"Master lacks columns: {sorted(missing)}")
    master["accepted_species"] = master["accepted_species"].map(_text)
    master = master.loc[master["accepted_species"].ne("")].drop_duplicates("accepted_species")
    scope = classify_scope(master, load_scope_config(scope_config_path))
    eligible_names = set(
        scope.loc[scope["angiosperm_analysis_eligible"], "accepted_species"].map(_text)
    )
    eligible_master = master.loc[master["accepted_species"].isin(eligible_names)].copy()
    master_by_name = {
        _key(record["accepted_species"]): record
        for record in eligible_master[["accepted_species", "genus", "family"]].to_dict("records")
    }
    output_dir.mkdir(parents=True, exist_ok=True)

    resolution_rows: list[dict[str, Any]] = []
    resolved: list[dict[str, Any]] = []
    unmatched_indices: list[int] = []
    for index, record in safe.iterrows():
        source_name = _text(record["source_scientific_name"])
        master_row = master_by_name.get(_key(source_name))
        if master_row is None:
            unmatched_indices.append(index)
            continue
        source_family = _text(record["source_family"])
        master_family = _text(master_row["family"])
        if source_family and master_family and _key(source_family) != _key(master_family):
            reason = "reject_family_conflict"
            accepted_species = ""
        else:
            reason = "accepted_name_exact"
            accepted_species = _text(master_row["accepted_species"])
            resolved.append(
                {
                    **record.to_dict(),
                    "accepted_species": accepted_species,
                    "match_method": "accepted_name_exact",
                }
            )
        resolution_rows.append(
            {
                "source_scientific_name": source_name,
                "source_family": source_family,
                "standardized_value": _text(record["standardized_value"]),
                "accepted_species": accepted_species,
                "name_match_method": "accepted_name_exact" if accepted_species else "unmatched",
                "match_status": "matched" if accepted_species else "unmatched",
                "reason": reason,
            }
        )

    gbif_cache_path = output_dir / "gbif_species_match_cache.jsonl.gz"
    gbif_records: dict[str, dict[str, Any]] = {}
    if resolve_gbif and unmatched_indices:
        unresolved = safe.loc[unmatched_indices].rename(
            columns={"source_scientific_name": "genus_species", "source_family": "Family"}
        )
        # Persist bounded batches so an interrupted run can resume without
        # repeating every completed name request.
        batch_size = 200
        for start in range(0, len(unresolved), batch_size):
            gbif_records.update(
                fetch_gbif_records(
                    unresolved.iloc[start : start + batch_size],
                    endpoint=_text(config["gbif_reconciliation"]["endpoint"]),
                    cache_path=gbif_cache_path,
                    workers=gbif_workers,
                )
            )

    for index in unmatched_indices:
        record = safe.loc[index].to_dict()
        source_name = _text(record["source_scientific_name"])
        cached = gbif_records.get(source_name, {})
        payload = cached.get("payload") if isinstance(cached, dict) else {}
        payload = payload if isinstance(payload, dict) else {}
        if resolve_gbif and not _text(cached.get("error")):
            accepted_species, reason = strict_gbif_resolution(
                source_name=source_name,
                source_family=_text(record["source_family"]),
                payload=payload,
                master_by_name=master_by_name,
                reconciliation=config["gbif_reconciliation"],
            )
        elif resolve_gbif:
            accepted_species, reason = "", "gbif_request_error"
        else:
            accepted_species, reason = "", "gbif_resolution_disabled"
        if accepted_species:
            resolved.append(
                {
                    **record,
                    "accepted_species": accepted_species,
                    "match_method": "gbif_exact_synonym_to_accepted",
                }
            )
        resolution_rows.append(
            {
                "source_scientific_name": source_name,
                "source_family": _text(record["source_family"]),
                "standardized_value": _text(record["standardized_value"]),
                "accepted_species": accepted_species,
                "name_match_method": (
                    "gbif_exact_synonym_to_accepted" if accepted_species else "unmatched"
                ),
                "match_status": "matched" if accepted_species else "unmatched",
                "reason": reason,
            }
        )

    resolved_frame = pd.DataFrame(resolved)
    conflict_species: set[str] = set()
    conflict_rows = pd.DataFrame()
    if not resolved_frame.empty:
        value_counts = resolved_frame.groupby("accepted_species")[
            "standardized_value"
        ].nunique()
        conflict_species = set(value_counts.loc[value_counts.gt(1)].index)
        conflict_rows = resolved_frame.loc[
            resolved_frame["accepted_species"].isin(conflict_species)
        ].copy()
        resolved_frame = resolved_frame.loc[
            ~resolved_frame["accepted_species"].isin(conflict_species)
        ].copy()
        resolved_frame["match_priority"] = resolved_frame["match_method"].map(
            {"accepted_name_exact": 0, "gbif_exact_synonym_to_accepted": 1}
        )
        resolved_frame = (
            resolved_frame.sort_values(
                ["accepted_species", "match_priority", "source_record_id"]
            )
            .drop_duplicates("accepted_species", keep="first")
            .reset_index(drop=True)
        )

    meyer_values = _meyer_values(meyer_candidates_csv)
    candidate_rows: list[dict[str, Any]] = []
    for record in resolved_frame.to_dict("records"):
        master_record = master_by_name[_key(record["accepted_species"])]
        candidate_rows.append(
            _candidate(
                record,
                accepted_species=_text(record["accepted_species"]),
                genus=_text(master_record["genus"]),
                family=_text(master_record["family"]),
                match_method=_text(record["match_method"]),
                config=config,
                meyer_values=meyer_values,
            )
        )
    direct = pd.DataFrame(candidate_rows, columns=CANDIDATE_COLUMNS).sort_values(
        ["accepted_species", "standardized_value"]
    )
    match_audit = pd.DataFrame(resolution_rows).sort_values(
        ["source_scientific_name", "standardized_value"]
    )

    direct_path = output_dir / "trait_candidates.csv.gz"
    match_path = output_dir / "name_match_audit.csv.gz"
    unmatched_path = output_dir / "unmatched_taxa.csv.gz"
    exclusions_path = output_dir / "source_exclusions.csv.gz"
    conflicts_path = output_dir / "source_conflicts.csv.gz"
    _write_csv_gzip(direct, direct_path)
    _write_csv_gzip(match_audit, match_path)
    _write_csv_gzip(
        match_audit.loc[match_audit["match_status"].eq("unmatched")], unmatched_path
    )
    _write_csv_gzip(exclusions, exclusions_path)
    _write_csv_gzip(conflict_rows, conflicts_path)

    report: dict[str, Any] = {
        "contract": "ferrer_2024_conservative_species_direct_si_v1",
        "source_name": config["source"]["name"],
        "source_path": str(source_xlsx),
        "source_sha256": observed_hash,
        "source_url": config["source"]["dataset_url"],
        "source_license": config["source"]["license"],
        "n_source_rows": len(source_rows),
        "n_safe_direct_rows_before_name_matching": len(safe),
        "safe_direct_rows_by_value": {
            str(key): int(value)
            for key, value in safe["standardized_value"].value_counts().items()
        },
        "n_direct_candidate_species": int(direct["accepted_species"].nunique()),
        "direct_candidates_by_value": {
            str(key): int(value)
            for key, value in direct["standardized_value"].value_counts().items()
        },
        "direct_candidates_by_name_match_method": {
            str(key): int(value)
            for key, value in direct["name_match_method"].value_counts().items()
        },
        "n_source_rows_excluded_from_direct": len(exclusions),
        "source_exclusions_by_reason": {
            str(key): int(value) for key, value in exclusions["reason"].value_counts().items()
        },
        "n_conflicting_master_mappings_quarantined": len(conflict_species),
        "conflicting_master_species": sorted(conflict_species),
        "n_unmatched_safe_rows": int(match_audit["match_status"].eq("unmatched").sum()),
        "n_lineages_bridged_to_meyer": int(
            direct["source_lineage"].eq("doi:10.5061/dryad.cc2fqz6hr").sum()
        ),
        "genus_inference_emitted": False,
        "outputs": {
            "trait_candidates": str(direct_path),
            "name_match_audit": str(match_path),
            "unmatched_taxa": str(unmatched_path),
            "source_exclusions": str(exclusions_path),
            "source_conflicts": str(conflicts_path),
            "gbif_match_cache": str(gbif_cache_path) if gbif_cache_path.exists() else "",
        },
        "policy": (
            "Only reference-supported SC and phenotype-confirmed SI/SC-SI are species-direct. "
            "SS, SC-SS and dataset-level genus-transferred SI are excluded. Conflicting values "
            "are quarantined; all genus Low is rebuilt downstream from all accepted direct "
            "evidence on genus x trait_name."
        ),
    }
    manifest_path = output_dir / "acquisition_manifest.json"
    manifest_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    report["manifest"] = str(manifest_path)
    return report


@app.command("download")
def download(
    output_xlsx: Annotated[Path, typer.Option()],
    config_path: Annotated[
        Path, typer.Option(exists=True)
    ] = Path("config/ferrer_2024_si_traits.yml"),
) -> None:
    """Download the exact hash-pinned supplementary workbook."""
    report = download_source_workbook(output_xlsx, config_path=config_path)
    typer.echo(json.dumps(report, indent=2))


@app.command("prepare")
def prepare(
    source_xlsx: Annotated[Path, typer.Option(exists=True)],
    master_csv: Annotated[Path, typer.Option(exists=True)] = Path(
        "data/v2/staging/gbif/collected/island_taxa.csv"
    ),
    output_dir: Annotated[Path, typer.Option()] = Path(
        "data/v2/staging/traits/bulk/ferrer_2024_si_traits"
    ),
    config_path: Annotated[Path, typer.Option(exists=True)] = Path(
        "config/ferrer_2024_si_traits.yml"
    ),
    scope_config_path: Annotated[Path, typer.Option(exists=True)] = Path(
        "config/angiosperm_scope.yml"
    ),
    meyer_candidates_csv: Annotated[Path | None, typer.Option()] = Path(
        "data/v2/staging/traits/bulk/meyer_2026_reproductive_traits/trait_candidates.csv.gz"
    ),
    resolve_gbif: Annotated[
        bool, typer.Option("--resolve-gbif/--no-resolve-gbif")
    ] = True,
    gbif_workers: Annotated[int, typer.Option(min=1, max=32)] = 8,
) -> None:
    """Build audited direct candidates; emit no independent genus inference."""
    report = prepare_ferrer_2024(
        source_xlsx=source_xlsx,
        master_csv=master_csv,
        output_dir=output_dir,
        config_path=config_path,
        scope_config_path=scope_config_path,
        meyer_candidates_csv=meyer_candidates_csv,
        resolve_gbif=resolve_gbif,
        gbif_workers=gbif_workers,
    )
    typer.echo(
        f"Prepared {report['n_direct_candidate_species']:,} direct species; "
        f"quarantined {report['n_conflicting_master_mappings_quarantined']:,} conflicts and "
        f"left {report['n_unmatched_safe_rows']:,} safe rows unmatched."
    )


if __name__ == "__main__":
    app()
