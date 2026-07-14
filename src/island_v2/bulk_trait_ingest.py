"""Ingest a downloaded plant-trait table into reviewable v2 candidate evidence.

This is the high-throughput acquisition path for the global species master. It
never calls an LLM, never guesses an unmapped source code, and never promotes a
candidate to curated evidence. A source-specific YAML profile declares column
names and any permissible value translations; all unmapped rows survive in a
separate audit file.
"""

from __future__ import annotations

import hashlib
import json
import re
import unicodedata
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.schemas import SourceReliability, SourceType

app = typer.Typer(add_completion=False, no_args_is_help=True)

MASTER_REQUIRED_COLUMNS = {"accepted_species", "genus", "family", "n_islands", "n_records"}
PROFILE_REQUIRED_COLUMNS = {"scientific_name", "trait_name", "trait_value"}
RANK_MARKERS = {"subsp", "subsp.", "ssp", "ssp.", "var", "var.", "f", "f.", "forma", "cv", "cv."}


def _normalise_text(value: Any) -> str:
    """Normalise whitespace and Unicode without attempting taxonomic synonymy."""
    text = unicodedata.normalize("NFKC", str(value or "")).strip()
    return re.sub(r"\s+", " ", text)


def _key(value: Any) -> str:
    return _normalise_text(value).casefold()


def _binomial_key(value: Any) -> str | None:
    """Return genus+epithet only when the input is not infra-specific.

    This permits a harmless match from ``Poa annua L.`` to ``Poa annua`` but
    deliberately refuses ``Poa annua subsp. annua``. The latter is a different
    assertion of taxonomic scope and must be reconciled explicitly upstream.
    """
    text = _normalise_text(value).replace("\u00d7", "x")
    tokens = [token.strip(".,;:()[]") for token in text.split()]
    lowered = [token.casefold() for token in tokens]
    if any(token in RANK_MARKERS for token in lowered):
        return None
    botanical = [token for token in tokens if re.fullmatch(r"[A-Za-z][A-Za-z-]*", token)]
    if len(botanical) < 2:
        return None
    genus, epithet = botanical[:2]
    if epithet.casefold() in {"sp", "spp", "cf", "aff"}:
        return None
    return f"{genus.casefold()} {epithet.casefold()}"


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise typer.BadParameter(f"File does not exist: {path}")
    parsed = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(parsed, dict):
        raise typer.BadParameter(f"Expected a mapping in YAML file: {path}")
    return parsed


def _profile(config: dict[str, Any], profile_name: str) -> dict[str, Any]:
    profiles = config.get("source_profiles") or {}
    profile = profiles.get(profile_name)
    if not isinstance(profile, dict):
        raise typer.BadParameter(
            f"Unknown source profile '{profile_name}'. Available: {', '.join(sorted(profiles)) or 'none'}"
        )
    columns = profile.get("columns") or {}
    missing = PROFILE_REQUIRED_COLUMNS.difference(columns)
    if missing:
        raise typer.BadParameter(
            f"Profile '{profile_name}' lacks required column mappings: {sorted(missing)}"
        )
    return profile


def _load_master(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise typer.BadParameter(f"Master CSV does not exist: {path}")
    master = pd.read_csv(path, dtype=str).fillna("")
    missing = MASTER_REQUIRED_COLUMNS.difference(master.columns)
    if missing:
        raise typer.BadParameter(f"Master CSV lacks columns: {sorted(missing)}")
    master["accepted_species"] = master["accepted_species"].map(_normalise_text)
    master = master.loc[master["accepted_species"].ne("")].copy()
    master["n_islands"] = pd.to_numeric(master["n_islands"], errors="coerce").fillna(0).astype(int)
    master["n_records"] = pd.to_numeric(master["n_records"], errors="coerce").fillna(0).astype(int)
    master = master.drop_duplicates(subset=["accepted_species"], keep="first").reset_index(
        drop=True
    )
    master["_exact_key"] = master["accepted_species"].map(_key)
    master["_binomial_key"] = master["accepted_species"].map(_binomial_key)
    return master


def _unique_index(table: pd.DataFrame, key_column: str) -> dict[str, int]:
    valid = table.loc[table[key_column].notna() & table[key_column].ne("")].copy()
    counts = valid.groupby(key_column).size()
    unique_keys = set(counts.loc[counts.eq(1)].index)
    return {
        str(row[key_column]): int(index)
        for index, row in valid.loc[valid[key_column].isin(unique_keys)].iterrows()
    }


def _match_species(
    source_name: Any, exact_index: dict[str, int], binomial_index: dict[str, int]
) -> tuple[int | None, str]:
    exact = _key(source_name)
    if exact in exact_index:
        return exact_index[exact], "accepted_name_exact"
    binomial = _binomial_key(source_name)
    if binomial and binomial in binomial_index:
        return binomial_index[binomial], "accepted_name_author_stripped"
    return None, "unmatched"


def _optional_value(row: pd.Series, source_column: Any) -> str:
    if not source_column or str(source_column) not in row.index:
        return ""
    return _normalise_text(row.get(str(source_column), ""))


def _trait_name(raw_name: str, profile: dict[str, Any], ontology: dict[str, Any]) -> str | None:
    key = _key(raw_name)
    mapping = {
        _key(source): target for source, target in (profile.get("trait_name_map") or {}).items()
    }
    candidate = mapping.get(key, raw_name)
    candidate = _normalise_text(candidate)
    return candidate if candidate in (ontology.get("traits") or {}) else None


def _trait_value(
    raw_value: str, trait_name: str, profile: dict[str, Any], ontology: dict[str, Any]
) -> str | None:
    trait_map = (profile.get("trait_value_map") or {}).get(trait_name) or {}
    mapping = {_key(source): _normalise_text(target) for source, target in trait_map.items()}
    candidate = mapping.get(_key(raw_value), _normalise_text(raw_value))
    allowed = set(
        ((ontology.get("traits") or {}).get(trait_name) or {}).get("allowed_values") or []
    )
    return candidate if candidate in allowed else None


def _role_by_trait(config: dict[str, Any]) -> dict[str, str]:
    output: dict[str, str] = {}
    for role, traits in (config.get("analysis_roles") or {}).items():
        for trait in traits or []:
            output[str(trait)] = str(role)
    return output


def _validate_source_defaults(profile: dict[str, Any]) -> dict[str, str]:
    defaults = {
        str(key): _normalise_text(value)
        for key, value in (profile.get("source_defaults") or {}).items()
    }
    source_type = defaults.get("source_type", "")
    reliability = defaults.get("source_reliability", "")
    if source_type not in {item.value for item in SourceType}:
        raise typer.BadParameter(f"Invalid profile source_type: {source_type!r}")
    if reliability not in {item.value for item in SourceReliability}:
        raise typer.BadParameter(f"Invalid profile source_reliability: {reliability!r}")
    confidence = defaults.get("confidence", "")
    if confidence not in {"high", "medium", "low"}:
        raise typer.BadParameter(f"Profile confidence must be high, medium, or low: {confidence!r}")
    return defaults


def _candidate_columns() -> list[str]:
    return [
        "evidence_id",
        "source_name",
        "source_record_id",
        "accepted_species",
        "genus",
        "family",
        "name_match_method",
        "trait_name",
        "raw_value",
        "standardized_value",
        "candidate_kind",
        "evidence_scope",
        "source_type",
        "source_reliability",
        "source_citation",
        "source_url",
        "evidence_excerpt",
        "supporting_taxa",
        "inference_rule",
        "confidence",
        "inference_note",
        "needs_human_review",
        "review_status",
        "analysis_role",
    ]


def ingest_bulk_traits(
    *,
    master_csv: Path,
    source_csv: Path,
    output_dir: Path,
    profiles_path: Path,
    profile_name: str,
    source_name: str,
    source_citation_override: str = "",
    source_url_override: str = "",
) -> dict[str, Any]:
    """Ingest one declared bulk source and return a compact run report."""
    config = _load_yaml(profiles_path)
    profile = _profile(config, profile_name)
    defaults = _validate_source_defaults(profile)
    ontology = _load_yaml(Path("config/trait_ontology.yml"))
    master = _load_master(master_csv)
    columns = profile["columns"]
    read = profile.get("read") or {}
    delimiter = str(read.get("delimiter", ","))
    encoding = str(read.get("encoding", "utf-8"))
    if not source_csv.exists():
        raise typer.BadParameter(f"Source CSV does not exist: {source_csv}")
    source = pd.read_csv(source_csv, dtype=str, sep=delimiter, encoding=encoding).fillna("")
    required_source_columns = {str(columns[name]) for name in PROFILE_REQUIRED_COLUMNS}
    missing_source_columns = required_source_columns.difference(source.columns)
    if missing_source_columns:
        raise typer.BadParameter(
            f"Source file lacks columns required by profile '{profile_name}': {sorted(missing_source_columns)}"
        )

    exact_index = _unique_index(master, "_exact_key")
    binomial_index = _unique_index(master, "_binomial_key")
    roles = _role_by_trait(config)
    output_dir.mkdir(parents=True, exist_ok=True)

    candidates: list[dict[str, Any]] = []
    matches: list[dict[str, Any]] = []
    unmapped: list[dict[str, Any]] = []
    source_citation_override = _normalise_text(source_citation_override)
    source_url_override = _normalise_text(source_url_override)

    for row_number, (_, source_row) in enumerate(source.iterrows(), start=2):
        source_scientific_name = _optional_value(source_row, columns["scientific_name"])
        source_trait_name = _optional_value(source_row, columns["trait_name"])
        source_trait_value = _optional_value(source_row, columns["trait_value"])
        master_index, match_method = _match_species(
            source_scientific_name, exact_index, binomial_index
        )
        record_id = _optional_value(source_row, columns.get("source_record_id")) or str(row_number)
        base_audit = {
            "source_name": source_name,
            "source_row_number": row_number,
            "source_record_id": record_id,
            "source_scientific_name": source_scientific_name,
            "source_trait_name": source_trait_name,
            "source_trait_value": source_trait_value,
            "name_match_method": match_method,
        }
        if master_index is None:
            matches.append({**base_audit, "accepted_species": "", "match_status": "unmatched"})
            continue

        master_row = master.loc[master_index]
        matches.append(
            {
                **base_audit,
                "accepted_species": master_row["accepted_species"],
                "family": master_row["family"],
                "match_status": "matched",
            }
        )
        canonical_trait = _trait_name(source_trait_name, profile, ontology)
        if canonical_trait is None:
            unmapped.append(
                {
                    **base_audit,
                    "accepted_species": master_row["accepted_species"],
                    "reason": "unmapped_trait_name",
                }
            )
            continue
        canonical_value = _trait_value(source_trait_value, canonical_trait, profile, ontology)
        if canonical_value is None:
            unmapped.append(
                {
                    **base_audit,
                    "accepted_species": master_row["accepted_species"],
                    "canonical_trait_name": canonical_trait,
                    "reason": "unmapped_or_invalid_trait_value",
                }
            )
            continue

        row_citation = (
            _optional_value(source_row, columns.get("source_citation"))
            or source_citation_override
            or source_name
        )
        row_url = _optional_value(source_row, columns.get("source_url")) or source_url_override
        excerpt = _optional_value(source_row, columns.get("evidence_excerpt"))
        if not excerpt:
            excerpt = f"Bulk source {source_name} record {record_id}: {source_trait_name} = {source_trait_value}."
        candidates.append(
            {
                "evidence_id": f"bulk:{source_name}:{record_id}",
                "source_name": source_name,
                "source_record_id": record_id,
                "accepted_species": master_row["accepted_species"],
                "genus": master_row["genus"],
                "family": master_row["family"],
                "name_match_method": match_method,
                "trait_name": canonical_trait,
                "raw_value": source_trait_value,
                "standardized_value": canonical_value,
                "candidate_kind": "source_backed",
                "evidence_scope": "species_direct",
                "source_type": defaults["source_type"],
                "source_reliability": defaults["source_reliability"],
                "source_citation": row_citation,
                "source_url": row_url,
                "evidence_excerpt": excerpt,
                "supporting_taxa": "[]",
                "inference_rule": "",
                "confidence": defaults["confidence"],
                "inference_note": "Direct species-level bulk-source observation; pending review.",
                "needs_human_review": True,
                "review_status": "pending",
                "analysis_role": roles.get(canonical_trait, "unclassified"),
            }
        )

    candidate_table = pd.DataFrame(candidates, columns=_candidate_columns())
    if not candidate_table.empty:
        candidate_table = candidate_table.drop_duplicates(
            subset=[
                "accepted_species",
                "trait_name",
                "standardized_value",
                "source_name",
                "source_record_id",
            ],
            keep="first",
        )
    match_table = pd.DataFrame(matches)
    unmatched_table = match_table.loc[
        match_table.get("match_status", pd.Series(dtype=str)).eq("unmatched")
    ].copy()
    unmapped_table = pd.DataFrame(unmapped)

    source_slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", source_name).strip("_") or "bulk_source"
    candidate_path = output_dir / f"bulk_trait_candidates_{source_slug}.csv"
    match_path = output_dir / f"bulk_trait_name_matches_{source_slug}.csv"
    unmatched_path = output_dir / f"bulk_trait_unmatched_taxa_{source_slug}.csv"
    unmapped_path = output_dir / f"bulk_trait_unmapped_values_{source_slug}.csv"
    candidate_table.to_csv(candidate_path, index=False)
    match_table.to_csv(match_path, index=False)
    unmatched_table.to_csv(unmatched_path, index=False)
    unmapped_table.to_csv(unmapped_path, index=False)

    total_species = len(master)
    total_pairs = int(master["n_islands"].sum())
    coverage_rows: list[dict[str, Any]] = []
    for trait_name in ontology.get("traits") or {}:
        selected = (
            candidate_table.loc[candidate_table["trait_name"].eq(trait_name)]
            if not candidate_table.empty
            else candidate_table
        )
        species = set(selected["accepted_species"].tolist()) if not selected.empty else set()
        covered_master = master.loc[master["accepted_species"].isin(species)]
        covered_pairs = int(covered_master["n_islands"].sum()) if not covered_master.empty else 0
        coverage_rows.append(
            {
                "source_name": source_name,
                "trait_name": trait_name,
                "analysis_role": roles.get(trait_name, "unclassified"),
                "n_candidate_rows": int(len(selected)),
                "n_species_with_direct_evidence": len(species),
                "n_master_species": total_species,
                "pct_master_species_covered": round(100 * len(species) / total_species, 6)
                if total_species
                else 0.0,
                "n_island_species_pairs_covered": covered_pairs,
                "n_master_island_species_pairs": total_pairs,
                "pct_island_species_pairs_covered": round(100 * covered_pairs / total_pairs, 6)
                if total_pairs
                else 0.0,
            }
        )
    coverage_table = pd.DataFrame(coverage_rows)
    coverage_path = output_dir / f"bulk_trait_coverage_by_trait_{source_slug}.csv"
    coverage_table.to_csv(coverage_path, index=False)

    family_rows: list[dict[str, Any]] = []
    if not candidate_table.empty:
        for (family, trait_name), group in candidate_table.groupby(
            ["family", "trait_name"], dropna=False
        ):
            species = set(group["accepted_species"])
            family_master = master.loc[master["family"].eq(family)]
            family_rows.append(
                {
                    "source_name": source_name,
                    "family": family,
                    "trait_name": trait_name,
                    "analysis_role": roles.get(trait_name, "unclassified"),
                    "n_species_with_direct_evidence": len(species),
                    "n_master_species_in_family": int(len(family_master)),
                    "pct_family_species_covered": round(100 * len(species) / len(family_master), 6)
                    if len(family_master)
                    else 0.0,
                    "n_island_species_pairs_covered": int(
                        family_master.loc[
                            family_master["accepted_species"].isin(species), "n_islands"
                        ].sum()
                    ),
                }
            )
    family_coverage_path = output_dir / f"bulk_trait_coverage_by_family_{source_slug}.csv"
    pd.DataFrame(family_rows).to_csv(family_coverage_path, index=False)

    report = {
        "source_name": source_name,
        "profile_name": profile_name,
        "source_path": str(source_csv),
        "source_sha256": _sha256_file(source_csv),
        "master_path": str(master_csv),
        "n_source_rows": int(len(source)),
        "n_name_matched_rows": int(
            match_table.get("match_status", pd.Series(dtype=str)).eq("matched").sum()
        ),
        "n_unmatched_taxon_rows": int(len(unmatched_table)),
        "n_unmapped_value_rows": int(len(unmapped_table)),
        "n_candidate_rows": int(len(candidate_table)),
        "n_species_with_any_direct_evidence": int(candidate_table["accepted_species"].nunique())
        if not candidate_table.empty
        else 0,
        "n_master_species": total_species,
        "candidate_csv": str(candidate_path),
        "name_match_csv": str(match_path),
        "unmatched_taxa_csv": str(unmatched_path),
        "unmapped_values_csv": str(unmapped_path),
        "coverage_by_trait_csv": str(coverage_path),
        "coverage_by_family_csv": str(family_coverage_path),
        "curation_rule": "All imported observations are source-backed pending candidates; no values are automatically curated.",
        "non_inference_rule": "Unmatched taxa and unmapped source labels/values are retained in audits and never guessed.",
    }
    report_path = output_dir / f"bulk_trait_ingest_manifest_{source_slug}.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    report["manifest_json"] = str(report_path)
    return report


@app.command("ingest")
def ingest(
    master_csv: Path = typer.Option(..., exists=True, help="Current island species master CSV."),
    source_csv: Path = typer.Option(..., exists=True, help="Downloaded long-format trait CSV/TSV."),
    output_dir: Path = typer.Option(
        ..., help="Directory for pending candidates and coverage audits."
    ),
    profiles_path: Path = typer.Option(Path("config/bulk_trait_sources.yml"), exists=True),
    profile_name: str = typer.Option("v2_canonical_long_csv", help="Named source profile in YAML."),
    source_name: str = typer.Option(..., help="Stable short name for this source version."),
    source_citation: str = typer.Option(
        "", help="Dataset-level citation used when a row does not provide one."
    ),
    source_url: str = typer.Option(
        "", help="Dataset-level URL used when a row does not provide one."
    ),
) -> None:
    """Match a bulk source to the global master and write pending evidence plus coverage."""
    report = ingest_bulk_traits(
        master_csv=master_csv,
        source_csv=source_csv,
        output_dir=output_dir,
        profiles_path=profiles_path,
        profile_name=profile_name,
        source_name=source_name,
        source_citation_override=source_citation,
        source_url_override=source_url,
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


@app.command("profiles")
def profiles(
    profiles_path: Path = typer.Option(Path("config/bulk_trait_sources.yml"), exists=True),
) -> None:
    """List source profile names without loading any external dataset."""
    config = _load_yaml(profiles_path)
    for name, profile in sorted((config.get("source_profiles") or {}).items()):
        description = _normalise_text((profile or {}).get("description", ""))
        typer.echo(f"{name}\t{description}")


@app.command("prepare-eol")
def prepare_eol(
    eol_archive: Path = typer.Option(
        ..., exists=True, help="EOL traits_all.zip or extracted archive directory."
    ),
    eol_content_dir: Path | None = typer.Option(
        None,
        exists=True,
        file_okay=False,
        help=(
            "Optional separately extracted archive directory. This avoids Python ZIP "
            "streaming for very large exports while retaining the original archive checksum."
        ),
    ),
    output_dir: Path = typer.Option(
        ..., help="Directory for EOL long CSV and unmapped-term audit."
    ),
    config_path: Path = typer.Option(Path("config/eol_traitbank_terms.yml"), exists=True),
    ontology_path: Path = typer.Option(Path("config/trait_ontology.yml"), exists=True),
    source_name: str = typer.Option(
        "eol_traitbank_all", help="Stable source/release name for manifests."
    ),
    chunksize: int = typer.Option(
        200_000, min=1, help="traits.csv rows streamed per pandas chunk."
    ),
    max_trait_rows: int | None = typer.Option(
        None, min=1, help="Optional smoke-test cap; omit for full archive."
    ),
) -> None:
    """Convert EOL TraitBank all-traits export to canonical v2 long rows."""
    from island_v2.eol_traitbank import prepare_eol_traitbank_long

    report = prepare_eol_traitbank_long(
        eol_archive=eol_archive,
        eol_content_dir=eol_content_dir,
        output_dir=output_dir,
        config_path=config_path,
        ontology_path=ontology_path,
        source_name=source_name,
        chunksize=chunksize,
        max_trait_rows=max_trait_rows,
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
