"""Acquire CC-BY self-compatibility and autofertility indices from Razanajatovo et al.

The source reports continuous indices.  This importer keeps the mapping
conservative: consistently high self-compatibility indices support SC; low
indices do not prove SI.  Autofertility indices are mapped separately to
autonomous-selfing capacity and are never used as mating-system labels.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import pandas as pd
import typer
import yaml

from island_v2.angiosperm_scope import classify_scope, load_config as load_scope_config
from island_v2.meyer_2026_reproductive_traits import (
    CANDIDATE_COLUMNS,
    _key,
    _sha256,
    _text,
    _write_csv_gzip,
    fetch_gbif_records,
    strict_gbif_resolution,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.callback()
def main() -> None:
    """Prepare audited index-derived reproductive candidates."""


def _load_yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise typer.BadParameter(f"Expected a mapping in {path}")
    return value


def _numeric_values(group: pd.DataFrame, columns: list[str]) -> list[float]:
    values: list[float] = []
    for column in columns:
        values.extend(pd.to_numeric(group[column], errors="coerce").dropna().astype(float).tolist())
    return sorted(values)


def classify_self_compatibility(
    values: list[float], *, lower: float, upper: float
) -> str:
    """Map measured compatibility indices without treating low seed set as SI."""
    if not values:
        return ""
    if all(value >= upper for value in values):
        return "SC"
    if any(value > lower for value in values):
        return "mixed_or_variable"
    return ""


def classify_autofertility(values: list[float], *, lower: float, upper: float) -> str:
    """Map measured autofertility independently from compatibility."""
    if not values:
        return ""
    if all(value <= lower for value in values):
        return "absent"
    if all(value >= upper for value in values):
        return "autonomous"
    return "mixed_or_variable"


def _family(value: Any, aliases: Mapping[str, str]) -> str:
    text = _text(value)
    return _text(aliases.get(text, text))


def _candidate(
    *,
    accepted_species: str,
    genus: str,
    family: str,
    trait_name: str,
    standardized_value: str,
    raw_values: list[float],
    source_names: list[str],
    source_references: list[str],
    name_match_method: str,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    source = config["source"]
    record_id = f"{accepted_species}:{trait_name}"
    values_text = ", ".join(f"{value:.6g}" for value in raw_values)
    if trait_name == "self_incompatibility":
        interpretation = (
            "consistently high indices support self-compatibility"
            if standardized_value == "SC"
            else "indices are intermediate or heterogeneous"
        )
    else:
        interpretation = {
            "absent": "indices consistently indicate little/no autonomous selfing",
            "autonomous": "indices consistently indicate high autonomous selfing",
            "mixed_or_variable": "indices are intermediate or heterogeneous",
        }[standardized_value]
    return {
        "evidence_id": f"razanajatovo2016:{record_id}",
        "source_name": source["name"],
        "source_record_id": record_id,
        "accepted_species": accepted_species,
        "genus": genus,
        "family": family,
        "name_match_method": name_match_method,
        "trait_name": trait_name,
        "raw_value": values_text,
        "standardized_value": standardized_value,
        "candidate_kind": "source_backed",
        "evidence_scope": "species_direct",
        "source_type": "curated_trait_database",
        "source_reliability": "B_curated_database_or_institution",
        "source_citation": source["citation"],
        "source_url": source["article_url"],
        "evidence_excerpt": (
            f"Reported index values for {' / '.join(source_names)}: [{values_text}]; "
            f"{interpretation}. Source reference keys: {', '.join(source_references)}."
        ),
        "supporting_taxa": "[]",
        "inference_rule": "",
        "confidence": "medium",
        "inference_note": (
            "Direct species-level experimental index synthesis. Compatibility and "
            "autofertility are mapped independently; low compatibility is not coded SI."
        ),
        "needs_human_review": True,
        "review_status": "pending",
        "analysis_role": "evidence_rich_subset",
    }


def prepare_razanajatovo_2016(
    *,
    workbook_path: Path,
    sheet_csv: Path,
    master_csv: Path,
    output_dir: Path,
    config_path: Path = Path("config/razanajatovo_2016_reproductive_traits.yml"),
    scope_config_path: Path = Path("config/angiosperm_scope.yml"),
    resolve_gbif: bool = True,
    gbif_workers: int = 8,
) -> dict[str, Any]:
    config = _load_yaml(config_path)
    source_config = config["source"]
    for path, expected_key in (
        (workbook_path, "expected_workbook_sha256"),
        (sheet_csv, "expected_sheet_csv_sha256"),
    ):
        observed = _sha256(path)
        expected = _text(source_config[expected_key])
        if observed != expected:
            raise typer.BadParameter(
                f"Source SHA-256 mismatch for {path}: expected {expected}, observed {observed}"
            )

    source = pd.read_csv(sheet_csv, dtype=str).fillna("")
    columns = config["columns"]
    required = {
        columns["scientific_name"],
        columns["family"],
        columns["source_reference"],
        *columns["self_compatibility_indices"],
        *columns["autofertility_indices"],
    }
    missing = required.difference(source.columns)
    if missing:
        raise typer.BadParameter(f"Source sheet CSV lacks columns: {sorted(missing)}")
    if len(source) != int(source_config["expected_rows"]):
        raise typer.BadParameter(
            f"Expected {source_config['expected_rows']} sheet rows, found {len(source)}"
        )
    source["_source_species"] = (
        source[columns["scientific_name"]].map(_text).str.replace("_", " ", regex=False)
    )

    aliases = {_text(key): _text(value) for key, value in config.get("family_aliases", {}).items()}
    lower = float(config["thresholds"]["lower"])
    upper = float(config["thresholds"]["upper"])
    source_species_rows: list[dict[str, Any]] = []
    for source_species, group in source.groupby("_source_species", sort=True):
        sc_values = _numeric_values(group, list(columns["self_compatibility_indices"]))
        af_values = _numeric_values(group, list(columns["autofertility_indices"]))
        sc_value = classify_self_compatibility(sc_values, lower=lower, upper=upper)
        af_value = classify_autofertility(af_values, lower=lower, upper=upper)
        if not sc_value and not af_value:
            continue
        families = sorted(
            {
                _family(value, aliases)
                for value in group[columns["family"]]
                if _text(value)
            }
        )
        references = sorted({_text(value) for value in group[columns["source_reference"]] if _text(value)})
        source_species_rows.append(
            {
                "source_species": source_species,
                "source_family": families[0] if len(families) == 1 else "",
                "family_conflict_in_source": len(families) > 1,
                "source_references": references,
                "sc_values": sc_values,
                "af_values": af_values,
            }
        )
    useful = pd.DataFrame(source_species_rows)

    master = pd.read_csv(master_csv, dtype=str).fillna("")
    required_master = {"accepted_species", "genus", "family"}
    if missing_master := required_master.difference(master.columns):
        raise typer.BadParameter(f"Master lacks columns: {sorted(missing_master)}")
    master = master.loc[master["accepted_species"].map(_text).ne("")].drop_duplicates(
        "accepted_species"
    )
    scope = classify_scope(master, load_scope_config(scope_config_path))
    eligible = set(
        scope.loc[scope["angiosperm_analysis_eligible"], "accepted_species"].map(_text)
    )
    master = master.loc[master["accepted_species"].isin(eligible)].copy()
    master_by_name = {
        _key(row["accepted_species"]): row
        for row in master[["accepted_species", "genus", "family"]].to_dict("records")
    }

    resolved_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    unmatched: list[dict[str, Any]] = []
    for record in useful.to_dict("records"):
        source_name = _text(record["source_species"])
        source_family = _text(record["source_family"])
        master_row = master_by_name.get(_key(source_name))
        if record["family_conflict_in_source"]:
            audit_rows.append(
                {
                    "source_scientific_name": source_name,
                    "source_family": source_family,
                    "accepted_species": "",
                    "match_status": "unmatched",
                    "name_match_method": "unmatched",
                    "reason": "source_family_conflict",
                }
            )
            continue
        if master_row is not None:
            master_family = _text(master_row["family"])
            if source_family and _key(source_family) != _key(master_family):
                audit_rows.append(
                    {
                        "source_scientific_name": source_name,
                        "source_family": source_family,
                        "accepted_species": "",
                        "match_status": "unmatched",
                        "name_match_method": "unmatched",
                        "reason": "accepted_name_family_conflict",
                    }
                )
                continue
            resolved_rows.append(
                {**record, "accepted_species": master_row["accepted_species"], "match_method": "accepted_name_exact"}
            )
            audit_rows.append(
                {
                    "source_scientific_name": source_name,
                    "source_family": source_family,
                    "accepted_species": master_row["accepted_species"],
                    "match_status": "matched",
                    "name_match_method": "accepted_name_exact",
                    "reason": "accepted_name_exact",
                }
            )
        else:
            unmatched.append(record)

    gbif_cache_path = output_dir / "gbif_species_match_cache.jsonl.gz"
    gbif_records: dict[str, dict[str, Any]] = {}
    if resolve_gbif and unmatched:
        query_frame = pd.DataFrame(
            [
                {"genus_species": row["source_species"], "Family": row["source_family"]}
                for row in unmatched
            ]
        )
        gbif_records = fetch_gbif_records(
            query_frame,
            endpoint=_text(config["gbif_reconciliation"]["endpoint"]),
            cache_path=gbif_cache_path,
            workers=gbif_workers,
        )
    for record in unmatched:
        source_name = _text(record["source_species"])
        cached = gbif_records.get(source_name, {})
        payload = cached.get("payload") if isinstance(cached, dict) else {}
        payload = payload if isinstance(payload, dict) else {}
        if resolve_gbif and not _text(cached.get("error")):
            accepted, reason = strict_gbif_resolution(
                source_name=source_name,
                source_family=_text(record["source_family"]),
                payload=payload,
                master_by_name=master_by_name,
                reconciliation=config["gbif_reconciliation"],
            )
        elif resolve_gbif:
            accepted, reason = "", "gbif_request_error"
        else:
            accepted, reason = "", "gbif_resolution_disabled"
        if accepted:
            resolved_rows.append(
                {**record, "accepted_species": accepted, "match_method": "gbif_exact_synonym_to_accepted"}
            )
        audit_rows.append(
            {
                "source_scientific_name": source_name,
                "source_family": record["source_family"],
                "accepted_species": accepted,
                "match_status": "matched" if accepted else "unmatched",
                "name_match_method": "gbif_exact_synonym_to_accepted" if accepted else "unmatched",
                "reason": reason,
                "gbif_match_type": _text(payload.get("matchType")),
                "gbif_status": _text(payload.get("status")),
                "gbif_confidence": _text(payload.get("confidence")),
                "gbif_usage_key": _text(payload.get("usageKey")),
                "gbif_accepted_usage_key": _text(payload.get("acceptedUsageKey")),
            }
        )

    candidates: list[dict[str, Any]] = []
    resolved = pd.DataFrame(resolved_rows)
    for accepted_species, group in resolved.groupby("accepted_species", sort=True):
        master_row = master_by_name[_key(accepted_species)]
        sc_values = sorted(value for values in group["sc_values"] for value in values)
        af_values = sorted(value for values in group["af_values"] for value in values)
        sc_value = classify_self_compatibility(sc_values, lower=lower, upper=upper)
        af_value = classify_autofertility(af_values, lower=lower, upper=upper)
        source_names = sorted(set(group["source_species"]))
        references = sorted({reference for values in group["source_references"] for reference in values})
        methods = sorted(set(group["match_method"]))
        match_method = methods[0] if len(methods) == 1 else "multiple_exact_name_reconciliations"
        if sc_value:
            candidates.append(
                _candidate(
                    accepted_species=accepted_species,
                    genus=_text(master_row["genus"]),
                    family=_text(master_row["family"]),
                    trait_name="self_incompatibility",
                    standardized_value=sc_value,
                    raw_values=sc_values,
                    source_names=source_names,
                    source_references=references,
                    name_match_method=match_method,
                    config=config,
                )
            )
        if af_value:
            candidates.append(
                _candidate(
                    accepted_species=accepted_species,
                    genus=_text(master_row["genus"]),
                    family=_text(master_row["family"]),
                    trait_name="autonomous_selfing_capacity",
                    standardized_value=af_value,
                    raw_values=af_values,
                    source_names=source_names,
                    source_references=references,
                    name_match_method=match_method,
                    config=config,
                )
            )

    candidate_table = pd.DataFrame(candidates, columns=CANDIDATE_COLUMNS).sort_values(
        ["accepted_species", "trait_name"]
    )
    audit = pd.DataFrame(audit_rows).fillna("").sort_values("source_scientific_name")
    output_dir.mkdir(parents=True, exist_ok=True)
    candidate_path = output_dir / "trait_candidates.csv.gz"
    audit_path = output_dir / "name_match_audit.csv.gz"
    unmatched_path = output_dir / "unmatched_taxa.csv"
    _write_csv_gzip(candidate_table, candidate_path)
    _write_csv_gzip(audit, audit_path)
    audit.loc[audit["match_status"].eq("unmatched")].to_csv(unmatched_path, index=False)

    report = {
        "version": "1.0",
        "source_name": source_config["name"],
        "source_license": source_config["license"],
        "workbook_sha256": _sha256(workbook_path),
        "sheet_csv_sha256": _sha256(sheet_csv),
        "n_source_observations": int(len(source)),
        "n_source_species": int(source["_source_species"].nunique()),
        "n_source_species_with_conservative_candidate": int(len(useful)),
        "n_direct_candidate_species": int(candidate_table["accepted_species"].nunique()),
        "n_direct_candidate_cells": int(len(candidate_table)),
        "candidate_counts_by_trait_value": {
            f"{trait}:{value}": int(count)
            for (trait, value), count in candidate_table.groupby(
                ["trait_name", "standardized_value"]
            ).size().items()
        },
        "name_match_counts": {
            str(key): int(value) for key, value in audit["name_match_method"].value_counts().items()
        },
        "n_unmatched_useful_source_species": int(audit["match_status"].eq("unmatched").sum()),
        "outputs": {
            "trait_candidates": str(candidate_path),
            "name_match_audit": str(audit_path),
            "unmatched_taxa": str(unmatched_path),
            "gbif_match_cache": str(gbif_cache_path) if gbif_cache_path.exists() else "",
        },
        "policy": (
            "High compatibility indices support SC; low compatibility indices never become SI. "
            "Autofertility maps only to autonomous_selfing_capacity. Exact synonyms are accepted "
            "only through the strict GBIF species/family contract."
        ),
    }
    manifest_path = output_dir / "acquisition_manifest.json"
    manifest_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    report["manifest"] = str(manifest_path)
    return report


@app.command("prepare")
def prepare(
    workbook_path: Path = typer.Option(..., exists=True),
    sheet_csv: Path = typer.Option(..., exists=True),
    master_csv: Path = typer.Option(
        Path("data/v2/staging/gbif/collected/island_taxa.csv"), exists=True
    ),
    output_dir: Path = typer.Option(
        Path("data/v2/staging/traits/bulk/razanajatovo_2016_reproductive_traits")
    ),
    config_path: Path = typer.Option(
        Path("config/razanajatovo_2016_reproductive_traits.yml"), exists=True
    ),
    scope_config_path: Path = typer.Option(Path("config/angiosperm_scope.yml"), exists=True),
    resolve_gbif: bool = typer.Option(True, "--resolve-gbif/--no-resolve-gbif"),
    gbif_workers: int = typer.Option(8, min=1, max=32),
) -> None:
    report = prepare_razanajatovo_2016(
        workbook_path=workbook_path,
        sheet_csv=sheet_csv,
        master_csv=master_csv,
        output_dir=output_dir,
        config_path=config_path,
        scope_config_path=scope_config_path,
        resolve_gbif=resolve_gbif,
        gbif_workers=gbif_workers,
    )
    typer.echo(
        f"Prepared {report['n_direct_candidate_cells']:,} reproductive cells for "
        f"{report['n_direct_candidate_species']:,} species; "
        f"{report['n_unmatched_useful_source_species']:,} useful source species remain unmatched."
    )


if __name__ == "__main__":
    app()
