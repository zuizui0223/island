"""Build species-direct floral-symmetry evidence from Zell et al.'s BSdb.

BSdb contains separate provenance columns for individually collected
species-level traits and for genus/family scores.  This adapter only accepts
rows with ``trait.source`` and never promotes ``genus_trait_source`` or
``family_trait_source``.  It also excludes the ``wind`` category because that
is a pollination mode, not a floral-symmetry state.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Annotated

import pandas as pd
import typer

from island_v2.bsdb_reproductive_checkpoint import (
    ARTICLE_DOI,
    SOURCE_ARTIFACT,
    SOURCE_COMMIT,
    SOURCE_RUN_ID,
    SOURCE_URL,
    _id,
    _lineage_key,
    _sha256,
    _text,
    apply_reviewed_audit,
)
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS

app = typer.Typer(add_completion=False)

AUDIT_SIZE = 200
VALUE_MAP = {
    "radial": "actinomorphic",
    "bilateral": "zygomorphic",
    "slightly zygomorphic": "zygomorphic",
    "asymmetric": "asymmetric",
}


def _trait_source_lineage_key(value: object) -> str:
    text = _text(value)
    delimiter = r"\s*\|\s*"
    if "http" not in text.casefold() and "," in text:
        delimiter = r"\s*(?:\||,)\s*"
    parts = sorted(
        {
            key
            for part in re.split(delimiter, text)
            if (key := _lineage_key(part))
        }
    )
    return "__".join(parts)


def _selection_reason(
    row: pd.Series,
    master_family: dict[str, str],
    current_pairs: set[tuple[str, str]],
) -> str:
    species = _text(row.get("tnrs_Sciname"))
    if _text(row.get("tnrs_infrasp")) or _text(row.get("Infrasp")):
        return "infraspecific_record"
    if not re.fullmatch(r"[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+", species):
        return "not_exact_species_rank"
    if species not in master_family:
        return "not_in_strict_island_universe"
    source_family = _text(row.get("tnrs_family"))
    if source_family and master_family[species] and source_family != master_family[species]:
        return "family_conflict"
    if not _text(row.get("trait.source")):
        return "missing_individually_collected_species_source"
    if _text(row.get("FlrSymCheck")).casefold() not in VALUE_MAP:
        return "unsupported_floral_symmetry_state"
    if (species, "floral_symmetry") in current_pairs:
        return "already_resolved_direct_pair"
    return "selected"


def _audit_sample(evidence: pd.DataFrame, audit_size: int) -> pd.DataFrame:
    if audit_size <= 0:
        raise ValueError("audit_size must be positive")
    per_value = max(1, audit_size // max(1, evidence["normalized_value"].nunique()))
    selected: list[pd.DataFrame] = []
    for _, group in evidence.groupby("normalized_value", sort=True):
        work = group.copy()
        work["_audit_order"] = work["source_record_id"].map(
            lambda item: _id("bsdb-symmetry-audit-v1", item)
        )
        selected.append(work.sort_values("_audit_order").head(per_value))
    sample = pd.concat(selected, ignore_index=False)
    if len(sample) < audit_size:
        remaining = evidence.loc[~evidence.index.isin(sample.index)].copy()
        remaining["_audit_order"] = remaining["source_record_id"].map(
            lambda item: _id("bsdb-symmetry-audit-v1", item)
        )
        sample = pd.concat(
            [sample, remaining.sort_values("_audit_order").head(audit_size - len(sample))]
        )
    sample = sample.sort_values(["normalized_value", "_audit_order"]).head(audit_size)
    if len(sample) != audit_size:
        raise ValueError(f"BSdb symmetry audit has fewer than {audit_size} rows")
    return sample


def build_checkpoint(
    source: pd.DataFrame,
    baseline_direct: pd.DataFrame,
    master: pd.DataFrame,
    *,
    audit_size: int = AUDIT_SIZE,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return common evidence, a reproducible audit sample, and exclusions."""

    required_source = {
        "Infrasp",
        "FlrSymCheck",
        "trait.source",
        "tnrs_status",
        "tnrs_family",
        "tnrs_Sciname",
        "tnrs_infrasp",
    }
    if missing := required_source.difference(source.columns):
        raise ValueError(f"BSdb source missing columns: {sorted(missing)}")
    if missing := {"accepted_species", "family"}.difference(master.columns):
        raise ValueError(f"master list missing columns: {sorted(missing)}")
    if missing := {"accepted_species", "trait_name"}.difference(
        baseline_direct.columns
    ):
        raise ValueError(f"direct ledger missing columns: {sorted(missing)}")

    master_family = {
        _text(row.accepted_species): _text(row.family)
        for row in master.itertuples(index=False)
    }
    current_pairs = {
        (_text(row.accepted_species), _text(row.trait_name))
        for row in baseline_direct.itertuples(index=False)
    }
    selection = source.copy().fillna("")
    selection["source_row"] = selection.index + 2
    selection["selection_reason"] = selection.apply(
        _selection_reason,
        axis=1,
        args=(master_family, current_pairs),
    )
    accepted = selection.loc[selection["selection_reason"].eq("selected")].copy()

    rows: list[dict[str, str]] = []
    for row in accepted.to_dict("records"):
        species = _text(row["tnrs_Sciname"])
        raw_value = _text(row["FlrSymCheck"])
        value = VALUE_MAP[raw_value.casefold()]
        source_key = _text(row["trait.source"])
        source_row = _text(row["source_row"])
        status = _text(row["tnrs_status"]).casefold()
        rows.append(
            {
                "accepted_species": species,
                "axis": "floral_structural_complexity",
                "trait_name": "floral_symmetry",
                "normalized_value": value,
                "quality": "high",
                "source_group": "zell_2025_bsdb",
                "source_provider": "Zell et al. 2025 Breeding System Database",
                "source_url": SOURCE_URL,
                "source_record_id": (
                    f"bsdb-symmetry:{_id(source_row, species, source_key, raw_value)}"
                ),
                "source_citation": (
                    f"Zell et al. (2025), New Phytologist, DOI {ARTICLE_DOI}; "
                    f"individually collected species-level trait source: {source_key}"
                ),
                "source_excerpt": (
                    f"tnrs_Sciname={species}; "
                    f"tnrs_status={_text(row['tnrs_status'])}; "
                    f"tnrs_family={_text(row['tnrs_family'])}; "
                    f"FlrSymCheck={raw_value}; trait.source={source_key}"
                ),
                "evidence_scope": (
                    "synonym_direct" if status == "synonym" else "species_direct"
                ),
                "name_match_method": (
                    "exact_synonym_to_accepted"
                    if status == "synonym"
                    else "accepted_name_exact"
                ),
                "source_lineage": (
                    f"bsdb-trait-original:{_trait_source_lineage_key(source_key)}"
                ),
                "lineage_method": (
                    "individually_collected_species_trait_source_not_family_or_genus"
                ),
                "source_run_id": SOURCE_RUN_ID,
                "source_artifact": SOURCE_ARTIFACT,
                "source_file": SOURCE_ARTIFACT,
                "acceptance_contract": (
                    "published_database_exact_species_floral_symmetry_trait_v1"
                ),
            }
        )
    evidence = pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)
    if evidence.empty:
        raise ValueError("BSdb morphology checkpoint selected no evidence")
    if evidence["source_record_id"].duplicated().any():
        raise ValueError("BSdb morphology checkpoint generated duplicate record IDs")

    sample = _audit_sample(evidence, audit_size)
    audit = pd.DataFrame(
        {
            "candidate_id": sample["source_record_id"],
            "trait_name": sample["trait_name"],
            "accepted_species": sample["accepted_species"],
            "normalized_value": sample["normalized_value"],
            "source_url": sample["source_url"],
            "source_excerpt": sample["source_excerpt"],
            "accepted_correct": "",
            "cultivar_status": "",
            "reviewer": "",
            "reviewed_at_utc": "",
            "audit_reason": "",
            "audit_hash": "",
        }
    ).reset_index(drop=True)
    return evidence, audit, selection


@app.command("build")
def build(
    source_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    baseline_direct_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    master_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
    reviewed_audit_csv: Annotated[
        Path | None, typer.Option(exists=True, dir_okay=False)
    ] = None,
) -> None:
    source = pd.read_csv(source_csv, dtype=str).fillna("")
    baseline = pd.read_csv(baseline_direct_csv, dtype=str).fillna("")
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    evidence, audit, selection = build_checkpoint(source, baseline, master)
    reviewed_input_hash = ""
    if reviewed_audit_csv is not None:
        reviewed = pd.read_csv(reviewed_audit_csv, dtype=str).fillna("")
        audit = apply_reviewed_audit(audit, reviewed)
        reviewed_input_hash = _sha256(reviewed_audit_csv)

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "bsdb_morphology_evidence_20260811.csv.gz"
    audit_path = output_dir / "bsdb_morphology_manual_audit_200_20260811.csv"
    selection_path = output_dir / "bsdb_morphology_selection_audit_20260811.csv.gz"
    evidence.to_csv(
        evidence_path,
        index=False,
        compression={"method": "gzip", "mtime": 0},
        lineterminator="\n",
    )
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    selection.to_csv(
        selection_path,
        index=False,
        compression={"method": "gzip", "mtime": 0},
        lineterminator="\n",
    )
    manifest = {
        "contract": "bsdb_morphology_checkpoint_v1",
        "source_article_doi": ARTICLE_DOI,
        "source_commit": SOURCE_COMMIT,
        "source_run_id": SOURCE_RUN_ID,
        "source_artifact": SOURCE_ARTIFACT,
        "source_url": SOURCE_URL,
        "input_hashes": {
            "source_csv_sha256": _sha256(source_csv),
            "baseline_direct_csv_sha256": _sha256(baseline_direct_csv),
            "master_csv_sha256": _sha256(master_csv),
            "reviewed_audit_csv_sha256": reviewed_input_hash,
        },
        "counts": {
            "source_rows": len(source),
            "selected_evidence_rows": len(evidence),
            "selected_species_trait": int(
                evidence[["accepted_species", "trait_name"]]
                .drop_duplicates()
                .shape[0]
            ),
            "selected_species": int(evidence["accepted_species"].nunique()),
            "source_lineages": int(evidence["source_lineage"].nunique()),
            "audit_rows": len(audit),
            "selection_reasons": {
                str(key): int(value)
                for key, value in selection["selection_reason"].value_counts().items()
            },
            "selected_by_value": {
                str(key): int(value)
                for key, value in evidence["normalized_value"].value_counts().items()
            },
        },
        "output_hashes": {
            evidence_path.name: _sha256(evidence_path),
            audit_path.name: _sha256(audit_path),
            selection_path.name: _sha256(selection_path),
        },
        "manual_audit_required_before_promotion": reviewed_audit_csv is None,
        "species_level_trait_source_required": True,
        "genus_trait_source_accepted": False,
        "family_trait_source_accepted": False,
        "wind_mapped_to_floral_symmetry": False,
        "family_inference": False,
        "global_fallback": False,
        "already_resolved_pairs_researched": False,
    }
    manifest_path = output_dir / "bsdb_morphology_checkpoint_manifest_20260811.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    typer.echo(json.dumps(manifest, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    app()
