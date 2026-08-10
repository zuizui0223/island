"""Build a conservative species-row checkpoint from Wang et al. (2023).

The published Data S2 table labels 279,877 angiosperm species as
actinomorphic or zygomorphic, but it does not map each row back to one of the
many upstream floras used by the compilation.  This adapter therefore treats
the table as Medium evidence, never as a second independent lineage when an
existing direct record already exists, and scales only families that pass an
independent species-direct validation gate.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS
from island_v2.open_web_common import reviewed_source_package_evidence

app = typer.Typer(add_completion=False, no_args_is_help=True)

ARTICLE_DOI = "10.1126/sciadv.adg2555"
DATASET_DOI = "10.5061/dryad.ghx3ffbrh"
SOURCE_URL = "https://doi.org/10.5061/dryad.ghx3ffbrh"
SOURCE_RUN_ID = "europepmc:PMC10599613:supplementaryFiles:20260811"
SOURCE_ARTIFACT = "sciadv.adg2555_data_s1_and_s2.zip/adg2555_Data_S2.rdata"
SOURCE_RDATA_SHA256 = "c381206d3f99825fc766e66425a123cea055682432e448358727ecf831d5c45b"
REVIEWER = "Codex independent direct-ledger cross-validation"
REVIEWED_AT_UTC = "2026-08-11T00:00:00Z"
DATASET_LINEAGE = "dataset:dryad.ghx3ffbrh:v7:upstream-row-lineage-unmapped"

MIN_FAMILY_AUDIT = 10
MIN_FAMILY_PRECISION = 0.95
CONTROVERSIAL_FAMILIES = {
    "Apiaceae",
    "Asteraceae",
    "Brassicaceae",
    "Cyperaceae",
    "Poaceae",
}
VALUE_MAP = {
    "BiSymmetry": "zygomorphic",
    "RadSymmetry": "actinomorphic",
}
SOURCE_COLUMNS = [
    "source_row",
    "order",
    "family",
    "genus",
    "accepted_species",
    "raw_symmetry",
    "note1",
    "note2",
]


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _id(*parts: object) -> str:
    value = "\n".join(_text(part) for part in parts)
    return hashlib.sha256(value.encode("utf-8")).hexdigest()[:24]


def _first_lineage(value: object) -> str:
    parts = sorted({_text(part) for part in _text(value).split("|") if _text(part)})
    return parts[0] if parts else DATASET_LINEAGE


def _source_excerpt(row: pd.Series) -> str:
    return (
        f"Data S2 row {_text(row['source_row'])}: "
        f"Species_E2={_text(row['accepted_species'])}; "
        f"Family_E={_text(row['family'])}; Genus_E={_text(row['genus'])}; "
        f"Symmtry={_text(row['raw_symmetry'])}; "
        f"Note1={_text(row['note1'])}; Note2={_text(row['note2'])}"
    )


def build_checkpoint(
    source: pd.DataFrame,
    baseline_direct: pd.DataFrame,
    master: pd.DataFrame,
    *,
    min_family_audit: int = MIN_FAMILY_AUDIT,
    min_family_precision: float = MIN_FAMILY_PRECISION,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return evidence, independent audit, selection audit, and family gate."""

    if missing := set(SOURCE_COLUMNS).difference(source.columns):
        raise ValueError(f"Wang source snapshot missing columns: {sorted(missing)}")
    required_direct = {
        "accepted_species",
        "trait_name",
        "normalized_value",
        "source_lineages",
    }
    if missing := required_direct.difference(baseline_direct.columns):
        raise ValueError(f"direct ledger missing columns: {sorted(missing)}")
    if missing := {"accepted_species", "family"}.difference(master.columns):
        raise ValueError(f"master list missing columns: {sorted(missing)}")
    if min_family_audit <= 0:
        raise ValueError("min_family_audit must be positive")
    if not 0 <= min_family_precision <= 1:
        raise ValueError("min_family_precision must be between zero and one")

    work = source.copy().fillna("")
    for column in SOURCE_COLUMNS:
        work[column] = work[column].map(_text)
    work["normalized_value"] = work["raw_symmetry"].map(VALUE_MAP).fillna("")
    duplicate_species = set(
        work.loc[work["accepted_species"].duplicated(keep=False), "accepted_species"]
    )

    master_work = master[["accepted_species", "family"]].copy().fillna("")
    master_work = master_work.map(_text).drop_duplicates("accepted_species")
    master_family = dict(zip(master_work["accepted_species"], master_work["family"], strict=True))

    all_direct = baseline_direct.loc[
        baseline_direct["trait_name"].map(_text).eq("floral_symmetry")
    ].copy()
    all_direct = all_direct.drop_duplicates("accepted_species")
    all_direct_values = dict(
        zip(all_direct["accepted_species"], all_direct["normalized_value"], strict=True)
    )
    direct_lineages = dict(
        zip(all_direct["accepted_species"], all_direct["source_lineages"], strict=True)
    )
    direct = all_direct.loc[
        all_direct["normalized_value"].map(_text).isin({"actinomorphic", "zygomorphic"})
    ].copy()
    direct_values = dict(zip(direct["accepted_species"], direct["normalized_value"], strict=True))

    work["master_family"] = work["accepted_species"].map(master_family).fillna("")
    work["existing_direct_value"] = work["accepted_species"].map(all_direct_values).fillna("")
    work["direct_value"] = work["accepted_species"].map(direct_values).fillna("")
    work["direct_agreement"] = work["direct_value"].ne("") & work["normalized_value"].eq(
        work["direct_value"]
    )
    comparable = work.loc[
        work["master_family"].ne("")
        & work["family"].eq(work["master_family"])
        & work["direct_value"].ne("")
        & ~work["accepted_species"].isin(duplicate_species)
        & work["normalized_value"].ne("")
    ].copy()

    family_rows: list[dict[str, Any]] = []
    for family in sorted(set(work["family"])):
        group = comparable.loc[comparable["family"].eq(family)]
        reviewed = len(group)
        correct = int(group["direct_agreement"].sum())
        precision = correct / reviewed if reviewed else None
        controversial = family in CONTROVERSIAL_FAMILIES
        approved = bool(
            not controversial
            and reviewed >= min_family_audit
            and precision is not None
            and precision >= min_family_precision
        )
        family_rows.append(
            {
                "family": family,
                "reviewed_species": reviewed,
                "accepted_correct": correct,
                "precision": precision,
                "controversial_definition": controversial,
                "production_approved": approved,
            }
        )
    family_audit = pd.DataFrame(family_rows)
    approved_families = set(family_audit.loc[family_audit["production_approved"], "family"])

    def reason(row: pd.Series) -> str:
        species = row["accepted_species"]
        if not re.fullmatch(r"[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+", species):
            return "not_exact_species_rank"
        if species in duplicate_species:
            return "duplicate_source_species"
        if not row["normalized_value"]:
            return "unsupported_symmetry_state"
        if not row["master_family"]:
            return "not_in_strict_island_universe"
        if row["family"] != row["master_family"]:
            return "family_conflict"
        if row["family"] in CONTROVERSIAL_FAMILIES:
            return "controversial_family_definition"
        if row["family"] not in approved_families:
            family_row = family_audit.loc[family_audit["family"].eq(row["family"])]
            if family_row.empty or int(family_row.iloc[0]["reviewed_species"]) < min_family_audit:
                return "insufficient_independent_family_audit"
            return "family_precision_below_threshold"
        return "selected_family_validated_species_row"

    work["selection_reason"] = work.apply(reason, axis=1)
    selected = work.loc[work["selection_reason"].eq("selected_family_validated_species_row")].copy()

    evidence_rows: list[dict[str, str]] = []
    for _, row in selected.iterrows():
        species = row["accepted_species"]
        overlapping = species in direct_lineages
        lineage = (
            _first_lineage(direct_lineages.get(species, "")) if overlapping else DATASET_LINEAGE
        )
        evidence_rows.append(
            {
                "accepted_species": species,
                "axis": "floral_structural_complexity",
                "trait_name": "floral_symmetry",
                "normalized_value": row["normalized_value"],
                "quality": "medium",
                "source_group": "wang_2023_floral_symmetry",
                "source_provider": "Wang et al. 2023 species-level floral symmetry dataset",
                "source_url": SOURCE_URL,
                "source_record_id": (
                    "wang2023-symmetry:" + _id(row["source_row"], species, row["raw_symmetry"])
                ),
                "source_citation": (
                    f"Wang et al. (2023), Science Advances, DOI {ARTICLE_DOI}; "
                    "supplementary Data S2"
                ),
                "source_excerpt": _source_excerpt(row),
                "evidence_scope": "species_direct",
                "name_match_method": "accepted_name_and_family_exact",
                "source_lineage": lineage,
                "lineage_method": (
                    "conservative_reuse_existing_direct_lineage_due_unknown_upstream_origin"
                    if overlapping
                    else "compiled_dataset_single_lineage_upstream_row_mapping_unavailable"
                ),
                "source_run_id": SOURCE_RUN_ID,
                "source_artifact": SOURCE_ARTIFACT,
                "source_file": "wang_2023_species_symmetry_snapshot.csv.gz",
                "acceptance_contract": ("peer_reviewed_species_row_family_validated_medium_v1"),
            }
        )
    evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS)
    if evidence.empty:
        raise ValueError("Wang symmetry checkpoint selected no evidence")
    if evidence["source_record_id"].duplicated().any():
        raise ValueError("Wang symmetry checkpoint generated duplicate record IDs")

    selected_by_species = evidence.set_index("accepted_species")
    audit_source = comparable.loc[comparable["family"].isin(approved_families)].copy()
    nonbinary_direct_conflicts = selected.loc[
        selected["existing_direct_value"].ne("")
        & selected["direct_value"].eq("")
        & selected["existing_direct_value"].ne(selected["normalized_value"])
    ].copy()
    audit_source = pd.concat(
        [audit_source, nonbinary_direct_conflicts], ignore_index=True
    ).drop_duplicates("accepted_species")
    audit_rows: list[dict[str, str]] = []
    for _, row in audit_source.sort_values(["family", "accepted_species"]).iterrows():
        evidence_row = selected_by_species.loc[row["accepted_species"]]
        accepted = bool(row["normalized_value"] == row["existing_direct_value"])
        excerpt = _source_excerpt(row)
        audit_rows.append(
            {
                "candidate_id": evidence_row["source_record_id"],
                "trait_name": "floral_symmetry",
                "accepted_species": row["accepted_species"],
                "normalized_value": row["normalized_value"],
                "source_url": SOURCE_URL,
                "source_excerpt": excerpt,
                "accepted_correct": str(accepted).lower(),
                "cultivar_status": "wild_or_unspecified_species_record",
                "reviewer": REVIEWER,
                "reviewed_at_utc": REVIEWED_AT_UTC,
                "audit_reason": (
                    "agrees_with_independent_resolved_species_direct_state"
                    if accepted
                    else "conflicts_with_independent_resolved_species_direct_state"
                ),
                "audit_hash": _id(
                    row["accepted_species"],
                    row["normalized_value"],
                    row["existing_direct_value"],
                    accepted,
                ),
                "provider": "Wang et al. 2023 Data S2",
                "source_lineage": evidence_row["source_lineage"],
                "source_citation": evidence_row["source_citation"],
                "exact_supporting_quote": excerpt,
                "name_match_method": "accepted_name_and_family_exact",
            }
        )
    audit = pd.DataFrame(audit_rows)
    return evidence, audit, work, family_audit


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    if path.name.endswith(".gz"):
        frame.to_csv(
            path,
            index=False,
            lineterminator="\n",
            compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
        )
    else:
        frame.to_csv(path, index=False, lineterminator="\n")


def write_checkpoint(
    source_csv: Path,
    baseline_direct_csv: Path,
    master_csv: Path,
    output_dir: Path,
) -> dict[str, Any]:
    source = pd.read_csv(source_csv, dtype=str).fillna("")
    baseline = pd.read_csv(baseline_direct_csv, dtype=str).fillna("")
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    evidence, audit, selection, family_audit = build_checkpoint(source, baseline, master)
    output_dir.mkdir(parents=True, exist_ok=True)
    files = {
        "evidence": output_dir / "wang_2023_symmetry_evidence_20260811.csv.gz",
        "audit": output_dir / "wang_2023_symmetry_independent_audit_20260811.csv.gz",
        "selection": output_dir / "wang_2023_symmetry_selection_audit_20260811.csv.gz",
        "family": output_dir / "wang_2023_symmetry_family_gate_20260811.csv",
    }
    _write_csv(evidence, files["evidence"])
    _write_csv(audit, files["audit"])
    _write_csv(selection, files["selection"])
    _write_csv(family_audit, files["family"])
    reviewed = len(audit)
    correct = int(audit["accepted_correct"].map(lambda value: value == "true").sum())
    manifest = {
        "contract": "wang_2023_family_validated_species_symmetry_checkpoint_v1",
        "source_article_doi": ARTICLE_DOI,
        "source_dataset_doi": DATASET_DOI,
        "source_rdata_sha256": SOURCE_RDATA_SHA256,
        "source_run_id": SOURCE_RUN_ID,
        "source_artifact": SOURCE_ARTIFACT,
        "input_hashes": {
            "source_snapshot": _sha256(source_csv),
            "baseline_direct": _sha256(baseline_direct_csv),
            "master": _sha256(master_csv),
        },
        "gate": {
            "minimum_family_audit": MIN_FAMILY_AUDIT,
            "minimum_family_precision": MIN_FAMILY_PRECISION,
            "controversial_families_excluded": sorted(CONTROVERSIAL_FAMILIES),
            "approved_families": int(family_audit["production_approved"].sum()),
            "reviewed_species": reviewed,
            "accepted_correct": correct,
            "precision": correct / reviewed,
        },
        "counts": {
            "source_rows": len(source),
            "selected_evidence_rows": len(evidence),
            "selected_species": int(evidence["accepted_species"].nunique()),
            "selected_by_value": {
                str(key): int(value)
                for key, value in evidence["normalized_value"].value_counts().items()
            },
            "selection_reasons": {
                str(key): int(value)
                for key, value in selection["selection_reason"].value_counts().items()
            },
        },
        "guardrails": {
            "quality": "medium",
            "family_inference": False,
            "global_fallback": False,
            "axis_only_join": False,
            "trait_specific_join": "genus_x_trait_name",
            "upstream_row_lineage_available": False,
            "overlap_counted_as_independent_lineage": False,
        },
        "output_hashes": {path.name: _sha256(path) for path in files.values()},
    }
    manifest_path = output_dir / "wang_2023_symmetry_checkpoint_manifest_20260811.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return manifest


def combine_source_packages(
    existing_evidence_csv: Path,
    existing_audit_csv: Path,
    wang_output_dir: Path,
    combined_output_dir: Path,
) -> dict[str, Any]:
    existing_evidence = pd.read_csv(existing_evidence_csv, dtype=str).fillna("")
    existing_audit = pd.read_csv(existing_audit_csv, dtype=str).fillna("")
    wang_evidence = pd.read_csv(
        wang_output_dir / "wang_2023_symmetry_evidence_20260811.csv.gz",
        dtype=str,
    ).fillna("")
    wang_audit = pd.read_csv(
        wang_output_dir / "wang_2023_symmetry_independent_audit_20260811.csv.gz",
        dtype=str,
    ).fillna("")
    evidence = pd.concat([existing_evidence, wang_evidence], ignore_index=True)
    audit = pd.concat([existing_audit, wang_audit], ignore_index=True).fillna("")
    if evidence["source_record_id"].duplicated().any():
        raise ValueError("combined source package has duplicate source_record_id")
    if audit["candidate_id"].duplicated().any():
        raise ValueError("combined source package has duplicate audit candidate_id")
    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)
    combined_output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = combined_output_dir / "bsdb_csiro_wang_symmetry_evidence_20260811.csv.gz"
    audit_path = combined_output_dir / "bsdb_csiro_wang_symmetry_independent_audit_20260811.csv.gz"
    trait_path = combined_output_dir / "source_package_trait_audit_20260811.csv"
    _write_csv(evidence, evidence_path)
    _write_csv(audit, audit_path)
    _write_csv(scopes, trait_path)
    manifest = {
        "contract": "bsdb_csiro_wang_combined_reviewed_source_package_v1",
        "candidate_evidence_rows": len(evidence),
        "reviewed_rows": len(audit),
        "selected_evidence_rows": len(selected),
        "source_package_gate": summary,
        "trait_audit": scopes.to_dict("records"),
        "component_hashes": {
            existing_evidence_csv.name: _sha256(existing_evidence_csv),
            existing_audit_csv.name: _sha256(existing_audit_csv),
            "wang_evidence": _sha256(
                wang_output_dir / "wang_2023_symmetry_evidence_20260811.csv.gz"
            ),
            "wang_audit": _sha256(
                wang_output_dir / "wang_2023_symmetry_independent_audit_20260811.csv.gz"
            ),
        },
        "output_hashes": {
            evidence_path.name: _sha256(evidence_path),
            audit_path.name: _sha256(audit_path),
            trait_path.name: _sha256(trait_path),
        },
    }
    manifest_path = (
        combined_output_dir / "bsdb_csiro_wang_symmetry_source_package_manifest_20260811.json"
    )
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return manifest


@app.command("build")
def build(
    source_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    baseline_direct_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    master_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
    existing_evidence_csv: Annotated[Path | None, typer.Option(exists=True, dir_okay=False)] = None,
    existing_audit_csv: Annotated[Path | None, typer.Option(exists=True, dir_okay=False)] = None,
    combined_output_dir: Annotated[Path | None, typer.Option(file_okay=False)] = None,
) -> None:
    manifest = write_checkpoint(
        source_csv,
        baseline_direct_csv,
        master_csv,
        output_dir,
    )
    if any(
        value is not None
        for value in (existing_evidence_csv, existing_audit_csv, combined_output_dir)
    ):
        if not all(
            value is not None
            for value in (existing_evidence_csv, existing_audit_csv, combined_output_dir)
        ):
            raise typer.BadParameter(
                "existing evidence, existing audit, and combined output must be supplied together"
            )
        manifest["combined"] = combine_source_packages(
            existing_evidence_csv,
            existing_audit_csv,
            output_dir,
            combined_output_dir,
        )
    typer.echo(json.dumps(manifest, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    app()
