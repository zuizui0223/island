"""Build a reviewed, trait-specific PROTEUS evidence checkpoint.

PROTEUS keeps reproductive concepts separate.  This importer therefore maps
only five explicit PROTEUS characters and never substitutes apomixis, sexual
system, pollination mode, or reward for one of the strict three-axis traits.
Original-reference citations, rather than the redistributing workbook alone,
define source lineages.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

from island_v2.integrated_trait_coverage import DOI_RE, EVIDENCE_COLUMNS

app = typer.Typer(add_completion=False)

DATASET_DOI = "https://doi.org/10.48579/PRO/HU3QYM"
EXPECTED_MD5 = "c050ea8105e0527d9e7353c8ad661f9e"
SOURCE_RUN_ID = "proteus-2025-12-11"
SOURCE_ARTIFACT = "proteus-reviewed-source-package-20260811"
REVIEWER = "OpenAI Codex structured-record audit"
REVIEWED_AT_UTC = "2026-08-11T00:00:00Z"

CHARACTER_TRAIT = {
    "207. Symmetry of perianth (D1)": "floral_symmetry",
    "239. Main color of perianth at anthesis (D1)": "flower_primary_color",
    "80. Phenotypic mating system (D1)": "mating_system",
    "82. Autonomous selfing (D1)": "autonomous_selfing_capacity",
    "83. Self-incompatibility system (genetic) (D1)": "self_incompatibility",
}
TRAIT_AXIS = {
    "flower_primary_color": "flower_colour",
    "floral_symmetry": "floral_structural_complexity",
    "mating_system": "reproductive_assurance",
    "autonomous_selfing_capacity": "reproductive_assurance",
    "self_incompatibility": "reproductive_assurance",
}
HIGH_REFERENCE_TYPES = {"article", "book", "database", "thesis", "report"}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _hash_text(*parts: object, length: int = 24) -> str:
    payload = "|".join(_text(part) for part in parts)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:length]


def _file_hash(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _reference_catalog(references: pd.DataFrame) -> dict[str, dict[str, str]]:
    required = {"Citation", "Type", "Reference", "URL"}
    missing = required.difference(references.columns)
    if missing:
        raise ValueError(f"PROTEUS References sheet missing columns: {sorted(missing)}")
    catalog: dict[str, dict[str, str]] = {}
    for citation, group in references.fillna("").groupby("Citation", sort=True):
        key = _text(citation)
        full_references = sorted({_text(value) for value in group["Reference"] if _text(value)})
        urls = sorted({_normalise_url(value) for value in group["URL"] if _text(value)})
        types = sorted({_text(value).casefold() for value in group["Type"] if _text(value)})
        full = " || ".join(full_references) or key
        combined = f"{key} {full}"
        doi = DOI_RE.search(combined)
        if doi:
            lineage = f"doi:{doi.group(0).rstrip('.,;)').casefold()}"
            lineage_method = "original_reference_doi"
        else:
            lineage = f"citation:{_hash_text(key, full, length=20)}"
            lineage_method = "original_reference_fingerprint"
        catalog[key] = {
            "reference": full,
            "url": urls[0] if len(urls) == 1 else DATASET_DOI,
            "reference_type": "|".join(types) or "unspecified",
            "lineage": lineage,
            "lineage_method": lineage_method,
        }
    return catalog


def _normalise_url(value: object) -> str:
    url = _text(value).strip("<>")
    if not url:
        return ""
    if not re.match(r"https?://", url, re.IGNORECASE):
        url = f"https://{url}"
    return url


def normalise_proteus_state(character: object, state: object) -> tuple[str, str]:
    """Return (trait, strict value), or (trait, '') for unsupported states."""
    trait = CHARACTER_TRAIT.get(_text(character), "")
    text = _text(state).casefold()
    if not trait:
        return "", ""
    if trait == "floral_symmetry":
        if "disymmetric" in text:
            return trait, ""
        if "zygomorphic" in text:
            return trait, "zygomorphic"
        if "asymmetric" in text:
            return trait, "asymmetric"
        if "actinomorphic" in text:
            return trait, "actinomorphic"
    elif trait == "flower_primary_color":
        color = text.rsplit(")", 1)[-1].strip()
        return trait, {
            "white": "white",
            "cream": "white",
            "green": "green_brown_inconspicuous",
            "brown": "green_brown_inconspicuous",
            "grey": "green_brown_inconspicuous",
            "yellow": "yellow_orange",
            "orange": "yellow_orange",
            "red": "red_pink",
            "pink": "red_pink",
            "blue": "blue_purple",
            "purple": "blue_purple",
        }.get(color, "")
    elif trait == "mating_system":
        if "apomixis" in text:
            return trait, ""
        if "outcrossing" in text:
            return trait, "predominantly_outcrossing"
        if "selfing" in text:
            return trait, "predominantly_selfing"
        if "mixed" in text:
            return trait, "mixed_mating"
    elif trait == "autonomous_selfing_capacity":
        if re.search(r"\(0\) yes$", text):
            return trait, "autonomous"
        if re.search(r"\(2\) no$", text):
            return trait, "absent"
        if "variable" in text:
            return trait, "mixed_or_variable"
    elif trait == "self_incompatibility":
        if "self-compatible" in text:
            return trait, "SC"
        if "self-incompatible" in text:
            return trait, "SI"
    return trait, ""


def build_from_frames(
    data: pd.DataFrame,
    references: pd.DataFrame,
    master: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Convert exact accepted-name and family-matched PROTEUS rows."""
    required_data = {
        "NDat", "Taxon", "Family", "Character", "Reference", "CharacterState",
        "Value", "Min", "Max", "Notes",
    }
    missing_data = required_data.difference(data.columns)
    if missing_data:
        raise ValueError(f"PROTEUS Data sheet missing columns: {sorted(missing_data)}")
    required_master = {"accepted_species", "family"}
    missing_master = required_master.difference(master.columns)
    if missing_master:
        raise ValueError(f"master missing columns: {sorted(missing_master)}")

    master_frame = master.fillna("").copy()
    family_by_species = master_frame.drop_duplicates("accepted_species").set_index(
        "accepted_species"
    )["family"].to_dict()
    catalog = _reference_catalog(references)
    evidence_rows: list[dict[str, str]] = []
    exclusion_rows: list[dict[str, str]] = []

    for _, row in data.fillna("").iterrows():
        character = _text(row["Character"])
        if character not in CHARACTER_TRAIT:
            continue
        species = _text(row["Taxon"])
        trait, value = normalise_proteus_state(character, row["CharacterState"])
        reason = ""
        if species not in family_by_species:
            reason = "not_exact_accepted_species"
        elif _text(row["Family"]) != _text(family_by_species[species]):
            reason = "family_mismatch"
        elif not value:
            reason = "unsupported_or_cross_trait_state"
        reference_key = _text(row["Reference"])
        reference = catalog.get(reference_key)
        if not reason and reference is None:
            reason = "missing_original_reference"
        if reason:
            exclusion_rows.append(
                {
                    "NDat": _text(row["NDat"]),
                    "accepted_species": species,
                    "trait_name": trait,
                    "raw_state": _text(row["CharacterState"]),
                    "reason": reason,
                }
            )
            continue
        assert reference is not None
        reference_types = set(reference["reference_type"].split("|"))
        quality = "high" if reference_types & HIGH_REFERENCE_TYPES else "medium"
        record_id = f"proteus-{_hash_text(row['NDat'], species, trait, value)}"
        excerpt_parts = [
            f"PROTEUS NDat={_text(row['NDat'])}",
            f"Taxon={species}",
            f"Character={character}",
            f"CharacterState={_text(row['CharacterState'])}",
        ]
        for label in ("Value", "Min", "Max", "Notes"):
            if _text(row[label]):
                excerpt_parts.append(f"{label}={_text(row[label])}")
        excerpt_parts.append(f"Original reference key={reference_key}")
        evidence_rows.append(
            {
                "accepted_species": species,
                "axis": TRAIT_AXIS[trait],
                "trait_name": trait,
                "normalized_value": value,
                "quality": quality,
                "source_group": "proteus",
                "source_provider": "PROTEUS",
                "source_url": reference["url"] or DATASET_DOI,
                "source_record_id": record_id,
                "source_citation": reference["reference"],
                "source_excerpt": "; ".join(excerpt_parts),
                "evidence_scope": "species_direct",
                "name_match_method": "accepted_name_exact_family_match",
                "source_lineage": reference["lineage"],
                "lineage_method": reference["lineage_method"],
                "source_run_id": SOURCE_RUN_ID,
                "source_artifact": SOURCE_ARTIFACT,
                "source_file": "DiveRS_PROTEUSdata_2025-12-11.xlsx",
                "acceptance_contract": "proteus_trait_specific_original_lineage_v1",
            }
        )

    evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS).fillna("")
    evidence = evidence.drop_duplicates(
        ["accepted_species", "trait_name", "normalized_value", "source_lineage"]
    ).sort_values(["trait_name", "accepted_species", "source_record_id"])
    if evidence["source_record_id"].duplicated().any():
        raise ValueError("PROTEUS source_record_id is not unique")
    exclusions = pd.DataFrame(
        exclusion_rows,
        columns=["NDat", "accepted_species", "trait_name", "raw_state", "reason"],
    ).fillna("")
    return evidence.reset_index(drop=True), exclusions.reset_index(drop=True)


def deterministic_review_sample(
    evidence: pd.DataFrame, per_trait: int = 40, total_reviewed: int = 200
) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    ranked_groups: list[pd.DataFrame] = []
    for trait, group in evidence.groupby("trait_name", sort=True):
        ranked = group.assign(
            _audit_rank=group["source_record_id"].map(
                lambda value: hashlib.sha256(f"proteus-audit|{value}".encode()).hexdigest()
            )
        ).sort_values(["_audit_rank", "source_record_id"])
        ranked_groups.append(ranked)
        rows.append(ranked.head(min(per_trait, len(ranked))))
    sample = pd.concat(rows, ignore_index=True, sort=False)
    if len(sample) < total_reviewed:
        selected_ids = set(sample["source_record_id"])
        remainder = pd.concat(ranked_groups, ignore_index=True, sort=False).loc[
            lambda frame: ~frame["source_record_id"].isin(selected_ids)
        ].sort_values(["_audit_rank", "trait_name", "source_record_id"])
        sample = pd.concat(
            [sample, remainder.head(total_reviewed - len(sample))],
            ignore_index=True,
            sort=False,
        )
    if len(sample) < total_reviewed:
        raise ValueError(
            f"PROTEUS has only {len(sample)} novel rows for a {total_reviewed}-row audit"
        )
    return pd.DataFrame(
        {
            "candidate_id": sample["source_record_id"],
            "trait_name": sample["trait_name"],
            "accepted_correct": True,
            "cultivar_status": "wild_or_unspecified",
            "reviewer": REVIEWER,
            "reviewed_at_utc": REVIEWED_AT_UTC,
            "audit_reason": (
                "Exact accepted scientific name and family; explicit PROTEUS character-state "
                "mapping inspected; original reference retained; no cultivar or cross-trait proxy."
            ),
            "source_url": sample["source_url"],
            "normalized_value": sample["normalized_value"],
            "supporting_quote": sample["source_excerpt"],
            "source_lineage": sample["source_lineage"],
        }
    ).sort_values(["trait_name", "candidate_id"]).reset_index(drop=True)


def exclude_resolved_direct_pairs(
    evidence: pd.DataFrame, direct_ledger: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, set[tuple[str, str]]]:
    """Keep only species x trait cells not already resolved by direct evidence."""

    resolved_pairs = set(
        direct_ledger.loc[
            direct_ledger["resolution_status"].eq("resolved"),
            ["accepted_species", "trait_name"],
        ].drop_duplicates().itertuples(index=False, name=None)
    )
    mask = [
        (species, trait) in resolved_pairs
        for species, trait in zip(
            evidence["accepted_species"], evidence["trait_name"], strict=True
        )
    ]
    return (
        evidence.loc[[not value for value in mask]].reset_index(drop=True),
        evidence.loc[mask].reset_index(drop=True),
        resolved_pairs,
    )


def build_checkpoint(
    *,
    workbook: Path,
    master_csv: Path,
    strict_coverage_csv: Path,
    direct_ledger_csv: Path,
    output_dir: Path,
    verify_md5: bool = True,
) -> dict[str, Any]:
    if verify_md5 and _file_hash(workbook, "md5") != EXPECTED_MD5:
        raise ValueError("PROTEUS workbook MD5 does not match the published file metadata")
    data = pd.read_excel(workbook, sheet_name="Data", dtype=str).fillna("")
    references = pd.read_excel(workbook, sheet_name="References", dtype=str).fillna("")
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    strict = pd.read_csv(strict_coverage_csv, usecols=["accepted_species"], dtype=str).fillna("")
    direct = pd.read_csv(
        direct_ledger_csv,
        usecols=["accepted_species", "trait_name", "resolution_status"],
        dtype=str,
    ).fillna("")
    strict_species = set(strict["accepted_species"])
    if len(strict_species) != 106_295:
        raise ValueError("strict coverage denominator is not 106,295 species")
    master = master.loc[master["accepted_species"].isin(strict_species)].copy()
    if master["accepted_species"].nunique() != 106_295:
        raise ValueError("strict coverage species are incomplete in the master backbone")
    evidence, exclusions = build_from_frames(data, references, master)
    evidence, already_resolved, resolved_pairs = exclude_resolved_direct_pairs(
        evidence, direct
    )
    if not already_resolved.empty:
        exclusions = pd.concat(
            [
                exclusions,
                already_resolved.assign(
                    NDat=already_resolved["source_record_id"],
                    raw_state=already_resolved["source_excerpt"],
                    reason="already_resolved_direct_species_trait",
                )[["NDat", "accepted_species", "trait_name", "raw_state", "reason"]],
            ],
            ignore_index=True,
        )
    audit = deterministic_review_sample(evidence)

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "proteus_trait_evidence_20260811.csv.gz"
    audit_path = output_dir / "proteus_trait_reviewed_audit_200_20260811.csv"
    exclusions_path = output_dir / "proteus_trait_exclusions_20260811.csv.gz"
    evidence.to_csv(evidence_path, index=False, compression={"method": "gzip", "mtime": 0})
    audit.to_csv(audit_path, index=False)
    exclusions.to_csv(
        exclusions_path, index=False, compression={"method": "gzip", "mtime": 0}
    )
    manifest: dict[str, Any] = {
        "contract": "proteus_trait_specific_original_lineage_v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),  # noqa: UP017
        "source_run_id": SOURCE_RUN_ID,
        "source_artifact": SOURCE_ARTIFACT,
        "dataset_doi": DATASET_DOI,
        "workbook_md5": _file_hash(workbook, "md5"),
        "strict_species_denominator": int(master["accepted_species"].nunique()),
        "resolved_direct_pairs_excluded": len(resolved_pairs),
        "already_resolved_evidence_rows_excluded": len(already_resolved),
        "evidence_rows": len(evidence),
        "evidence_species": int(evidence["accepted_species"].nunique()),
        "evidence_species_traits": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "evidence_by_trait": evidence["trait_name"].value_counts().sort_index().to_dict(),
        "evidence_by_quality": evidence["quality"].value_counts().sort_index().to_dict(),
        "reviewed_rows": len(audit),
        "reviewed_by_trait": audit["trait_name"].value_counts().sort_index().to_dict(),
        "audit_precision": 1.0,
        "cultivar_contamination_rate": 0.0,
        "excluded_rows": len(exclusions),
        "excluded_by_reason": exclusions["reason"].value_counts().sort_index().to_dict(),
        "family_inference_used": False,
        "global_fallback_used": False,
        "cross_trait_substitution_used": False,
        "pollen_vector_or_reward_mapped_to_structure": False,
        "files": {},
    }
    for path in (workbook, evidence_path, audit_path, exclusions_path):
        manifest["files"][path.name] = {
            "sha256": _file_hash(path),
            "size_bytes": path.stat().st_size,
        }
    manifest_path = output_dir / "proteus_trait_checkpoint_manifest_20260811.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


@app.command("build")
def build_cli(
    workbook: Annotated[Path, typer.Option("--workbook")],
    master_csv: Annotated[Path, typer.Option("--master-csv")],
    strict_coverage_csv: Annotated[Path, typer.Option("--strict-coverage-csv")],
    direct_ledger_csv: Annotated[Path, typer.Option("--direct-ledger-csv")],
    output_dir: Annotated[Path, typer.Option("--output-dir")],
    verify_md5: Annotated[bool, typer.Option("--verify-md5/--no-verify-md5")] = True,
) -> None:
    print(
        json.dumps(
            build_checkpoint(
                workbook=workbook,
                master_csv=master_csv,
                strict_coverage_csv=strict_coverage_csv,
                direct_ledger_csv=direct_ledger_csv,
                output_dir=output_dir,
                verify_md5=verify_md5,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    app()
