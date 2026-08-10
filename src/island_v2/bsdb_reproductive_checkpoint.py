"""Build a strict reproductive checkpoint from Zell et al.'s BSdb release.

The source is a published species-level breeding-system database.  This
adapter keeps only SC/SI records that can be reconciled exactly between the
source's WCVP-backed TNRS name and the island master list, rejects
infraspecific and family-conflicting rows, and skips every species x trait
pair already present in the prior formal direct ledger.  Original study keys
remain the source lineage so the BSdb redistribution is never counted as an
independent study.
"""

from __future__ import annotations

import hashlib
import json
import re
import unicodedata
from pathlib import Path
from typing import Annotated

import pandas as pd
import typer

from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS

app = typer.Typer(add_completion=False)

SOURCE_COMMIT = "9e87946d1e3121d39e657b702cf9f92ccc10936e"
SOURCE_RUN_ID = f"dirtyplants-BSdb@{SOURCE_COMMIT}"
SOURCE_ARTIFACT = "Zell_df_12_29_23.csv"
SOURCE_URL = (
    "https://raw.githubusercontent.com/dirtyplants/BSdb/"
    f"{SOURCE_COMMIT}/{SOURCE_ARTIFACT}"
)
ARTICLE_DOI = "10.1111/nph.20234"
AUDIT_PER_VALUE = 100
AUDIT_REVIEWED_AT = "2026-08-11T00:00:00Z"


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _id(*parts: object) -> str:
    payload = "|".join(_text(part) for part in parts)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


def _lineage_key(value: object) -> str:
    text = unicodedata.normalize("NFKD", _text(value))
    text = "".join(char for char in text if not unicodedata.combining(char))
    text = text.encode("ascii", errors="ignore").decode("ascii").casefold()
    return re.sub(r"[^a-z0-9]+", "_", text).strip("_")


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
    if _text(row.get("BreedingSys")) not in {"SC", "SI"}:
        return "unsupported_or_missing_breeding_system"
    if not _text(row.get("bs.Source")):
        return "missing_original_source_lineage"
    source_family = _text(row.get("tnrs_family"))
    if source_family and master_family[species] and source_family != master_family[species]:
        return "family_conflict"
    if (species, "self_incompatibility") in current_pairs:
        return "already_resolved_direct_pair"
    return "selected"


def build_checkpoint(
    source: pd.DataFrame,
    baseline_direct: pd.DataFrame,
    master: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return common evidence, a reproducible audit sample, and exclusions."""

    required_source = {
        "Infrasp",
        "BreedingSys",
        "ISI_value",
        "ISItype",
        "bs.Source",
        "tnrs_status",
        "tnrs_family",
        "tnrs_Sciname",
        "tnrs_infrasp",
    }
    missing = required_source.difference(source.columns)
    if missing:
        raise ValueError(f"BSdb source missing columns: {sorted(missing)}")
    required_master = {"accepted_species", "family"}
    if missing := required_master.difference(master.columns):
        raise ValueError(f"master list missing columns: {sorted(missing)}")
    required_direct = {"accepted_species", "trait_name"}
    if missing := required_direct.difference(baseline_direct.columns):
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

    evidence_rows: list[dict[str, str]] = []
    for row in accepted.to_dict("records"):
        species = _text(row["tnrs_Sciname"])
        value = _text(row["BreedingSys"])
        source_key = _text(row["bs.Source"])
        source_row = _text(row["source_row"])
        status = _text(row["tnrs_status"]).casefold()
        lineage_key = _lineage_key(source_key)
        evidence_rows.append(
            {
                "accepted_species": species,
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "normalized_value": value,
                "quality": "high",
                "source_group": "zell_2025_bsdb",
                "source_provider": "Zell et al. 2025 Breeding System Database",
                "source_url": SOURCE_URL,
                "source_record_id": f"bsdb:{_id(source_row, species, source_key, value)}",
                "source_citation": (
                    f"Zell et al. (2025), New Phytologist, DOI {ARTICLE_DOI}; "
                    f"underlying source key: {source_key}"
                ),
                "source_excerpt": (
                    f"tnrs_Sciname={species}; "
                    f"tnrs_status={_text(row['tnrs_status'])}; "
                    f"tnrs_family={_text(row['tnrs_family'])}; "
                    f"BreedingSys={value}; "
                    f"ISI_value={_text(row['ISI_value']) or 'NA'}; "
                    f"ISItype={_text(row['ISItype']) or 'NA'}; "
                    f"bs.Source={source_key}"
                ),
                "evidence_scope": (
                    "synonym_direct" if status == "synonym" else "species_direct"
                ),
                "name_match_method": (
                    "exact_synonym_to_accepted"
                    if status == "synonym"
                    else "accepted_name_exact"
                ),
                "source_lineage": f"bsdb-original:{lineage_key}",
                "lineage_method": "underlying_study_key_not_bsdb_redistributor",
                "source_run_id": SOURCE_RUN_ID,
                "source_artifact": SOURCE_ARTIFACT,
                "source_file": SOURCE_ARTIFACT,
                "acceptance_contract": (
                    "published_database_exact_species_reproductive_trait_v1"
                ),
            }
        )
    evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS)
    if evidence.empty:
        raise ValueError("BSdb checkpoint selected no evidence")
    if evidence["source_record_id"].duplicated().any():
        raise ValueError("BSdb checkpoint generated duplicate source_record_id values")

    sample_parts: list[pd.DataFrame] = []
    for value in ("SC", "SI"):
        group = evidence.loc[evidence["normalized_value"].eq(value)].copy()
        group["_audit_order"] = group["source_record_id"].map(
            lambda item: _id("bsdb-audit-v1", item)
        )
        group = group.sort_values("_audit_order").head(AUDIT_PER_VALUE)
        if len(group) != AUDIT_PER_VALUE:
            raise ValueError(f"BSdb audit has fewer than {AUDIT_PER_VALUE} {value} rows")
        sample_parts.append(group)
    sample = pd.concat(sample_parts, ignore_index=True)
    audit_rows: list[dict[str, str]] = []
    for row in sample.to_dict("records"):
        audit_rows.append(
            {
                "candidate_id": row["source_record_id"],
                "trait_name": row["trait_name"],
                "accepted_species": row["accepted_species"],
                "normalized_value": row["normalized_value"],
                "source_url": row["source_url"],
                "source_excerpt": row["source_excerpt"],
                "accepted_correct": "",
                "cultivar_status": "",
                "reviewer": "",
                "reviewed_at_utc": "",
                "audit_reason": "",
                "audit_hash": "",
            }
        )
    return evidence, pd.DataFrame(audit_rows), selection


def apply_reviewed_audit(template: pd.DataFrame, reviewed: pd.DataFrame) -> pd.DataFrame:
    """Attach review decisions without permitting evidence fields to change."""

    immutable = [
        "candidate_id",
        "trait_name",
        "accepted_species",
        "normalized_value",
        "source_url",
        "source_excerpt",
    ]
    decisions = [
        "accepted_correct",
        "cultivar_status",
        "reviewer",
        "reviewed_at_utc",
        "audit_reason",
        "audit_hash",
    ]
    expected = set(immutable + decisions)
    if missing := expected.difference(reviewed.columns):
        raise ValueError(f"reviewed audit missing columns: {sorted(missing)}")
    if reviewed["candidate_id"].duplicated().any():
        raise ValueError("reviewed audit contains duplicate candidate_id values")
    aligned = template[["candidate_id"]].merge(
        reviewed[list(expected)],
        on="candidate_id",
        how="left",
        validate="one_to_one",
    )
    for column in immutable[1:]:
        if not template[column].map(_text).eq(aligned[column].map(_text)).all():
            raise ValueError(f"reviewed audit changed immutable column: {column}")
    if aligned[decisions].apply(lambda column: column.map(_text).eq("").any()).any():
        raise ValueError("reviewed audit has incomplete decisions or provenance")
    result = template.copy()
    result[decisions] = aligned[decisions].to_numpy()
    return result


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
    evidence_path = output_dir / "bsdb_reproductive_evidence_20260811.csv.gz"
    audit_path = output_dir / "bsdb_reproductive_manual_audit_200_20260811.csv"
    selection_path = output_dir / "bsdb_reproductive_selection_audit_20260811.csv.gz"
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
        "contract": "bsdb_reproductive_checkpoint_v1",
        "source_article_doi": ARTICLE_DOI,
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
                evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
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
        "family_inference": False,
        "global_fallback": False,
        "already_resolved_pairs_researched": False,
        "trait_substitution": False,
    }
    manifest_path = output_dir / "bsdb_reproductive_checkpoint_manifest_20260811.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    typer.echo(json.dumps(manifest, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    app()
