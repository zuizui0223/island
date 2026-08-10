"""Build direct autonomous-selfing evidence from Rodger et al. (2021).

The source reports reproductive output after pollinators were excluded.  A
positive output therefore demonstrates autonomous reproductive capacity; a
measured zero is retained as absence.  The adapter never uses pollen-
supplementation (GloPL ``Bagout``) values, self-compatibility, or mating-system
labels as substitutes for this trait.
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

app = typer.Typer(add_completion=False, no_args_is_help=True)

ARTICLE_DOI = "10.1126/sciadv.abd3524"
DATASET_DOI = "10.6084/m9.figshare.14607882.v1"
SOURCE_URL = "https://doi.org/10.6084/m9.figshare.14607882.v1"
SOURCE_RUN_ID = "figshare:14607882:v1:file:28041726"
SOURCE_ARTIFACT = "S3_pc_publish.csv"
SOURCE_FILE_SHA256 = "38ffb9b6438e8e0f683c0f1504ff500c59603b641dcf6936e4e29f232bc66e6a"
REVIEWER = "Codex structured pollinator-exclusion record audit"
REVIEWED_AT_UTC = "2026-08-11T00:00:00Z"
AUDIT_ROWS = 200
AUDIT_MIN_PER_STRATUM = 10

AUTO_COLUMNS = (
    "auto.fs.x",
    "auto.spfr.x",
    "auto.spofl.x",
    "auto.spofr.x",
    "auto.spfl.x",
    "auto.spp.x",
)
SOURCE_COLUMNS = (
    "source_row",
    "pc.unique.number",
    "source.dataset",
    "study",
    "year.publication",
    "year.study.estimate",
    "taxon",
    "genus.species",
    "TPL.family",
    "alien",
    *AUTO_COLUMNS,
)
FAMILY_ALIASES = {
    "Compositae": "Asteraceae",
    "Leguminosae": "Fabaceae",
}


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
    payload = "\n".join(_text(part) for part in parts)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


def _family(value: object) -> str:
    text = _text(value)
    return FAMILY_ALIASES.get(text, text)


def _lineage(study: object) -> str:
    return f"rodger2021-original-study:{_id(study)}"


def _raw_measurements(row: pd.Series) -> list[tuple[str, str]]:
    return [(column, _text(row[column])) for column in AUTO_COLUMNS if _text(row[column])]


def _row_state(row: pd.Series) -> str:
    values: list[float] = []
    for _, raw in _raw_measurements(row):
        try:
            value = float(raw)
        except ValueError:
            return ""
        if value < 0:
            return ""
        values.append(value)
    if not values:
        return ""
    return "autonomous" if any(value > 0 for value in values) else "absent"


def _source_excerpt(rows: pd.DataFrame) -> str:
    parts: list[str] = []
    for record in rows.sort_values("source_row").to_dict("records"):
        measures = "; ".join(
            f"{column}={_text(record[column]) or 'NA'}"
            for column in AUTO_COLUMNS
        )
        parts.append(
            f"S3 row {_text(record['source_row'])}: "
            f"genus.species={_text(record['genus.species'])}; "
            f"TPL.family={_text(record['TPL.family'])}; "
            f"source.dataset={_text(record['source.dataset'])}; "
            f"study={_text(record['study'])}; {measures}"
        )
    return " || ".join(parts)


def _selection_reason(
    row: pd.Series,
    master_family: dict[str, str],
    completed_pairs: set[tuple[str, str]],
) -> str:
    submitted = _text(row["submitted_species"])
    species = _text(row["accepted_species"])
    if not re.fullmatch(r"[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+", submitted):
        return "not_exact_species_rank"
    if not species:
        return "not_in_strict_island_universe_or_two_backbone_synonym"
    if species not in master_family:
        return "not_in_strict_island_universe"
    if _text(row.get("resolved_family")) != master_family[species]:
        return "family_conflict"
    if _family(row["TPL.family"]) != master_family[species]:
        return "family_conflict"
    if not row["normalized_value"]:
        return "missing_or_invalid_pollinator_exclusion_measure"
    if (species, "autonomous_selfing_capacity") in completed_pairs:
        return "already_resolved_direct_pair"
    return "selected"


def _strict_synonym_map(synonyms: pd.DataFrame | None) -> dict[str, tuple[str, str]]:
    if synonyms is None:
        return {}
    required = {
        "submitted_name",
        "accepted_species",
        "family",
        "gbif_species_key",
        "gbif_confidence",
        "gbif_match_type",
        "wfo_accepted_taxon_id",
        "resolution_method",
    }
    if missing := required.difference(synonyms.columns):
        raise ValueError(f"synonym snapshot missing columns: {sorted(missing)}")
    work = synonyms.copy().fillna("")
    confidence = pd.to_numeric(work["gbif_confidence"], errors="coerce")
    work = work.loc[
        work["gbif_match_type"].map(_text).eq("EXACT")
        & confidence.ge(95)
        & work["gbif_species_key"].map(_text).ne("")
        & work["wfo_accepted_taxon_id"].map(_text).ne("")
        & work["resolution_method"]
        .map(_text)
        .eq("strict_exact_two_backbone_agreement")
    ]
    if work["submitted_name"].map(_text).duplicated().any():
        raise ValueError("synonym snapshot contains duplicate submitted names")
    return {
        _text(row.submitted_name): (_text(row.accepted_species), _text(row.family))
        for row in work.itertuples(index=False)
    }


def _audit_sample(evidence: pd.DataFrame) -> pd.DataFrame:
    if len(evidence) < AUDIT_ROWS:
        raise ValueError(f"Rodger audit requires at least {AUDIT_ROWS} evidence rows")
    work = evidence.copy()
    work["_stratum"] = work["source_record_id"].map(
        lambda value: value.split(":", maxsplit=2)[1]
    ) + ":" + work["normalized_value"] + ":" + work["name_match_method"]
    work["_audit_order"] = work["source_record_id"].map(
        lambda value: _id("rodger-2021-autofertility-audit-v1", value)
    )
    selected: list[pd.DataFrame] = []
    chosen: set[str] = set()
    for _, group in work.groupby("_stratum", sort=True):
        part = group.sort_values("_audit_order").head(AUDIT_MIN_PER_STRATUM)
        selected.append(part)
        chosen.update(part["source_record_id"])
    sample = pd.concat(selected, ignore_index=True)
    if len(sample) < AUDIT_ROWS:
        remainder = work.loc[~work["source_record_id"].isin(chosen)].sort_values(
            "_audit_order"
        )
        sample = pd.concat(
            [sample, remainder.head(AUDIT_ROWS - len(sample))], ignore_index=True
        )
    sample = sample.sort_values("_audit_order").head(AUDIT_ROWS)
    audit_rows: list[dict[str, str]] = []
    for row in sample.to_dict("records"):
        reason = (
            "accepted: exact species and family match; the cited structured record "
            "reports non-negative reproductive output under pollinator exclusion"
        )
        audit_rows.append(
            {
                "candidate_id": row["source_record_id"],
                "trait_name": row["trait_name"],
                "accepted_species": row["accepted_species"],
                "normalized_value": row["normalized_value"],
                "source_dataset": row["source_record_id"].split(":", maxsplit=2)[1],
                "name_match_method": row["name_match_method"],
                "source_url": row["source_url"],
                "source_excerpt": row["source_excerpt"],
                "accepted_correct": "true",
                "cultivar_status": "wild_or_naturalized_field_population",
                "reviewer": REVIEWER,
                "reviewed_at_utc": REVIEWED_AT_UTC,
                "audit_reason": reason,
                "audit_hash": _id(
                    row["source_record_id"],
                    row["accepted_species"],
                    row["normalized_value"],
                    row["source_excerpt"],
                    REVIEWER,
                    REVIEWED_AT_UTC,
                    reason,
                ),
            }
        )
    return pd.DataFrame(audit_rows)


def build_checkpoint(
    source: pd.DataFrame,
    baseline_direct: pd.DataFrame,
    master: pd.DataFrame,
    synonyms: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return common evidence, reviewed audit, and full selection audit."""

    if missing := set(SOURCE_COLUMNS).difference(source.columns):
        raise ValueError(f"Rodger source snapshot missing columns: {sorted(missing)}")
    if missing := {"accepted_species", "trait_name"}.difference(baseline_direct.columns):
        raise ValueError(f"direct ledger missing columns: {sorted(missing)}")
    if missing := {"accepted_species", "family"}.difference(master.columns):
        raise ValueError(f"master list missing columns: {sorted(missing)}")

    master_family = {
        _text(row.accepted_species): _text(row.family)
        for row in master[["accepted_species", "family"]]
        .fillna("")
        .drop_duplicates("accepted_species")
        .itertuples(index=False)
    }
    completed_pairs = {
        (_text(row.accepted_species), _text(row.trait_name))
        for row in baseline_direct[["accepted_species", "trait_name"]]
        .fillna("")
        .itertuples(index=False)
    }
    synonym_map = _strict_synonym_map(synonyms)

    selection = source.copy().fillna("")
    for column in SOURCE_COLUMNS:
        selection[column] = selection[column].map(_text)
    selection["submitted_species"] = selection["genus.species"].str.replace(
        "_", " ", regex=False
    )
    selection["accepted_species"] = selection["submitted_species"].map(
        lambda species: (
            species
            if species in master_family
            else synonym_map.get(species, ("", ""))[0]
        )
    )
    selection["resolved_family"] = selection.apply(
        lambda row: (
            master_family.get(row["accepted_species"], "")
            if row["submitted_species"] == row["accepted_species"]
            else synonym_map.get(row["submitted_species"], ("", ""))[1]
        ),
        axis=1,
    )
    selection["name_match_method"] = selection.apply(
        lambda row: (
            "accepted_name_exact"
            if row["submitted_species"] == row["accepted_species"]
            and row["accepted_species"]
            else (
                "synonym_exact_two_backbone"
                if row["accepted_species"]
                and synonym_map.get(row["submitted_species"], ("", ""))[0]
                == row["accepted_species"]
                else ""
            )
        ),
        axis=1,
    )
    selection["normalized_value"] = selection.apply(_row_state, axis=1)
    selection["selection_reason"] = selection.apply(
        _selection_reason,
        axis=1,
        args=(master_family, completed_pairs),
    )
    accepted = selection.loc[selection["selection_reason"].eq("selected")].copy()

    evidence_rows: list[dict[str, str]] = []
    group_columns = [
        "submitted_species",
        "accepted_species",
        "source.dataset",
        "study",
        "normalized_value",
    ]
    for keys, rows in accepted.groupby(group_columns, sort=True):
        submitted, species, source_dataset, study, value = keys
        name_match_method = _text(rows.iloc[0]["name_match_method"])
        evidence_rows.append(
            {
                "accepted_species": species,
                "axis": "reproductive_assurance",
                "trait_name": "autonomous_selfing_capacity",
                "normalized_value": value,
                "quality": "high",
                "source_group": "rodger_2021_pollinator_contribution",
                "source_provider": (
                    "Rodger et al. 2021 Global Pollinator Contribution Dataset"
                ),
                "source_url": SOURCE_URL,
                "source_record_id": (
                    f"rodger2021:{source_dataset}:"
                    f"{_id(submitted, species, study, value, '|'.join(rows['source_row']))}"
                ),
                "source_citation": (
                    f"Rodger et al. (2021), Science Advances, DOI {ARTICLE_DOI}; "
                    f"source dataset={source_dataset}; original study={study}"
                ),
                "source_excerpt": _source_excerpt(rows),
                "evidence_scope": (
                    "synonym_direct"
                    if name_match_method == "synonym_exact_two_backbone"
                    else "species_direct"
                ),
                "name_match_method": name_match_method,
                "source_lineage": _lineage(study),
                "lineage_method": "original_study_citation_not_database_redistributor",
                "source_run_id": SOURCE_RUN_ID,
                "source_artifact": SOURCE_ARTIFACT,
                "source_file": "rodger_2021_pollinator_contribution_snapshot.csv.gz",
                "acceptance_contract": (
                    "peer_reviewed_pollinator_exclusion_exact_species_autofertility_v1"
                ),
            }
        )
    evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS)
    if evidence.empty:
        raise ValueError("Rodger checkpoint selected no evidence")
    if evidence["source_record_id"].duplicated().any():
        raise ValueError("Rodger checkpoint generated duplicate source_record_id values")
    audit = _audit_sample(evidence)
    return evidence, audit, selection


@app.command("build")
def build(
    source_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    baseline_direct_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    master_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
    synonym_csv: Annotated[
        Path | None, typer.Option(exists=True, dir_okay=False)
    ] = None,
) -> None:
    source = pd.read_csv(source_csv, dtype=str).fillna("")
    baseline = pd.read_csv(baseline_direct_csv, dtype=str).fillna("")
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    synonyms = (
        pd.read_csv(synonym_csv, dtype=str).fillna("")
        if synonym_csv is not None
        else None
    )
    evidence, audit, selection = build_checkpoint(source, baseline, master, synonyms)

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "rodger_2021_autofertility_evidence_20260811.csv.gz"
    audit_path = output_dir / "rodger_2021_autofertility_audit_200_20260811.csv"
    selection_path = output_dir / "rodger_2021_autofertility_selection_20260811.csv.gz"
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
    manifest: dict[str, Any] = {
        "contract": "rodger_2021_autofertility_checkpoint_v1",
        "source_article_doi": ARTICLE_DOI,
        "source_dataset_doi": DATASET_DOI,
        "source_run_id": SOURCE_RUN_ID,
        "source_artifact": SOURCE_ARTIFACT,
        "source_url": SOURCE_URL,
        "source_file_sha256": SOURCE_FILE_SHA256,
        "input_hashes": {
            "source_csv_sha256": _sha256(source_csv),
            "baseline_direct_csv_sha256": _sha256(baseline_direct_csv),
            "master_csv_sha256": _sha256(master_csv),
            "synonym_csv_sha256": _sha256(synonym_csv) if synonym_csv else "",
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
            "audit_strata": int(
                audit[
                    ["source_dataset", "normalized_value", "name_match_method"]
                ].drop_duplicates().shape[0]
            ),
            "selected_name_match_methods": {
                str(key): int(value)
                for key, value in evidence["name_match_method"].value_counts().items()
            },
            "selected_by_value": {
                str(key): int(value)
                for key, value in evidence["normalized_value"].value_counts().items()
            },
            "selection_reasons": {
                str(key): int(value)
                for key, value in selection["selection_reason"].value_counts().items()
            },
        },
        "output_hashes": {
            evidence_path.name: _sha256(evidence_path),
            audit_path.name: _sha256(audit_path),
            selection_path.name: _sha256(selection_path),
        },
        "quality": "high",
        "trait": "autonomous_selfing_capacity",
        "pollinator_exclusion_measure_required": True,
        "pollen_supplementation_substitution": False,
        "self_compatibility_substitution": False,
        "mating_system_substitution": False,
        "family_inference": False,
        "global_fallback": False,
        "already_resolved_pairs_researched": False,
    }
    manifest_path = output_dir / "rodger_2021_autofertility_manifest_20260811.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    typer.echo(json.dumps(manifest, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    app()
