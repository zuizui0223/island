"""Build strict external-congener reproductive evidence for genus rules.

The checkpoint uses species that are outside the fixed 106,295-species island
denominator.  Their source records remain species-direct evidence, but the
records are exported with an explicit external-congener scope so they can only
train trait-specific genus rules.  WFO and GBIF must agree exactly on accepted
species and family; fuzzy, ambiguous, infraspecific, and target-universe rows
fail closed.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import unicodedata
from collections.abc import Iterable
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

from island_v2.angiosperm_scope import classify_scope
from island_v2.angiosperm_scope import load_config as load_scope_config
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS

app = typer.Typer(add_completion=False, no_args_is_help=True)

BSDB_COMMIT = "9e87946d1e3121d39e657b702cf9f92ccc10936e"
BSDB_ARTIFACT = "Zell_df_12_29_23.csv"
BSDB_SHA256 = "f935c8150b3d719aba5f62f14335c1a185304155403bf50db1e3ef1393fc55f4"
BSDB_URL = (
    "https://raw.githubusercontent.com/dirtyplants/BSdb/"
    f"{BSDB_COMMIT}/{BSDB_ARTIFACT}"
)
BSDB_ARTICLE_DOI = "10.1111/nph.20234"

RODGER_ARTICLE_DOI = "10.1126/sciadv.abd3524"
RODGER_DATASET_DOI = "10.6084/m9.figshare.14607882.v1"
RODGER_URL = f"https://doi.org/{RODGER_DATASET_DOI}"
RODGER_RUN_ID = "figshare:14607882:v1:file:28041726"
RODGER_ARTIFACT = "S3_pc_publish.csv"
RODGER_SHA256 = "fa745c578f3537933fafedc1d36b4ea266348cd83d7f6cbb231c253b0f348d3f"

CONTRACT = "external_congener_species_direct_strict_two_backbone_v1"
AUTO_COLUMNS = (
    "auto.fs.x",
    "auto.spfr.x",
    "auto.spofl.x",
    "auto.spofr.x",
    "auto.spfl.x",
    "auto.spp.x",
)
FAMILY_ALIASES = {"Compositae": "Asteraceae", "Leguminosae": "Fabaceae"}
BINOMIAL = re.compile(r"[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+")


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _family(value: object) -> str:
    text = _text(value)
    return FAMILY_ALIASES.get(text, text)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _stable_id(*parts: object) -> str:
    payload = "\x1f".join(_text(part) for part in parts)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


def _lineage_key(value: object) -> str:
    text = unicodedata.normalize("NFKD", _text(value))
    text = "".join(char for char in text if not unicodedata.combining(char))
    text = text.encode("ascii", errors="ignore").decode("ascii").casefold()
    return re.sub(r"[^a-z0-9]+", "_", text).strip("_")


def _read_jsonl(path: Path) -> list[dict[str, Any]]:
    if path.suffix == ".gz":
        with gzip.open(path, mode="rt", encoding="utf-8") as handle:
            return [json.loads(line) for line in handle if line.strip()]
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _wfo_accepted_species(record: dict[str, Any]) -> str:
    match = (record.get("wfo") or {}).get("match") or {}
    placement = _text(match.get("placement"))
    if not placement.startswith("Code/Plantae/"):
        return ""
    parts = [part for part in placement.split("$", 1)[0].split("/") if part]
    return " ".join(parts[-2:]) if len(parts) >= 2 else ""


def _expected_families(
    names_and_families: Iterable[tuple[str, str]],
) -> dict[str, str]:
    grouped: dict[str, set[str]] = {}
    for name, family in names_and_families:
        if len(name.split()) != 2:
            continue
        grouped.setdefault(name, set()).add(family)
    ambiguous = {name: values for name, values in grouped.items() if len(values) != 1}
    if ambiguous:
        raise ValueError(f"source names have ambiguous families: {sorted(ambiguous)[:5]}")
    return {name: next(iter(values)) for name, values in grouped.items()}


def _mapping_reason(
    record: dict[str, Any],
    *,
    source_name: str,
    source_family: str,
    target_species: set[str],
) -> tuple[str, str]:
    wfo = record.get("wfo") or {}
    gbif = record.get("gbif") or {}
    match = wfo.get("match") or {}
    accepted = _wfo_accepted_species(record)
    if not BINOMIAL.fullmatch(source_name):
        return "source_name_not_species_binomial", accepted
    if record.get("wfo_status") != 200 or not match:
        return "wfo_no_unambiguous_exact_match", accepted
    if (wfo.get("parsedName") or {}).get("rank") != "species":
        return "wfo_not_species_rank", accepted
    if _text((wfo.get("parsedName") or {}).get("canonical_form")) != source_name:
        return "wfo_canonical_input_mismatch", accepted
    if (wfo.get("params") or {}).get("fuzzyNameParts") != 0:
        return "wfo_fuzzy_match_rejected", accepted
    if not (wfo.get("params") or {}).get("checkHomonyms"):
        return "wfo_homonym_check_not_enabled", accepted
    if not (wfo.get("params") or {}).get("checkRank"):
        return "wfo_rank_check_not_enabled", accepted
    if record.get("gbif_status") != 200:
        return "gbif_no_response", accepted
    if gbif.get("matchType") != "EXACT":
        return "gbif_not_exact", accepted
    if gbif.get("rank") != "SPECIES" or gbif.get("kingdom") != "Plantae":
        return "gbif_not_plant_species", accepted
    if accepted != _text(gbif.get("species")):
        return "backbones_disagree", accepted
    placement_parts = set(_text(match.get("placement")).split("/"))
    if source_family not in placement_parts or _text(gbif.get("family")) != source_family:
        return "family_conflict", accepted
    if not BINOMIAL.fullmatch(accepted):
        return "accepted_not_species_binomial", accepted
    if accepted in target_species:
        return "mapped_into_fixed_target_universe", accepted
    return "accepted_strict_two_backbone", accepted


def build_mapping_audit(
    *,
    records: list[dict[str, Any]],
    expected_families: dict[str, str],
    target_species: set[str],
    provider: str,
) -> pd.DataFrame:
    response_names = [_text(record.get("source_name")) for record in records]
    if len(response_names) != len(set(response_names)):
        raise ValueError(f"{provider} response checkpoint has duplicate source names")
    if set(response_names) != set(expected_families):
        missing = sorted(set(expected_families).difference(response_names))
        extra = sorted(set(response_names).difference(expected_families))
        raise ValueError(
            f"{provider} response checkpoint does not cover source names exactly; "
            f"missing={missing[:5]}, extra={extra[:5]}"
        )
    rows: list[dict[str, object]] = []
    for record in records:
        source_name = _text(record.get("source_name"))
        source_family = expected_families[source_name]
        reason, accepted = _mapping_reason(
            record,
            source_name=source_name,
            source_family=source_family,
            target_species=target_species,
        )
        wfo = record.get("wfo") or {}
        gbif = record.get("gbif") or {}
        match = wfo.get("match") or {}
        rows.append(
            {
                "provider": provider,
                "source_name": source_name,
                "source_family": source_family,
                "accepted_species": (
                    accepted if reason == "accepted_strict_two_backbone" else ""
                ),
                "mapping_reason": reason,
                "wfo_accepted_species": accepted,
                "wfo_match_id": _text(match.get("wfo_id")),
                "wfo_classification_version": _text(
                    (wfo.get("params") or {}).get("classificationVersion")
                ),
                "gbif_accepted_species": _text(gbif.get("species")),
                "gbif_family": _text(gbif.get("family")),
                "gbif_usage_key": _text(
                    gbif.get("acceptedUsageKey") or gbif.get("usageKey")
                ),
                "retrieved_at_utc": _text(record.get("retrieved_at_utc")),
            }
        )
    return pd.DataFrame(rows).sort_values(["provider", "source_name"]).reset_index(drop=True)


def _mapping_dict(audit: pd.DataFrame) -> dict[str, str]:
    accepted = audit.loc[audit["mapping_reason"].eq("accepted_strict_two_backbone")]
    return dict(zip(accepted["source_name"], accepted["accepted_species"], strict=True))


def build_bsdb_evidence(
    source: pd.DataFrame,
    mapping_audit: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    mapping = _mapping_dict(mapping_audit)
    selection: list[dict[str, object]] = []
    evidence: list[dict[str, str]] = []
    for index, row in source.fillna("").iterrows():
        source_name = " ".join(_text(row.get("tnrs_Sciname")).split()[:2])
        value = _text(row.get("BreedingSys"))
        source_key = _text(row.get("bs.Source"))
        reason = "selected"
        if value not in {"SC", "SI"}:
            reason = "unsupported_or_missing_breeding_system"
        elif _text(row.get("tnrs_infrasp")) or _text(row.get("Infrasp")):
            reason = "infraspecific_record"
        elif not source_key:
            reason = "missing_original_source_lineage"
        elif source_name not in mapping:
            reason = "not_strict_external_two_backbone"
        accepted = mapping.get(source_name, "")
        selection.append(
            {
                "provider": "bsdb",
                "source_row": index + 2,
                "source_name": source_name,
                "accepted_species": accepted,
                "trait_name": "self_incompatibility",
                "normalized_value": value,
                "selection_reason": reason,
            }
        )
        if reason != "selected":
            continue
        record = {column: "" for column in EVIDENCE_COLUMNS}
        record.update(
            {
                "accepted_species": accepted,
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "normalized_value": value,
                "quality": "high",
                "source_group": "zell_2025_bsdb_external_congener",
                "source_provider": "Zell et al. 2025 Breeding System Database",
                "source_url": BSDB_URL,
                "source_record_id": (
                    "bsdb-external:"
                    + _stable_id(index + 2, source_name, accepted, source_key, value)
                ),
                "source_citation": (
                    f"Zell et al. (2025), New Phytologist, DOI {BSDB_ARTICLE_DOI}; "
                    f"underlying source key: {source_key}"
                ),
                "source_excerpt": (
                    f"source_row={index + 2}; tnrs_Sciname={source_name}; "
                    f"strict_accepted_species={accepted}; "
                    f"tnrs_family={_text(row.get('tnrs_family'))}; "
                    f"BreedingSys={value}; "
                    f"ISI_value={_text(row.get('ISI_value')) or 'NA'}; "
                    f"ISItype={_text(row.get('ISItype')) or 'NA'}; "
                    f"bs.Source={source_key}"
                ),
                "evidence_scope": "external_congener_species_direct",
                "name_match_method": "strict_wfo_gbif_two_backbone",
                "source_lineage": f"bsdb-original:{_lineage_key(source_key)}",
                "lineage_method": "underlying_study_key_not_bsdb_redistributor",
                "source_run_id": f"dirtyplants-BSdb@{BSDB_COMMIT}",
                "source_artifact": BSDB_ARTIFACT,
                "source_file": BSDB_ARTIFACT,
                "acceptance_contract": CONTRACT,
            }
        )
        evidence.append(record)
    return (
        pd.DataFrame(evidence, columns=EVIDENCE_COLUMNS),
        pd.DataFrame(selection),
    )


def _rodger_state(row: pd.Series) -> str:
    values: list[float] = []
    for column in AUTO_COLUMNS:
        raw = _text(row.get(column))
        if not raw:
            continue
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


def _rodger_excerpt(rows: pd.DataFrame) -> str:
    excerpts: list[str] = []
    for record in rows.sort_values("source_row").to_dict("records"):
        measures = "; ".join(
            f"{column}={_text(record.get(column)) or 'NA'}" for column in AUTO_COLUMNS
        )
        excerpts.append(
            f"S3 row {_text(record.get('source_row'))}: "
            f"genus.species={_text(record.get('genus.species'))}; "
            f"TPL.family={_text(record.get('TPL.family'))}; "
            f"source.dataset={_text(record.get('source.dataset'))}; "
            f"study={_text(record.get('study'))}; {measures}"
        )
    return " || ".join(excerpts)


def build_rodger_evidence(
    source: pd.DataFrame,
    mapping_audit: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    mapping = _mapping_dict(mapping_audit)
    work = source.fillna("").copy()
    work["submitted_species"] = work["genus.species"].str.replace(
        "_", " ", regex=False
    ).str.strip()
    work["accepted_species"] = work["submitted_species"].map(mapping).fillna("")
    work["normalized_value"] = work.apply(_rodger_state, axis=1)
    work["selection_reason"] = "selected"
    work.loc[
        ~work["submitted_species"].map(lambda name: bool(BINOMIAL.fullmatch(name))),
        "selection_reason",
    ] = "not_exact_species_rank"
    work.loc[
        work["accepted_species"].eq(""), "selection_reason"
    ] = "not_strict_external_two_backbone"
    work.loc[
        work["normalized_value"].eq(""), "selection_reason"
    ] = "missing_or_invalid_pollinator_exclusion_measure"
    accepted = work.loc[work["selection_reason"].eq("selected")].copy()
    evidence: list[dict[str, str]] = []
    groups = [
        "submitted_species",
        "accepted_species",
        "source.dataset",
        "study",
        "normalized_value",
    ]
    for keys, rows in accepted.groupby(groups, sort=True):
        submitted, species, source_dataset, study, value = keys
        record = {column: "" for column in EVIDENCE_COLUMNS}
        record.update(
            {
                "accepted_species": species,
                "axis": "reproductive_assurance",
                "trait_name": "autonomous_selfing_capacity",
                "normalized_value": value,
                "quality": "high",
                "source_group": "rodger_2021_external_congener",
                "source_provider": (
                    "Rodger et al. 2021 Global Pollinator Contribution Dataset"
                ),
                "source_url": RODGER_URL,
                "source_record_id": (
                    f"rodger2021-external:{source_dataset}:"
                    + _stable_id(
                        submitted,
                        species,
                        study,
                        value,
                        "|".join(rows["source_row"]),
                    )
                ),
                "source_citation": (
                    f"Rodger et al. (2021), Science Advances, DOI "
                    f"{RODGER_ARTICLE_DOI}; source dataset={source_dataset}; "
                    f"original study={study}"
                ),
                "source_excerpt": _rodger_excerpt(rows),
                "evidence_scope": "external_congener_species_direct",
                "name_match_method": "strict_wfo_gbif_two_backbone",
                "source_lineage": f"rodger2021-original-study:{_stable_id(study)}",
                "lineage_method": (
                    "original_study_citation_not_database_redistributor"
                ),
                "source_run_id": RODGER_RUN_ID,
                "source_artifact": RODGER_ARTIFACT,
                "source_file": "rodger_2021_pollinator_contribution_snapshot.csv.gz",
                "acceptance_contract": CONTRACT,
            }
        )
        evidence.append(record)
    selection = work[
        [
            "source_row",
            "submitted_species",
            "accepted_species",
            "normalized_value",
            "selection_reason",
        ]
    ].rename(columns={"submitted_species": "source_name"})
    selection.insert(0, "provider", "rodger_2021")
    selection.insert(4, "trait_name", "autonomous_selfing_capacity")
    return pd.DataFrame(evidence, columns=EVIDENCE_COLUMNS), selection


def build_checkpoint(
    *,
    bsdb_source: pd.DataFrame,
    bsdb_records: list[dict[str, Any]],
    rodger_source: pd.DataFrame,
    rodger_records: list[dict[str, Any]],
    master: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    target_species = set(master["accepted_species"].map(_text))
    bsdb_names = _expected_families(
        (
            " ".join(_text(row.get("tnrs_Sciname")).split()[:2]),
            _family(row.get("tnrs_family")),
        )
        for row in bsdb_source.to_dict("records")
        if _text(row.get("BreedingSys")) in {"SC", "SI"}
        and " ".join(_text(row.get("tnrs_Sciname")).split()[:2]) not in target_species
    )
    rodger_names = _expected_families(
        (
            _text(row.get("genus.species")).replace("_", " "),
            _family(row.get("TPL.family")),
        )
        for row in rodger_source.to_dict("records")
        if _text(row.get("genus.species")).replace("_", " ") not in target_species
    )
    bsdb_mapping = build_mapping_audit(
        records=bsdb_records,
        expected_families=bsdb_names,
        target_species=target_species,
        provider="bsdb",
    )
    rodger_mapping = build_mapping_audit(
        records=rodger_records,
        expected_families=rodger_names,
        target_species=target_species,
        provider="rodger_2021",
    )
    bsdb_evidence, bsdb_selection = build_bsdb_evidence(bsdb_source, bsdb_mapping)
    rodger_evidence, rodger_selection = build_rodger_evidence(
        rodger_source, rodger_mapping
    )
    evidence = pd.concat([bsdb_evidence, rodger_evidence], ignore_index=True)
    evidence = evidence.sort_values(
        ["accepted_species", "trait_name", "source_lineage", "source_record_id"]
    ).reset_index(drop=True)
    if evidence["source_record_id"].duplicated().any():
        raise ValueError("external congener checkpoint has duplicate source record IDs")
    if set(evidence["accepted_species"]).intersection(target_species):
        raise ValueError("external congener evidence entered the fixed target universe")
    mapping = pd.concat([bsdb_mapping, rodger_mapping], ignore_index=True).sort_values(
        ["provider", "source_name"]
    )
    selection = pd.concat(
        [bsdb_selection, rodger_selection], ignore_index=True, sort=False
    ).fillna("")
    return evidence, mapping.reset_index(drop=True), selection.reset_index(drop=True)


@app.command("build")
def build(
    bsdb_source_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    bsdb_responses_jsonl: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    rodger_source_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    rodger_responses_jsonl: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    master_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    angiosperm_scope_yaml: Annotated[
        Path, typer.Option(exists=True, dir_okay=False)
    ],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
    expected_species: Annotated[int, typer.Option(min=1)] = 106_295,
) -> None:
    if _sha256(bsdb_source_csv) != BSDB_SHA256:
        raise ValueError("BSDB source hash does not match pinned commit artifact")
    if _sha256(rodger_source_csv) != RODGER_SHA256:
        raise ValueError("Rodger source snapshot hash mismatch")
    bsdb_source = pd.read_csv(bsdb_source_csv, dtype=str).fillna("")
    rodger_source = pd.read_csv(rodger_source_csv, dtype=str).fillna("")
    raw_master = pd.read_csv(master_csv, dtype=str).fillna("")
    scope = classify_scope(raw_master, load_scope_config(angiosperm_scope_yaml))
    eligible_species = set(
        scope.loc[scope["angiosperm_analysis_eligible"], "accepted_species"]
    )
    master = raw_master.loc[raw_master["accepted_species"].isin(eligible_species)].copy()
    if master["accepted_species"].nunique() != expected_species:
        raise ValueError(
            f"fixed target has {master['accepted_species'].nunique()} species, "
            f"expected {expected_species}"
        )
    evidence, mapping, selection = build_checkpoint(
        bsdb_source=bsdb_source,
        bsdb_records=_read_jsonl(bsdb_responses_jsonl),
        rodger_source=rodger_source,
        rodger_records=_read_jsonl(rodger_responses_jsonl),
        master=master,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "wave41_external_congener_reproductive_evidence.csv.gz"
    mapping_path = output_dir / "wave41_external_congener_mapping_audit.csv.gz"
    selection_path = output_dir / "wave41_external_congener_selection_audit.csv.gz"
    for frame, path in (
        (evidence, evidence_path),
        (mapping, mapping_path),
        (selection, selection_path),
    ):
        frame.to_csv(
            path,
            index=False,
            compression={"method": "gzip", "mtime": 0},
            lineterminator="\n",
        )
    summary = {
        "contract": "wave41_external_congener_reproduction_checkpoint_v1",
        "fixed_target_species": int(master["accepted_species"].nunique()),
        "query_cost_usd": 0,
        "queries": {
            "wfo": len(mapping),
            "gbif": len(mapping),
            "total": len(mapping) * 2,
        },
        "mapping_reasons": {
            provider: {
                str(key): int(value)
                for key, value in group["mapping_reason"].value_counts().items()
            }
            for provider, group in mapping.groupby("provider", sort=True)
        },
        "evidence": {
            "rows": len(evidence),
            "species": int(evidence["accepted_species"].nunique()),
            "species_trait": int(
                evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
            ),
            "source_lineages": int(evidence["source_lineage"].nunique()),
            "by_trait": {
                str(key): int(value)
                for key, value in evidence["trait_name"].value_counts().items()
            },
            "by_value": {
                str(key): int(value)
                for key, value in evidence["normalized_value"].value_counts().items()
            },
            "entered_fixed_target_direct_coverage": 0,
        },
        "selection_reasons": {
            provider: {
                str(key): int(value)
                for key, value in group["selection_reason"].value_counts().items()
            }
            for provider, group in selection.groupby("provider", sort=True)
        },
        "input_sha256": {
            "bsdb_source_csv": _sha256(bsdb_source_csv),
            "bsdb_responses_jsonl": _sha256(bsdb_responses_jsonl),
            "rodger_source_csv": _sha256(rodger_source_csv),
            "rodger_responses_jsonl": _sha256(rodger_responses_jsonl),
            "master_csv": _sha256(master_csv),
            "angiosperm_scope_yaml": _sha256(angiosperm_scope_yaml),
        },
        "output_sha256": {
            path.name: _sha256(path)
            for path in (evidence_path, mapping_path, selection_path)
        },
        "checks": {
            "fixed_target_excluded_from_external_evidence": True,
            "wfo_gbif_species_and_family_agreement_required": True,
            "external_evidence_not_confirmatory_direct": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
            "trait_substitution_absent": True,
        },
    }
    summary_path = output_dir / "wave41_external_congener_checkpoint_summary.json"
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    typer.echo(json.dumps(summary, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    app()
