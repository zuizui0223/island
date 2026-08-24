"""Recover BSdb reproductive evidence hidden behind strict synonyms.

The acquisition queries every BSdb species name that failed the fixed-master
join against the public WFO matching service and the GBIF species matcher.
This checkpoint is deliberately fail-closed: a synonym is usable only when
both backbones independently resolve an exact species-rank name to the same
fixed-master species and all three family assignments agree.  Frozen response
records make the decision reproducible without repeating the Web requests.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import unicodedata
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.targeted_support2_wave15_checkpoint import (
    AUDIT_COLUMNS,
    EVIDENCE_COLUMNS,
    _candidate_id,
)
from island_v2.targeted_support2_wave15_checkpoint import (
    build_audit as _base_build_audit,
)

SOURCE_COMMIT = "9e87946d1e3121d39e657b702cf9f92ccc10936e"
SOURCE_URL = (
    "https://raw.githubusercontent.com/dirtyplants/BSdb/"
    f"{SOURCE_COMMIT}/Zell_df_12_29_23.csv"
)
SOURCE_SHA256 = "f935c8150b3d719aba5f62f14335c1a185304155403bf50db1e3ef1393fc55f4"
ARTICLE_DOI = "10.1111/nph.20234"
SOURCE_GROUP = "bsdb_two_backbone_synonym_wave17_20260820"
RETRIEVED_AT = "2026-08-20T08:00:00Z"


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _lineage_key(value: object) -> str:
    text = unicodedata.normalize("NFKD", _text(value))
    text = "".join(char for char in text if not unicodedata.combining(char))
    text = text.encode("ascii", errors="ignore").decode("ascii").casefold()
    return re.sub(r"[^a-z0-9]+", "_", text).strip("_")


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_response_records(path: Path) -> list[dict[str, Any]]:
    if path.suffix == ".gz":
        with gzip.open(path, mode="rt", encoding="utf-8") as handle:
            return [json.loads(line) for line in handle if line.strip()]
    with path.open(mode="rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _wfo_accepted(record: dict[str, Any]) -> str:
    match = (record.get("wfo") or {}).get("match") or {}
    placement = _text(match.get("placement"))
    if not placement:
        return ""
    path = placement.split("$", 1)[0].split("/")
    return " ".join(path[-2:]) if len(path) >= 2 else ""


def _mapping_reason(
    record: dict[str, Any], master_family: dict[str, str]
) -> tuple[str, str]:
    wfo = record.get("wfo") or {}
    gbif = record.get("gbif") or {}
    match = wfo.get("match") or {}
    accepted = _wfo_accepted(record)
    if record.get("wfo_status") != 200 or not match:
        return "wfo_no_exact_match", ""
    if wfo.get("parsedName", {}).get("rank") != "species":
        return "wfo_not_species_rank", ""
    if wfo.get("params", {}).get("fuzzyNameParts") != 0:
        return "wfo_fuzzy_match_rejected", ""
    if record.get("gbif_status") != 200:
        return "gbif_no_response", ""
    if gbif.get("matchType") != "EXACT":
        return "gbif_not_exact", ""
    if gbif.get("rank") != "SPECIES" or gbif.get("kingdom") != "Plantae":
        return "gbif_not_plant_species", ""
    if accepted != _text(gbif.get("species")):
        return "backbones_disagree", ""
    if accepted not in master_family:
        return "agreed_name_not_in_fixed_master", ""
    families = {
        _text(record.get("source_family")),
        _text(gbif.get("family")),
        master_family[accepted],
    }
    if "" in families or len(families) != 1:
        return "family_conflict", ""
    return "accepted_strict_two_backbone", accepted


def build_mapping_audit(
    records: Iterable[dict[str, Any]], master: pd.DataFrame
) -> pd.DataFrame:
    required_master = {"accepted_species", "family"}
    if missing := required_master.difference(master.columns):
        raise ValueError(f"master missing columns: {sorted(missing)}")
    master_family = {
        _text(row.accepted_species): _text(row.family)
        for row in master.itertuples(index=False)
    }
    rows: list[dict[str, str]] = []
    for record in records:
        wfo = record.get("wfo") or {}
        gbif = record.get("gbif") or {}
        reason, accepted = _mapping_reason(record, master_family)
        canonical = json.dumps(record, ensure_ascii=False, sort_keys=True)
        rows.append(
            {
                "source_name": _text(record.get("source_name")),
                "source_family": _text(record.get("source_family")),
                "accepted_species": accepted,
                "mapping_reason": reason,
                "wfo_classification_version": _text(
                    wfo.get("params", {}).get("classificationVersion")
                ),
                "wfo_match_id": _text((wfo.get("match") or {}).get("wfo_id")),
                "wfo_accepted_species": _wfo_accepted(record),
                "gbif_match_type": _text(gbif.get("matchType")),
                "gbif_rank": _text(gbif.get("rank")),
                "gbif_kingdom": _text(gbif.get("kingdom")),
                "gbif_usage_key": _text(gbif.get("usageKey")),
                "gbif_accepted_usage_key": _text(
                    gbif.get("acceptedUsageKey") or gbif.get("usageKey")
                ),
                "gbif_accepted_species": _text(gbif.get("species")),
                "gbif_family": _text(gbif.get("family")),
                "reviewer": "Codex strict mapping audit",
                "reviewed_at_utc": RETRIEVED_AT,
                "response_record_sha256": hashlib.sha256(
                    canonical.encode("utf-8")
                ).hexdigest(),
            }
        )
    return pd.DataFrame(rows)


def build_bsdb_evidence(
    source: pd.DataFrame,
    mapping_audit: pd.DataFrame,
    baseline_direct: pd.DataFrame,
    master: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return novel direct evidence and a complete BSdb-row selection audit."""
    required_source = {
        "Infrasp",
        "BreedingSys",
        "ISI_value",
        "ISItype",
        "bs.Source",
        "tnrs_family",
        "tnrs_Sciname",
        "tnrs_infrasp",
    }
    if missing := required_source.difference(source.columns):
        raise ValueError(f"BSdb source missing columns: {sorted(missing)}")
    if missing := {"accepted_species", "trait_name"}.difference(
        baseline_direct.columns
    ):
        raise ValueError(f"baseline direct missing columns: {sorted(missing)}")
    if missing := {"accepted_species", "family"}.difference(master.columns):
        raise ValueError(f"master missing columns: {sorted(missing)}")

    accepted_mapping = mapping_audit.loc[
        mapping_audit["mapping_reason"].eq("accepted_strict_two_backbone")
    ].copy()
    if accepted_mapping["source_name"].duplicated().any():
        raise ValueError("strict mapping audit contains duplicate source names")
    mapping_by_name = accepted_mapping.set_index("source_name").to_dict("index")
    current_pairs = {
        (_text(row.accepted_species), _text(row.trait_name))
        for row in baseline_direct.itertuples(index=False)
    }
    master_family = dict(zip(master["accepted_species"], master["family"], strict=False))

    selection_rows: list[dict[str, str]] = []
    evidence_rows: list[dict[str, str]] = []
    for source_index, row in enumerate(source.fillna("").to_dict("records"), start=2):
        source_name = _text(row["tnrs_Sciname"])
        if source_name not in mapping_by_name:
            continue
        mapping = mapping_by_name[source_name]
        accepted = _text(mapping["accepted_species"])
        value = _text(row["BreedingSys"])
        reason = "selected"
        if _text(row["Infrasp"]) or _text(row["tnrs_infrasp"]):
            reason = "infraspecific_record"
        elif value not in {"SC", "SI"}:
            reason = "unsupported_breeding_system"
        elif not _text(row["bs.Source"]):
            reason = "missing_original_source_lineage"
        elif _text(row["tnrs_family"]) != _text(master_family.get(accepted)):
            reason = "family_conflict"
        elif (accepted, "self_incompatibility") in current_pairs:
            reason = "already_resolved_direct_pair"
        selection_rows.append(
            {
                "source_row": str(source_index),
                "source_name": source_name,
                "accepted_species": accepted,
                "BreedingSys": value,
                "bs.Source": _text(row["bs.Source"]),
                "selection_reason": reason,
            }
        )
        if reason != "selected":
            continue

        source_key = _text(row["bs.Source"])
        lineage = f"bsdb-original:{_lineage_key(source_key)}"
        excerpt = (
            f"tnrs_Sciname={source_name}; accepted_species={accepted}; "
            f"BreedingSys={value}; ISI_value={_text(row['ISI_value']) or 'NA'}; "
            f"ISItype={_text(row['ISItype']) or 'NA'}; bs.Source={source_key}"
        )
        evidence_rows.append(
            {
                "candidate_id": _candidate_id(
                    accepted, "self_incompatibility", lineage, value
                ),
                "accepted_species": accepted,
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "raw_value": value,
                "normalized_value": value,
                "evidence_quality": "high",
                "evidence_scope": "synonym_direct",
                "source_group": SOURCE_GROUP,
                "source_provider": "Zell et al. 2025 Breeding System Database",
                "source_url": SOURCE_URL,
                "page_title": "BSdb Zell_df_12_29_23.csv",
                "source_citation": (
                    f"Zell et al. (2025), New Phytologist, DOI {ARTICLE_DOI}; "
                    f"underlying source key: {source_key}"
                ),
                "source_excerpt": excerpt,
                "source_record_id": (
                    f"bsdb-synonym:{source_index}:{accepted}:self_incompatibility"
                ),
                "source_lineage": lineage,
                "lineage_method": "underlying_study_key_not_bsdb_redistributor",
                "name_resolution_lineage": (
                    f"wfo:{mapping['wfo_classification_version']}:"
                    f"{mapping['wfo_match_id']};gbif:usage:"
                    f"{mapping['gbif_usage_key']}:accepted:"
                    f"{mapping['gbif_accepted_usage_key']};"
                    "fixed_master_species_family_exact"
                ),
                "name_match_method": "synonym_exact_two_backbone",
                "matched_page_name": source_name,
                "source_tier": "A",
                "source_type": "peer_reviewed_species_trait_database",
                "domain": "raw.githubusercontent.com",
                "language": "en",
                "wild_cultivated_cultivar_status": (
                    "published_species_level_database_record_not_cultivar_limited"
                ),
                "evidence_status": "accepted_individual_source_backed_audit",
                "content_sha256": SOURCE_SHA256,
                "content_sha256_basis": "retrieved_pinned_bsdb_release_csv_bytes",
                "retrieved_at_utc": RETRIEVED_AT,
                "query": "wfo_2026_06_and_gbif_exact_synonym_recovery",
                "search_rank": "",
                "inference_rule": "",
            }
        )
    evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS)
    selection = pd.DataFrame(selection_rows)
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("BSdb synonym evidence contains duplicate candidate IDs")
    return evidence, selection


def build_review_audit(evidence: pd.DataFrame) -> pd.DataFrame:
    audit = _base_build_audit(evidence)
    audit["reviewer"] = "Codex Wave17 BSdb and regional-flora audit"
    audit["reviewed_at_utc"] = RETRIEVED_AT
    audit["decision_reason"] = (
        "Accepted after fixed-master family check, strict WFO and GBIF exact synonym "
        "agreement or exact regional-flora species treatment, direct-trait ontology "
        "review, original-lineage preservation, and cultivar screening."
    )
    return audit.loc[:, AUDIT_COLUMNS]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-csv", required=True, type=Path)
    parser.add_argument("--responses-jsonl", required=True, type=Path)
    parser.add_argument("--baseline-direct-csv", required=True, type=Path)
    parser.add_argument("--master-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    source = pd.read_csv(args.source_csv, dtype=str).fillna("")
    baseline = pd.read_csv(args.baseline_direct_csv, dtype=str).fillna("")
    master = pd.read_csv(args.master_csv, dtype=str).fillna("")
    records = read_response_records(args.responses_jsonl)
    mapping = build_mapping_audit(records, master)
    evidence, selection = build_bsdb_evidence(source, mapping, baseline, master)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(
        args.output_dir / "bsdb_two_backbone_mapping_audit.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
        lineterminator="\n",
    )
    evidence.to_csv(
        args.output_dir / "bsdb_two_backbone_evidence.csv",
        index=False,
        lineterminator="\n",
    )
    selection.to_csv(
        args.output_dir / "bsdb_two_backbone_selection_audit.csv",
        index=False,
        lineterminator="\n",
    )
    print(
        json.dumps(
            {
                "source_rows": len(source),
                "response_records": len(records),
                "strict_mappings": int(
                    mapping["mapping_reason"].eq("accepted_strict_two_backbone").sum()
                ),
                "selected_evidence_rows": len(evidence),
                "responses_sha256": _file_sha256(args.responses_jsonl),
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
