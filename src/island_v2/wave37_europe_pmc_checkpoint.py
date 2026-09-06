"""Fail-closed Wave 37 Europe PMC review checkpoint.

The acquisition lane searches high-value ``genus x reproductive trait`` tasks,
but this checkpoint only promotes individually reviewed, exact-species evidence.
Genus inference remains the responsibility of the shared trait-specific
all-evidence audit; this module never writes Validated Low rows itself.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.all_evidence_trait_audit import (
    apply_genus_rules,
    build_rule_audit,
    dedupe_direct_lineages,
    load_ontology,
    load_reviewed_direct_supplements,
    resolve_direct_cells,
)
from island_v2.wave35_trait_overlay import (
    AXES,
    QUALITY_RANK,
    _aggregate_direct,
    _aggregate_low,
    _parse_composition,
    _serialize_composition,
    _split_pipe,
    _validate_direct,
    _validate_low,
)

EXPECTED_QUEUE_ROWS = 500
EXPECTED_MACHINE_CANDIDATES = 44
EXPECTED_SPECIES = 106_295
SOURCE_GROUP = "europe_pmc_reproduction_wave37"
RETRIEVED_AT = "2026-08-28T03:58:12Z"
EVIDENCE_COLUMNS = [
    "accepted_species",
    "axis",
    "trait_name",
    "normalized_value",
    "quality",
    "source_group",
    "source_provider",
    "source_url",
    "source_record_id",
    "source_citation",
    "source_excerpt",
    "evidence_scope",
    "name_match_method",
    "source_lineage",
    "lineage_method",
    "source_run_id",
    "source_artifact",
    "source_file",
    "acceptance_contract",
]

MANUAL_RECORDS = (
    {
        "candidate_id": "MANUAL-Saussurea-polylepis-SI",
        "accepted_species": "Saussurea polylepis",
        "trait_name": "self_incompatibility",
        "normalized_value": "SI",
        "quality": "high",
        "pmcid": "PMC8031399",
        "doi": "10.1371/journal.pone.0249752",
        "source_excerpt": (
            "Although S . obvallata within Saussurea was reported as a "
            "self-compatible species [ 72 ], the fact that self-pollination could "
            "not produce mature fruit directly supports self-incompatibility of "
            "S . polylepis ."
        ),
        "source_citation": (
            "Genetic diversity and structure of Saussurea polylepis. PLOS ONE. "
            "DOI 10.1371/journal.pone.0249752."
        ),
        "name_match_method": "accepted_name_exact_manual_subject_correction",
        "acceptance_contract": (
            "primary_article_exact_species_controlled_self_pollination_high_v1"
        ),
        "review_reason": (
            "The parser attached the verb to 'Saussurea was'; manual review of the "
            "same full-text sentence identifies S. polylepis as the experimental "
            "subject and the failed self-pollination as direct SI evidence."
        ),
    },
    {
        "candidate_id": "MANUAL-Vitex-doniana-mating-system",
        "accepted_species": "Vitex doniana",
        "trait_name": "mating_system",
        "normalized_value": "mixed_mating",
        "quality": "medium",
        "pmcid": "PMC4972488",
        "doi": "10.1093/aobpla/plw051",
        "source_excerpt": (
            "Such a floral morphology, combined with flowering asynchrony could "
            "suggest a mixed mating system, combining both geitonogamy and "
            "outcrossing."
        ),
        "source_citation": (
            "Breeding biology of Vitex doniana. AoB PLANTS. DOI 10.1093/aobpla/plw051."
        ),
        "name_match_method": "accepted_name_exact",
        "acceptance_contract": (
            "primary_article_exact_species_interpretive_mating_statement_medium_v1"
        ),
        "review_reason": (
            "The primary article explicitly proposes mixed mating for V. doniana; "
            "the qualified wording is retained as Medium rather than High."
        ),
    },
    {
        "candidate_id": "MANUAL-Vitex-doniana-autonomous-selfing",
        "accepted_species": "Vitex doniana",
        "trait_name": "autonomous_selfing_capacity",
        "normalized_value": "absent",
        "quality": "high",
        "pmcid": "PMC4972488",
        "doi": "10.1093/aobpla/plw051",
        "source_excerpt": (
            "Despite the observed self-compatibility, V. doniana did not produce "
            "any fruits and seeds after spontaneous selfing ( SFI = 0)."
        ),
        "source_citation": (
            "Breeding biology of Vitex doniana. AoB PLANTS. DOI 10.1093/aobpla/plw051."
        ),
        "name_match_method": "accepted_name_exact",
        "acceptance_contract": (
            "primary_article_exact_species_spontaneous_selfing_experiment_high_v1"
        ),
        "review_reason": (
            "The primary article reports zero fruit and seed after spontaneous "
            "selfing; this is autonomous-selfing evidence and is not substituted "
            "for self-compatibility."
        ),
    },
)


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", mtime=0) as handle:
        handle.write(frame.to_csv(index=False, lineterminator="\n").encode("utf-8"))
    path.write_bytes(buffer.getvalue())


def _read_xml_text(path: Path) -> str:
    with gzip.open(path, "rb") as handle:
        root = ET.parse(handle).getroot()
    return " ".join(" ".join(root.itertext()).split())


def _verify_queue(queue_path: Path) -> pd.DataFrame:
    queue = pd.read_csv(queue_path, dtype=str).fillna("")
    required = {"genus", "trait_name", "current_support_actual", "state_count_actual"}
    missing = required.difference(queue.columns)
    if missing:
        raise ValueError(f"frozen queue missing columns: {sorted(missing)}")
    if len(queue) != EXPECTED_QUEUE_ROWS:
        raise ValueError(f"frozen queue has {len(queue)} rows, expected 500")
    if queue.duplicated(["genus", "trait_name"]).any():
        raise ValueError("frozen queue contains duplicate genus x trait tasks")
    if not queue["state_count_actual"].eq("1").all():
        raise ValueError("frozen queue contains a direct-state conflict")
    if not queue["current_support_actual"].isin({"1", "2"}).all():
        raise ValueError("frozen queue contains an unsupported support level")
    return queue


def _verify_backbone_snapshot(path: Path) -> dict[str, Any]:
    snapshot = json.loads(path.read_text(encoding="utf-8"))
    expected = "Samanea saman"
    placements = snapshot.get("wfo", {}).get("candidates", [])
    if not placements or any(
        "/Samanea/saman$Albizia/saman" not in _text(row.get("placement")) for row in placements
    ):
        raise ValueError("WFO candidates do not unanimously resolve Albizia saman")
    gbif = snapshot.get("gbif", {})
    if (
        gbif.get("match_type") != "EXACT"
        or gbif.get("rank") != "SPECIES"
        or gbif.get("status") != "SYNONYM"
        or gbif.get("accepted_species") != expected
        or gbif.get("family") != "Fabaceae"
    ):
        raise ValueError("GBIF does not exactly resolve Albizia saman to Samanea saman")
    if snapshot.get("accepted_species") != expected:
        raise ValueError("backbone snapshot accepted species changed")
    return snapshot


def _candidate_evidence(
    reviewed: pd.DataFrame,
    *,
    source_run_id: str,
    source_artifact: str,
    source_file: str,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    promoted = reviewed.loc[reviewed["promote_to_ledger"].str.casefold().eq("true")]
    for record in promoted.to_dict("records"):
        species = _text(record.get("accepted_species_override")) or _text(
            record.get("accepted_species")
        )
        value = _text(record.get("normalized_value_override")) or _text(
            record.get("normalized_value")
        )
        doi = _text(record.get("doi"))
        pmcid = _text(record.get("pmcid"))
        lineage = _text(record.get("source_lineage_override")) or _text(
            record.get("source_lineage")
        )
        name_method = _text(record.get("name_match_method"))
        if _text(record.get("accepted_species_override")):
            name_method = "strict_two_backbone_exact_synonym"
        rows.append(
            {
                "accepted_species": species,
                "axis": "reproductive_assurance",
                "trait_name": _text(record.get("trait_name")),
                "normalized_value": value,
                "quality": _text(record.get("quality")).casefold(),
                "source_group": SOURCE_GROUP,
                "source_provider": "Europe PMC full-text XML",
                "source_url": _text(record.get("source_url")),
                "source_record_id": (
                    f"{_text(record.get('candidate_id'))}:{species}:"
                    f"{_text(record.get('trait_name'))}"
                ),
                "source_citation": (
                    f"{_text(record.get('article_title'))}. DOI {doi or 'not assigned'}; {pmcid}."
                ),
                "source_excerpt": _text(record.get("exact_supporting_quote")),
                "evidence_scope": _text(record.get("evidence_scope")),
                "name_match_method": name_method,
                "source_lineage": lineage,
                "lineage_method": (
                    "cited_original_source_doi"
                    if lineage != _text(record.get("source_lineage"))
                    else "article_doi"
                ),
                "source_run_id": source_run_id,
                "source_artifact": source_artifact,
                "source_file": source_file,
                "acceptance_contract": ("reviewed_europe_pmc_exact_species_reproductive_trait_v1"),
            }
        )
    return rows


def _manual_evidence(
    raw_xml_dir: Path,
    *,
    source_run_id: str,
    source_artifact: str,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    rows: list[dict[str, str]] = []
    audits: list[dict[str, str]] = []
    xml_cache: dict[str, str] = {}
    for record in MANUAL_RECORDS:
        pmcid = record["pmcid"]
        xml_path = raw_xml_dir / f"{pmcid}.xml.gz"
        if not xml_path.exists():
            raise ValueError(f"manual-review XML missing: {xml_path.name}")
        xml_text = xml_cache.setdefault(pmcid, _read_xml_text(xml_path))
        quote = _text(record["source_excerpt"])
        if quote not in xml_text:
            raise ValueError(f"manual quote absent from fetched {pmcid} full text")
        lineage = f"doi:{record['doi']}"
        rows.append(
            {
                "accepted_species": record["accepted_species"],
                "axis": "reproductive_assurance",
                "trait_name": record["trait_name"],
                "normalized_value": record["normalized_value"],
                "quality": record["quality"],
                "source_group": SOURCE_GROUP,
                "source_provider": "Europe PMC full-text XML",
                "source_url": f"https://europepmc.org/articles/{pmcid}",
                "source_record_id": (
                    f"{record['candidate_id']}:{record['accepted_species']}:{record['trait_name']}"
                ),
                "source_citation": record["source_citation"],
                "source_excerpt": quote,
                "evidence_scope": "species_direct",
                "name_match_method": record["name_match_method"],
                "source_lineage": lineage,
                "lineage_method": "article_doi",
                "source_run_id": source_run_id,
                "source_artifact": source_artifact,
                "source_file": str(xml_path),
                "acceptance_contract": record["acceptance_contract"],
            }
        )
        audits.append(
            {
                "candidate_id": record["candidate_id"],
                "review_decision": "manual_correction_accepted",
                "promote_to_ledger": "true",
                "accepted_species": record["accepted_species"],
                "trait_name": record["trait_name"],
                "normalized_value": record["normalized_value"],
                "quality": record["quality"],
                "source_lineage": lineage,
                "source_url": f"https://europepmc.org/articles/{pmcid}",
                "exact_supporting_quote": quote,
                "review_reason": record["review_reason"],
                "reviewer": "Codex Wave37 source-backed reproduction audit",
                "reviewed_at_utc": RETRIEVED_AT,
            }
        )
    return rows, audits


def build_checkpoint(
    *,
    candidate_csv: Path,
    raw_xml_dir: Path,
    decisions_csv: Path,
    backbone_json: Path,
    frozen_queue_csv: Path,
    master_csv: Path,
    source_run_id: str,
    source_artifact: str,
    output_dir: Path,
    acquisition_manifest_json: Path | None = None,
    expected_species: int = EXPECTED_SPECIES,
) -> dict[str, Any]:
    queue = _verify_queue(frozen_queue_csv)
    _verify_backbone_snapshot(backbone_json)
    candidates = pd.read_csv(candidate_csv, dtype=str).fillna("")
    decisions = pd.read_csv(decisions_csv, dtype=str).fillna("")
    if len(candidates) != EXPECTED_MACHINE_CANDIDATES:
        raise ValueError(
            f"machine candidate count changed: {len(candidates)} != {EXPECTED_MACHINE_CANDIDATES}"
        )
    if (
        candidates["candidate_id"].duplicated().any()
        or decisions["candidate_id"].duplicated().any()
    ):
        raise ValueError("candidate and decision IDs must be unique")
    if set(candidates["candidate_id"]) != set(decisions["candidate_id"]):
        raise ValueError("every machine candidate must have exactly one frozen decision")

    reviewed = decisions.merge(
        candidates,
        on="candidate_id",
        how="left",
        validate="one_to_one",
        suffixes=("_decision", ""),
    )
    promoted_count = int(reviewed["promote_to_ledger"].str.casefold().eq("true").sum())
    if promoted_count != 8:
        raise ValueError(f"expected 8 promoted machine candidates, got {promoted_count}")

    master = pd.read_csv(master_csv, dtype=str).fillna("")
    master_species = set(master["accepted_species"])
    if len(master_species) < expected_species:
        raise ValueError(
            f"master has only {len(master_species)} species; fixed denominator is "
            f"{expected_species}"
        )
    candidate_rows = _candidate_evidence(
        reviewed,
        source_run_id=source_run_id,
        source_artifact=source_artifact,
        source_file=str(candidate_csv),
    )
    manual_rows, manual_audit = _manual_evidence(
        raw_xml_dir,
        source_run_id=source_run_id,
        source_artifact=source_artifact,
    )
    evidence = pd.DataFrame(candidate_rows + manual_rows, columns=EVIDENCE_COLUMNS)
    if len(evidence) != 11:
        raise ValueError(f"expected 11 reviewed evidence rows, got {len(evidence)}")
    if evidence.duplicated(["accepted_species", "trait_name", "source_lineage"]).any():
        raise ValueError("reviewed evidence contains duplicate species x trait x lineage")
    outside_master = sorted(set(evidence["accepted_species"]) - master_species)
    if outside_master:
        raise ValueError(f"promoted evidence is outside fixed master: {outside_master}")
    if not evidence["quality"].isin({"high", "medium"}).all():
        raise ValueError("reviewed evidence contains an invalid quality")
    if not evidence["evidence_scope"].isin({"species_direct", "synonym_direct"}).all():
        raise ValueError("reviewed evidence contains a non-direct scope")
    if (
        evidence[EVIDENCE_COLUMNS]
        .apply(lambda column: column.astype(str).str.strip().eq("").any())
        .any()
    ):
        raise ValueError("reviewed evidence contains incomplete provenance")

    machine_audit = reviewed.copy()
    machine_audit["reviewer"] = "Codex Wave37 source-backed reproduction audit"
    machine_audit["reviewed_at_utc"] = RETRIEVED_AT
    audit = pd.concat(
        [machine_audit, pd.DataFrame(manual_audit)],
        ignore_index=True,
        sort=False,
    ).fillna("")
    audit = audit.sort_values("candidate_id", kind="stable").reset_index(drop=True)
    evidence = evidence.sort_values(
        ["accepted_species", "trait_name", "source_lineage"], kind="stable"
    ).reset_index(drop=True)

    acquisition_manifest = (
        json.loads(acquisition_manifest_json.read_text(encoding="utf-8"))
        if acquisition_manifest_json is not None
        else {}
    )
    if acquisition_manifest and (
        int(acquisition_manifest.get("queries", -1)) != EXPECTED_QUEUE_ROWS
        or int(acquisition_manifest.get("candidate_rows", -1)) != EXPECTED_MACHINE_CANDIDATES
    ):
        raise ValueError("acquisition manifest does not match frozen Wave37 contract")

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "wave37_reviewed_direct_evidence.csv.gz"
    audit_path = output_dir / "wave37_candidate_review_audit.csv.gz"
    _write_gzip_csv(evidence, evidence_path)
    _write_gzip_csv(audit, audit_path)
    summary = {
        "contract": "wave37_europe_pmc_fail_closed_checkpoint_v1",
        "fixed_species_denominator": expected_species,
        "fixed_species_axis_denominator": expected_species * 3,
        "frozen_queue_rows": len(queue),
        "frozen_queue_sha256": _sha256(frozen_queue_csv),
        "machine_candidates": len(candidates),
        "machine_candidates_promoted": promoted_count,
        "manual_corrections_promoted": len(manual_rows),
        "reviewed_direct_evidence_rows": len(evidence),
        "reviewed_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "quality_counts": evidence["quality"].value_counts().sort_index().to_dict(),
        "trait_counts": evidence["trait_name"].value_counts().sort_index().to_dict(),
        "decision_counts": decisions["review_decision"].value_counts().sort_index().to_dict(),
        "source_run_id": source_run_id,
        "source_artifact": source_artifact,
        "acquisition": acquisition_manifest,
        "safeguards": {
            "search_snippets_as_evidence": False,
            "fulltext_required": True,
            "candidate_decision_coverage_required": True,
            "two_backbone_synonym_consensus_required": True,
            "trait_substitution": False,
            "validated_low_created_here": False,
            "family_inference": False,
            "global_fallback": False,
        },
        "input_sha256": {
            "candidates": _sha256(candidate_csv),
            "decisions": _sha256(decisions_csv),
            "backbone_snapshot": _sha256(backbone_json),
            "frozen_queue": _sha256(frozen_queue_csv),
        },
        "output_sha256": {
            evidence_path.name: _sha256(evidence_path),
            audit_path.name: _sha256(audit_path),
        },
    }
    (output_dir / "wave37_checkpoint_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def _coverage_summary(frame: pd.DataFrame, expected_species: int) -> dict[str, Any]:
    resolved = frame["quality"].ne("")
    per_species = frame.assign(_resolved=resolved).groupby("accepted_species")["_resolved"].sum()
    by_axis: dict[str, Any] = {}
    for axis in AXES:
        part = frame.loc[frame["axis"].eq(axis)]
        filled = int(part["quality"].ne("").sum())
        by_axis[axis] = {
            "denominator": expected_species,
            "filled_species": filled,
            "fill_rate": filled / expected_species,
            "unresolved_species": expected_species - filled,
            "quality_counts": {
                quality: int(part["quality"].eq(quality).sum())
                for quality in ("high", "medium", "low")
            },
        }
    return {
        "filled_species_axis": int(resolved.sum()),
        "unresolved_species_axis": int((~resolved).sum()),
        "quality_counts": {
            quality: int(frame["quality"].eq(quality).sum())
            for quality in ("high", "medium", "low")
        },
        "by_axis": by_axis,
        "species_by_filled_axis_count": {
            str(count): int(per_species.eq(count).sum()) for count in range(4)
        },
    }


def _validate_coverage(frame: pd.DataFrame, expected_species: int) -> None:
    required = {
        "accepted_species",
        "axis",
        "trait_composition",
        "trait_names",
        "source_groups",
        "source_lineages",
        "quality",
    }
    if missing := required.difference(frame.columns):
        raise ValueError(f"coverage lacks columns: {sorted(missing)}")
    expected_cells = expected_species * len(AXES)
    if len(frame) != expected_cells:
        raise ValueError(f"coverage has {len(frame)} rows, expected {expected_cells}")
    if frame[["accepted_species", "axis"]].duplicated().any():
        raise ValueError("coverage contains duplicate species x axis cells")
    if frame["accepted_species"].nunique() != expected_species:
        raise ValueError("coverage species denominator changed")
    if set(frame["axis"]) != set(AXES):
        raise ValueError("coverage axis set changed")
    if set(frame["quality"]) - set(QUALITY_RANK):
        raise ValueError("coverage contains an unknown quality")


def _eligible(frame: pd.DataFrame) -> pd.Series:
    return frame["eligible"].astype(str).str.casefold().eq("true")


def build_touched_rule_rebuild(
    *,
    formal_resolved_direct_csv: Path,
    checkpoint_dir: Path,
    supplemental_direct_evidence_csvs: tuple[Path, ...],
    master_csv: Path,
    ontology_yaml: Path,
    output_dir: Path,
) -> dict[str, Any]:
    """Rebuild only Wave37-touched genus x trait rules with common logic.

    The formal resolved ledger is expanded back to one row per retained source
    lineage. All lineages for a resolved species x trait necessarily support
    the selected state set, so this is sufficient for the common lineage LOO
    audit while avoiding a second import of the 118k-row public-web artifact.
    """
    evidence_path = checkpoint_dir / "wave37_reviewed_direct_evidence.csv.gz"
    reviewed = pd.read_csv(evidence_path, dtype=str).fillna("")
    touched_keys = {
        (str(row.accepted_species).split()[0], str(row.trait_name))
        for row in reviewed.itertuples(index=False)
    }
    formal = pd.read_csv(formal_resolved_direct_csv, dtype=str).fillna("")
    formal["genus"] = formal["accepted_species"].str.split().str[0]
    formal = formal.loc[
        formal.apply(
            lambda row: (str(row["genus"]), str(row["trait_name"])) in touched_keys,
            axis=1,
        )
        & formal["resolution_status"].eq("resolved")
    ].copy()
    reconstructed: list[dict[str, str]] = []
    for row in formal.to_dict("records"):
        lineages = _split_pipe(row.get("source_lineages", ""))
        if not lineages:
            raise ValueError("formal resolved direct cell lacks source lineage")
        for lineage in sorted(lineages):
            reconstructed.append(
                {
                    "accepted_species": _text(row.get("accepted_species")),
                    "axis": _text(row.get("axis")),
                    "trait_name": _text(row.get("trait_name")),
                    "normalized_value": _text(row.get("state_set")),
                    "quality": _text(row.get("quality")).casefold(),
                    "source_group": _text(row.get("source_groups")),
                    "source_provider": "formal_resolved_direct_lineage",
                    "source_url": "artifact://formal-wave33-resolved-direct",
                    "source_record_id": lineage,
                    "source_citation": "formal resolved direct lineage",
                    "source_excerpt": "",
                    "evidence_scope": "species_direct",
                    "name_match_method": "formal_resolved_species_identity",
                    "source_lineage": lineage,
                    "lineage_method": "expanded_from_formal_resolved_cell",
                    "source_run_id": "32932103226",
                    "source_artifact": "integrated-trait-coverage-32932103226",
                    "source_file": str(formal_resolved_direct_csv),
                    "acceptance_contract": "formal_resolved_direct_lineage_reuse_v1",
                }
            )
    formal_lineages = pd.DataFrame(reconstructed, columns=EVIDENCE_COLUMNS)
    supplemental = load_reviewed_direct_supplements(supplemental_direct_evidence_csvs)
    supplemental["genus"] = supplemental["accepted_species"].str.split().str[0]
    supplemental = supplemental.loc[
        supplemental.apply(
            lambda row: (str(row["genus"]), str(row["trait_name"])) in touched_keys,
            axis=1,
        )
    ].drop(columns="genus")
    raw = pd.concat([formal_lineages, supplemental], ignore_index=True).fillna("")
    ontology = load_ontology(ontology_yaml)
    lineages, duplicates = dedupe_direct_lineages(raw, ontology)
    direct_cells, conflicts = resolve_direct_cells(lineages)
    direct_cells["genus"] = direct_cells["accepted_species"].str.split().str[0]
    lineages["genus"] = lineages["accepted_species"].str.split().str[0]
    old_low = pd.DataFrame(columns=["genus", "trait_name", "state_set"])
    rules = build_rule_audit(direct_cells, lineages, old_low)
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    rebuilt_low = apply_genus_rules(master, direct_cells, rules, "current_min3")
    source_touched = direct_cells["source_groups"].map(
        lambda value: SOURCE_GROUP in _split_pipe(value)
    )
    wave37_direct = direct_cells.loc[source_touched].copy()
    if len(wave37_direct) != 7:
        raise ValueError(f"expected 7 resolved Wave37 direct cells, got {len(wave37_direct)}")
    vitex_rule = rules.loc[
        rules["setting"].eq("current_min3")
        & _eligible(rules)
        & rules["genus"].eq("Vitex")
        & rules["trait_name"].eq("mating_system")
    ]
    if len(vitex_rule) != 1:
        raise ValueError("touched rebuild did not validate Vitex x mating_system")
    vitex_low = rebuilt_low.loc[
        rebuilt_low["genus"].eq("Vitex") & rebuilt_low["trait_name"].eq("mating_system")
    ].copy()
    if len(vitex_low) != 116:
        raise ValueError(f"expected 116 Vitex Low candidates, got {len(vitex_low)}")

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "resolved_direct_species_trait.csv.gz": wave37_direct,
        "rebuilt_all_evidence_validated_low.csv.gz": vitex_low,
        "trait_specific_genus_rule_audit.csv.gz": rules,
        "wave37_touched_source_lineages.csv.gz": lineages,
        "wave37_touched_source_lineage_duplicates.csv.gz": duplicates,
        "wave37_touched_source_lineage_conflicts.csv.gz": conflicts,
    }
    for name, frame in outputs.items():
        _write_gzip_csv(frame, output_dir / name)
    summary = {
        "contract": "wave37_common_trait_specific_touched_rebuild_v1",
        "formal_resolved_cells_reused": len(formal),
        "formal_source_lineages_reconstructed": len(formal_lineages),
        "supplemental_touched_rows": len(supplemental),
        "touched_genus_trait": len(touched_keys),
        "resolved_wave37_direct_species_trait": len(wave37_direct),
        "validated_vitex_mating_rule": True,
        "vitex_low_candidates_before_baseline_filter": len(vitex_low),
        "checks": {
            "shared_trait_specific_rule_builder": True,
            "genus_axis_trait_join": True,
            "lineage_leave_one_out_required": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
        },
        "input_sha256": {
            formal_resolved_direct_csv.name: _sha256(formal_resolved_direct_csv),
            evidence_path.name: _sha256(evidence_path),
            ontology_yaml.name: _sha256(ontology_yaml),
            **{path.name: _sha256(path) for path in supplemental_direct_evidence_csvs},
        },
        "artifact_sha256": {name: _sha256(output_dir / name) for name in outputs},
    }
    (output_dir / "wave37_touched_rebuild_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def build_formal_overlay(
    *,
    baseline_csv: Path,
    all_evidence_dir: Path,
    previous_rule_audit_csv: Path,
    checkpoint_dir: Path,
    output_dir: Path,
    expected_species: int = EXPECTED_SPECIES,
    formal_wave33_run_id: int = 32932103226,
    formal_wave36_run_id: int = 33137984367,
) -> dict[str, Any]:
    """Build the lossless Wave37 overlay from shared trait-specific outputs."""
    evidence_path = checkpoint_dir / "wave37_reviewed_direct_evidence.csv.gz"
    checkpoint_summary_path = checkpoint_dir / "wave37_checkpoint_summary.json"
    resolved_path = all_evidence_dir / "resolved_direct_species_trait.csv.gz"
    rebuilt_low_path = all_evidence_dir / "rebuilt_all_evidence_validated_low.csv.gz"
    rules_path = all_evidence_dir / "trait_specific_genus_rule_audit.csv.gz"
    required = (
        baseline_csv,
        evidence_path,
        checkpoint_summary_path,
        resolved_path,
        rebuilt_low_path,
        rules_path,
        previous_rule_audit_csv,
    )
    if missing := [str(path) for path in required if not path.is_file()]:
        raise ValueError(f"Wave37 overlay inputs missing: {missing}")

    baseline = pd.read_csv(baseline_csv, dtype=str).fillna("")
    _validate_coverage(baseline, expected_species)
    reviewed = pd.read_csv(evidence_path, dtype=str).fillna("")
    resolved = pd.read_csv(resolved_path, dtype=str).fillna("")
    rebuilt_low = pd.read_csv(rebuilt_low_path, dtype=str).fillna("")
    rules = pd.read_csv(rules_path, dtype=str).fillna("")
    previous_rules = pd.read_csv(previous_rule_audit_csv, dtype=str).fillna("")

    touched_keys = {
        (str(row.accepted_species).split()[0], str(row.trait_name))
        for row in reviewed.itertuples(index=False)
    }
    source_touched = resolved["source_groups"].map(lambda value: SOURCE_GROUP in _split_pipe(value))
    direct = _validate_direct(resolved.loc[source_touched].copy())
    if len(direct) != 7:
        raise ValueError(f"expected 7 resolved Wave37 direct cells, got {len(direct)}")

    current_rules = rules.loc[rules["setting"].eq("current_min3") & _eligible(rules)].copy()
    prior_keys = {
        (str(row.genus), str(row.trait_name))
        for row in previous_rules.loc[
            previous_rules["setting"].eq("current_min3") & _eligible(previous_rules)
        ].itertuples(index=False)
    }
    new_rules = current_rules.loc[
        current_rules.apply(
            lambda row: (
                (str(row["genus"]), str(row["trait_name"])) in touched_keys
                and (str(row["genus"]), str(row["trait_name"])) not in prior_keys
            ),
            axis=1,
        )
    ].copy()
    if set(zip(new_rules["genus"], new_rules["trait_name"], strict=True)) != {
        ("Vitex", "mating_system")
    }:
        raise ValueError("Wave37 must create exactly the Vitex x mating_system rule")

    low = rebuilt_low.merge(
        new_rules[["genus", "axis", "trait_name"]],
        on=["genus", "axis", "trait_name"],
        how="inner",
        validate="many_to_one",
    )
    baseline_keys = baseline[["accepted_species", "axis", "quality"]]
    low = low.merge(
        baseline_keys,
        on=["accepted_species", "axis"],
        how="left",
        validate="many_to_one",
        suffixes=("", "_baseline"),
    )
    low = low.loc[low["quality_baseline"].eq("")].drop(columns="quality_baseline")
    _validate_low(low, new_rules)
    if len(low) != 114:
        raise ValueError(f"expected 114 new Vitex Low rows, got {len(low)}")

    result = baseline.copy().set_index(["accepted_species", "axis"])
    changes: list[dict[str, str]] = []
    for row in _aggregate_direct(direct).itertuples(index=False):
        key = (row.accepted_species, row.axis)
        if key not in result.index:
            raise ValueError(f"direct evidence outside denominator: {key}")
        before_quality = str(result.loc[key, "quality"])
        if before_quality in {"", "low"}:
            composition = _parse_composition(row.trait_composition, axis=row.axis)
            groups = _split_pipe(row.source_groups)
            lineages = _split_pipe(row.source_lineages)
        else:
            composition = _parse_composition(result.loc[key, "trait_composition"], axis=row.axis)
            composition.update(_parse_composition(row.trait_composition, axis=row.axis))
            groups = _split_pipe(result.loc[key, "source_groups"]) | _split_pipe(row.source_groups)
            lineages = _split_pipe(result.loc[key, "source_lineages"]) | _split_pipe(
                row.source_lineages
            )
        after_quality = max((before_quality, row.quality), key=QUALITY_RANK.__getitem__)
        result.loc[key, "trait_composition"] = _serialize_composition(composition)
        result.loc[key, "trait_names"] = "|".join(sorted(composition))
        result.loc[key, "source_groups"] = "|".join(sorted(groups))
        result.loc[key, "source_lineages"] = "|".join(sorted(lineages))
        result.loc[key, "quality"] = after_quality
        changes.append(
            {
                "accepted_species": row.accepted_species,
                "axis": row.axis,
                "action": (
                    "direct_fill"
                    if before_quality == ""
                    else (
                        "direct_upgrade"
                        if QUALITY_RANK[after_quality] > QUALITY_RANK[before_quality]
                        else "direct_enrichment"
                    )
                ),
                "quality_before": before_quality,
                "quality_after": after_quality,
                "trait_names": row.trait_names,
                "source_groups": row.source_groups,
            }
        )

    for row in _aggregate_low(low).itertuples(index=False):
        key = (row.accepted_species, row.axis)
        if str(result.loc[key, "quality"]):
            raise ValueError(f"Low attempts to replace a resolved cell: {key}")
        for column in (
            "trait_composition",
            "trait_names",
            "source_groups",
            "source_lineages",
            "quality",
        ):
            result.loc[key, column] = getattr(row, column)
        changes.append(
            {
                "accepted_species": row.accepted_species,
                "axis": row.axis,
                "action": "validated_low_fill",
                "quality_before": "",
                "quality_after": "low",
                "trait_names": row.trait_names,
                "source_groups": row.source_groups,
            }
        )

    result = result.reset_index().sort_values(["accepted_species", "axis"])
    _validate_coverage(result, expected_species)
    before = _coverage_summary(baseline, expected_species)
    after = _coverage_summary(result, expected_species)
    changes_frame = pd.DataFrame(changes).sort_values(["accepted_species", "axis", "action"])
    comparison = baseline[["accepted_species", "axis", "quality"]].merge(
        result[["accepted_species", "axis", "quality"]],
        on=["accepted_species", "axis"],
        suffixes=("_before", "_after"),
        validate="one_to_one",
    )
    was_filled = comparison["quality_before"].ne("")
    now_filled = comparison["quality_after"].ne("")
    loss = int((was_filled & ~now_filled).sum())
    gain = int((~was_filled & now_filled).sum())
    if loss:
        raise ValueError(f"Wave37 coverage loss must be zero, observed {loss}")
    if (
        comparison["quality_after"].map(QUALITY_RANK)
        < comparison["quality_before"].map(QUALITY_RANK)
    ).any():
        raise ValueError("Wave37 downgraded an existing quality")
    if gain != 116:
        raise ValueError(f"Wave37 expected exact gain 116, observed {gain}")

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "wave37_species_axis_coverage.csv.gz": result,
        "wave37_resolved_direct_species_trait.csv.gz": direct,
        "wave37_candidate_validated_low_species_trait.csv.gz": low,
        "wave37_provider_touched_new_rule_audit.csv.gz": new_rules,
        "wave37_change_audit.csv.gz": changes_frame,
    }
    for name, frame in outputs.items():
        _write_gzip_csv(frame, output_dir / name)
    summary = {
        "contract": "wave37_lossless_europe_pmc_reproduction_overlay_v1",
        "formal_wave33_run_id": formal_wave33_run_id,
        "formal_wave36_run_id": formal_wave36_run_id,
        "wave36_before": before,
        "wave37_after": after,
        "delta": {
            "gross_gain_species_axis": gain,
            "loss_species_axis": loss,
            "net_gain_species_axis": gain,
            "by_axis_net_gain": {
                axis: after["by_axis"][axis]["filled_species"]
                - before["by_axis"][axis]["filled_species"]
                for axis in AXES
            },
            "action_counts": {
                str(key): int(value)
                for key, value in changes_frame["action"].value_counts().items()
            },
            "reviewed_direct_species_trait": len(reviewed),
            "resolved_direct_species_trait": len(direct),
            "new_validated_low_species_trait": len(low),
        },
        "new_eligible_rules": [
            f"{row.genus} x {row.trait_name}" for row in new_rules.itertuples(index=False)
        ],
        "checks": {
            "fixed_denominator": True,
            "quality_precedence_high_medium_low": True,
            "direct_conflicts_excluded": True,
            "trait_specific_genus_join": True,
            "lineage_leave_one_out_required": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
            "baseline_loss_zero": loss == 0,
        },
        "checkpoint": json.loads(checkpoint_summary_path.read_text(encoding="utf-8")),
        "input_sha256": {path.name: _sha256(path) for path in required},
    }
    summary["artifact_sha256"] = {name: _sha256(output_dir / name) for name in outputs}
    (output_dir / "wave37_coverage_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-csv", required=True, type=Path)
    parser.add_argument("--raw-xml-dir", required=True, type=Path)
    parser.add_argument("--decisions-csv", required=True, type=Path)
    parser.add_argument("--backbone-json", required=True, type=Path)
    parser.add_argument("--frozen-queue-csv", required=True, type=Path)
    parser.add_argument("--master-csv", required=True, type=Path)
    parser.add_argument("--source-run-id", required=True)
    parser.add_argument("--source-artifact", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--acquisition-manifest-json", type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    args = parser.parse_args()
    summary = build_checkpoint(
        candidate_csv=args.candidate_csv,
        raw_xml_dir=args.raw_xml_dir,
        decisions_csv=args.decisions_csv,
        backbone_json=args.backbone_json,
        frozen_queue_csv=args.frozen_queue_csv,
        master_csv=args.master_csv,
        source_run_id=args.source_run_id,
        source_artifact=args.source_artifact,
        output_dir=args.output_dir,
        acquisition_manifest_json=args.acquisition_manifest_json,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


def overlay_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-csv", required=True, type=Path)
    parser.add_argument("--all-evidence-dir", required=True, type=Path)
    parser.add_argument("--previous-rule-audit-csv", required=True, type=Path)
    parser.add_argument("--checkpoint-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    parser.add_argument("--formal-wave33-run-id", type=int, default=32932103226)
    parser.add_argument("--formal-wave36-run-id", type=int, default=33137984367)
    args = parser.parse_args()
    summary = build_formal_overlay(
        baseline_csv=args.baseline_csv,
        all_evidence_dir=args.all_evidence_dir,
        previous_rule_audit_csv=args.previous_rule_audit_csv,
        checkpoint_dir=args.checkpoint_dir,
        output_dir=args.output_dir,
        expected_species=args.expected_species,
        formal_wave33_run_id=args.formal_wave33_run_id,
        formal_wave36_run_id=args.formal_wave36_run_id,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


def touched_rebuild_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--formal-resolved-direct-csv", required=True, type=Path)
    parser.add_argument("--checkpoint-dir", required=True, type=Path)
    parser.add_argument(
        "--supplemental-direct-evidence-csv",
        action="append",
        default=[],
        type=Path,
    )
    parser.add_argument("--master-csv", required=True, type=Path)
    parser.add_argument("--ontology-yaml", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    summary = build_touched_rule_rebuild(
        formal_resolved_direct_csv=args.formal_resolved_direct_csv,
        checkpoint_dir=args.checkpoint_dir,
        supplemental_direct_evidence_csvs=tuple(args.supplemental_direct_evidence_csv),
        master_csv=args.master_csv,
        ontology_yaml=args.ontology_yaml,
        output_dir=args.output_dir,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
