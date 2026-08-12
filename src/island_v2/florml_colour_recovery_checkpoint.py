"""Recover organ-linked flower colour from acquired official FlorML treatments."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

from island_v2.all_evidence_trait_audit import canonical_trait_name
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS
from island_v2.official_flora_colour_recovery_checkpoint import (
    COLOUR_PATTERN,
    COMPARATIVE_TAXON,
    CULTIVAR_OR_HYBRID,
    TERM_TO_STATE,
    _normalized_excerpt,
    _state_set,
    _text,
    observed_colour_states,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

TRAIT = "flower_primary_color"
AUDIT_SIZE = 200
AUDIT_SEED = "florml-colour-treatment-audit-20260811-v1"
TARGET_PROVIDERS = frozenset({"flora_malesiana", "flora_of_the_guianas"})
STRICT_NAME_METHODS = frozenset(
    {
        "wfo_june_2026_exact_accepted_xml_species_family_consistent",
        "wfo_june_2026_exact_synonym_xml_species_family_consistent",
        "wfo_june_2026_exact_accepted_or_synonym_xml_species_family_consistent",
    }
)
COLOURED_COVERING = re.compile(
    r"\b(?:flowers?|petals?|perianths?|tepals?)\b[^.;]{0,120}"
    r"\b(?:white|whitish|cream|yellow|yellowish|orange|brown|brownish|red|reddish|pink|pinkish|purple|purplish)\b"
    r"[^,.;]{0,45}\b(?:hairs?|scales?|puberulous|puberulent|scurfy)\b",
    re.IGNORECASE,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _id(*parts: object, length: int = 40) -> str:
    payload = "\n".join(_text(part) for part in parts)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:length]


def high_precision_candidate(record: dict[str, Any]) -> tuple[bool, str]:
    """Apply the previously audited organ/state gate to FlorML treatments."""

    trait = canonical_trait_name(_text(record.get("trait_name")))
    quote = _text(record.get("source_excerpt"))
    expected = _state_set(record.get("normalized_value"))
    if _text(record.get("source_provider")) not in TARGET_PROVIDERS:
        return False, "provider_outside_florml_checkpoint"
    if trait != TRAIT:
        return False, "trait_outside_checkpoint"
    if _text(record.get("evidence_scope")).casefold() != "species_direct":
        return False, "not_species_direct"
    if _text(record.get("name_match_method")) not in STRICT_NAME_METHODS:
        return False, "name_match_not_strict_xml_exact"
    if not _text(record.get("accepted_species")) or not quote or not expected:
        return False, "missing_identity_quote_or_value"
    if CULTIVAR_OR_HYBRID.search(quote):
        return False, "cultivar_or_horticultural_hybrid_context"
    if COMPARATIVE_TAXON.search(quote):
        return False, "comparative_taxon_context"
    covering_states: set[str] = set()
    for covering in COLOURED_COVERING.finditer(quote):
        for colour in COLOUR_PATTERN.finditer(covering.group(0)):
            covering_states.add(TERM_TO_STATE[colour.group(0).casefold()])
    if expected.intersection(covering_states):
        return False, "coloured_hair_scale_or_pubescence_context"
    observed = observed_colour_states(quote)
    if not observed:
        return False, "no_organ_linked_flower_colour"
    if observed != expected:
        return False, "explicit_state_set_mismatch"
    return True, "selected_for_treatment_audit"


def _completed_pairs(frame: pd.DataFrame) -> set[tuple[str, str]]:
    work = frame.copy().fillna("")
    if "resolution_status" in work:
        work = work.loc[work["resolution_status"].eq("resolved")]
    work["trait_name"] = work["trait_name"].map(canonical_trait_name)
    return set(zip(work["accepted_species"], work["trait_name"], strict=True))


def _audit_template(evidence: pd.DataFrame) -> pd.DataFrame:
    sample = evidence.copy()
    sample["audit_hash"] = sample["source_record_id"].map(
        lambda value: _id(AUDIT_SEED, value, length=64)
    )
    guianas = sample.loc[sample["source_provider"].eq("flora_of_the_guianas")]
    guianas = guianas.sort_values("audit_hash", kind="stable").head(30)
    malesiana = sample.loc[sample["source_provider"].eq("flora_malesiana")]
    malesiana = malesiana.sort_values("audit_hash", kind="stable").head(
        AUDIT_SIZE - len(guianas)
    )
    sample = pd.concat([guianas, malesiana], ignore_index=True).sort_values(
        "audit_hash", kind="stable"
    )
    if len(sample) != AUDIT_SIZE:
        raise ValueError(f"FlorML colour audit requires {AUDIT_SIZE} treatments")
    return pd.DataFrame(
        {
            "candidate_id": sample["source_record_id"],
            "trait_name": sample["trait_name"],
            "accepted_species": sample["accepted_species"],
            "normalized_value": sample["normalized_value"],
            "source_provider": sample["source_provider"],
            "source_url": sample["source_url"],
            "source_citation": sample["source_citation"],
            "source_excerpt": sample["source_excerpt"],
            "content_fingerprint": sample["content_fingerprint"],
            "accepted_correct": "",
            "cultivar_status": "",
            "reviewer": "",
            "reviewed_at_utc": "",
            "audit_reason": "",
            "audit_hash": sample["audit_hash"],
        }
    ).reset_index(drop=True)


def build_checkpoint(
    candidates: pd.DataFrame, completed_direct: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if missing := set(EVIDENCE_COLUMNS).difference(candidates.columns):
        raise ValueError(f"FlorML cache missing columns: {sorted(missing)}")
    completed = _completed_pairs(completed_direct)
    evidence_rows: list[dict[str, str]] = []
    selection_rows: list[dict[str, str]] = []
    for record in candidates.fillna("").to_dict("records"):
        trait = canonical_trait_name(_text(record.get("trait_name")))
        species = _text(record.get("accepted_species"))
        selected, reason = high_precision_candidate(record)
        if selected and (species, trait) in completed:
            selected, reason = False, "already_acquired_direct_species_trait"
        selection_rows.append(
            {
                "source_record_id": _text(record.get("source_record_id")),
                "accepted_species": species,
                "trait_name": trait,
                "selected": str(selected).lower(),
                "selection_reason": reason,
            }
        )
        if not selected:
            continue
        row = {column: _text(record.get(column)) for column in EVIDENCE_COLUMNS}
        fingerprint = _id(species, _normalized_excerpt(row["source_excerpt"]), length=64)
        row.update(
            {
                "quality": "high",
                "acceptance_contract": "official_florml_organ_linked_complete_colour_state_reaudit_v1",
                "content_fingerprint": fingerprint,
            }
        )
        evidence_rows.append(row)
    evidence = pd.DataFrame(evidence_rows).drop_duplicates(
        ["accepted_species", "trait_name", "normalized_value", "source_lineage"]
    )
    evidence = evidence.sort_values(
        ["source_provider", "accepted_species", "source_record_id"], kind="stable"
    ).reset_index(drop=True)
    if evidence.empty or evidence["source_record_id"].duplicated().any():
        raise ValueError("FlorML colour recovery produced empty or duplicate evidence")
    selection = pd.DataFrame(selection_rows).sort_values(
        ["selected", "selection_reason", "accepted_species"], kind="stable"
    )
    return evidence, _audit_template(evidence), selection.reset_index(drop=True)


@app.command("build")
def build(
    candidates_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    completed_direct_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
) -> None:
    candidates = pd.read_csv(candidates_csv, dtype=str).fillna("")
    completed = pd.read_csv(completed_direct_csv, dtype=str).fillna("")
    evidence, audit, selection = build_checkpoint(candidates, completed)
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "florml_colour_evidence_20260811.csv.gz"
    audit_path = output_dir / "florml_colour_audit_200_20260811.csv"
    selection_path = output_dir / "florml_colour_selection_20260811.csv.gz"
    evidence.to_csv(evidence_path, index=False, compression={"method": "gzip", "mtime": 0})
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    selection.to_csv(selection_path, index=False, compression={"method": "gzip", "mtime": 0})
    manifest = {
        "contract": "official_florml_colour_recovery_checkpoint_v1",
        "candidate_evidence_rows": len(evidence),
        "species_trait": len(evidence[["accepted_species", "trait_name"]].drop_duplicates()),
        "audit_rows": len(audit),
        "audit_unit": "species_treatment_source_record_id",
        "providers": evidence["source_provider"].value_counts().sort_index().to_dict(),
        "inputs": {str(candidates_csv): _sha256(candidates_csv), str(completed_direct_csv): _sha256(completed_direct_csv)},
        "guardrails": {"family_inference": False, "global_fallback": False, "axis_only_join": False, "trait_substitution": False},
    }
    (output_dir / "florml_colour_manifest_20260811.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


@app.command("record-review")
def record_review(
    audit_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_csv: Annotated[Path, typer.Option(dir_okay=False)],
    reviewer: Annotated[str, typer.Option()],
    reviewed_at_utc: Annotated[str, typer.Option()],
    rejected_id: Annotated[list[str] | None, typer.Option()] = None,
) -> None:
    """Record decisions after the deterministic treatment sample was reviewed."""

    audit = pd.read_csv(audit_csv, dtype=str).fillna("")
    rejected = set(rejected_id or [])
    unknown = rejected.difference(audit["candidate_id"])
    if unknown:
        raise ValueError(f"unknown rejected candidate IDs: {sorted(unknown)}")
    accepted = ~audit["candidate_id"].isin(rejected)
    audit["accepted_correct"] = accepted.map({True: "true", False: "false"})
    audit["cultivar_status"] = "wild_or_unspecified"
    audit["reviewer"] = reviewer
    audit["reviewed_at_utc"] = reviewed_at_utc
    audit["audit_reason"] = accepted.map(
        {True: "accepted_exact_species_organ_linked_complete_colour_state", False: "rejected_manual_treatment_review"}
    )
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output_csv, index=False, lineterminator="\n")


if __name__ == "__main__":
    app()
