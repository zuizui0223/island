"""Recover high-precision whole-flower size from acquired FlorML treatments."""

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
from island_v2.official_flora_colour_recovery_checkpoint import _text

app = typer.Typer(add_completion=False, no_args_is_help=True)

TRAIT = "flower_size_class"
AUDIT_SIZE = 200
AUDIT_SEED = "florml-size-strict-v2-audit"
TARGET_PROVIDERS = frozenset({"flora_malesiana", "flora_of_the_guianas"})
STRICT_NAME_METHODS = frozenset(
    {
        "wfo_june_2026_exact_accepted_xml_species_family_consistent",
        "wfo_june_2026_exact_synonym_xml_species_family_consistent",
        "wfo_june_2026_exact_accepted_or_synonym_xml_species_family_consistent",
    }
)
MEASUREMENT = (
    r"(?:c\.|ca\.|about|up to|at least|to)?\s*"
    r"\(?\d+(?:\.\d+)?(?:\s*[-–]\s*\d+(?:\.\d+)?)?"
    r"(?:\s*\([^)]*\))?\s*(?:mm|cm)\s*"
    r"(?:long|wide|diam(?:\.|eter)?|across|deep|in all)"
)
WHOLE_FLOWER_MEASUREMENT = re.compile(
    r"\b(flowers?(?![-‑–\s]*bearing)|corolla)\b"
    rf"([^.;]{{0,65}}?)({MEASUREMENT})",
    re.IGNORECASE,
)
OTHER_ORGAN = re.compile(
    r"\b(?:flowering\s+stems?|stems?|axes?|axis|peduncles?|pedicels?|"
    r"pedicellate|hairs?|receptacles?|bracts?|bracteoles?|sepals?|petals?|"
    r"perianths?|perigone|glumes?|bristles?|stamens?|anthers?|lobes?|tubes?|"
    r"disks?|bands?|belts?|buds?|inflorescences?|branches?|rays?|umbels?|"
    r"spikes?|cymes?|racemes?|panicles?|thyrses?|heads?|mouths?|throats?|"
    r"racemose|spicate|clustered|surrounded|subtended)\b",
    re.IGNORECASE,
)
SEX_OR_INCOMPLETE_STAGE = re.compile(
    r"\b(?:male|female|staminate|pistillate|young|immature|not quite mature|"
    r"mature bud|in bud|prior to opening)\b|[\u2640\u2642]",
    re.IGNORECASE,
)


def _id(*parts: object, length: int = 64) -> str:
    payload = "\n".join(_text(part) for part in parts)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:length]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strict_size_candidate(record: dict[str, Any]) -> tuple[bool, str]:
    """Require an explicit whole-flower/corolla measurement and reject parts."""

    quote = " ".join(_text(record.get("source_excerpt")).split())
    if _text(record.get("source_provider")) not in TARGET_PROVIDERS:
        return False, "provider_outside_florml_checkpoint"
    if canonical_trait_name(_text(record.get("trait_name"))) != TRAIT:
        return False, "trait_outside_checkpoint"
    if _text(record.get("evidence_scope")).casefold() != "species_direct":
        return False, "not_species_direct"
    if _text(record.get("name_match_method")) not in STRICT_NAME_METHODS:
        return False, "name_match_not_strict_xml_exact"
    if not _text(record.get("accepted_species")) or not quote:
        return False, "missing_identity_or_quote"
    if SEX_OR_INCOMPLETE_STAGE.search(quote):
        return False, "sex_specific_or_incomplete_stage_measurement"
    if re.search(r"\d+h-\d", quote, re.IGNORECASE):
        return False, "ambiguous_ocr_measurement"

    for match in WHOLE_FLOWER_MEASUREMENT.finditer(quote):
        start = max(
            quote.rfind(";", 0, match.start()),
            quote.rfind(". ", 0, match.start()),
            0,
        )
        preceding_clause = quote[start : match.start()]
        following_clause = quote[match.end() : match.end() + 65]
        following_clause = following_clause.split(";", 1)[0].split(". ", 1)[0]
        if any(
            OTHER_ORGAN.search(value)
            for value in (preceding_clause, match.group(2), following_clause)
        ):
            continue
        if re.search(r"\b(?:on|forming|with)\b", match.group(2), re.IGNORECASE):
            continue
        return True, "selected_explicit_whole_flower_or_corolla_measurement"
    return False, "no_unambiguous_whole_flower_measurement"


def _completed_pairs(frame: pd.DataFrame) -> set[tuple[str, str]]:
    work = frame.copy().fillna("")
    if "resolution_status" in work:
        work = work.loc[work["resolution_status"].eq("resolved")]
    work["trait_name"] = work["trait_name"].map(canonical_trait_name)
    return set(zip(work["accepted_species"], work["trait_name"], strict=True))


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
        selected, reason = strict_size_candidate(record)
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
        row.update(
            {
                "quality": "high",
                "acceptance_contract": "official_florml_explicit_whole_flower_size_v1",
                "content_fingerprint": _id(
                    species, " ".join(row["source_excerpt"].casefold().split())
                ),
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
        raise ValueError("FlorML size recovery produced empty or duplicate evidence")
    selection = pd.DataFrame(selection_rows).sort_values(
        ["selected", "selection_reason", "accepted_species"], kind="stable"
    )
    return evidence, _audit_template(evidence), selection.reset_index(drop=True)


def _audit_template(evidence: pd.DataFrame) -> pd.DataFrame:
    sample = evidence.copy()
    sample["audit_hash"] = sample["source_record_id"].map(
        lambda value: _id(f"{AUDIT_SEED}|{value}")
    )
    guianas = sample.loc[sample["source_provider"].eq("flora_of_the_guianas")]
    guianas = guianas.sort_values("audit_hash", kind="stable").head(20)
    malesiana = sample.loc[sample["source_provider"].eq("flora_malesiana")]
    malesiana = malesiana.sort_values("audit_hash", kind="stable").head(
        AUDIT_SIZE - len(guianas)
    )
    sample = pd.concat([guianas, malesiana], ignore_index=True).sort_values(
        "audit_hash", kind="stable"
    )
    if len(sample) != AUDIT_SIZE:
        raise ValueError(f"FlorML size audit requires {AUDIT_SIZE} treatments")
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
    evidence_path = output_dir / "florml_size_evidence_20260812.csv.gz"
    audit_path = output_dir / "florml_size_audit_200_20260812.csv"
    selection_path = output_dir / "florml_size_selection_20260812.csv.gz"
    evidence.to_csv(evidence_path, index=False, compression={"method": "gzip", "mtime": 0})
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    selection.to_csv(
        selection_path, index=False, compression={"method": "gzip", "mtime": 0}
    )
    manifest = {
        "contract": "official_florml_whole_flower_size_recovery_v1",
        "candidate_evidence_rows": len(evidence),
        "species_trait": len(evidence[["accepted_species", "trait_name"]].drop_duplicates()),
        "audit_rows": len(audit),
        "audit_seed": AUDIT_SEED,
        "providers": evidence["source_provider"].value_counts().sort_index().to_dict(),
        "inputs": {
            str(candidates_csv): _sha256(candidates_csv),
            str(completed_direct_csv): _sha256(completed_direct_csv),
        },
        "guardrails": {
            "family_inference": False,
            "global_fallback": False,
            "axis_only_join": False,
            "trait_substitution": False,
            "whole_flower_measurement_only": True,
        },
    }
    (output_dir / "florml_size_manifest_20260812.json").write_text(
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
        {
            True: "accepted_explicit_whole_flower_or_corolla_measurement",
            False: "rejected_ambiguous_question_mark_flower_measurement",
        }
    )
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output_csv, index=False, lineterminator="\n")


if __name__ == "__main__":
    app()
