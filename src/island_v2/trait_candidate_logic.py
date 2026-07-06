"""Logic audit for broad, reviewable trait-extraction candidates.

The global extraction stage optimizes coverage. A weak but explicitly labelled
source is therefore still useful. A biologically invalid shortcut is different:
for example, coding autonomous selfing from self-compatibility alone, or calling
a genus-level claim species-direct. This module quarantines those logical errors
without deleting the raw candidate table or filtering by source reliability.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

_REQUIRED_COLUMNS = {
    "trait_name",
    "standardized_value",
    "candidate_kind",
    "evidence_scope",
    "source_type",
    "source_reliability",
    "source_citation",
    "source_url",
    "evidence_excerpt",
    "supporting_taxa",
    "inference_rule",
    "confidence",
}

_SOURCE_BACKED_SCOPES = {"species_direct", "species_indirect"}
_INFERENCE_SCOPES = {"genus_inference", "family_inference"}
_AUTONOMY_TERMS = re.compile(
    r"autonom|autogam|automatic self|spontaneous self|without pollinator|unassisted self",
    flags=re.IGNORECASE,
)
_COMPATIBILITY_TERMS = re.compile(r"self[ -]?compatib|\bSC\b", flags=re.IGNORECASE)
_POLLINATION_TERMS = re.compile(
    r"pollinat|visitor|forag|bee|bumble|fly|bird|moth|butterfl|bat|wind|water|entomophil|ornithophil",
    flags=re.IGNORECASE,
)
_MORPHOLOGY_TERMS = re.compile(
    r"colour|color|corolla|tubular|tube|symmetr|flower shape|red|pink|white|yellow|blue|purple",
    flags=re.IGNORECASE,
)


def _text(row: pd.Series) -> str:
    return " ".join(
        str(row.get(column, "") or "")
        for column in ("raw_value", "evidence_excerpt", "inference_rule", "inference_note")
    )


def _parsed_supporting_taxa(value: Any) -> list[str]:
    if isinstance(value, list):
        return [str(item).strip() for item in value if str(item).strip()]
    if not value or str(value).strip() in {"", "[]", "null", "None"}:
        return []
    try:
        parsed = json.loads(str(value))
    except json.JSONDecodeError:
        return [str(value).strip()]
    if isinstance(parsed, list):
        return [str(item).strip() for item in parsed if str(item).strip()]
    return [str(value).strip()]


def audit_candidate(row: pd.Series, ontology: dict[str, Any]) -> tuple[list[str], list[str]]:
    """Return (hard_errors, review_warnings) for one candidate row.

    Hard errors quarantine the row from coverage-ready candidates. Warnings retain
    the row because low confidence and incomplete evidence are legitimate first-pass
    results in the broad acquisition track.
    """
    errors: list[str] = []
    warnings: list[str] = []
    trait_name = str(row.get("trait_name", "") or "")
    value = str(row.get("standardized_value", "") or "")
    kind = str(row.get("candidate_kind", "") or "")
    scope = str(row.get("evidence_scope", "") or "")
    source_type = str(row.get("source_type", "") or "")
    source_reliability = str(row.get("source_reliability", "") or "")
    citation = str(row.get("source_citation", "") or "").strip()
    confidence = str(row.get("confidence", "") or "")
    source_url = str(row.get("source_url", "") or "").strip()
    text = _text(row)

    traits = ontology.get("traits") or {}
    if trait_name not in traits:
        errors.append("unknown_trait_name")
    else:
        allowed = set(traits[trait_name].get("allowed_values") or [])
        if value not in allowed:
            errors.append("standardized_value_outside_ontology")

    if kind == "source_backed":
        if scope not in _SOURCE_BACKED_SCOPES:
            errors.append("source_backed_candidate_has_non_species_scope")
        if source_type == "none_found" or source_reliability == "none" or citation in {"", "none_found"}:
            errors.append("source_backed_candidate_missing_usable_source")
        if not source_url:
            warnings.append("source_backed_candidate_without_stable_url")
    elif kind == "hierarchical_inference":
        if scope not in _INFERENCE_SCOPES:
            errors.append("hierarchical_inference_missing_genus_or_family_scope")
        if not _parsed_supporting_taxa(row.get("supporting_taxa")):
            errors.append("hierarchical_inference_missing_supporting_taxa")
        if not str(row.get("inference_rule", "") or "").strip():
            errors.append("hierarchical_inference_missing_transfer_rule")
        if confidence not in {"medium", "low"}:
            errors.append("hierarchical_inference_has_invalid_confidence")
    elif kind == "unresolved":
        if scope != "unresolved" or value != "unresolved":
            errors.append("unresolved_candidate_not_coded_as_unresolved")
        if source_type != "none_found" or source_reliability != "none" or confidence != "unresolved":
            errors.append("unresolved_candidate_has_nonmissing_evidence_metadata")
    else:
        errors.append("unknown_candidate_kind")

    if trait_name == "autonomous_selfing_capacity" and value != "unresolved":
        if _COMPATIBILITY_TERMS.search(text) and not _AUTONOMY_TERMS.search(text):
            errors.append("autonomous_selfing_inferred_from_compatibility_only")

    if trait_name == "mating_system" and value != "unresolved":
        if _COMPATIBILITY_TERMS.search(text) and not re.search(
            r"outcross|selfing|selfed|mixed mating|mating system|outcrossing", text, re.IGNORECASE
        ):
            warnings.append("mating_system_may_be_based_on_compatibility_only")

    if trait_name == "pollination_functional_guild" and value != "unresolved":
        has_pollination = bool(_POLLINATION_TERMS.search(text))
        has_morphology = bool(_MORPHOLOGY_TERMS.search(text))
        if not has_pollination:
            errors.append("pollination_guild_without_pollination_evidence")
        elif has_morphology and not re.search(r"pollinat|visitor|forag|entomophil|ornithophil", text, re.IGNORECASE):
            errors.append("pollination_guild_inferred_from_flower_morphology")

    return errors, warnings


def audit_table(table: pd.DataFrame, ontology: dict[str, Any]) -> pd.DataFrame:
    missing = _REQUIRED_COLUMNS.difference(table.columns)
    if missing:
        raise typer.BadParameter(f"Candidate table missing required columns: {sorted(missing)}")
    audited = table.copy()
    findings = [audit_candidate(row, ontology) for _, row in audited.iterrows()]
    audited["logic_errors_json"] = [json.dumps(errors, ensure_ascii=False) for errors, _ in findings]
    audited["logic_warnings_json"] = [json.dumps(warnings, ensure_ascii=False) for _, warnings in findings]
    audited["logic_status"] = ["quarantine" if errors else "pass" for errors, _ in findings]
    return audited


@app.command("audit")
def audit(
    candidates_csv: Path = typer.Option(..., exists=True),
    ontology_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    """Audit broad candidates; quarantine logical errors but retain all raw rows."""
    table = pd.read_csv(candidates_csv, dtype=str).fillna("")
    ontology = yaml.safe_load(ontology_path.read_text(encoding="utf-8"))
    audited = audit_table(table, ontology)
    output_dir.mkdir(parents=True, exist_ok=True)
    audited.to_csv(output_dir / "candidate_logic_audit.csv", index=False)
    audited.loc[audited["logic_status"].eq("pass")].to_csv(
        output_dir / "candidate_logic_pass.csv", index=False
    )
    audited.loc[audited["logic_status"].eq("quarantine")].to_csv(
        output_dir / "candidate_logic_quarantine.csv", index=False
    )
    report = {
        "n_input_rows": int(len(audited)),
        "n_logic_pass": int(audited["logic_status"].eq("pass").sum()),
        "n_quarantined": int(audited["logic_status"].eq("quarantine").sum()),
        "n_low_confidence_pass": int(
            audited.loc[audited["logic_status"].eq("pass"), "confidence"].eq("low").sum()
        ),
        "n_source_reliability_D_pass": int(
            audited.loc[audited["logic_status"].eq("pass"), "source_reliability"].eq("D_unvetted_web").sum()
        ),
        "policy": "Low-confidence and D-tier candidates are retained; only logical contradictions are quarantined.",
    }
    (output_dir / "candidate_logic_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
