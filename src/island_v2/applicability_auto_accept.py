"""Apply human-preregistered rules to auto-accept Phase 0.5 applicability drafts.

The applicability freeze normally requires a human to flip each agent-drafted
source-region row from ``agent_drafted_pending_review`` to ``accepted`` after
checking its sources. This module keeps that accountability while removing the
per-item click: a human registers a decision rule once (config/applicability_
auto_accept_rules.yml) and owns it, and this command auto-accepts every drafted
row the rule matches, stamping ``reviewer=auto_rule:<rule_id>`` so the provenance
is transparent -- it is never recorded as if a human inspected the source.

It is outcome-blind: rules reference only the biogeographic region, documented
native Bombus status, evidence count, citation, and locator. Rows that do not
match a rule (or whose evidence is a placeholder) stay pending and are reported
as human exceptions. Accepting a region also accepts its linked evidence, so the
downstream ``bombus_applicability`` validator (accepted region -> accepted
evidence) stays satisfied.
"""

from __future__ import annotations

import json
import re
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

PENDING_STATES = {"", "pending", "agent_drafted_pending_review", "agent_drafted", "draft"}


@app.callback()
def main() -> None:
    """Auto-accept applicability drafts under human-preregistered rules."""


def load_rules(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or "source_region_rules" not in config:
        raise typer.BadParameter("auto-accept rules must contain a source_region_rules mapping")
    return config


def _is_pending(status: str) -> bool:
    return str(status).strip().lower() in PENDING_STATES


def _is_placeholder(value: str, markers: list[str]) -> bool:
    text = str(value or "").strip().lower()
    if not text:
        return True
    return any(marker.lower() in text for marker in markers)


def auto_accept(
    registry: pd.DataFrame,
    region_evidence: pd.DataFrame,
    assignments: pd.DataFrame,
    assignment_evidence: pd.DataFrame,
    rules: dict[str, Any],
    decision_date: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Return the four curation tables with rule-matched rows accepted, plus a report."""
    registry = registry.copy()
    region_evidence = region_evidence.copy()
    assignments = assignments.copy()
    assignment_evidence = assignment_evidence.copy()
    markers = list(rules.get("placeholder_locator_markers") or [])

    sr = rules["source_region_rules"]
    min_ev = int(sr.get("min_independent_evidence_records", 1))
    region_rules = sr.get("rules") or []

    accepted_regions: list[dict[str, str]] = []
    pending_regions: list[dict[str, str]] = []
    accepted_evidence_ids: set[str] = set()

    for idx, row in registry.iterrows():
        region_id = str(row["source_region_id"]).strip()
        if not _is_pending(row.get("review_status", "")):
            continue
        bio = str(row.get("biogeographic_region", "")).strip()
        native = str(row.get("native_bombus_status", "")).strip()
        applicability = str(row.get("applicability", "")).strip()
        match = next(
            (
                r
                for r in region_rules
                if re.search(str(r["match_biogeographic_region_regex"]), bio)
                and native == str(r["require_native_bombus_status"])
                and applicability == str(r["accept_applicability"])
            ),
            None,
        )
        if match is None:
            pending_regions.append({"source_region_id": region_id, "reason": "no matching rule"})
            continue
        ev = region_evidence.loc[
            region_evidence["source_region_id"].astype(str).str.strip() == region_id
        ]
        if len(ev) < min_ev:
            pending_regions.append({"source_region_id": region_id, "reason": "insufficient evidence records"})
            continue
        if sr.get("require_source_citation", True) and ev["source_citation"].astype(str).str.strip().eq("").any():
            pending_regions.append({"source_region_id": region_id, "reason": "missing source citation"})
            continue
        if sr.get("require_real_source_locator", False):
            if ev["source_locator"].apply(lambda v: _is_placeholder(v, markers)).any():
                pending_regions.append({"source_region_id": region_id, "reason": "placeholder source locator"})
                continue
        stamp = f"auto_rule:{match['rule_id']}"
        registry.loc[idx, "review_status"] = "accepted"
        registry.loc[idx, "reviewer"] = stamp
        registry.loc[idx, "decision_date"] = decision_date
        accepted_evidence_ids.update(ev["source_region_evidence_id"].astype(str).str.strip())
        accepted_regions.append({"source_region_id": region_id, "rule_id": match["rule_id"]})

    # Accept the linked source-region evidence so the applicability validator
    # (accepted region -> accepted evidence) is satisfied.
    for idx, row in region_evidence.iterrows():
        if str(row["source_region_evidence_id"]).strip() in accepted_evidence_ids and _is_pending(
            row.get("review_status", "")
        ):
            region_evidence.loc[idx, "review_status"] = "accepted"
            region_evidence.loc[idx, "reviewer"] = "auto_rule:source_region"
            region_evidence.loc[idx, "decision_date"] = decision_date

    accepted_region_ids = {r["source_region_id"] for r in accepted_regions} | set(
        registry.loc[registry["review_status"].astype(str).str.strip().str.lower() == "accepted", "source_region_id"]
        .astype(str)
        .str.strip()
    )

    ar = rules.get("assignment_rules") or {}
    accepted_assignments: list[dict[str, str]] = []
    pending_assignments: list[dict[str, str]] = []
    accepted_assignment_ev_ids: set[str] = set()
    if ar.get("auto_accept", False):
        for idx, row in assignments.iterrows():
            island_id = str(row["island_id"]).strip()
            if not _is_pending(row.get("review_status", "")):
                continue
            region_id = str(row.get("source_region_id", "")).strip()
            if ar.get("require_source_region_accepted", True) and region_id not in accepted_region_ids:
                pending_assignments.append({"island_id": island_id, "reason": "source region not accepted"})
                continue
            ev_id = str(row.get("assignment_evidence_id", "")).strip()
            ev = assignment_evidence.loc[
                assignment_evidence["assignment_evidence_id"].astype(str).str.strip() == ev_id
            ]
            if ev.empty:
                pending_assignments.append({"island_id": island_id, "reason": "assignment evidence missing"})
                continue
            if ar.get("require_assignment_citation", True) and ev["source_citation"].astype(str).str.strip().eq("").any():
                pending_assignments.append({"island_id": island_id, "reason": "missing assignment citation"})
                continue
            if ar.get("require_real_assignment_locator", True):
                if ev["source_locator"].apply(lambda v: _is_placeholder(v, markers)).any():
                    pending_assignments.append({"island_id": island_id, "reason": "placeholder assignment locator"})
                    continue
            assignments.loc[idx, "review_status"] = "accepted"
            assignments.loc[idx, "reviewer"] = "auto_rule:assignment"
            assignments.loc[idx, "decision_date"] = decision_date
            accepted_assignment_ev_ids.add(ev_id)
            accepted_assignments.append({"island_id": island_id, "source_region_id": region_id})

    for idx, row in assignment_evidence.iterrows():
        if str(row["assignment_evidence_id"]).strip() in accepted_assignment_ev_ids and _is_pending(
            row.get("review_status", "")
        ):
            assignment_evidence.loc[idx, "review_status"] = "accepted"
            assignment_evidence.loc[idx, "reviewer"] = "auto_rule:assignment"
            assignment_evidence.loc[idx, "decision_date"] = decision_date

    report = {
        "decision_date": decision_date,
        "registered_by": rules.get("registered_by", ""),
        "n_regions_auto_accepted": len(accepted_regions),
        "regions_auto_accepted": accepted_regions,
        "regions_left_pending": pending_regions,
        "n_assignments_auto_accepted": len(accepted_assignments),
        "assignments_auto_accepted": accepted_assignments,
        "assignments_left_pending": pending_assignments,
        "note": (
            "Auto-accepted rows are stamped reviewer=auto_rule:<rule_id>; the "
            "registered rule owner is accountable. Nothing here decides trait "
            "values, establishment, or analysis inclusion."
        ),
    }
    return registry, region_evidence, assignments, assignment_evidence, report


@app.command("apply")
def apply(
    rules_path: Path = typer.Option(Path("config/applicability_auto_accept_rules.yml"), exists=True),
    registry_csv: Path = typer.Option(..., exists=True),
    region_evidence_csv: Path = typer.Option(..., exists=True),
    assignments_csv: Path = typer.Option(..., exists=True),
    assignment_evidence_csv: Path = typer.Option(..., exists=True),
    decision_date: str = typer.Option(date.today().isoformat()),
    write: bool = typer.Option(False, help="Write accepted rows back to the CSVs (default: dry-run)."),
    report_path: Path = typer.Option(None, help="Optional path to write the JSON report."),
) -> None:
    """Auto-accept rule-matched applicability drafts; dry-run unless --write."""
    rules = load_rules(rules_path)
    registry = pd.read_csv(registry_csv, dtype=str).fillna("")
    region_evidence = pd.read_csv(region_evidence_csv, dtype=str).fillna("")
    assignments = pd.read_csv(assignments_csv, dtype=str).fillna("")
    assignment_evidence = pd.read_csv(assignment_evidence_csv, dtype=str).fillna("")

    registry, region_evidence, assignments, assignment_evidence, report = auto_accept(
        registry, region_evidence, assignments, assignment_evidence, rules, decision_date
    )
    if write:
        registry.to_csv(registry_csv, index=False)
        region_evidence.to_csv(region_evidence_csv, index=False)
        assignments.to_csv(assignments_csv, index=False)
        assignment_evidence.to_csv(assignment_evidence_csv, index=False)
    if report_path:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    mode = "WROTE" if write else "DRY-RUN"
    typer.echo(
        f"[{mode}] auto-accepted {report['n_regions_auto_accepted']} source region(s), "
        f"{report['n_assignments_auto_accepted']} assignment(s); "
        f"{len(report['regions_left_pending'])} region(s) and "
        f"{len(report['assignments_left_pending'])} assignment(s) left pending."
    )


if __name__ == "__main__":
    app()
