"""Bounded retrieval of explicit flower-visit and pollination interaction claims.

The output is an unreviewed evidence table. It does not infer effective
pollination, functional guilds, or biological absence from missing records.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import time
from datetime import UTC, datetime
from pathlib import Path

import httpx
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.callback()
def main() -> None:
    """Collect bounded explicit flower-visit and pollination interaction claims."""


API_URL = "https://api.globalbioticinteractions.org/interaction.csv"
MAX_TAXA = 100
RELATIONS = {
    "pollinates": "pollination_claim",
    "pollinated by": "pollination_claim",
    "visits flowers of": "flower_visit_claim",
    "flowers visited by": "flower_visit_claim",
}


def _text(value: object) -> str:
    return " ".join(str(value or "").split())


def _get(row: dict[str, str], *keys: str) -> str:
    for key in keys:
        if _text(row.get(key, "")):
            return _text(row[key])
    return ""


def _query(client: httpx.Client, taxon: str, orientation: str) -> list[dict[str, str]]:
    parameter = "sourceTaxon" if orientation == "source" else "targetTaxon"
    response = client.get(API_URL, params={parameter: taxon, "includeObservations": "true"})
    response.raise_for_status()
    return list(csv.DictReader(io.StringIO(response.text)))


def _record(taxon: str, orientation: str, row: dict[str, str]) -> dict[str, str] | None:
    relation = _get(row, "interaction_type", "interactionTypeName").casefold()
    claim = RELATIONS.get(relation)
    if claim is None:
        return None
    source = _get(row, "source_taxon_name", "sourceTaxonName")
    target = _get(row, "target_taxon_name", "targetTaxonName")
    focal = _text(taxon).casefold()
    source_match = source.casefold() == focal
    target_match = target.casefold() == focal
    partner = target if source_match and not target_match else source if target_match and not source_match else ""
    identity = "|".join([taxon, orientation, relation, source, target, _get(row, "study_source_id", "studySourceId"), _get(row, "event_date", "eventDate")])
    return {
        "evidence_id": "globi:" + hashlib.sha256(identity.encode("utf-8")).hexdigest()[:20],
        "query_taxon": taxon,
        "query_orientation": orientation,
        "query_taxon_exact_in_record": str(source_match or target_match).lower(),
        "interaction_type": relation,
        "interaction_claim_class": claim,
        "source_taxon_name": source,
        "target_taxon_name": target,
        "partner_taxon_name": partner,
        "study_title": _get(row, "study_title", "studyTitle"),
        "study_url": _get(row, "study_url", "studyUrl"),
        "study_doi": _get(row, "study_doi", "studyDoi"),
        "study_citation": _get(row, "study_citation", "studyCitation"),
        "study_source_citation": _get(row, "study_source_citation", "studySourceCitation"),
        "study_source_id": _get(row, "study_source_id", "studySourceId"),
        "locality": _get(row, "locality"),
        "latitude": _get(row, "latitude", "decimalLatitude"),
        "longitude": _get(row, "longitude", "decimalLongitude"),
        "event_date": _get(row, "event_date", "eventDate"),
        "effectiveness_class": "not_assessed",
        "review_status": "unreviewed_source_backed",
        "decision_boundary": "Explicit interaction claim only; no automatic guild, effectiveness, absence, or inclusion decision.",
    }


def collect(taxa_csv: Path, output_dir: Path, max_taxa: int = MAX_TAXA, pause_seconds: float = 0.2) -> dict[str, object]:
    if not 1 <= max_taxa <= MAX_TAXA:
        raise typer.BadParameter(f"max_taxa must be 1..{MAX_TAXA}")
    taxa = pd.read_csv(taxa_csv, dtype=str).fillna("")
    if "accepted_species" not in taxa.columns:
        raise typer.BadParameter("taxa CSV must contain accepted_species")
    species = taxa["accepted_species"].map(_text)
    selected = species.loc[species.ne("")].drop_duplicates().head(max_taxa).tolist()
    records: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    with httpx.Client(timeout=45.0, follow_redirects=True, headers={"User-Agent": "island-floral-v2/0.1"}) as client:
        for taxon in selected:
            for orientation in ("source", "target"):
                try:
                    for row in _query(client, taxon, orientation):
                        record = _record(taxon, orientation, row)
                        if record is not None:
                            records.append(record)
                except Exception as exc:  # noqa: BLE001
                    errors.append({"query_taxon": taxon, "query_orientation": orientation, "error": str(exc)})
                if pause_seconds > 0:
                    time.sleep(pause_seconds)
    columns = [
        "evidence_id", "query_taxon", "query_orientation", "query_taxon_exact_in_record", "interaction_type", "interaction_claim_class", "source_taxon_name", "target_taxon_name", "partner_taxon_name", "study_title", "study_url", "study_doi", "study_citation", "study_source_citation", "study_source_id", "locality", "latitude", "longitude", "event_date", "effectiveness_class", "review_status", "decision_boundary",
    ]
    evidence = pd.DataFrame(records, columns=columns).drop_duplicates(subset=["evidence_id"])
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(output_dir / "globi_interaction_evidence.csv", index=False)
    pd.DataFrame(errors, columns=["query_taxon", "query_orientation", "error"]).to_csv(output_dir / "globi_query_errors.csv", index=False)
    summary = {
        "source": "Global Biotic Interactions (GloBI)",
        "retrieved_at_utc": datetime.now(UTC).isoformat(),
        "n_input_taxa": len(selected),
        "n_query_attempts": len(selected) * 2,
        "n_query_errors": len(errors),
        "n_explicit_records": int(len(evidence)),
        "evidence_status": "unreviewed_source_backed",
        "non_inference_rule": "Missing records are not absence; visitor records are not effectiveness; no guild is inferred automatically.",
    }
    (output_dir / "globi_interaction_evidence_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


@app.command("collect")
def collect_command(
    taxa_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_taxa: int = typer.Option(MAX_TAXA, min=1, max=MAX_TAXA),
    pause_seconds: float = typer.Option(0.2, min=0.0),
) -> None:
    typer.echo(json.dumps(collect(taxa_csv, output_dir, max_taxa, pause_seconds), ensure_ascii=False))


if __name__ == "__main__":
    app()
