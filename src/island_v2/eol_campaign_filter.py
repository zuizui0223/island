"""Filter EOL TraitBank candidates to the exact global campaign species master.

The generic bulk-trait importer may be run against the broader island taxon table.
This module performs the analysis-facing gate: only species present in the frozen
105,612-species global campaign ledger are retained. Excluded rows are preserved
for audit and are never interpreted as biological absences.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str).fillna("")


def _write_gzip(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(
        path,
        index=False,
        compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
    )


def _ontology_traits(path: Path) -> list[str]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    traits = payload.get("traits") or {}
    return [str(name) for name in traits]


def filter_eol_candidates_to_campaign(
    *,
    campaign_ledger: Path,
    candidate_csv: Path,
    output_dir: Path,
    ontology_path: Path = Path("config/trait_ontology.yml"),
    source_summary: Path | None = None,
) -> dict[str, Any]:
    """Retain direct EOL candidates whose species occur in the campaign ledger."""
    for path, label in (
        (campaign_ledger, "campaign ledger"),
        (candidate_csv, "EOL candidate CSV"),
        (ontology_path, "trait ontology"),
    ):
        if not path.exists():
            raise typer.BadParameter(f"Missing {label}: {path}")

    ledger = _read_csv(campaign_ledger)
    candidates = _read_csv(candidate_csv)
    if "accepted_species" not in ledger.columns:
        raise typer.BadParameter("Campaign ledger must contain accepted_species.")
    required = {"accepted_species", "trait_name"}
    missing = sorted(required.difference(candidates.columns))
    if missing:
        raise typer.BadParameter(f"EOL candidates lack required columns: {missing}")

    campaign_species = {
        value.strip()
        for value in ledger["accepted_species"].astype(str)
        if value.strip()
    }
    if not campaign_species:
        raise typer.BadParameter("Campaign ledger contains no accepted species.")

    accepted = candidates["accepted_species"].astype(str).str.strip().isin(campaign_species)
    retained = candidates.loc[accepted].copy()
    excluded = candidates.loc[~accepted].copy()

    # Preserve input ordering while removing exact duplicate evidence rows.
    retained = retained.drop_duplicates().reset_index(drop=True)
    excluded = excluded.drop_duplicates().reset_index(drop=True)

    output_dir.mkdir(parents=True, exist_ok=True)
    retained_path = output_dir / "trait_candidates_campaign.csv.gz"
    excluded_path = output_dir / "excluded_noncampaign_candidates.csv.gz"
    coverage_path = output_dir / "coverage_by_trait_campaign.csv"
    summary_path = output_dir / "coverage_summary_campaign.json"
    _write_gzip(retained, retained_path)
    _write_gzip(excluded, excluded_path)

    coverage_rows: list[dict[str, Any]] = []
    for trait_name in _ontology_traits(ontology_path):
        selected = retained.loc[retained["trait_name"].eq(trait_name)]
        n_species = int(selected["accepted_species"].nunique()) if not selected.empty else 0
        coverage_rows.append(
            {
                "trait_name": trait_name,
                "n_candidate_rows": int(len(selected)),
                "n_species_with_direct_evidence": n_species,
                "n_campaign_species": len(campaign_species),
                "pct_campaign_species_covered": round(
                    100 * n_species / len(campaign_species), 6
                ),
            }
        )
    coverage = pd.DataFrame(coverage_rows)
    coverage.to_csv(coverage_path, index=False)

    upstream: dict[str, Any] = {}
    if source_summary and source_summary.exists():
        upstream = json.loads(source_summary.read_text(encoding="utf-8"))

    summary = {
        "version": "1.0",
        "filter_policy": (
            "Direct EOL evidence candidates are retained only when accepted_species "
            "occurs in the frozen global campaign ledger. Excluded rows remain an audit, "
            "not biological absences."
        ),
        "campaign_ledger": str(campaign_ledger),
        "candidate_csv": str(candidate_csv),
        "n_campaign_species": len(campaign_species),
        "n_input_candidate_rows": int(len(candidates)),
        "n_retained_candidate_rows": int(len(retained)),
        "n_retained_species": int(retained["accepted_species"].nunique())
        if not retained.empty
        else 0,
        "n_excluded_candidate_rows": int(len(excluded)),
        "n_excluded_species": int(excluded["accepted_species"].nunique())
        if not excluded.empty
        else 0,
        "retained_candidates_csv": str(retained_path),
        "excluded_candidates_csv": str(excluded_path),
        "coverage_by_trait_csv": str(coverage_path),
        "upstream_source_name": upstream.get("source_name", ""),
        "upstream_archive_sha256": upstream.get("archive_sha256", ""),
        "human_review_required": True,
    }
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    summary["summary_json"] = str(summary_path)
    return summary


@app.command("filter")
def filter_command(
    campaign_ledger: Path = typer.Option(..., exists=True),
    candidate_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    ontology_path: Path = typer.Option(Path("config/trait_ontology.yml"), exists=True),
    source_summary: Path | None = typer.Option(None),
) -> None:
    report = filter_eol_candidates_to_campaign(
        campaign_ledger=campaign_ledger,
        candidate_csv=candidate_csv,
        output_dir=output_dir,
        ontology_path=ontology_path,
        source_summary=source_summary,
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
