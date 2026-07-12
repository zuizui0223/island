"""Consolidate global machine trait candidates into auditable snapshot artifacts."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

REPRODUCTIVE_TRAITS = {
    "pollen_vector_mode",
    "self_incompatibility",
    "autonomous_selfing_capacity",
    "mating_system",
    "cleistogamy",
}
FLORAL_ACCESS_TRAITS = {
    "floral_symmetry",
    "floral_form",
    "tube_depth_class",
    "flower_size_class",
    "flower_primary_color",
}
ALTERNATIVE_GUILD_TRAITS = {"pollination_functional_guild"}


def _read_candidate_waves(campaign_dir: Path) -> pd.DataFrame:
    paths = sorted((campaign_dir / "waves").glob("wave_*/machine_candidates.csv.gz"))
    if not paths:
        return pd.DataFrame()
    frames = []
    for path in paths:
        frame = pd.read_csv(path, dtype=str).fillna("")
        frame["wave_id"] = path.parent.name
        frames.append(frame)
    return pd.concat(frames, ignore_index=True, sort=False).drop_duplicates()


def build_snapshot(campaign_dir: Path, output_dir: Path) -> dict[str, object]:
    status_path = campaign_dir / "campaign_status.json"
    ledger_path = campaign_dir / "campaign_ledger.csv.gz"
    if not status_path.exists() or not ledger_path.exists():
        raise typer.BadParameter("campaign_status.json and campaign_ledger.csv.gz are required")

    status = json.loads(status_path.read_text(encoding="utf-8"))
    ledger = pd.read_csv(ledger_path, dtype=str).fillna("")
    candidates = _read_candidate_waves(campaign_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if candidates.empty:
        candidates = pd.DataFrame(
            columns=[
                "accepted_species",
                "trait_name",
                "candidate_value",
                "source_url",
                "source_citation",
                "source_excerpt",
                "source_type",
                "evidence_scope",
                "matched_term",
                "confidence",
                "campaign_task",
                "campaign_phase",
                "target_for_task",
                "wave_id",
            ]
        )
    candidates["review_status"] = "machine_candidate_unreviewed"
    candidates["analysis_eligible"] = False
    candidates["zero_hit_is_absence"] = False

    target_mask = candidates.get("target_for_task", "").astype(str).str.lower().isin(
        {"true", "1", "yes"}
    )
    target = candidates.loc[target_mask].copy()
    target.to_csv(output_dir / "global_trait_machine_candidates.csv.gz", index=False, compression="gzip")

    reproductive = target.loc[target["trait_name"].isin(REPRODUCTIVE_TRAITS)].copy()
    floral = target.loc[target["trait_name"].isin(FLORAL_ACCESS_TRAITS)].copy()
    guild = target.loc[target["trait_name"].isin(ALTERNATIVE_GUILD_TRAITS)].copy()
    reproductive.to_csv(output_dir / "reproductive_and_pollen_vector_candidates.csv.gz", index=False, compression="gzip")
    floral.to_csv(output_dir / "floral_access_candidates.csv.gz", index=False, compression="gzip")
    guild.to_csv(output_dir / "alternative_pollinator_guild_candidates.csv.gz", index=False, compression="gzip")

    trait_counts = (
        target.groupby("trait_name", dropna=False)
        .agg(n_candidate_rows=("accepted_species", "size"), n_species=("accepted_species", "nunique"))
        .reset_index()
        .sort_values("trait_name")
    )
    trait_counts.to_csv(output_dir / "trait_candidate_coverage.csv", index=False)

    ledger.to_csv(output_dir / "campaign_ledger_snapshot.csv.gz", index=False, compression="gzip")
    summary = {
        "scope": "machine_candidate_snapshot_not_curated_trait_values",
        "n_global_species": int(status.get("n_global_species", len(ledger))),
        "n_candidate_rows_all": int(len(candidates)),
        "n_target_candidate_rows": int(len(target)),
        "n_candidate_species": int(target["accepted_species"].nunique()) if not target.empty else 0,
        "n_reproductive_rows": int(len(reproductive)),
        "n_floral_access_rows": int(len(floral)),
        "n_alternative_guild_rows": int(len(guild)),
        "active_task": status.get("active_task", ""),
        "tasks": status.get("tasks", {}),
        "interpretation": (
            "These rows are source-backed machine candidates only. They are not accepted trait "
            "values and are not eligible for confirmatory analysis until review/adjudication."
        ),
    }
    (output_dir / "trait_artifact_snapshot_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return summary


@app.command("build")
def build(
    campaign_dir: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    typer.echo(json.dumps(build_snapshot(campaign_dir, output_dir), ensure_ascii=False))


if __name__ == "__main__":
    app()
