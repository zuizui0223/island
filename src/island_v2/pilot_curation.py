"""Build island-specific curation worklists from taxon-intake artifacts.

The raw GBIF collector correctly retains occurrence candidates, but occurrence
presence is not sufficient for either an island flora denominator or floral-trait
analysis. This module combines three explicit inputs:

* the frozen island manifest;
* non-curated taxonomy and establishment queues; and
* a small reviewed registry that identifies the first pilot islands and names
  their best available flora anchor or source gap.

It never accepts a taxon, infers island nativeness, or releases a taxon to the
confirmatory angiosperm universe. It produces a curation worklist and a separate
*trait-candidate staging list* whose rows remain blocked until establishment and
taxonomic review have been accepted.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_REGISTRY_COLUMNS = {
    "island_id",
    "canonical_island_name",
    "country",
    "country_iso3",
    "source_region_id",
    "flora_anchor_status",
    "flora_anchor_citation",
    "flora_anchor_url",
    "flora_anchor_scope",
    "curation_route",
    "identity_review_status",
}
REQUIRED_TAXON_COLUMNS = {
    "reported_name",
    "reported_genus",
    "reported_family",
    "n_islands",
    "n_records",
    "name_rank_triage",
    "provisional_accepted_species",
    "provisional_taxonomic_group",
    "provisional_match_status",
    "taxon_scope_review_status",
    "establishment_review_status",
    "review_priority",
    "next_review_action",
    "gbif_match_url",
}
REQUIRED_ESTABLISHMENT_COLUMNS = {
    "island_id",
    "species",
    "provisional_accepted_species",
    "provisional_taxonomic_group",
    "provisional_match_status",
    "n_records",
    "n_unique_gbif_ids",
    "basis_of_record_set",
    "island_establishment_status",
    "review_status",
}

WORKLIST_COLUMNS = [
    "island_id",
    "canonical_island_name",
    "country",
    "country_iso3",
    "source_region_id",
    "area_km2",
    "reported_species",
    "provisional_accepted_species",
    "reported_genus",
    "reported_family",
    "n_records",
    "n_unique_gbif_ids",
    "basis_of_record_set",
    "name_rank_triage",
    "provisional_taxonomic_group",
    "provisional_match_status",
    "taxon_scope_review_status",
    "island_establishment_status",
    "establishment_review_status",
    "flora_anchor_status",
    "flora_anchor_citation",
    "flora_anchor_url",
    "flora_anchor_scope",
    "curation_route",
    "review_stage",
    "review_priority",
    "selection_rationale",
    "required_evidence",
    "gbif_match_url",
]

TRAIT_STAGING_COLUMNS = [
    "accepted_species",
    "genus",
    "family",
    "island_id",
    "canonical_island_name",
    "country_iso3",
    "source_region_id",
    "n_records",
    "flora_anchor_status",
    "trait_candidate_status",
    "release_gate",
    "selection_rationale",
]


@app.callback()
def main() -> None:
    """Build review worklists for initial island-flora curation pilots."""


def _text(values: pd.Series) -> pd.Series:
    return values.fillna("").astype(str).str.strip()


def _require(table: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = columns.difference(table.columns)
    if missing:
        raise typer.BadParameter(f"{label} missing columns: {', '.join(sorted(missing))}")


def _unique_nonblank(table: pd.DataFrame, column: str, label: str) -> None:
    values = _text(table[column])
    if values.eq("").any():
        raise typer.BadParameter(f"{label} has blank {column} values")
    duplicated = values.loc[values.duplicated()].unique().tolist()
    if duplicated:
        raise typer.BadParameter(f"{label} has duplicate {column} values: {duplicated[:5]}")


def _priority_rank(value: str) -> int:
    return {"high": 0, "medium": 1, "normal": 2}.get(str(value).strip().lower(), 3)


def classify_review_stage(row: dict[str, Any]) -> tuple[str, str, str]:
    """Return stage, rationale, and minimum evidence without making decisions."""
    rank = str(row.get("name_rank_triage", ""))
    group = str(row.get("provisional_taxonomic_group", ""))
    match = str(row.get("provisional_match_status", ""))
    anchor = str(row.get("flora_anchor_status", ""))
    if rank not in {"binomial_candidate", "infra_or_hybrid_candidate"}:
        return (
            "resolve_taxon_name",
            "Name is not a resolved species-level candidate; it cannot enter flora or trait work.",
            "Accepted species-level name from a taxonomic authority; record the source and reviewer.",
        )
    if group in {"non_seed_vascular", "gymnosperm"}:
        return (
            "confirm_non_angiosperm_scope",
            "Candidate is outside the floral/pollen-vector denominator once its broad group is reviewed.",
            "Accepted broad-group decision with a taxonomic source; retain in raw audit but exclude from angiosperm denominator.",
        )
    if group != "angiosperm" or match != "exact_backbone_candidate":
        return (
            "resolve_taxonomy_before_membership",
            "Backbone evidence is not an exact species-level angiosperm candidate.",
            "Accepted taxonomic name and broad group, with source and reviewer.",
        )
    if anchor == "source_gap":
        return (
            "resolve_island_flora_source_gap",
            "No island-level flora anchor is registered, so native/introduced/cultivated status cannot be accepted yet.",
            "Locate and cite a Robben Island flora, vegetation survey, herbarium checklist, or other island-specific authoritative source.",
        )
    return (
        "review_island_establishment",
        "Species-level angiosperm candidate has an island flora anchor but remains uncurated.",
        "Island or regional flora evidence that states presence and, where available, native/introduced/cultivated status; record page/table and reviewer.",
    )


def build_worklist(
    island_manifest: pd.DataFrame,
    taxon_intake: pd.DataFrame,
    establishment_queue: pd.DataFrame,
    registry: pd.DataFrame,
    restrict_to_registry: bool = False,
) -> pd.DataFrame:
    """Merge review artifacts with pilot-island context and route each row.

    By default the pilot registry must cover every island in the establishment
    queue (fail-closed: never silently drop an island you meant to curate). When
    ``restrict_to_registry`` is set, the establishment queue is first scoped to
    the registry's islands, so a globally growing collection can be pilot-curated
    for just the declared pilot islands without tripping that guard.
    """
    _require(island_manifest, {"island_id", "area_km2"}, "island manifest")
    _require(taxon_intake, REQUIRED_TAXON_COLUMNS, "taxon intake")
    _require(establishment_queue, REQUIRED_ESTABLISHMENT_COLUMNS, "establishment queue")
    _require(registry, REQUIRED_REGISTRY_COLUMNS, "pilot island registry")
    _unique_nonblank(registry, "island_id", "pilot island registry")

    manifest = island_manifest.copy()
    manifest["island_id"] = _text(manifest["island_id"])
    manifest["area_km2"] = pd.to_numeric(manifest["area_km2"], errors="coerce")
    registry = registry.copy()
    for column in registry.columns:
        registry[column] = _text(registry[column])
    taxon = taxon_intake.copy()
    for column in taxon.columns:
        taxon[column] = _text(taxon[column])
    taxon = taxon.drop_duplicates("reported_name")
    establishment = establishment_queue.copy()
    for column in establishment.columns:
        establishment[column] = _text(establishment[column])
    if restrict_to_registry:
        registry_islands = set(registry["island_id"])
        establishment = establishment.loc[establishment["island_id"].isin(registry_islands)].copy()
    establishment = establishment.rename(columns={"species": "reported_species"})

    work = establishment.merge(
        registry,
        on="island_id",
        how="left",
        validate="many_to_one",
        suffixes=("", "_registry"),
    )
    missing_registry = _text(work["canonical_island_name"]).eq("")
    if missing_registry.any():
        ids = work.loc[missing_registry, "island_id"].drop_duplicates().tolist()
        raise typer.BadParameter(
            "pilot registry does not cover island IDs in establishment queue: " f"{ids[:5]}"
        )
    work = work.merge(manifest[["island_id", "area_km2"]], on="island_id", how="left", validate="many_to_one")
    work = work.merge(
        taxon,
        left_on="reported_species",
        right_on="reported_name",
        how="left",
        validate="many_to_one",
        suffixes=("", "_taxon"),
    )
    missing_taxon = _text(work["reported_name"]).eq("")
    if missing_taxon.any():
        names = work.loc[missing_taxon, "reported_species"].drop_duplicates().tolist()
        raise typer.BadParameter("Taxon intake does not cover establishment species: " f"{names[:5]}")

    work["n_records"] = pd.to_numeric(work["n_records"], errors="coerce").fillna(0).astype(int)
    work["n_unique_gbif_ids"] = pd.to_numeric(work["n_unique_gbif_ids"], errors="coerce").fillna(0).astype(int)
    stages = work.apply(lambda row: classify_review_stage(row.to_dict()), axis=1)
    work[["review_stage", "selection_rationale", "required_evidence"]] = pd.DataFrame(
        stages.tolist(), index=work.index
    )
    work["review_priority"] = work["review_priority"].replace("", "normal")
    work = work.sort_values(
        by=["canonical_island_name", "review_priority", "n_records", "reported_species"],
        key=lambda column: (
            column.map(_priority_rank)
            if column.name == "review_priority"
            else (-column if column.name == "n_records" else column)
        ),
    )
    for column in WORKLIST_COLUMNS:
        if column not in work.columns:
            work[column] = ""
    return work[WORKLIST_COLUMNS].reset_index(drop=True)


def build_trait_staging(worklist: pd.DataFrame, max_per_island: int) -> pd.DataFrame:
    """Select a bounded, outcome-blind trait candidate pilot without releasing rows."""
    if max_per_island < 1:
        raise typer.BadParameter("max_per_island must be at least 1")
    eligible = worklist.loc[
        worklist["review_stage"].eq("review_island_establishment")
        & worklist["provisional_taxonomic_group"].eq("angiosperm")
        & worklist["provisional_match_status"].eq("exact_backbone_candidate")
    ].copy()
    if eligible.empty:
        return pd.DataFrame(columns=TRAIT_STAGING_COLUMNS)
    eligible = eligible.sort_values(["canonical_island_name", "n_records", "reported_species"], ascending=[True, False, True])
    selected = eligible.groupby("island_id", group_keys=False).head(max_per_island).copy()
    selected = selected.rename(
        columns={
            "provisional_accepted_species": "accepted_species",
            "reported_genus": "genus",
            "reported_family": "family",
        }
    )
    selected["trait_candidate_status"] = "staged_not_curated"
    selected["release_gate"] = (
        "Release only after accepted taxon scope and accepted island establishment; "
        "candidate traits never define inclusion."
    )
    selected["selection_rationale"] = (
        "Pilot selection uses only exact taxonomic candidate status, island context, and GBIF record count; "
        "it does not use floral traits, pollinator state, or model outcome."
    )
    return selected[TRAIT_STAGING_COLUMNS].drop_duplicates("accepted_species").reset_index(drop=True)


def summary(worklist: pd.DataFrame, trait_staging: pd.DataFrame) -> dict[str, Any]:
    """Summarise work without treating review queues as accepted data."""
    return {
        "n_island_species_worklist_rows": int(len(worklist)),
        "n_pilot_islands": int(worklist["island_id"].nunique()),
        "n_candidate_angiosperm_rows": int(worklist["provisional_taxonomic_group"].eq("angiosperm").sum()),
        "n_rows_requiring_name_resolution": int(worklist["review_stage"].eq("resolve_taxon_name").sum()),
        "n_rows_requiring_non_angiosperm_scope_review": int(
            worklist["review_stage"].eq("confirm_non_angiosperm_scope").sum()
        ),
        "n_rows_blocked_by_island_flora_source_gap": int(
            worklist["review_stage"].eq("resolve_island_flora_source_gap").sum()
        ),
        "n_rows_ready_for_island_establishment_review": int(
            worklist["review_stage"].eq("review_island_establishment").sum()
        ),
        "n_staged_trait_candidates": int(len(trait_staging)),
        "trait_staging_note": "All staged traits remain non-curated until both taxon scope and island establishment are accepted.",
        "by_island": {
            island: {
                "canonical_island_name": group["canonical_island_name"].iloc[0],
                "n_rows": int(len(group)),
                "flora_anchor_status": group["flora_anchor_status"].iloc[0],
                "curation_route": group["curation_route"].iloc[0],
                "review_stage_counts": {
                    str(stage): int(count)
                    for stage, count in group["review_stage"].value_counts().sort_index().items()
                },
            }
            for island, group in worklist.groupby("island_id", sort=True)
        },
    }


@app.command("build")
def build(
    island_manifest_csv: Path = typer.Option(..., exists=True),
    taxon_intake_csv: Path = typer.Option(..., exists=True),
    establishment_queue_csv: Path = typer.Option(..., exists=True),
    pilot_registry_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_trait_candidates_per_island: int = typer.Option(25, min=1),
    restrict_to_pilot_registry: bool = typer.Option(
        False,
        help="Scope the establishment queue to the pilot-registry islands instead of "
        "requiring the registry to cover every collected island.",
    ),
) -> None:
    """Write island-specific curation worklists and blocked trait staging rows."""
    worklist = build_worklist(
        pd.read_csv(island_manifest_csv, dtype=str).fillna(""),
        pd.read_csv(taxon_intake_csv, dtype=str).fillna(""),
        pd.read_csv(establishment_queue_csv, dtype=str).fillna(""),
        pd.read_csv(pilot_registry_csv, dtype=str).fillna(""),
        restrict_to_registry=restrict_to_pilot_registry,
    )
    trait_staging = build_trait_staging(worklist, max_trait_candidates_per_island)
    result = summary(worklist, trait_staging)
    output_dir.mkdir(parents=True, exist_ok=True)
    worklist.to_csv(output_dir / "island_curation_worklist.csv", index=False)
    trait_staging.to_csv(output_dir / "trait_candidate_pilot.csv", index=False)
    (output_dir / "pilot_curation_summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"Built {result['n_island_species_worklist_rows']} island-species curation rows and "
        f"staged {result['n_staged_trait_candidates']} blocked trait candidates."
    )


if __name__ == "__main__":
    app()
