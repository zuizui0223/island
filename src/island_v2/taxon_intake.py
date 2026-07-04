"""Build a reviewable taxonomy and establishment intake from raw GBIF candidates.

The collection stage deliberately retains a broad Tracheophyta occurrence universe.
Before any taxon can enter a whole-angiosperm denominator or trait workflow, its
name rank, backbone match, broad taxonomic group, and island establishment status
must be reviewed. This module creates *candidates*, never automatic curation.

GBIF's public backbone match endpoint is used only as a transparent name-matching
source. It does not decide whether a taxon is a native island flora member, and
it does not automatically set the accepted taxonomic scope used by downstream
analyses.
"""

from __future__ import annotations

import json
import re
import time
from pathlib import Path
from typing import Any

import httpx
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

GBIF_MATCH_URL = "https://api.gbif.org/v1/species/match"
REQUIRED_TAXA_COLUMNS = {"accepted_species", "genus", "family", "n_islands", "n_records"}
REQUIRED_ISLAND_SPECIES_COLUMNS = {
    "island_id",
    "species",
    "n_records",
    "n_unique_gbif_ids",
    "basis_of_record_set",
}

NON_SEED_CLASSES = {
    "Equisetopsida",
    "Lycopodiopsida",
    "Marattiopsida",
    "Ophioglossopsida",
    "Polypodiopsida",
    "Psilotopsida",
}
GYMNOSPERM_CLASSES = {"Cycadopsida", "Ginkgoopsida", "Gnetopsida", "Pinopsida"}
ANGIOSPERM_CLASSES = {"Dicotyledonae", "Liliopsida", "Magnoliopsida", "Monocotyledonae"}

QUEUE_COLUMNS = [
    "reported_name",
    "reported_genus",
    "reported_family",
    "n_islands",
    "n_records",
    "name_rank_triage",
    "gbif_usage_key",
    "gbif_scientific_name",
    "gbif_canonical_name",
    "gbif_rank",
    "gbif_status",
    "gbif_match_type",
    "gbif_confidence",
    "gbif_synonym",
    "gbif_species",
    "gbif_genus",
    "gbif_family",
    "gbif_class",
    "gbif_phylum",
    "provisional_accepted_species",
    "provisional_taxonomic_group",
    "provisional_match_status",
    "taxon_scope_review_status",
    "establishment_review_status",
    "review_priority",
    "next_review_action",
    "gbif_match_url",
]


@app.callback()
def main() -> None:
    """Create a non-curated intake queue before floral-trait extraction."""


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def classify_name_rank(name: str) -> str:
    """Conservative lexical triage; it does not establish accepted taxonomy."""
    value = " ".join(str(name or "").strip().split())
    lower = value.lower()
    if not value:
        return "blank"
    if any(marker in lower for marker in (" sp.", " spp.", " cf.", " aff.", " indet", " unknown")):
        return "indeterminate_or_placeholder"
    tokens = value.replace("×", " x ").split()
    if len(tokens) == 1:
        return "genus_or_higher_rank"
    if len(tokens) >= 2 and tokens[0][:1].isupper() and tokens[1][:1].islower():
        if len(tokens) == 2:
            return "binomial_candidate"
        return "infra_or_hybrid_candidate"
    return "unresolved_name_rank"


def provisional_group_from_backbone(match: dict[str, Any]) -> str:
    """Return a conservative broad group from a GBIF backbone candidate."""
    taxon_class = str(match.get("class") or "").strip()
    if taxon_class in NON_SEED_CLASSES:
        return "non_seed_vascular"
    if taxon_class in GYMNOSPERM_CLASSES:
        return "gymnosperm"
    if taxon_class in ANGIOSPERM_CLASSES:
        return "angiosperm"
    return "unresolved"


def provisional_species_from_backbone(match: dict[str, Any]) -> str:
    """Use the backbone species field only as a candidate name for review."""
    species = str(match.get("species") or "").strip()
    canonical = str(match.get("canonicalName") or "").strip()
    rank = str(match.get("rank") or "").strip().upper()
    if species:
        return species
    if rank == "SPECIES" and canonical:
        return canonical
    return ""


def provisional_match_status(name_rank: str, match: dict[str, Any]) -> str:
    match_type = str(match.get("matchType") or "").upper()
    confidence = pd.to_numeric(pd.Series([match.get("confidence")]), errors="coerce").iloc[0]
    rank = str(match.get("rank") or "").upper()
    if name_rank not in {"binomial_candidate", "infra_or_hybrid_candidate"}:
        return "rank_review_required"
    if match_type == "EXACT" and rank in {"SPECIES", "SUBSPECIES", "VARIETY", "FORM"}:
        return "exact_backbone_candidate"
    if match_type in {"FUZZY", "HIGHERRANK"} or pd.notna(confidence):
        return "manual_backbone_review_required"
    return "unmatched"


def _match_url(name: str) -> str:
    return f"{GBIF_MATCH_URL}?name={httpx.QueryParams({'name': name})['name']}"


def normalize_backbone_match(reported_name: str, response: dict[str, Any] | None) -> dict[str, Any]:
    """Create one normalized staging row from a raw GBIF backbone response."""
    match = response or {}
    name_rank = classify_name_rank(reported_name)
    provisional_group = provisional_group_from_backbone(match)
    provisional_species = provisional_species_from_backbone(match)
    match_status = provisional_match_status(name_rank, match)
    return {
        "name_rank_triage": name_rank,
        "gbif_usage_key": str(match.get("usageKey") or ""),
        "gbif_scientific_name": str(match.get("scientificName") or ""),
        "gbif_canonical_name": str(match.get("canonicalName") or ""),
        "gbif_rank": str(match.get("rank") or ""),
        "gbif_status": str(match.get("status") or ""),
        "gbif_match_type": str(match.get("matchType") or ""),
        "gbif_confidence": match.get("confidence") if match.get("confidence") is not None else "",
        "gbif_synonym": bool(match.get("synonym", False)),
        "gbif_species": str(match.get("species") or ""),
        "gbif_genus": str(match.get("genus") or ""),
        "gbif_family": str(match.get("family") or ""),
        "gbif_class": str(match.get("class") or ""),
        "gbif_phylum": str(match.get("phylum") or ""),
        "provisional_accepted_species": provisional_species,
        "provisional_taxonomic_group": provisional_group,
        "provisional_match_status": match_status,
        "taxon_scope_review_status": "pending",
        "establishment_review_status": "pending",
        "gbif_match_url": _match_url(reported_name),
    }


def review_priority(row: dict[str, Any]) -> tuple[str, str]:
    """Prioritize high-leverage candidates without accepting any of them."""
    rank = str(row["name_rank_triage"])
    group = str(row["provisional_taxonomic_group"])
    status = str(row["provisional_match_status"])
    n_records = int(pd.to_numeric(pd.Series([row.get("n_records", 0)]), errors="coerce").fillna(0).iloc[0])
    if rank not in {"binomial_candidate", "infra_or_hybrid_candidate"}:
        return "high", "resolve_name_rank_before_any_trait_or_scope_decision"
    if group in {"non_seed_vascular", "gymnosperm"}:
        return "medium", "review_broad_group_then_exclude_from_angiosperm_denominator_if_confirmed"
    if status != "exact_backbone_candidate":
        return "high", "review_backbone_match_and_accepted_name"
    if group != "angiosperm":
        return "medium", "review_broad_group_before_trait_work"
    if n_records >= 5:
        return "high", "review_taxonomy_and_island_establishment_before_trait_search"
    return "normal", "review_taxonomy_and_island_establishment_before_trait_search"


def build_taxon_intake(taxa: pd.DataFrame, matches: dict[str, dict[str, Any] | None]) -> pd.DataFrame:
    """Build a review queue from collector taxa plus externally retrieved matches."""
    missing = REQUIRED_TAXA_COLUMNS.difference(taxa.columns)
    if missing:
        raise typer.BadParameter(f"taxa table missing columns: {', '.join(sorted(missing))}")
    work = taxa.copy()
    for column in REQUIRED_TAXA_COLUMNS:
        work[column] = _text(work[column])
    work = work.loc[work["accepted_species"].ne("")].drop_duplicates("accepted_species")
    rows: list[dict[str, Any]] = []
    for item in work.sort_values("accepted_species").to_dict("records"):
        reported_name = str(item["accepted_species"])
        row = {
            "reported_name": reported_name,
            "reported_genus": str(item["genus"]),
            "reported_family": str(item["family"]),
            "n_islands": int(pd.to_numeric(pd.Series([item["n_islands"]]), errors="coerce").fillna(0).iloc[0]),
            "n_records": int(pd.to_numeric(pd.Series([item["n_records"]]), errors="coerce").fillna(0).iloc[0]),
            **normalize_backbone_match(reported_name, matches.get(reported_name)),
        }
        priority, action = review_priority(row)
        row["review_priority"] = priority
        row["next_review_action"] = action
        rows.append(row)
    return pd.DataFrame(rows, columns=QUEUE_COLUMNS)


def build_island_establishment_queue(island_species: pd.DataFrame, intake: pd.DataFrame) -> pd.DataFrame:
    """Create a per-island establishment review queue from retained occurrence evidence.

    GBIF record flags are evidence flags only. A blank or native-looking GBIF
    field is not a native-island classification; every row remains pending.
    """
    missing = REQUIRED_ISLAND_SPECIES_COLUMNS.difference(island_species.columns)
    if missing:
        raise typer.BadParameter(
            f"island species table missing columns: {', '.join(sorted(missing))}"
        )
    work = island_species.copy()
    for column in REQUIRED_ISLAND_SPECIES_COLUMNS:
        work[column] = _text(work[column])
    candidate = intake[[
        "reported_name",
        "provisional_accepted_species",
        "provisional_taxonomic_group",
        "provisional_match_status",
        "review_priority",
    ]].rename(columns={"reported_name": "species"})
    out = work.merge(candidate, on="species", how="left", validate="many_to_one")
    out["record_evidence_status"] = "GBIF_occurrence_only"
    out["island_establishment_status"] = "unresolved"
    out["review_status"] = "pending"
    out["next_review_action"] = "check_island_or_regional_flora_for_native_introduced_or_cultivated_status"
    columns = [
        "island_id",
        "species",
        "provisional_accepted_species",
        "provisional_taxonomic_group",
        "provisional_match_status",
        "n_records",
        "n_unique_gbif_ids",
        "basis_of_record_set",
        "record_evidence_status",
        "island_establishment_status",
        "review_priority",
        "review_status",
        "next_review_action",
    ]
    for column in columns:
        if column not in out.columns:
            out[column] = ""
    return out[columns].sort_values(["review_priority", "island_id", "species"]).reset_index(drop=True)


def fetch_matches(
    names: list[str],
    cache_dir: Path,
    pause_seconds: float = 0.05,
    timeout_seconds: float = 30.0,
) -> dict[str, dict[str, Any] | None]:
    """Retrieve and cache GBIF backbone candidates by name.

    Cached JSON is always preferred. Network or HTTP errors yield ``None`` for a
    candidate and remain reviewable as unmatched rather than aborting the batch.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    matches: dict[str, dict[str, Any] | None] = {}
    with httpx.Client(timeout=timeout_seconds, follow_redirects=True) as client:
        for name in names:
            safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "blank"
            cache_path = cache_dir / f"{safe}.json"
            if cache_path.exists():
                try:
                    matches[name] = json.loads(cache_path.read_text(encoding="utf-8"))
                except json.JSONDecodeError:
                    matches[name] = None
                continue
            try:
                response = client.get(GBIF_MATCH_URL, params={"name": name})
                response.raise_for_status()
                payload = response.json()
                cache_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
                matches[name] = payload
            except (httpx.HTTPError, ValueError):
                matches[name] = None
            if pause_seconds > 0:
                time.sleep(pause_seconds)
    return matches


@app.command("build")
def build(
    taxa_csv: Path = typer.Option(..., exists=True, help="Collector island_taxa.csv."),
    island_species_csv: Path = typer.Option(..., exists=True, help="Collector island_species_occurrences.csv."),
    output_dir: Path = typer.Option(..., help="Directory for reviewable taxon-intake products."),
    offline: bool = typer.Option(False, help="Do not call GBIF; use only cached backbone responses."),
    pause_seconds: float = typer.Option(0.05, min=0.0, help="Pause between uncached GBIF calls."),
) -> None:
    """Write taxonomy and island-establishment review queues; never curate automatically."""
    taxa = pd.read_csv(taxa_csv, dtype=str).fillna("")
    island_species = pd.read_csv(island_species_csv, dtype=str).fillna("")
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = output_dir / "gbif_backbone_raw"
    names = _text(taxa["accepted_species"]).loc[lambda values: values.ne("")].drop_duplicates().tolist()
    matches = (
        {name: None for name in names}
        if offline
        else fetch_matches(names, cache_dir=cache_dir, pause_seconds=pause_seconds)
    )
    intake = build_taxon_intake(taxa, matches)
    establishment = build_island_establishment_queue(island_species, intake)
    intake.to_csv(output_dir / "taxon_intake_review_queue.csv", index=False)
    establishment.to_csv(output_dir / "island_establishment_review_queue.csv", index=False)
    scope_template = intake.loc[
        intake["provisional_accepted_species"].astype(str).str.strip().ne(""),
        ["provisional_accepted_species", "provisional_taxonomic_group"],
    ].drop_duplicates("provisional_accepted_species")
    scope_template = scope_template.rename(columns={"provisional_accepted_species": "accepted_species"})
    scope_template["decision_basis"] = "GBIF backbone candidate; human review required"
    scope_template["source_citation"] = "GBIF Backbone Taxonomy match endpoint"
    scope_template["source_url"] = "https://api.gbif.org/v1/species/match"
    scope_template["review_status"] = "pending"
    scope_template["reviewer"] = ""
    scope_template["decision_date"] = ""
    scope_template.to_csv(output_dir / "taxon_scope_decisions_prefilled.csv", index=False)
    summary = {
        "n_input_taxa": int(len(taxa)),
        "n_queue_rows": int(len(intake)),
        "n_exact_backbone_candidates": int(intake["provisional_match_status"].eq("exact_backbone_candidate").sum()),
        "n_provisional_angiosperm_candidates": int(intake["provisional_taxonomic_group"].eq("angiosperm").sum()),
        "n_non_seed_vascular_candidates": int(intake["provisional_taxonomic_group"].eq("non_seed_vascular").sum()),
        "n_gymnosperm_candidates": int(intake["provisional_taxonomic_group"].eq("gymnosperm").sum()),
        "n_unresolved_group_candidates": int(intake["provisional_taxonomic_group"].eq("unresolved").sum()),
        "n_island_establishment_queue_rows": int(len(establishment)),
        "offline": bool(offline),
    }
    (output_dir / "taxon_intake_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"Built {summary['n_queue_rows']} taxon review rows and "
        f"{summary['n_island_establishment_queue_rows']} island-establishment review rows."
    )


if __name__ == "__main__":
    app()
