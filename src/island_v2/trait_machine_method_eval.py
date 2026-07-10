"""Evaluate trait-acquisition methods for no-human-review machine outputs.

The outputs here are deliberately named machine candidates/selections. They are
not human-adjudicated curated truth. The goal is to compare source lanes and
produce usable, provenance-rich machine layers when the workflow must avoid
manual review.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2 import trait_candidate_review_queue as queue_tools

app = typer.Typer(add_completion=False, no_args_is_help=True)

COLOUR_VALUE_HINTS = {
    "white": "white",
    "cream": "white",
    "creamy": "white",
    "yellow": "yellow_orange",
    "yellowish": "yellow_orange",
    "orange": "yellow_orange",
    "red": "red_pink",
    "reddish": "red_pink",
    "pink": "red_pink",
    "purple": "blue_purple",
    "blue": "blue_purple",
    "green": "green_brown_inconspicuous",
    "brown": "green_brown_inconspicuous",
    "brownish": "green_brown_inconspicuous",
    "black": "green_brown_inconspicuous",
}

MACHINE_COLUMNS = [
    "machine_candidate_id",
    "accepted_species",
    "trait_layer",
    "trait_name",
    "candidate_value",
    "machine_final_value",
    "candidate_class",
    "machine_method",
    "machine_use_tier",
    "evidence_level",
    "source_lane",
    "evidence_scope",
    "source",
    "source_type",
    "source_url",
    "source_citation",
    "basis_traits",
    "basis_family",
    "raw_description",
    "no_human_review_policy",
]


@app.callback()
def main() -> None:
    """Compare source lanes and export no-human-review machine layers."""


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _read_csv(path: Path | None) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, dtype=str).fillna("")


def _candidate_id(parts: list[str]) -> str:
    digest = hashlib.sha1("\x1f".join(parts).encode("utf-8")).hexdigest()[:16]
    return f"machine_trait:{digest}"


def load_ontology(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    return data if isinstance(data, dict) else {}


def machine_final_value(trait_name: str, candidate_value: str, ontology: dict[str, Any]) -> str:
    """Map a raw candidate value to the controlled vocabulary when possible."""
    trait = _text(trait_name)
    value = _text(candidate_value)
    allowed = queue_tools._allowed_final_values(trait, ontology)  # noqa: SLF001
    if allowed is None:
        return value
    if value in allowed:
        return value
    if trait == "flower_primary_color":
        mapped = COLOUR_VALUE_HINTS.get(value.casefold())
        if mapped in allowed:
            return str(mapped)
    return ""


def _method_for(source_lane: str, evidence_scope: str, candidate_class: str) -> tuple[str, str, str]:
    if source_lane == "web_reported_scout" and evidence_scope == "species_direct":
        return ("web_reported_species_direct", "main_text_reported", "high")
    if source_lane == "web_reported_scout":
        return ("web_reported_species_indirect", "sensitivity_text_reported", "medium")
    if source_lane == "m0_descriptive_scout":
        return ("gbif_descriptive_scout", "sensitivity_text_reported", "low")
    if source_lane == "trait_proxy_generator" or candidate_class == "proxy":
        return ("rule_based_trait_proxy", "proxy_sensitivity_only", "low")
    if source_lane == "llm_reported_ecology" and evidence_scope == "species_direct":
        return ("llm_reported_ecology_excerpt", "main_text_reported", "high")
    if source_lane == "llm_reported_ecology":
        return ("llm_reported_ecology_excerpt", "sensitivity_text_reported", "medium")
    if source_lane == "visual_trait_candidate":
        return ("visual_image_candidate", "visual_sensitivity_only", "medium")
    return (source_lane or "unknown", "sensitivity_only", "low")


def _row_from_queue(row: pd.Series, ontology: dict[str, Any]) -> dict[str, str]:
    source_lane = _text(row.get("source_lane"))
    evidence_scope = _text(row.get("evidence_scope"))
    candidate_class = _text(row.get("candidate_class"))
    trait_name = _text(row.get("trait_name"))
    candidate_value = _text(row.get("candidate_value"))
    final_value = machine_final_value(trait_name, candidate_value, ontology)
    method, tier, evidence_level = _method_for(source_lane, evidence_scope, candidate_class)
    parts = [
        source_lane,
        _text(row.get("accepted_species")),
        trait_name,
        candidate_value,
        _text(row.get("source_url")),
        _text(row.get("raw_description")),
    ]
    return {
        "machine_candidate_id": _candidate_id(parts),
        "accepted_species": _text(row.get("accepted_species")),
        "trait_layer": _text(row.get("trait_layer")),
        "trait_name": trait_name,
        "candidate_value": candidate_value,
        "machine_final_value": final_value,
        "candidate_class": candidate_class,
        "machine_method": method,
        "machine_use_tier": tier if final_value else "unmapped_value_holdout",
        "evidence_level": evidence_level if final_value else "unmapped",
        "source_lane": source_lane,
        "evidence_scope": evidence_scope,
        "source": _text(row.get("source")),
        "source_type": _text(row.get("source_type")),
        "source_url": _text(row.get("source_url")),
        "source_citation": _text(row.get("source_citation")),
        "basis_traits": _text(row.get("basis_traits")),
        "basis_family": _text(row.get("basis_family")),
        "raw_description": _text(row.get("raw_description")),
        "no_human_review_policy": (
            "Machine-only layer. Use main_text_reported for primary no-review analyses; "
            "keep proxy and sensitivity tiers separate."
        ),
    }


def machine_candidates_from_queue(queue: pd.DataFrame, ontology: dict[str, Any]) -> pd.DataFrame:
    if queue.empty:
        return pd.DataFrame(columns=MACHINE_COLUMNS)
    rows = [_row_from_queue(row, ontology) for _, row in queue.iterrows()]
    return pd.DataFrame(rows, columns=MACHINE_COLUMNS).drop_duplicates("machine_candidate_id")


def select_machine_traits(candidates: pd.DataFrame) -> pd.DataFrame:
    """Aggregate high-confidence machine candidates into one row per species/trait."""
    columns = [
        "accepted_species",
        "trait_name",
        "machine_final_value",
        "machine_selection_rule",
        "n_supporting_rows",
        "supporting_candidate_ids",
        "supporting_methods",
        "supporting_sources",
        "raw_descriptions",
        "machine_use_tier",
    ]
    if candidates.empty:
        return pd.DataFrame(columns=columns)
    usable = candidates.loc[
        candidates["machine_use_tier"].eq("main_text_reported")
        & candidates["machine_final_value"].astype(str).str.strip().ne("")
    ].copy()
    rows: list[dict[str, str]] = []
    for (species, trait_name), group in usable.groupby(["accepted_species", "trait_name"], sort=True):
        values = sorted(v for v in group["machine_final_value"].drop_duplicates() if _text(v))
        if not values:
            continue
        if len(values) == 1:
            final_value = values[0]
            rule = "single_explicit_machine_value"
        elif trait_name == "flower_primary_color":
            final_value = "multicolored_variable"
            rule = "multiple_explicit_colour_bins"
        else:
            final_value = ""
            rule = "ambiguous_multiple_machine_values"
        if not final_value:
            continue
        rows.append(
            {
                "accepted_species": species,
                "trait_name": trait_name,
                "machine_final_value": final_value,
                "machine_selection_rule": rule,
                "n_supporting_rows": str(len(group)),
                "supporting_candidate_ids": "|".join(group["machine_candidate_id"]),
                "supporting_methods": "|".join(sorted(set(group["machine_method"]))),
                "supporting_sources": "|".join(sorted(set(group["source_url"]))),
                "raw_descriptions": " || ".join(group["raw_description"].astype(str).head(5)),
                "machine_use_tier": "main_text_reported",
            }
        )
    return pd.DataFrame(rows, columns=columns)


def sensitivity_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    if candidates.empty:
        return candidates
    return candidates.loc[
        candidates["machine_use_tier"].isin(
            ["sensitivity_text_reported", "visual_sensitivity_only", "proxy_sensitivity_only"]
        )
    ].reset_index(drop=True)


def method_comparison(
    species: pd.DataFrame,
    candidates: pd.DataFrame,
    globi: pd.DataFrame,
    visual_rows: int = 0,
    llm_rows: int = 0,
) -> pd.DataFrame:
    total_species = int(species["accepted_species"].astype(str).str.strip().ne("").sum()) if not species.empty else 0
    rows: list[dict[str, Any]] = []

    def add_candidate_method(method: str, recommended_role: str, reason: str) -> None:
        subset = candidates.loc[candidates["machine_method"].eq(method)] if not candidates.empty else candidates
        valid = subset.loc[subset["machine_final_value"].astype(str).str.strip().ne("")] if not subset.empty else subset
        conflict_groups = 0
        if not valid.empty:
            conflict_groups = int(
                valid.groupby(["accepted_species", "trait_name"])["machine_final_value"]
                .nunique()
                .gt(1)
                .sum()
            )
        rows.append(
            {
                "method": method,
                "n_rows": int(len(subset)),
                "n_species": int(subset["accepted_species"].nunique()) if not subset.empty else 0,
                "species_coverage_rate": (
                    round(float(subset["accepted_species"].nunique()) / total_species, 4)
                    if total_species and not subset.empty
                    else 0.0
                ),
                "n_valid_machine_values": int(len(valid)),
                "n_conflict_species_trait_groups": conflict_groups,
                "recommended_role": recommended_role,
                "reason": reason,
            }
        )

    add_candidate_method(
        "web_reported_species_direct",
        "best_current_no_review_trait_lane",
        "Explicit species-page text with source URLs/excerpts; highest precision among current trait methods.",
    )
    add_candidate_method(
        "web_reported_species_indirect",
        "sensitivity_only",
        "Useful lead, but genus/indirect pages are weaker without manual checking.",
    )
    add_candidate_method(
        "gbif_descriptive_scout",
        "not_recommended_as_main",
        "Very sparse in the validation pilot and prone to description-context noise.",
    )
    add_candidate_method(
        "rule_based_trait_proxy",
        "proxy_sensitivity_only",
        "Expands coverage but is phenotype/family inference, not source-stated trait evidence.",
    )
    add_candidate_method(
        "llm_reported_ecology_excerpt",
        "next_main_text_lane_for_reproduction_and_visitors",
        "Source-excerpted reproductive, selfing, pollen-vector, or visitor statements; usable as main no-review text evidence when species-direct.",
    )
    add_candidate_method(
        "visual_image_candidate",
        "later_auxiliary_visible_traits",
        "Useful for visible M0 traits only; heavier than text and not suitable for reproductive or visitor states.",
    )
    rows.append(
        {
            "method": "globi_interaction_claims",
            "n_rows": int(len(globi)),
            "n_species": int(globi["query_taxon"].nunique()) if not globi.empty else 0,
            "species_coverage_rate": (
                round(float(globi["query_taxon"].nunique()) / total_species, 4)
                if total_species and not globi.empty
                else 0.0
            ),
            "n_valid_machine_values": int(len(globi)),
            "n_conflict_species_trait_groups": 0,
            "recommended_role": "best_current_no_review_interaction_lane",
            "reason": "Source-backed pollination/flower-visit claims; sparse and not effectiveness or guild by itself.",
        }
    )
    return pd.DataFrame(rows)


def interaction_claims(globi: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "evidence_id",
        "query_taxon",
        "interaction_type",
        "interaction_claim_class",
        "partner_taxon_name",
        "study_url",
        "study_doi",
        "study_citation",
        "effectiveness_class",
        "machine_use_tier",
        "no_human_review_policy",
    ]
    if globi.empty:
        return pd.DataFrame(columns=columns)
    table = globi.copy()
    table["machine_use_tier"] = "source_backed_interaction_claim"
    table["no_human_review_policy"] = (
        "Machine-only interaction layer. Treat as an explicit claim, not proof of "
        "effectiveness or a pollinator guild."
    )
    for column in columns:
        if column not in table.columns:
            table[column] = ""
    return table[columns].drop_duplicates("evidence_id").reset_index(drop=True)


def _markdown_table(table: pd.DataFrame) -> str:
    if table.empty:
        return "(no rows)"
    columns = list(table.columns)
    rows = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for record in table.astype(str).to_dict("records"):
        values = [record[column].replace("|", "\\|") for column in columns]
        rows.append("| " + " | ".join(values) + " |")
    return "\n".join(rows)


def write_markdown(
    path: Path,
    comparison: pd.DataFrame,
    selected: pd.DataFrame,
    sensitivity: pd.DataFrame,
    interactions: pd.DataFrame,
) -> None:
    best_traits = comparison.loc[
        comparison["method"].eq("web_reported_species_direct")
    ].to_dict("records")[0]
    best_interactions = comparison.loc[
        comparison["method"].eq("globi_interaction_claims")
    ].to_dict("records")[0]
    lines = [
        "# Machine Method Evaluation",
        "",
        "Human review was not required for this run. Outputs are machine layers, not human-adjudicated truth.",
        "",
        "## Recommendation",
        "",
        (
            "- Use `web_reported_species_direct` as the current main no-review trait lane: "
            f"{best_traits['n_rows']} rows across {best_traits['n_species']} species."
        ),
        (
            "- Use `globi_interaction_claims` as the current no-review interaction lane: "
            f"{best_interactions['n_rows']} records across {best_interactions['n_species']} taxa."
        ),
        "- Keep `rule_based_trait_proxy` and indirect/descriptive text as sensitivity layers.",
        "- Build the next LLM extractor for reproductive mode, selfing, and visitor statements; do not use images for those.",
        "- Use images later only for visible floral traits, because image acquisition/inference is heavier.",
        "",
        "## Output Counts",
        "",
        f"- `machine_trait_selected.csv`: {len(selected)} selected species-trait rows",
        f"- `machine_trait_sensitivity.csv`: {len(sensitivity)} sensitivity/proxy rows",
        f"- `machine_interaction_claims.csv`: {len(interactions)} interaction-claim rows",
        "",
        "## Method Table",
        "",
        _markdown_table(comparison),
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


@app.command("evaluate")
def evaluate(
    species_csv: Path = typer.Option(..., exists=True, help="Validation species CSV."),
    review_queue_csv: Path = typer.Option(..., exists=True, help="Candidate review queue CSV."),
    globi_csv: Path | None = typer.Option(None, exists=True, help="Optional GloBI interaction evidence CSV."),
    output_dir: Path = typer.Option(..., help="Directory for machine evaluation outputs."),
    ontology_path: Path = typer.Option(
        Path("config/trait_ontology.yml"),
        exists=True,
        help="Trait ontology for machine value mapping.",
    ),
) -> None:
    """Evaluate method performance and write no-human-review machine outputs."""
    species = _read_csv(species_csv)
    queue = _read_csv(review_queue_csv)
    globi = _read_csv(globi_csv)
    ontology = load_ontology(ontology_path)

    candidates = machine_candidates_from_queue(queue, ontology)
    selected = select_machine_traits(candidates)
    sensitivity = sensitivity_candidates(candidates)
    interactions = interaction_claims(globi)
    comparison = method_comparison(species, candidates, globi)

    output_dir.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(output_dir / "machine_trait_candidates.csv", index=False)
    selected.to_csv(output_dir / "machine_trait_selected.csv", index=False)
    sensitivity.to_csv(output_dir / "machine_trait_sensitivity.csv", index=False)
    interactions.to_csv(output_dir / "machine_interaction_claims.csv", index=False)
    comparison.to_csv(output_dir / "method_comparison.csv", index=False)
    summary = {
        "n_species": int(species["accepted_species"].nunique()) if not species.empty else 0,
        "n_machine_candidates": int(len(candidates)),
        "n_machine_trait_selected": int(len(selected)),
        "n_machine_trait_sensitivity": int(len(sensitivity)),
        "n_machine_interaction_claims": int(len(interactions)),
        "recommended_main_trait_method": "web_reported_species_direct",
        "recommended_main_interaction_method": "globi_interaction_claims",
        "next_build": "llm_reported_ecology_excerpt for reproductive mode, selfing capacity, and visitor statements.",
        "image_policy": "Use later for visible M0 traits only; do not infer reproductive or visitor states from images.",
        "human_review_required": False,
        "machine_layer_warning": "No human review was required, so outputs remain machine layers with provenance and method tiers.",
        "methods": comparison.to_dict("records"),
    }
    (output_dir / "method_comparison.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_markdown(
        output_dir / "method_comparison.md",
        comparison,
        selected,
        sensitivity,
        interactions,
    )
    typer.echo(
        f"Wrote machine evaluation: {len(selected)} selected trait row(s), "
        f"{len(sensitivity)} sensitivity row(s), {len(interactions)} interaction claim(s)."
    )


if __name__ == "__main__":
    app()
