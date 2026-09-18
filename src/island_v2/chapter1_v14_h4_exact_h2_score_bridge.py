"""Post-hoc functional bridge using the exact Chapter 1 v14 H2 species scores."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.chapter1_h5_glopl_global_distance import _truthy
from island_v2.chapter1_v14_h4_family_bridge import (
    aggregate_family_cells,
    fit_global_family_score,
    unique_effect_rows,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)
CONTRACT = "chapter1_v14_h4_exact_h2_score_bridge_v1"


def load_config(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict) or value.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected exact H2-score bridge contract")
    if value.get("inferential_role") != "posthoc_exact_H2_score_functional_triangulation":
        raise typer.BadParameter("exact H2-score bridge must remain post-hoc")
    return value


def _name_key(value: object) -> str:
    return " ".join(str(value).replace("_", " ").split()).casefold()


def exact_h2_scores(
    species_scores: pd.DataFrame,
    *,
    syndrome: str,
    score_name: str,
) -> pd.DataFrame:
    required = {
        "accepted_species",
        "syndrome",
        "soft_membership",
    }
    if missing := required - set(species_scores.columns):
        raise typer.BadParameter(f"H2 species score table missing columns: {sorted(missing)}")

    work = species_scores.loc[
        species_scores["syndrome"].astype(str).eq(str(syndrome)),
        ["accepted_species", "soft_membership"],
    ].copy()
    work[score_name] = pd.to_numeric(work["soft_membership"], errors="coerce")
    work = work.dropna(subset=["accepted_species", score_name])
    work["species_key"] = work["accepted_species"].map(_name_key)

    conflicts = (
        work.groupby("species_key")[score_name]
        .nunique()
        .gt(1)
    )
    if bool(conflicts.any()):
        raise typer.BadParameter("exact H2 score table has conflicting species scores")

    return (
        work[["species_key", score_name]]
        .drop_duplicates("species_key")
        .reset_index(drop=True)
    )


def analyse_exact_family(
    species_scores: pd.DataFrame,
    matched_effect_rows: pd.DataFrame,
    *,
    family: str,
    family_config: dict[str, Any],
    publication_total_weight: float,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    score_name = str(family_config["score_name"])
    scores = exact_h2_scores(
        species_scores,
        syndrome=str(family_config["H2_syndrome"]),
        score_name=score_name,
    )
    outcomes = unique_effect_rows(matched_effect_rows).copy()
    outcomes["species_key"] = outcomes["species_key"].map(_name_key)

    cells = aggregate_family_cells(
        outcomes,
        scores,
        score_name=score_name,
        publication_total_weight=publication_total_weight,
    )
    rows: list[dict[str, Any]] = []
    for analysis, part in {
        "primary": cells,
        "supplemental_only": cells.loc[
            cells["PL_Effect_Size_Type2"].astype(str).eq("Sup")
        ].copy(),
        "no_zero_constant": cells.loc[~_truthy(cells["Constant_added"])].copy(),
    }.items():
        rows.append(
            {
                "family": family,
                "H2_syndrome": str(family_config["H2_syndrome"]),
                "score_name": score_name,
                "analysis": analysis,
                "inferential_role": "posthoc_exact_H2_score_functional_triangulation",
                **fit_global_family_score(part, score_name=score_name),
            }
        )

    support = {
        "family": family,
        "H2_syndrome": str(family_config["H2_syndrome"]),
        "score_name": score_name,
        "n_H2_score_species": int(len(scores)),
        "n_matched_species": int(cells["species_key"].nunique()),
        "n_matched_publications": int(cells["study_key"].nunique()),
        "n_matched_sites": int(cells["site_key"].nunique()),
    }
    return rows, support


def run_exact_h2_bridge(
    species_scores: pd.DataFrame,
    reproductive_effect_rows: pd.DataFrame,
    architecture_effect_rows: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    parents = {
        "reproductive_assurance": reproductive_effect_rows,
        "floral_architecture": architecture_effect_rows,
    }
    rows: list[dict[str, Any]] = []
    support: dict[str, Any] = {}
    for family, family_config in config["families"].items():
        family_rows, family_support = analyse_exact_family(
            species_scores,
            parents[str(family_config["parent"])],
            family=str(family),
            family_config=family_config,
            publication_total_weight=float(config["analysis"]["publication_total_weight"]),
        )
        rows.extend(family_rows)
        support[str(family)] = family_support

    results = pd.DataFrame(rows)
    manifest = {
        "contract": CONTRACT,
        "inferential_role": "posthoc_exact_H2_score_functional_triangulation",
        "support": support,
        "all_primary_estimates_negative": bool(
            results.loc[results["analysis"].eq("primary"), "estimate"].lt(0).all()
        ),
        "historical_selection_identified": False,
        "mediation_identified": False,
        "prospective_replication": False,
    }
    return results, manifest


@app.command("analyse")
def analyse(
    h2_species_scores_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    reproductive_effect_rows_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    architecture_effect_rows_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    results, manifest = run_exact_h2_bridge(
        pd.read_csv(h2_species_scores_csv),
        pd.read_csv(reproductive_effect_rows_csv),
        pd.read_csv(architecture_effect_rows_csv),
        config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "h4_exact_h2_score_bridge_results.csv", index=False)
    (output_dir / "h4_exact_h2_score_bridge_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
