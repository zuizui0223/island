"""Outcome-blind discovery waves for N2 lineage functional dependency evidence.

The target universe is derived only from the frozen PR138 mainland source assignments
and native mainland source flora.  Island realised flora, genus-entry outcomes, N1/N2
effect directions, focal floral traits, and pollination guild labels are not inputs.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

ASSIGNMENT_REQUIRED = {"island_id", "source_mode", "source_rank", "entity_ID"}
FLORA_REQUIRED = {"entity_ID", "work_species"}
RANKED_COLUMNS = [
    "accepted_genus",
    "n_islands_with_source_available_genus",
    "n_selected_mainland_source_entities_containing_genus",
    "global_rank",
    "wave",
    "wave_rank",
    "source_mode",
]
HOLDOUT_COLUMNS = [
    "raw_genus_token",
    "reason",
    "n_species_rows",
    "n_selected_mainland_source_entities",
]
GENUS_TOKEN = re.compile(r"^[A-Z][A-Za-z-]+$")


def load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("N2 dependency discovery config must be a mapping")
    if payload.get("contract") != "chapter1_nee_n2_dependency_discovery_v1":
        raise typer.BadParameter("unexpected N2 dependency discovery contract")
    return payload


def _require(table: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = columns.difference(table.columns)
    if missing:
        raise ValueError(f"{label} missing columns: {sorted(missing)}")


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def source_genus_token(work_species: object) -> tuple[str, str | None]:
    """Return the literal first taxonomic token and any primary-search holdout reason."""
    value = "" if work_species is None else " ".join(str(work_species).split())
    if not value:
        return "", "missing_work_species"
    token = value.split(" ", 1)[0]
    if token.startswith("×"):
        return token, "hybrid_or_nothogenus_marker"
    if not GENUS_TOKEN.fullmatch(token):
        return token, "nonstandard_genus_token"
    return token, None


def _validate_primary_assignments(assignments: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    _require(assignments, ASSIGNMENT_REQUIRED, "source assignments")
    work = assignments.copy()
    for column in ASSIGNMENT_REQUIRED:
        work[column] = _text(work[column])

    primary_mode = str(config["ranking"]["primary_source_mode"])
    expected_ranks = {int(value) for value in config["ranking"]["source_rank_set"]}
    primary = work.loc[work["source_mode"].eq(primary_mode)].copy()
    if primary.empty:
        raise ValueError(f"no source assignments for primary mode {primary_mode}")
    primary["source_rank"] = pd.to_numeric(primary["source_rank"], errors="raise").astype(int)
    if primary[["island_id", "source_rank"]].duplicated().any():
        raise ValueError("duplicate island x source_rank rows in primary source assignments")
    invalid_ranks = sorted(set(primary["source_rank"]).difference(expected_ranks))
    if invalid_ranks:
        raise ValueError(f"unexpected source ranks in primary mode: {invalid_ranks}")

    observed = primary.groupby("island_id", sort=False)["source_rank"].agg(lambda x: set(x.tolist()))
    bad = observed.loc[observed.map(lambda ranks: ranks != expected_ranks)]
    if not bad.empty:
        raise ValueError(
            f"primary source assignments require exactly ranks {sorted(expected_ranks)} for every island; "
            f"n_bad_islands={len(bad)}"
        )
    if primary["entity_ID"].eq("").any():
        raise ValueError("primary source assignments contain blank entity_ID")
    return primary


def build_ranked_universe(
    assignments: pd.DataFrame,
    source_flora: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Rank genera only by prespecified source-side opportunity."""
    primary = _validate_primary_assignments(assignments, config)
    _require(source_flora, FLORA_REQUIRED, "source flora")
    flora = source_flora.copy()
    flora["entity_ID"] = _text(flora["entity_ID"])
    flora["work_species"] = _text(flora["work_species"])

    selected_entities = set(primary["entity_ID"])
    flora = flora.loc[flora["entity_ID"].isin(selected_entities)].copy()
    if flora.empty:
        raise ValueError("source flora has no rows for selected primary source entities")

    parsed = flora["work_species"].map(source_genus_token)
    flora["raw_genus_token"] = parsed.map(lambda item: item[0])
    flora["holdout_reason"] = parsed.map(lambda item: item[1] or "")

    holdout_source = flora.loc[flora["holdout_reason"].ne("")].copy()
    if holdout_source.empty:
        holdouts = pd.DataFrame(columns=HOLDOUT_COLUMNS)
    else:
        holdouts = (
            holdout_source.groupby(["raw_genus_token", "holdout_reason"], dropna=False)
            .agg(
                n_species_rows=("work_species", "size"),
                n_selected_mainland_source_entities=("entity_ID", "nunique"),
            )
            .reset_index()
            .rename(columns={"holdout_reason": "reason"})
            .sort_values(["reason", "raw_genus_token"], kind="mergesort")
            .reset_index(drop=True)
        )

    eligible = flora.loc[flora["holdout_reason"].eq("")].copy()
    eligible = eligible.rename(columns={"raw_genus_token": "accepted_genus"})
    entity_genus = eligible[["entity_ID", "accepted_genus"]].drop_duplicates()
    memberships = primary[["island_id", "entity_ID"]].drop_duplicates()
    island_genus = memberships.merge(entity_genus, on="entity_ID", how="inner", validate="many_to_many")

    ranked = (
        island_genus.groupby("accepted_genus", sort=False)
        .agg(
            n_islands_with_source_available_genus=("island_id", "nunique"),
            n_selected_mainland_source_entities_containing_genus=("entity_ID", "nunique"),
        )
        .reset_index()
        .sort_values(
            [
                "n_islands_with_source_available_genus",
                "n_selected_mainland_source_entities_containing_genus",
                "accepted_genus",
            ],
            ascending=[False, False, True],
            kind="mergesort",
        )
        .reset_index(drop=True)
    )
    wave_size = int(config["waves"]["genera_per_wave"])
    if wave_size < 1:
        raise ValueError("genera_per_wave must be >= 1")
    ranked["global_rank"] = range(1, len(ranked) + 1)
    ranked["wave"] = ((ranked["global_rank"] - 1) // wave_size) + 1
    ranked["wave_rank"] = ((ranked["global_rank"] - 1) % wave_size) + 1
    ranked["source_mode"] = str(config["ranking"]["primary_source_mode"])
    ranked = ranked[RANKED_COLUMNS]

    receipt = {
        "contract": config["contract"],
        "status": "ranked_source_opportunity_universe_built_pre_dependency_search",
        "primary_source_mode": str(config["ranking"]["primary_source_mode"]),
        "n_islands": int(primary["island_id"].nunique()),
        "n_selected_mainland_source_entities": int(primary["entity_ID"].nunique()),
        "n_search_eligible_genera": int(len(ranked)),
        "n_holdout_genus_tokens": int(len(holdouts)),
        "n_waves": int(math.ceil(len(ranked) / wave_size)) if len(ranked) else 0,
        "genera_per_wave": wave_size,
        "uses_island_realised_flora": False,
        "uses_island_genus_entry": False,
        "uses_N1_effect_direction": False,
        "uses_N2_effect_direction": False,
        "uses_focal_flower_traits": False,
        "uses_pollination_guild": False,
        "dependency_evidence_opened": False,
        "N2_fitted": False,
    }
    return ranked, holdouts, receipt


def select_wave(ranked: pd.DataFrame, wave_number: int) -> pd.DataFrame:
    if wave_number < 1:
        raise ValueError("wave_number must be >= 1")
    if "wave" not in ranked.columns:
        raise ValueError("ranked universe missing wave column")
    return ranked.loc[pd.to_numeric(ranked["wave"], errors="coerce").eq(wave_number)].copy().reset_index(drop=True)


@app.command("build")
def build_command(
    source_assignments_csv: Path = typer.Option(..., exists=True),
    source_flora_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    wave_number: int = typer.Option(1, min=1),
    config_path: Path = typer.Option(Path("config/chapter1_nee_n2_dependency_discovery.yml")),
) -> None:
    config = load_config(config_path)
    assignments = pd.read_csv(source_assignments_csv, dtype=str).fillna("")
    flora = pd.read_csv(source_flora_csv, dtype=str).fillna("")
    ranked, holdouts, receipt = build_ranked_universe(assignments, flora, config)
    wave = select_wave(ranked, wave_number)
    output_dir.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(output_dir / "n2_dependency_discovery_ranked_universe.csv", index=False)
    holdouts.to_csv(output_dir / "n2_dependency_discovery_genus_holdouts.csv", index=False)
    wave.to_csv(output_dir / f"n2_dependency_discovery_wave_{wave_number:03d}.csv", index=False)
    receipt = {
        **receipt,
        "selected_wave": int(wave_number),
        "n_genera_in_selected_wave": int(len(wave)),
        "selected_wave_first_global_rank": int(wave["global_rank"].min()) if not wave.empty else None,
        "selected_wave_last_global_rank": int(wave["global_rank"].max()) if not wave.empty else None,
    }
    (output_dir / "n2_dependency_discovery_receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(receipt))


if __name__ == "__main__":
    app()
