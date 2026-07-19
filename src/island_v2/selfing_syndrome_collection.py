"""Plan and summarize the three-axis direct-evidence campaign.

Every accepted species and every axis has equal scheduling priority. Direct
collection continues independently of whether another axis is already covered;
missing values remain unresolved and global fallback is prohibited.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

DIRECT_TIERS = {"species_direct", "synonym_direct"}
INFERENCE_TIERS = {"genus_inference", "family_inference"}
ALLOWED_TIERS = DIRECT_TIERS | INFERENCE_TIERS | {"unresolved_no_evidence"}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    axes = config.get("axis_definitions") or {}
    if set(axes) != {"colour_vividness", "floral_complexity", "reproductive_assurance"}:
        raise typer.BadParameter("axis_definitions must contain exactly the three syndrome axes")
    for axis, rule in axes.items():
        traits = [str(value) for value in rule.get("target_traits") or []]
        if not traits:
            raise typer.BadParameter(f"axis {axis} has no target_traits")
    return config


def read_fills(root: Path) -> pd.DataFrame:
    paths = sorted(root.glob("shards/shard_*/fills.csv.gz"))
    if not paths:
        paths = sorted(root.glob("**/fills.csv*"))
    if not paths:
        raise typer.BadParameter(f"no fills.csv files found below {root}")
    frame = pd.concat((pd.read_csv(path, dtype=str).fillna("") for path in paths), ignore_index=True)
    required = {"accepted_species", "trait_name", "fill_tier", "filled_value"}
    missing = required.difference(frame.columns)
    if missing:
        raise typer.BadParameter(f"fill ledger missing columns: {sorted(missing)}")
    frame["accepted_species"] = frame["accepted_species"].map(_text)
    frame["trait_name"] = frame["trait_name"].map(_text)
    frame["fill_tier"] = frame["fill_tier"].map(_text)
    unknown_tiers = set(frame["fill_tier"]) - ALLOWED_TIERS
    if unknown_tiers:
        raise typer.BadParameter(
            "global or unsupported fallback tiers are prohibited: " + ", ".join(sorted(unknown_tiers))
        )
    return frame


def axis_status(frame: pd.DataFrame, config: dict[str, Any], round_index: int) -> pd.DataFrame:
    completion = config.get("primary_completion") or {}
    max_rounds = int(completion.get("max_direct_rounds", 6))
    species = sorted(set(frame["accepted_species"]))
    rows: list[dict[str, object]] = []

    for accepted_species in species:
        subset = frame.loc[frame["accepted_species"].eq(accepted_species)]
        row: dict[str, object] = {"accepted_species": accepted_species}
        missing_direct: list[str] = []
        unresolved_final: list[str] = []

        for axis, rule in config["axis_definitions"].items():
            traits = {str(value) for value in rule.get("target_traits") or []}
            axis_rows = subset.loc[subset["trait_name"].isin(traits)]
            direct = axis_rows.loc[
                axis_rows["fill_tier"].isin(DIRECT_TIERS)
                & axis_rows["filled_value"].ne("unresolved")
                & axis_rows["filled_value"].ne("")
            ]
            inferred = axis_rows.loc[
                axis_rows["fill_tier"].isin(INFERENCE_TIERS)
                & axis_rows["filled_value"].ne("unresolved")
                & axis_rows["filled_value"].ne("")
            ]
            row[f"{axis}_direct"] = not direct.empty
            row[f"{axis}_resolved"] = not direct.empty or not inferred.empty
            row[f"{axis}_resolution_tier"] = (
                "direct" if not direct.empty else "genus_or_family" if not inferred.empty else "unresolved"
            )
            row[f"{axis}_evidence_traits"] = "|".join(
                sorted(set((direct if not direct.empty else inferred)["trait_name"]))
            )
            if direct.empty:
                missing_direct.append(axis)
            if direct.empty and inferred.empty:
                unresolved_final.append(axis)

        row["all_three_axes_direct"] = not missing_direct
        row["all_three_axes_resolved"] = not unresolved_final
        row["missing_direct_axes"] = "|".join(missing_direct)
        row["unresolved_axes"] = "|".join(unresolved_final)
        row["direct_round"] = round_index
        row["continue_direct_acquisition"] = bool(missing_direct) and round_index < max_rounds
        row["ready_for_taxonomic_fallback"] = bool(missing_direct) and round_index >= max_rounds
        # Equal priority is the campaign contract: no species subset and no axis
        # is advanced ahead of another.
        row["collection_priority"] = 0
        rows.append(row)

    return pd.DataFrame(rows).sort_values(
        ["continue_direct_acquisition", "accepted_species"], ascending=[False, True]
    ).reset_index(drop=True)


def write_outputs(status: pd.DataFrame, output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    status.to_csv(output_dir / "selfing_syndrome_axis_status.csv.gz", index=False)
    queue = status.loc[status["continue_direct_acquisition"].astype(bool)].copy()
    queue.to_csv(output_dir / "direct_acquisition_queue.csv", index=False)
    fallback = status.loc[status["ready_for_taxonomic_fallback"].astype(bool)].copy()
    fallback.to_csv(output_dir / "taxonomic_fallback_queue.csv", index=False)

    summary = {
        "n_species": int(len(status)),
        "n_all_three_axes_direct": int(status["all_three_axes_direct"].sum()),
        "n_all_three_axes_resolved": int(status["all_three_axes_resolved"].sum()),
        "n_continue_direct_acquisition": int(status["continue_direct_acquisition"].sum()),
        "n_ready_for_taxonomic_fallback": int(status["ready_for_taxonomic_fallback"].sum()),
        "n_unresolved_after_family": int(status["unresolved_axes"].ne("").sum()),
        "axis_priority_policy": "equal",
        "global_fallback_used": False,
        "axis_direct_counts": {
            axis: int(status[f"{axis}_direct"].sum())
            for axis in ("colour_vividness", "floral_complexity", "reproductive_assurance")
        },
        "axis_resolved_counts": {
            axis: int(status[f"{axis}_resolved"].sum())
            for axis in ("colour_vividness", "floral_complexity", "reproductive_assurance")
        },
    }
    (output_dir / "selfing_syndrome_collection_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return summary


@app.command()
def plan(
    fill_root: Path = typer.Option(..., exists=True, file_okay=False),
    config_path: Path = typer.Option(Path("config/trait_fill_cascade.yml"), exists=True),
    output_dir: Path = typer.Option(Path("data/v2/staging/traits/selfing_syndrome_collection")),
    round_index: int = typer.Option(0, min=0),
) -> None:
    config = load_config(config_path)
    fills = read_fills(fill_root)
    status = axis_status(fills, config, round_index)
    summary = write_outputs(status, output_dir)
    typer.echo(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    app()
