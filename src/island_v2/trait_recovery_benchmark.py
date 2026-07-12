"""Benchmark competing trait-recovery lanes on a fixed species denominator.

The benchmark compares coverage without treating missing evidence as biological
absence. Each lane is a CSV of source-backed candidate rows. The evaluator
reports per-trait species coverage, union coverage, and the marginal gain from
adding each lane in a specified order.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

SPECIES_COLUMN = "accepted_species"
TRAIT_COLUMN = "trait_name"
VALUE_COLUMNS = (
    "provisional_candidate_value",
    "candidate_value",
    "standardized_value",
    "value",
)


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _value_column(table: pd.DataFrame) -> str:
    for column in VALUE_COLUMNS:
        if column in table.columns:
            return column
    raise typer.BadParameter(
        f"candidate table must contain one of value columns: {', '.join(VALUE_COLUMNS)}"
    )


def read_species(path: Path, species_column: str = SPECIES_COLUMN) -> list[str]:
    table = pd.read_csv(path, dtype=str).fillna("")
    if species_column not in table.columns:
        raise typer.BadParameter(f"species column {species_column!r} not found in {path}")
    seen: set[str] = set()
    ordered: list[str] = []
    for value in table[species_column]:
        species = _text(value)
        if species and species not in seen:
            seen.add(species)
            ordered.append(species)
    return ordered


def read_lane(path: Path, denominator: set[str]) -> set[tuple[str, str]]:
    table = pd.read_csv(path, dtype=str).fillna("")
    for column in (SPECIES_COLUMN, TRAIT_COLUMN):
        if column not in table.columns:
            raise typer.BadParameter(f"{column!r} not found in {path}")
    value_column = _value_column(table)
    recovered: set[tuple[str, str]] = set()
    for row in table.to_dict("records"):
        species = _text(row.get(SPECIES_COLUMN))
        trait = _text(row.get(TRAIT_COLUMN))
        value = _text(row.get(value_column))
        if species not in denominator or not trait or not value or value == "unresolved":
            continue
        recovered.add((species, trait))
    return recovered


def _coverage_rows(
    recovered: set[tuple[str, str]],
    denominator: set[str],
    lane: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    traits = sorted({trait for _, trait in recovered})
    for trait in traits:
        species = {sp for sp, tr in recovered if tr == trait}
        rows.append(
            {
                "lane": lane,
                "trait_name": trait,
                "n_species_recovered": len(species),
                "n_species_denominator": len(denominator),
                "coverage_fraction": len(species) / len(denominator) if denominator else 0.0,
            }
        )
    any_species = {species for species, _ in recovered}
    rows.append(
        {
            "lane": lane,
            "trait_name": "__any_trait__",
            "n_species_recovered": len(any_species),
            "n_species_denominator": len(denominator),
            "coverage_fraction": len(any_species) / len(denominator) if denominator else 0.0,
        }
    )
    return rows


def benchmark_lanes(
    *,
    species: list[str],
    lanes: dict[str, set[tuple[str, str]]],
    order: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    denominator = set(species)
    coverage_rows: list[dict[str, Any]] = []
    for lane, recovered in lanes.items():
        coverage_rows.extend(_coverage_rows(recovered, denominator, lane))

    union_all: set[tuple[str, str]] = set().union(*lanes.values()) if lanes else set()
    coverage_rows.extend(_coverage_rows(union_all, denominator, "__union_all__"))
    coverage = pd.DataFrame(coverage_rows)

    sequence = order or list(lanes)
    unknown = [lane for lane in sequence if lane not in lanes]
    if unknown:
        raise typer.BadParameter(f"unknown lane(s) in order: {unknown}")

    cumulative: set[tuple[str, str]] = set()
    marginal_rows: list[dict[str, Any]] = []
    for rank, lane in enumerate(sequence, start=1):
        before = set(cumulative)
        cumulative |= lanes[lane]
        new_pairs = cumulative - before
        new_species = {species for species, _ in new_pairs}
        marginal_rows.append(
            {
                "rank": rank,
                "lane": lane,
                "new_species_trait_pairs": len(new_pairs),
                "new_species_any_trait": len(new_species),
                "cumulative_species_trait_pairs": len(cumulative),
                "cumulative_species_any_trait": len({species for species, _ in cumulative}),
            }
        )
    marginal = pd.DataFrame(marginal_rows)

    summary = {
        "n_species_denominator": len(denominator),
        "n_lanes": len(lanes),
        "lane_order": sequence,
        "union_species_any_trait": len({species for species, _ in union_all}),
        "union_species_trait_pairs": len(union_all),
        "policy": (
            "Coverage is evidence recovery only. Missing rows remain unresolved and do not imply "
            "biological absence. Marginal gain counts only newly recovered species-trait pairs."
        ),
    }
    return coverage, marginal, summary


def parse_lane_specs(specs: list[str]) -> dict[str, Path]:
    lanes: dict[str, Path] = {}
    for spec in specs:
        if "=" not in spec:
            raise typer.BadParameter("lane must be NAME=PATH")
        name, raw_path = spec.split("=", 1)
        name = name.strip()
        path = Path(raw_path.strip())
        if not name or not path.exists():
            raise typer.BadParameter(f"invalid lane spec: {spec}")
        if name in lanes:
            raise typer.BadParameter(f"duplicate lane name: {name}")
        lanes[name] = path
    return lanes


@app.command("run")
def run_command(
    species_csv: Path = typer.Option(..., exists=True),
    lane: list[str] = typer.Option(..., "--lane", help="Repeat NAME=PATH for each method lane."),
    output_dir: Path = typer.Option(...),
    order: str = typer.Option("", help="Comma-separated lane order for marginal-gain analysis."),
    species_column: str = typer.Option(SPECIES_COLUMN),
) -> None:
    species = read_species(species_csv, species_column)
    lane_paths = parse_lane_specs(lane)
    denominator = set(species)
    lanes = {name: read_lane(path, denominator) for name, path in lane_paths.items()}
    sequence = [part.strip() for part in order.split(",") if part.strip()] or list(lanes)
    coverage, marginal, summary = benchmark_lanes(species=species, lanes=lanes, order=sequence)
    output_dir.mkdir(parents=True, exist_ok=True)
    coverage.to_csv(output_dir / "trait_recovery_coverage.csv", index=False)
    marginal.to_csv(output_dir / "trait_recovery_marginal_gain.csv", index=False)
    (output_dir / "trait_recovery_benchmark.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    app()
