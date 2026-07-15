"""Build an exact-island public-data snapshot for the six frozen Izu pilot islands."""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

IZU_ISLANDS = {
    "gshhg_2.3.7_h_3f48d601bf60ade348ae": "Izu Oshima",
    "gshhg_2.3.7_h_80bd057071c9043cebcc": "Niijima",
    "gshhg_2.3.7_h_ab48ddfcdcd4ccc62870": "Kozushima",
    "gshhg_2.3.7_h_cc2a3ac665156bf28969": "Miyakejima",
    "gshhg_2.3.7_h_78b9099b78b447bab1fd": "Mikurajima",
    "gshhg_2.3.7_h_79143661e5762bf8fa25": "Hachijojima",
}


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str).fillna("")


def _integer(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(0).astype(int)


def build(species_snapshot: Path, effort_csv: Path, output_dir: Path) -> dict[str, object]:
    species = _read(species_snapshot)
    effort = _read(effort_csv)
    required_species = {"island_id", "species", "n_records", "n_unique_gbif_ids"}
    required_effort = {"island_id", "n_records", "n_species", "n_datasets"}
    if not required_species.issubset(species.columns):
        raise ValueError(f"species snapshot missing {sorted(required_species - set(species.columns))}")
    if not required_effort.issubset(effort.columns):
        raise ValueError(f"effort table missing {sorted(required_effort - set(effort.columns))}")

    subset = species.loc[species["island_id"].isin(IZU_ISLANDS)].copy()
    subset["island_name"] = subset["island_id"].map(IZU_ISLANDS)
    subset["n_records"] = _integer(subset["n_records"])
    subset["n_unique_gbif_ids"] = _integer(subset["n_unique_gbif_ids"])
    subset = subset.sort_values(["island_name", "species"]).reset_index(drop=True)

    island_effort = effort.loc[effort["island_id"].isin(IZU_ISLANDS)].copy()
    island_effort["island_name"] = island_effort["island_id"].map(IZU_ISLANDS)
    for column in ["n_records", "n_species", "n_datasets"]:
        island_effort[column] = _integer(island_effort[column])
    island_effort = island_effort.sort_values("island_name").reset_index(drop=True)

    incidence_rows = []
    for species_name, group in subset.groupby("species", sort=True):
        islands = sorted(group["island_name"].unique())
        incidence_rows.append(
            {
                "species": species_name,
                "n_islands": len(islands),
                "islands": "|".join(islands),
                "n_records": int(group["n_records"].sum()),
                "n_unique_gbif_ids": int(group["n_unique_gbif_ids"].sum()),
            }
        )
    incidence = pd.DataFrame(incidence_rows)

    island_sets = {
        island: set(group["species"])
        for island, group in subset.groupby("island_name")
    }
    pair_rows = []
    for first, second in combinations(sorted(IZU_ISLANDS.values()), 2):
        a = island_sets.get(first, set())
        b = island_sets.get(second, set())
        union = a | b
        pair_rows.append(
            {
                "island_a": first,
                "island_b": second,
                "n_a": len(a),
                "n_b": len(b),
                "n_shared": len(a & b),
                "jaccard": (len(a & b) / len(union)) if union else 0.0,
            }
        )
    pairwise = pd.DataFrame(pair_rows)

    output_dir.mkdir(parents=True, exist_ok=True)
    subset.to_csv(output_dir / "izu_island_species.csv.gz", index=False, compression="gzip")
    island_effort.to_csv(output_dir / "izu_island_effort.csv", index=False)
    incidence.to_csv(output_dir / "izu_species_incidence.csv", index=False)
    pairwise.to_csv(output_dir / "izu_pairwise_jaccard.csv", index=False)

    missing = sorted(set(IZU_ISLANDS) - set(island_effort["island_id"]))
    summary = {
        "scope": "six frozen exact-island Izu pilot polygons",
        "n_islands_expected": len(IZU_ISLANDS),
        "n_islands_with_effort_rows": int(island_effort["island_id"].nunique()),
        "missing_island_ids": missing,
        "n_island_species_pairs": int(len(subset)),
        "n_unique_species_labels": int(subset["species"].nunique()),
        "n_total_records": int(subset["n_records"].sum()),
        "n_species_on_all_six_islands": int((incidence["n_islands"] == 6).sum()) if not incidence.empty else 0,
        "n_single_island_species_labels": int((incidence["n_islands"] == 1).sum()) if not incidence.empty else 0,
        "island_rows": island_effort[["island_name", "n_records", "n_species", "n_datasets"]].to_dict("records"),
        "limitations": [
            "Occurrence-based species labels are not a reviewed native flora.",
            "Cultivated, introduced, transient, and taxonomically unresolved records remain possible.",
            "Toshima, Shikinejima, and Aogashima are outside this six-polygon pilot output and require a supplemental exact-polygon acquisition.",
        ],
    }
    (output_dir / "izu_public_data_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return summary


@app.command("build")
def build_command(
    species_snapshot: Path = typer.Option(..., exists=True),
    effort_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    summary = build(species_snapshot, effort_csv, output_dir)
    typer.echo(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    app()
