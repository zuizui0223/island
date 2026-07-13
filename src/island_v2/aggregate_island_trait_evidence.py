"""Aggregate conflict-safe species trait states to island-level evidence tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import duckdb


COLORS = ("plain", "yellow_orange", "red_pink", "blue_purple")
FORMS = (
    "open_generalized",
    "tubular_trumpet",
    "zygomorphic_specialized",
    "composite_brush",
)


def aggregate_island_trait_evidence(
    species_states_parquet: Path,
    output_long_parquet: Path,
    output_wide_parquet: Path,
) -> None:
    output_long_parquet.parent.mkdir(parents=True, exist_ok=True)
    output_wide_parquet.parent.mkdir(parents=True, exist_ok=True)
    src = str(species_states_parquet).replace("'", "''")
    out_long = str(output_long_parquet).replace("'", "''")
    out_wide = str(output_wide_parquet).replace("'", "''")

    color_columns = ",\n          ".join(
        f"sum(CASE WHEN is_animal=1 AND color_cat='{cat}' THEN 1 ELSE 0 END) AS color_{cat}"
        for cat in COLORS
    )
    form_columns = ",\n          ".join(
        f"sum(CASE WHEN is_animal=1 AND form_cat='{cat}' THEN 1 ELSE 0 END) AS form_{cat}"
        for cat in FORMS
    )
    joint_columns = ",\n          ".join(
        f"sum(CASE WHEN is_animal=1 AND color_cat='{color}' AND form_cat='{form}' "
        f"THEN 1 ELSE 0 END) AS joint_{color}__{form}"
        for color in COLORS
        for form in FORMS
    )

    con = duckdb.connect()
    con.execute(
        f"""
        CREATE TEMP TABLE island_wide AS
        SELECT
          island_id,
          count(*) AS n_species_total,
          count(is_animal) AS animal_status_observed,
          sum(CASE WHEN is_animal=1 THEN 1 ELSE 0 END) AS n_animal_species,
          sum(CASE WHEN is_wind=1 THEN 1 ELSE 0 END) AS n_wind_species,
          sum(CASE WHEN is_mixed=1 THEN 1 ELSE 0 END) AS n_mixed_species,
          sum(CASE WHEN is_wind=1 OR is_mixed=1 THEN 1 ELSE 0 END) AS n_wind_mixed_species,
          sum(CASE WHEN any_trait_conflict THEN 1 ELSE 0 END) AS n_species_with_any_trait_conflict,
          sum(pollen_vector_conflict) AS n_pollen_vector_conflicts,
          sum(flower_color_conflict) AS n_flower_color_conflicts,
          sum(floral_form_conflict) AS n_floral_form_conflicts,
          sum(self_incompatibility_conflict) AS n_self_incompatibility_conflicts,
          sum(pollinator_guild_conflict) AS n_pollinator_guild_conflicts,
          sum(CASE WHEN is_animal=1 AND color_cat IS NOT NULL THEN 1 ELSE 0 END) AS color_trials,
          {color_columns},
          sum(CASE WHEN is_animal=1 AND form_cat IS NOT NULL THEN 1 ELSE 0 END) AS form_trials,
          {form_columns},
          sum(CASE WHEN is_animal=1 AND color_cat IS NOT NULL AND form_cat IS NOT NULL
                   THEN 1 ELSE 0 END) AS joint_color_form_trials,
          {joint_columns},
          sum(CASE WHEN is_animal=1 AND sc_binary IS NOT NULL THEN 1 ELSE 0 END) AS sc_trials,
          sum(CASE WHEN is_animal=1 AND sc_binary=1 THEN 1 ELSE 0 END) AS sc_successes,
          avg(CASE WHEN is_animal=1 THEN showy_alt END) AS showy_alt_guild_share,
          avg(CASE WHEN is_animal=1 THEN other_bee END) AS other_bee_guild_share,
          avg(CASE WHEN is_animal=1 THEN generalist_insect END) AS generalist_insect_guild_share
        FROM read_parquet('{src}')
        GROUP BY island_id
        """
    )
    con.execute(f"COPY island_wide TO '{out_wide}' (FORMAT PARQUET, COMPRESSION ZSTD)")

    long_parts: list[str] = []
    for category in COLORS:
        long_parts.append(
            f"SELECT island_id, 'flower_color' AS endpoint, '{category}' AS category, "
            f"color_{category} AS successes, color_trials AS trials FROM island_wide"
        )
    for category in FORMS:
        long_parts.append(
            f"SELECT island_id, 'floral_form', '{category}', "
            f"form_{category}, form_trials FROM island_wide"
        )
    for color in COLORS:
        for form in FORMS:
            category = f"{color}__{form}"
            long_parts.append(
                f"SELECT island_id, 'joint_color_form', '{category}', "
                f"joint_{category}, joint_color_form_trials FROM island_wide"
            )
    long_parts.append(
        "SELECT island_id, 'self_compatibility', 'self_compatible', "
        "sc_successes, sc_trials FROM island_wide"
    )
    union_sql = "\nUNION ALL\n".join(long_parts)
    con.execute(
        f"COPY ({union_sql}) TO '{out_long}' (FORMAT PARQUET, COMPRESSION ZSTD)"
    )
    con.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-states-parquet", type=Path, required=True)
    parser.add_argument("--output-long-parquet", type=Path, required=True)
    parser.add_argument("--output-wide-parquet", type=Path, required=True)
    args = parser.parse_args()
    aggregate_island_trait_evidence(
        args.species_states_parquet,
        args.output_long_parquet,
        args.output_wide_parquet,
    )


if __name__ == "__main__":
    main()
