"""Aggregate species-level trait states to island-level evidence tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import duckdb


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

    con = duckdb.connect()
    con.execute(
        f"""
        CREATE TEMP TABLE island_wide AS
        SELECT
          island_id,
          count(*) AS n_species_total,
          count(is_animal) AS animal_status_observed,
          sum(CASE WHEN is_animal=1 THEN 1 ELSE 0 END) AS n_animal_species,
          sum(CASE WHEN pollen_vector_mode='abiotic_wind' THEN 1 ELSE 0 END) AS n_wind_species,
          sum(CASE WHEN is_animal=1 AND color_cat IS NOT NULL THEN 1 ELSE 0 END) AS color_trials,
          sum(CASE WHEN is_animal=1 AND color_cat='plain' THEN 1 ELSE 0 END) AS color_plain,
          sum(CASE WHEN is_animal=1 AND color_cat='yellow_orange' THEN 1 ELSE 0 END) AS color_yellow_orange,
          sum(CASE WHEN is_animal=1 AND color_cat='red_pink' THEN 1 ELSE 0 END) AS color_red_pink,
          sum(CASE WHEN is_animal=1 AND color_cat='blue_purple' THEN 1 ELSE 0 END) AS color_blue_purple,
          sum(CASE WHEN is_animal=1 AND form_cat IS NOT NULL THEN 1 ELSE 0 END) AS form_trials,
          sum(CASE WHEN is_animal=1 AND form_cat='open_generalized' THEN 1 ELSE 0 END) AS form_open_generalized,
          sum(CASE WHEN is_animal=1 AND form_cat='tubular_trumpet' THEN 1 ELSE 0 END) AS form_tubular_trumpet,
          sum(CASE WHEN is_animal=1 AND form_cat='zygomorphic_specialized' THEN 1 ELSE 0 END) AS form_zygomorphic_specialized,
          sum(CASE WHEN is_animal=1 AND form_cat='composite_brush' THEN 1 ELSE 0 END) AS form_composite_brush,
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
    con.execute(
        f"""
        COPY (
          SELECT island_id, 'flower_color' AS endpoint, 'plain' AS category,
                 color_plain AS successes, color_trials AS trials FROM island_wide
          UNION ALL
          SELECT island_id, 'flower_color', 'yellow_orange', color_yellow_orange, color_trials FROM island_wide
          UNION ALL
          SELECT island_id, 'flower_color', 'red_pink', color_red_pink, color_trials FROM island_wide
          UNION ALL
          SELECT island_id, 'flower_color', 'blue_purple', color_blue_purple, color_trials FROM island_wide
          UNION ALL
          SELECT island_id, 'floral_form', 'open_generalized', form_open_generalized, form_trials FROM island_wide
          UNION ALL
          SELECT island_id, 'floral_form', 'tubular_trumpet', form_tubular_trumpet, form_trials FROM island_wide
          UNION ALL
          SELECT island_id, 'floral_form', 'zygomorphic_specialized', form_zygomorphic_specialized, form_trials FROM island_wide
          UNION ALL
          SELECT island_id, 'floral_form', 'composite_brush', form_composite_brush, form_trials FROM island_wide
          UNION ALL
          SELECT island_id, 'self_compatibility', 'self_compatible', sc_successes, sc_trials FROM island_wide
        ) TO '{out_long}' (FORMAT PARQUET, COMPRESSION ZSTD)
        """
    )


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
