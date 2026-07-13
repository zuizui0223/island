#!/usr/bin/env python3
"""Audit fine-grained colour/form categories before promoting them to models."""

from __future__ import annotations

import argparse
from pathlib import Path

import duckdb


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-states-parquet", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    src = str(args.species_states_parquet).replace("'", "''")
    out = str(args.output_dir).replace("'", "''")
    con = duckdb.connect()

    con.execute(
        f"""
        COPY (
          SELECT 'flower_color_fine' AS endpoint, color_fine_cat AS category,
                 count(*) AS island_species_records,
                 count(DISTINCT accepted_species) AS n_species,
                 count(DISTINCT island_id) AS n_islands
          FROM read_parquet('{src}')
          WHERE is_animal = 1 AND color_fine_cat IS NOT NULL
          GROUP BY color_fine_cat
          UNION ALL
          SELECT 'floral_form_fine', form_fine_cat,
                 count(*), count(DISTINCT accepted_species), count(DISTINCT island_id)
          FROM read_parquet('{src}')
          WHERE is_animal = 1 AND form_fine_cat IS NOT NULL
          GROUP BY form_fine_cat
          ORDER BY endpoint, island_species_records DESC
        ) TO '{out}/fine_category_support.csv' (HEADER, DELIMITER ',')
        """
    )

    con.execute(
        f"""
        COPY (
          SELECT color_fine_cat, form_fine_cat,
                 count(*) AS island_species_records,
                 count(DISTINCT accepted_species) AS n_species,
                 count(DISTINCT island_id) AS n_islands
          FROM read_parquet('{src}')
          WHERE is_animal = 1
            AND color_fine_cat IS NOT NULL
            AND form_fine_cat IS NOT NULL
          GROUP BY color_fine_cat, form_fine_cat
          ORDER BY island_species_records DESC
        ) TO '{out}/fine_joint_color_form_support.csv' (HEADER, DELIMITER ',')
        """
    )

    con.execute(
        f"""
        COPY (
          SELECT
            count(*) AS island_species_records,
            sum(flower_color_conflict) AS flower_color_conflicts,
            sum(floral_form_conflict) AS floral_form_conflicts,
            sum(pollen_vector_conflict) AS pollen_vector_conflicts,
            sum(self_incompatibility_conflict) AS self_incompatibility_conflicts,
            sum(pollinator_guild_conflict) AS pollinator_guild_conflicts,
            sum(CASE WHEN color_fine_cat IS NULL AND flower_primary_color IS NOT NULL THEN 1 ELSE 0 END)
              AS unresolved_nonnull_color_labels,
            sum(CASE WHEN form_fine_cat IS NULL AND floral_form IS NOT NULL THEN 1 ELSE 0 END)
              AS unresolved_nonnull_form_labels
          FROM read_parquet('{src}')
        ) TO '{out}/fine_category_conflict_audit.csv' (HEADER, DELIMITER ',')
        """
    )
    con.close()


if __name__ == "__main__":
    main()
