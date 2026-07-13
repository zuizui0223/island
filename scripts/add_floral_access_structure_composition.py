#!/usr/bin/env python3
"""Add floral access-structure compositions for all flora and animal-pollinated flora.

The access ontology collapses visually heterogeneous fine-form labels onto a single
functional axis: how a visitor gains access to floral rewards/reproductive organs.
It preserves unresolved/rare states as ``other`` rather than silently discarding them.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import duckdb

ACCESS_CLASSES = (
    "open_accessible",
    "intermediate_bell_funnel",
    "restricted_tubular",
    "specialized_restricted",
    "multi_flower_display",
    "other",
)

# Mapping from conflict-safe fine floral-form states to a functional access structure.
# Fine labels not explicitly mapped are retained as ``other``.
ACCESS_CASE = """
CASE
  WHEN form_fine_cat IN ('open_radial', 'bowl_cup', 'star') THEN 'open_accessible'
  WHEN form_fine_cat IN ('bell_campanulate', 'trumpet_funnel', 'salver') THEN 'intermediate_bell_funnel'
  WHEN form_fine_cat IN ('tubular', 'urn_urceolate') THEN 'restricted_tubular'
  WHEN form_fine_cat IN ('zygomorphic_unspecified', 'orchid', 'papilionoid', 'spurred', 'bilabiate') THEN 'specialized_restricted'
  WHEN form_fine_cat IN ('composite_head', 'brush_puff') THEN 'multi_flower_display'
  WHEN form_fine_cat IS NOT NULL THEN 'other'
  ELSE NULL
END
""".strip()


def add_access_composition(input_csv: Path, species_states_parquet: Path, output_csv: Path) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    inp = str(input_csv).replace("'", "''")
    states = str(species_states_parquet).replace("'", "''")
    out = str(output_csv).replace("'", "''")

    all_cols = ",\n        ".join(
        f"sum(CASE WHEN access_class='{c}' THEN 1 ELSE 0 END) AS all_access_{c}" for c in ACCESS_CLASSES
    )
    animal_cols = ",\n        ".join(
        f"sum(CASE WHEN is_animal=1 AND access_class='{c}' THEN 1 ELSE 0 END) AS animal_access_{c}"
        for c in ACCESS_CLASSES
    )

    con = duckdb.connect()
    con.execute(
        f"""
        CREATE TEMP TABLE access_species AS
        SELECT *, {ACCESS_CASE} AS access_class
        FROM read_parquet('{states}')
        """
    )
    con.execute(
        f"""
        CREATE TEMP TABLE access_by_island AS
        SELECT
          island_id,
          sum(CASE WHEN access_class IS NOT NULL THEN 1 ELSE 0 END) AS all_access_trials,
          {all_cols},
          sum(CASE WHEN is_animal=1 AND access_class IS NOT NULL THEN 1 ELSE 0 END) AS animal_access_trials,
          {animal_cols}
        FROM access_species
        GROUP BY island_id
        """
    )
    con.execute(
        f"""
        COPY (
          SELECT b.*, a.* EXCLUDE (island_id)
          FROM read_csv_auto('{inp}') b
          LEFT JOIN access_by_island a USING (island_id)
        ) TO '{out}' (HEADER, DELIMITER ',')
        """
    )
    con.close()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-csv", type=Path, required=True)
    p.add_argument("--species-states-parquet", type=Path, required=True)
    p.add_argument("--output-csv", type=Path, required=True)
    args = p.parse_args()
    add_access_composition(args.input_csv, args.species_states_parquet, args.output_csv)


if __name__ == "__main__":
    main()
