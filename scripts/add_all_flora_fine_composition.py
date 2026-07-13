#!/usr/bin/env python3
"""Add all-flora fine colour/form composition counts to the model table.

This deliberately differs from the existing animal-pollinated fine-composition
columns, which are retained for the northern Bombus mechanism analysis.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import duckdb

FINE_COLORS = ("white", "blue", "yellow", "red", "green", "cream", "rare_other")
FINE_FORMS = (
    "tubular",
    "composite_head",
    "star",
    "bell_campanulate",
    "open_radial",
    "trumpet_funnel",
    "rare_other",
)


def sql_quote(path: Path) -> str:
    return str(path).replace("'", "''")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-csv", type=Path, required=True)
    parser.add_argument("--species-states-parquet", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    args = parser.parse_args()

    model = sql_quote(args.model_csv)
    states = sql_quote(args.species_states_parquet)
    out = sql_quote(args.output_csv)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)

    color_cols = ",\n      ".join(
        f"sum(CASE WHEN color_fine_cat='{c}' THEN 1 ELSE 0 END) AS all_fine_color_{c}"
        for c in FINE_COLORS[:-1]
    )
    form_cols = ",\n      ".join(
        f"sum(CASE WHEN form_fine_cat='{f}' THEN 1 ELSE 0 END) AS all_fine_form_{f}"
        for f in FINE_FORMS[:-1]
    )
    known_colors = ",".join(f"'{c}'" for c in FINE_COLORS[:-1])
    known_forms = ",".join(f"'{f}'" for f in FINE_FORMS[:-1])

    con = duckdb.connect()
    con.execute(
        f"""
        CREATE TEMP TABLE all_flora AS
        SELECT
          island_id,
          sum(CASE WHEN color_fine_cat IS NOT NULL THEN 1 ELSE 0 END) AS all_fine_color_trials,
          {color_cols},
          sum(CASE WHEN color_fine_cat IS NOT NULL AND color_fine_cat NOT IN ({known_colors}) THEN 1 ELSE 0 END)
            AS all_fine_color_rare_other,
          sum(CASE WHEN form_fine_cat IS NOT NULL THEN 1 ELSE 0 END) AS all_fine_form_trials,
          {form_cols},
          sum(CASE WHEN form_fine_cat IS NOT NULL AND form_fine_cat NOT IN ({known_forms}) THEN 1 ELSE 0 END)
            AS all_fine_form_rare_other
        FROM read_parquet('{states}')
        GROUP BY island_id
        """
    )
    con.execute(
        f"""
        COPY (
          SELECT m.*, a.* EXCLUDE (island_id)
          FROM read_csv_auto('{model}', union_by_name=true) m
          LEFT JOIN all_flora a USING (island_id)
        ) TO '{out}' (HEADER, DELIMITER ',')
        """
    )
    con.close()


if __name__ == "__main__":
    main()
