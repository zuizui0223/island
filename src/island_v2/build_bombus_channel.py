"""Build an auditable Bombus-channel table while preserving available components."""

from __future__ import annotations

import argparse
from pathlib import Path

import duckdb


def build_bombus_channel(covariates_csv: Path, output_parquet: Path) -> None:
    output_parquet.parent.mkdir(parents=True, exist_ok=True)
    src = str(covariates_csv).replace("'", "''")
    out = str(output_parquet).replace("'", "''")

    con = duckdb.connect()
    columns = {
        row[0]
        for row in con.execute(f"DESCRIBE SELECT * FROM read_csv_auto('{src}')").fetchall()
    }
    if "bombus_channel_state" in columns:
        channel_expr = "bombus_channel_state"
    elif "environmental_compatibility_max" in columns:
        channel_expr = "environmental_compatibility_max"
    else:
        raise ValueError(
            "Covariates must contain bombus_channel_state or environmental_compatibility_max"
        )

    environmental_expr = (
        "environmental_compatibility_max"
        if "environmental_compatibility_max" in columns
        else "NULL::DOUBLE"
    )
    con.execute(
        f"""
        COPY (
          SELECT
            *,
            {environmental_expr} AS bombus_environmental_compatibility,
            least(1.0, greatest(0.0, {channel_expr})) AS bombus_channel_state_audited,
            1.0 - least(1.0, greatest(0.0, {channel_expr})) AS bombus_deficit
          FROM read_csv_auto('{src}')
        ) TO '{out}' (FORMAT PARQUET, COMPRESSION ZSTD)
        """
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--island-covariates-csv", type=Path, required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    args = parser.parse_args()
    build_bombus_channel(args.island_covariates_csv, args.output_parquet)


if __name__ == "__main__":
    main()
