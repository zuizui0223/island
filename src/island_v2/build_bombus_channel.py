"""Build an auditable Bombus-channel table while preserving component meaning.

The analysis-ready ``bombus_deficit`` column is retained for compatibility, but its
provenance is now explicit. A true occurrence-informed channel is preferred when
available; environmental compatibility alone is labelled as an environmental
mismatch proxy and must not be interpreted as observed Bombus absence.
"""

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

    environmental_expr = (
        "environmental_compatibility_max"
        if "environmental_compatibility_max" in columns
        else "NULL::DOUBLE"
    )

    # Precedence is deliberately fail-safe:
    # 1. use a precomputed occurrence-informed combined deficit when present;
    # 2. otherwise combine environmental compatibility and observation availability;
    # 3. otherwise use a declared primary channel state only on rows actually marked
    #    as primary Bombus-channel analyses;
    # 4. fall back to environmental mismatch, explicitly labelled as a proxy.
    if "bombus_channel_deficit" in columns:
        availability_expr = "1.0 - least(1.0, greatest(0.0, bombus_channel_deficit))"
        kind_expr = "'occurrence_informed_combined_channel'"
    elif {
        "observation_availability",
        "environmental_compatibility_max",
    }.issubset(columns):
        availability_expr = (
            "sqrt(least(1.0, greatest(0.0, environmental_compatibility_max)) * "
            "least(1.0, greatest(0.0, observation_availability)))"
        )
        kind_expr = "'occurrence_informed_combined_channel'"
    elif {
        "bombus_channel_state",
        "primary_bombus_channel_analysis",
        "environmental_compatibility_max",
    }.issubset(columns):
        availability_expr = (
            "CASE WHEN coalesce(primary_bombus_channel_analysis, FALSE) "
            "THEN least(1.0, greatest(0.0, bombus_channel_state)) "
            "ELSE least(1.0, greatest(0.0, environmental_compatibility_max)) END"
        )
        kind_expr = (
            "CASE WHEN coalesce(primary_bombus_channel_analysis, FALSE) "
            "THEN 'declared_primary_channel' ELSE 'environmental_mismatch_proxy' END"
        )
    elif "environmental_compatibility_max" in columns:
        availability_expr = "least(1.0, greatest(0.0, environmental_compatibility_max))"
        kind_expr = "'environmental_mismatch_proxy'"
    else:
        raise ValueError(
            "Covariates must contain an occurrence-informed Bombus channel or "
            "environmental_compatibility_max"
        )

    con.execute(
        f"""
        COPY (
          SELECT
            *,
            {environmental_expr} AS bombus_environmental_compatibility,
            {availability_expr} AS bombus_channel_state_audited,
            1.0 - ({availability_expr}) AS bombus_deficit,
            {kind_expr} AS bombus_predictor_kind,
            CASE
              WHEN ({kind_expr}) = 'environmental_mismatch_proxy'
                THEN 'Bombus environmental mismatch proxy; not observed absence'
              ELSE 'Occurrence-informed or declared Bombus channel deficit'
            END AS bombus_deficit_interpretation
          FROM read_csv_auto('{src}')
        ) TO '{out}' (FORMAT PARQUET, COMPRESSION ZSTD)
        """
    )
    con.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--island-covariates-csv", type=Path, required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    args = parser.parse_args()
    build_bombus_channel(args.island_covariates_csv, args.output_parquet)


if __name__ == "__main__":
    main()
