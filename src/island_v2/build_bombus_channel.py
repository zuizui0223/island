"""Build an auditable Bombus-channel table while preserving component meaning.

The analysis-ready ``bombus_deficit`` is fail-closed. A global all-species maximum
is retained only as a diagnostic environmental screen; it is never promoted to a
Bombus channel predictor unless the environmental score is explicitly constrained
to an island-specific biogeographic source pool. Occurrence-informed channel scores
remain preferred when available.
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

    # A source-pool environmental score is accepted only when the upstream product
    # explicitly says so, or when candidate-count columns prove that the scored pool
    # is a strict subset of the global Bombus set. This prevents the former 169-species
    # global maximum from masquerading as island-specific Bombus availability.
    if "environmental_source_pool_constrained" in columns:
        constrained_expr = "coalesce(environmental_source_pool_constrained, FALSE)"
    elif {
        "n_source_pool_bombus_species",
        "n_bombus_species_total",
    }.issubset(columns):
        constrained_expr = (
            "coalesce(n_source_pool_bombus_species, 0) > 0 AND "
            "coalesce(n_source_pool_bombus_species, 0) < coalesce(n_bombus_species_total, 0)"
        )
    else:
        constrained_expr = "FALSE"

    source_pool_env_available = (
        "environmental_compatibility_max" in columns
        and constrained_expr != "FALSE"
    )

    # Precedence:
    # 1. precomputed occurrence-informed combined channel;
    # 2. occurrence evidence + source-pool-constrained environmental compatibility;
    # 3. source-pool-constrained declared primary channel state;
    # 4. source-pool-constrained environmental compatibility alone.
    # Unconstrained global-max compatibility is diagnostic only and yields NULL for
    # the analysis-ready deficit.
    if "bombus_channel_deficit" in columns:
        availability_expr = (
            "1.0 - least(1.0, greatest(0.0, bombus_channel_deficit))"
        )
        kind_expr = "'occurrence_informed_combined_channel'"
    elif {
        "observation_availability",
        "environmental_compatibility_max",
    }.issubset(columns):
        availability_expr = (
            f"CASE WHEN ({constrained_expr}) THEN "
            "sqrt(least(1.0, greatest(0.0, environmental_compatibility_max)) * "
            "least(1.0, greatest(0.0, observation_availability))) "
            "ELSE least(1.0, greatest(0.0, observation_availability)) END"
        )
        kind_expr = (
            f"CASE WHEN ({constrained_expr}) "
            "THEN 'occurrence_informed_source_pool_channel' "
            "ELSE 'occurrence_only_channel' END"
        )
    elif {
        "bombus_channel_state",
        "primary_bombus_channel_analysis",
    }.issubset(columns):
        availability_expr = (
            f"CASE WHEN coalesce(primary_bombus_channel_analysis, FALSE) "
            f"AND ({constrained_expr}) "
            "THEN least(1.0, greatest(0.0, bombus_channel_state)) "
            "ELSE NULL::DOUBLE END"
        )
        kind_expr = (
            f"CASE WHEN coalesce(primary_bombus_channel_analysis, FALSE) "
            f"AND ({constrained_expr}) "
            "THEN 'declared_source_pool_channel' "
            "ELSE 'global_environmental_screen_only' END"
        )
    elif source_pool_env_available:
        availability_expr = (
            f"CASE WHEN ({constrained_expr}) "
            "THEN least(1.0, greatest(0.0, environmental_compatibility_max)) "
            "ELSE NULL::DOUBLE END"
        )
        kind_expr = (
            f"CASE WHEN ({constrained_expr}) "
            "THEN 'source_pool_environmental_compatibility' "
            "ELSE 'global_environmental_screen_only' END"
        )
    elif "environmental_compatibility_max" in columns:
        availability_expr = "NULL::DOUBLE"
        kind_expr = "'global_environmental_screen_only'"
    else:
        raise ValueError(
            "Covariates must contain an occurrence-informed Bombus channel or "
            "environmental compatibility diagnostics"
        )

    con.execute(
        f"""
        COPY (
          SELECT
            *,
            {environmental_expr} AS bombus_environmental_compatibility,
            ({constrained_expr}) AS bombus_environmental_source_pool_constrained,
            {availability_expr} AS bombus_channel_state_audited,
            CASE
              WHEN ({availability_expr}) IS NULL THEN NULL::DOUBLE
              ELSE 1.0 - ({availability_expr})
            END AS bombus_deficit,
            {kind_expr} AS bombus_predictor_kind,
            CASE
              WHEN ({kind_expr}) = 'global_environmental_screen_only'
                THEN 'Global all-species environmental maximum retained for diagnostics only; not a Bombus absence or channel predictor'
              WHEN ({kind_expr}) = 'source_pool_environmental_compatibility'
                THEN 'Island-specific source-pool environmental mismatch; not direct observed absence'
              ELSE 'Occurrence-informed or source-pool-constrained Bombus channel deficit'
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
