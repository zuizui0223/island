"""Prepare auditable island-level inputs for the v2 Bayesian M0-M4 main analysis.

This compatibility entry point orchestrates four explicit layers:
1. species trait states with unknown preserved as unknown;
2. island-level long and wide evidence aggregation;
3. Bombus-channel construction with components retained;
4. final model table plus a machine-readable analysis contract.

The frozen island universe is never reduced here. Islands with
``distance_to_continent_km == 0`` remain in the analysis-ready table and are
explicitly flagged. Canonical primary models exclude them at model-support
selection time; sensitivity analyses may retain them.
"""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

import duckdb

from island_v2.aggregate_island_trait_evidence import aggregate_island_trait_evidence
from island_v2.build_analysis_contract import build_analysis_contract
from island_v2.build_bombus_channel import build_bombus_channel
from island_v2.build_species_trait_states import build_species_trait_states
from island_v2.validate_analysis_inputs import validate_analysis_inputs


def prepare_main_input(
    master_csv_gz: Path,
    island_covariates_csv: Path,
    output_csv: Path,
    analysis_tier: str,
) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    audit_dir = output_csv.parent / "input_audit"
    audit_dir.mkdir(parents=True, exist_ok=True)

    # Persist species-level states so ontology coverage, conflicts, and fine
    # category support can be audited after the preparation job completes.
    species_states = audit_dir / "species_trait_states.parquet"

    with tempfile.TemporaryDirectory(prefix="island-v2-main-") as tmp:
        tmpdir = Path(tmp)
        evidence_wide = tmpdir / "island_trait_evidence_wide.parquet"
        evidence_long = audit_dir / "island_trait_evidence_long.parquet"
        bombus_channel = audit_dir / "bombus_channel_components.parquet"

        build_species_trait_states(master_csv_gz, species_states)
        aggregate_island_trait_evidence(species_states, evidence_long, evidence_wide)
        build_bombus_channel(island_covariates_csv, bombus_channel)

        out = str(output_csv).replace("'", "''")
        traits = str(evidence_wide).replace("'", "''")
        cov = str(bombus_channel).replace("'", "''")
        tier = analysis_tier.replace("'", "''")
        con = duckdb.connect()
        con.execute(
            f"""
            COPY (
              SELECT
                c.*,
                t.* EXCLUDE (island_id),
                '{tier}' AS analysis_tier,
                CASE
                  WHEN c.distance_to_continent_km = 0 THEN TRUE
                  ELSE FALSE
                END AS distance_zero_flag,
                CASE
                  WHEN c.distance_to_continent_km > 0 THEN TRUE
                  ELSE FALSE
                END AS primary_distance_support,
                CASE WHEN t.color_trials > 0
                  THEN 1.0*t.color_plain/t.color_trials END AS plain_share,
                CASE WHEN t.form_trials > 0
                  THEN 1.0*t.form_open_generalized/t.form_trials END AS generalized_share,
                CASE WHEN t.sc_trials > 0
                  THEN 1.0*t.sc_successes/t.sc_trials END AS sc_share
              FROM read_parquet('{cov}') c
              LEFT JOIN read_parquet('{traits}') t USING (island_id)
            ) TO '{out}' (HEADER, DELIMITER ',')
            """
        )

        zero_audit = str(audit_dir / "distance_zero_islands.csv").replace("'", "''")
        con.execute(
            f"""
            COPY (
              SELECT *
              FROM read_csv_auto('{out}')
              WHERE distance_to_continent_km = 0
              ORDER BY island_id
            ) TO '{zero_audit}' (HEADER, DELIMITER ',')
            """
        )
        con.close()

    contract = Path("config/v2_analysis_contract.yml")
    if contract.exists():
        build_analysis_contract(contract, audit_dir / "analysis_contract.json")
    validate_analysis_inputs(output_csv)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv-gz", type=Path, required=True)
    parser.add_argument("--island-covariates-csv", type=Path, required=True)
    parser.add_argument("--analysis-tier", default="primary")
    parser.add_argument("--output-csv", type=Path, required=True)
    args = parser.parse_args()
    prepare_main_input(
        args.master_csv_gz,
        args.island_covariates_csv,
        args.output_csv,
        args.analysis_tier,
    )


if __name__ == "__main__":
    main()
