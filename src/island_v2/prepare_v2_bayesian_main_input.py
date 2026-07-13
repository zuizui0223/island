"""Prepare compact island-level inputs for the preregistered v2 M0-M4 SEM ladder.

The integrated island x species x trait master is large. DuckDB pivots only the
traits required by config/m0_m4_sem.yml and aggregates them to one row per island.
The script is evidence-tier agnostic: run it separately for primary, broad, and
sensitivity_all inputs.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import duckdb


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv-gz", type=Path, required=True)
    parser.add_argument("--island-covariates-csv", type=Path, required=True)
    parser.add_argument("--analysis-tier", required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    args = parser.parse_args()
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect()
    master = str(args.master_csv_gz).replace("'", "''")
    cov = str(args.island_covariates_csv).replace("'", "''")
    out = str(args.output_csv).replace("'", "''")
    tier = args.analysis_tier.replace("'", "''")

    con.execute(
        f"""
        CREATE TEMP TABLE species_traits AS
        SELECT
          island_id,
          accepted_species,
          max(CASE WHEN trait_name='pollen_vector_mode' THEN filled_value END) AS pollen_vector_mode,
          max(CASE WHEN trait_name='pollination_functional_guild' THEN filled_value END) AS pollination_functional_guild,
          max(CASE WHEN trait_name='self_incompatibility' THEN filled_value END) AS self_incompatibility,
          max(CASE WHEN trait_name='autonomous_selfing_capacity' THEN filled_value END) AS autonomous_selfing_capacity
        FROM read_csv_auto('{master}', union_by_name=true)
        WHERE trait_name IN (
          'pollen_vector_mode',
          'pollination_functional_guild',
          'self_incompatibility',
          'autonomous_selfing_capacity'
        )
        GROUP BY island_id, accepted_species;
        """
    )

    con.execute(
        """
        CREATE TEMP TABLE classified AS
        SELECT *,
          CASE
            WHEN pollen_vector_mode='biotic' THEN 1
            WHEN pollen_vector_mode IS NULL THEN NULL
            ELSE 0
          END AS animal_binary,
          CASE
            WHEN self_incompatibility IN ('SC','likely_SC') THEN 1
            WHEN self_incompatibility IN ('SI','likely_SI','obligate_SI') THEN 0
            ELSE NULL
          END AS sc_binary,
          CASE
            WHEN autonomous_selfing_capacity='autonomous' THEN 1
            WHEN autonomous_selfing_capacity IS NULL THEN NULL
            ELSE 0
          END AS autonomous_binary,
          CASE
            WHEN pollination_functional_guild='bumblebees' THEN 1
            WHEN pollination_functional_guild IS NULL THEN NULL
            ELSE 0
          END AS floral_proxy_binary,
          CASE
            WHEN pollination_functional_guild IN ('other_bees','flies','birds','bats','wind','mixed') THEN 1
            WHEN pollination_functional_guild IS NULL THEN NULL
            ELSE 0
          END AS alternative_binary
        FROM species_traits;
        """
    )

    con.execute(
        """
        CREATE TEMP TABLE island_traits AS
        SELECT
          island_id,
          count(*) AS n_species_total,

          count(animal_binary) AS animal_trials,
          sum(CASE WHEN animal_binary=1 THEN 1 ELSE 0 END) AS animal_successes,

          count(sc_binary) AS sc_trials,
          sum(CASE WHEN sc_binary=1 THEN 1 ELSE 0 END) AS sc_successes,

          count(autonomous_binary) AS autonomous_trials,
          sum(CASE WHEN autonomous_binary=1 THEN 1 ELSE 0 END) AS autonomous_successes,

          count(floral_proxy_binary) AS floral_phenotype_trials,
          sum(CASE WHEN floral_proxy_binary=1 THEN 1 ELSE 0 END) AS floral_phenotype_successes,

          count(alternative_binary) AS alternative_trials,
          sum(CASE WHEN alternative_binary=1 THEN 1 ELSE 0 END) AS alternative_successes,

          count(sc_binary) + count(autonomous_binary) AS reproductive_assurance_trials,
          sum(CASE WHEN sc_binary=1 THEN 1 ELSE 0 END)
            + sum(CASE WHEN autonomous_binary=1 THEN 1 ELSE 0 END)
            AS reproductive_assurance_successes
        FROM classified
        GROUP BY island_id;
        """
    )

    con.execute(
        f"""
        COPY (
          SELECT
            c.*,
            t.* EXCLUDE (island_id),
            '{tier}' AS analysis_tier,
            c.environmental_compatibility_max AS bombus_channel_state,
            CASE
              WHEN t.animal_trials > 0
              THEN 1.0*t.animal_successes/t.animal_trials
            END AS animal_mediation_proxy,
            CASE
              WHEN t.sc_trials > 0
              THEN 1.0*t.sc_successes/t.sc_trials
            END AS self_compatibility_proxy,
            CASE
              WHEN t.autonomous_trials > 0
              THEN 1.0*t.autonomous_successes/t.autonomous_trials
            END AS autonomous_selfing_proxy,
            CASE
              WHEN t.floral_phenotype_trials > 0
              THEN 1.0*t.floral_phenotype_successes/t.floral_phenotype_trials
            END AS floral_phenotype_proxy,
            CASE
              WHEN t.reproductive_assurance_trials > 0
              THEN 1.0*t.reproductive_assurance_successes/t.reproductive_assurance_trials
            END AS reproductive_assurance_proxy,
            CASE
              WHEN t.alternative_trials > 0
              THEN least(1.0, 1.0*t.alternative_successes/t.alternative_trials)
            END AS alternative_functional_replacement_proxy
          FROM read_csv_auto('{cov}') c
          LEFT JOIN island_traits t USING (island_id)
        ) TO '{out}' (HEADER, DELIMITER ',');
        """
    )


if __name__ == "__main__":
    main()
