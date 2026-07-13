"""Prepare compact island-level inputs for the v2 Bayesian M0-M4 main analysis.

The source master is ~8.7M rows. DuckDB performs the species-level pivot and
island aggregation out-of-core so the Bayesian stage reads one compact table.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import duckdb


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv-gz", type=Path, required=True)
    parser.add_argument("--island-covariates-csv", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    args = parser.parse_args()
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect()
    master = str(args.master_csv_gz).replace("'", "''")
    cov = str(args.island_covariates_csv).replace("'", "''")
    out = str(args.output_csv).replace("'", "''")

    con.execute(
        f"""
        CREATE TEMP TABLE species_traits AS
        SELECT
          island_id,
          accepted_species,
          max(CASE WHEN trait_name='pollen_vector_mode' THEN filled_value END) AS pollen_vector_mode,
          max(CASE WHEN trait_name='flower_primary_color' THEN filled_value END) AS flower_primary_color,
          max(CASE WHEN trait_name='floral_form' THEN filled_value END) AS floral_form,
          max(CASE WHEN trait_name='self_incompatibility' THEN filled_value END) AS self_incompatibility,
          max(CASE WHEN trait_name='pollination_functional_guild' THEN filled_value END) AS pollination_functional_guild
        FROM read_csv_auto('{master}', union_by_name=true)
        WHERE trait_name IN (
          'pollen_vector_mode','flower_primary_color','floral_form',
          'self_incompatibility','pollination_functional_guild'
        )
        GROUP BY island_id, accepted_species;
        """
    )

    con.execute(
        """
        CREATE TEMP TABLE classified AS
        SELECT *,
          CASE
            WHEN pollen_vector_mode='biotic' THEN 1 ELSE 0
          END AS is_animal,
          CASE
            WHEN lower(coalesce(flower_primary_color,'')) ~ 'blue|purple|violet|lilac|lavender' THEN 'blue_purple'
            WHEN lower(coalesce(flower_primary_color,'')) ~ 'red|pink|magenta|scarlet|crimson' THEN 'red_pink'
            WHEN lower(coalesce(flower_primary_color,'')) ~ 'yellow|orange|gold' THEN 'yellow_orange'
            WHEN lower(coalesce(flower_primary_color,'')) ~ 'white|cream|pale|green|brown|inconspicuous|dull' THEN 'plain'
            ELSE NULL
          END AS color_cat,
          CASE
            WHEN lower(coalesce(floral_form,'')) ~ 'zygomorphic|orchid|papilion|spurr|bilabiate' THEN 'zygomorphic_specialized'
            WHEN lower(coalesce(floral_form,'')) ~ 'tubular|trumpet|funnel|salver|urn|urceolate|campanulate|bell' THEN 'tubular_trumpet'
            WHEN lower(coalesce(floral_form,'')) ~ 'composite|head|brush|puff' THEN 'composite_brush'
            WHEN lower(coalesce(floral_form,'')) ~ 'open|radial|actinomorphic|bowl|cup|star' THEN 'open_generalized'
            ELSE NULL
          END AS form_cat,
          CASE WHEN self_incompatibility IN ('SC','likely_SC') THEN 1
               WHEN self_incompatibility IN ('SI','likely_SI','obligate_SI') THEN 0
               ELSE NULL END AS sc_binary,
          CASE WHEN pollination_functional_guild IN ('birds','butterflies','moths','bats') THEN 1 ELSE 0 END AS showy_alt,
          CASE WHEN pollination_functional_guild IN ('other_bees','bees') THEN 1 ELSE 0 END AS other_bee,
          CASE WHEN pollination_functional_guild IN ('flies','beetles_wasps_ants','mixed','mixed_or_generalist') THEN 1 ELSE 0 END AS generalist_insect
        FROM species_traits;
        """
    )

    con.execute(
        """
        CREATE TEMP TABLE island_traits AS
        SELECT
          island_id,
          count(*) AS n_species_total,
          sum(is_animal) AS n_animal_species,
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
            1.0 - c.bombus_channel_state AS bombus_deficit,
            CASE WHEN t.color_trials > 0 THEN 1.0*t.color_plain/t.color_trials END AS plain_share,
            CASE WHEN t.form_trials > 0 THEN 1.0*t.form_open_generalized/t.form_trials END AS generalized_share,
            CASE WHEN t.sc_trials > 0 THEN 1.0*t.sc_successes/t.sc_trials END AS sc_share
          FROM read_csv_auto('{cov}') c
          LEFT JOIN island_traits t USING (island_id)
        ) TO '{out}' (HEADER, DELIMITER ',');
        """
    )


if __name__ == "__main__":
    main()
