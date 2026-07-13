"""Build species-level trait states without conflating unknown with observed absence."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import duckdb
import yaml


ONTOLOGY_DIR = Path(__file__).with_name("ontology")


def _load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _category_case(column: str, ontology: dict) -> str:
    clauses: list[str] = []
    for category, terms in ontology["categories"].items():
        pattern = "|".join(re.escape(str(term).lower()) for term in terms)
        clauses.append(
            f"WHEN regexp_matches(lower(coalesce({column}, '')), '(?:{pattern})') "
            f"THEN '{category}'"
        )
    return "CASE " + " ".join(clauses) + " ELSE NULL END"


def _guild_case(column: str, values: list[str]) -> str:
    quoted = ",".join("'" + str(value).replace("'", "''") + "'" for value in values)
    return (
        "CASE "
        f"WHEN {column} IS NULL OR trim({column}) = '' THEN NULL "
        f"WHEN {column} IN ({quoted}) THEN 1 "
        "ELSE 0 END"
    )


def build_species_trait_states(master_csv_gz: Path, output_parquet: Path) -> None:
    color = _load_yaml(ONTOLOGY_DIR / "flower_color.yml")
    form = _load_yaml(ONTOLOGY_DIR / "floral_form.yml")
    guild = _load_yaml(ONTOLOGY_DIR / "pollinator_guild.yml")
    output_parquet.parent.mkdir(parents=True, exist_ok=True)

    master = str(master_csv_gz).replace("'", "''")
    out = str(output_parquet).replace("'", "''")
    color_case = _category_case("flower_primary_color", color)
    form_case = _category_case("floral_form", form)
    showy_alt_case = _guild_case(
        "pollination_functional_guild", guild["groups"]["showy_alt"]["positive"]
    )
    other_bee_case = _guild_case(
        "pollination_functional_guild", guild["groups"]["other_bee"]["positive"]
    )
    generalist_case = _guild_case(
        "pollination_functional_guild",
        guild["groups"]["generalist_insect"]["positive"],
    )

    con = duckdb.connect()
    con.execute(
        f"""
        COPY (
          WITH species_traits AS (
            SELECT
              island_id,
              accepted_species,
              max(CASE WHEN trait_name='pollen_vector_mode' THEN filled_value END)
                AS pollen_vector_mode,
              max(CASE WHEN trait_name='flower_primary_color' THEN filled_value END)
                AS flower_primary_color,
              max(CASE WHEN trait_name='floral_form' THEN filled_value END)
                AS floral_form,
              max(CASE WHEN trait_name='self_incompatibility' THEN filled_value END)
                AS self_incompatibility,
              max(CASE WHEN trait_name='pollination_functional_guild' THEN filled_value END)
                AS pollination_functional_guild
            FROM read_csv_auto('{master}', union_by_name=true)
            WHERE trait_name IN (
              'pollen_vector_mode','flower_primary_color','floral_form',
              'self_incompatibility','pollination_functional_guild'
            )
            GROUP BY island_id, accepted_species
          )
          SELECT
            *,
            CASE
              WHEN pollen_vector_mode IS NULL OR trim(pollen_vector_mode) = '' THEN NULL
              WHEN pollen_vector_mode = 'biotic' THEN 1
              ELSE 0
            END AS is_animal,
            {color_case} AS color_cat,
            {form_case} AS form_cat,
            CASE
              WHEN self_incompatibility IN ('SC','likely_SC') THEN 1
              WHEN self_incompatibility IN ('SI','likely_SI','obligate_SI') THEN 0
              ELSE NULL
            END AS sc_binary,
            {showy_alt_case} AS showy_alt,
            {other_bee_case} AS other_bee,
            {generalist_case} AS generalist_insect
          FROM species_traits
        ) TO '{out}' (FORMAT PARQUET, COMPRESSION ZSTD)
        """
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv-gz", type=Path, required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    args = parser.parse_args()
    build_species_trait_states(args.master_csv_gz, args.output_parquet)


if __name__ == "__main__":
    main()
