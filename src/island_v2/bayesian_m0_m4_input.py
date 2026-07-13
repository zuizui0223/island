"""Build the minimal v2 Bayesian M0-M4 input from the locked island×species trait master.

The main analysis keeps all original trait categories in a long audit table, while planned
binary contrasts are created only for the Bayesian path models. Floral responses are
computed within animal-pollinated species, matching v1's biological scope while retaining
v2's finer categories for later multinomial models.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

SC_VALUES = {"SC", "likely_SC", "self_compatible", "likely_self_compatible"}
SI_VALUES = {"SI", "likely_SI", "obligate_SI", "self_incompatible", "likely_self_incompatible"}
ANIMAL_VECTOR_VALUES = {"biotic", "mixed"}
SHOWY_ALT_GUILDS = {"birds", "butterflies", "moths", "bats"}
ALL_ALT_GUILDS = {"bees", "other_bees", "flies", "butterflies", "moths", "birds", "bats", "mixed"}

PLAIN_COLOR_RE = re.compile(r"white|cream|pale|green|brown|inconspicuous|dull", re.I)
CONSPICUOUS_COLOR_RE = re.compile(
    r"yellow|orange|gold|red|pink|magenta|scarlet|crimson|blue|purple|violet|lilac|lavender",
    re.I,
)
GENERALIZED_FORM_RE = re.compile(
    r"open_radial|brush_puff|composite_head|actinomorphic|radial|rotate|open|bowl|cup|"
    r"star|campanulate|bell|composite|head|brush|puff|globular|spherical|daisy|radiate",
    re.I,
)
SPECIALIZED_FORM_RE = re.compile(
    r"tubular|trumpet|funnel|salver|urceolate|urn|zygomorphic|bilabiate|bilateral|"
    r"papilion|orchid|labiate|spurred",
    re.I,
)


def _pivot_species(master_csv: Path) -> pd.DataFrame:
    usecols = [
        "island_id",
        "accepted_species",
        "trait_name",
        "filled_value",
        "environmental_compatibility_max",
    ]
    long = pd.read_csv(master_csv, usecols=usecols)
    long["island_id"] = long["island_id"].astype(str)
    long["accepted_species"] = long["accepted_species"].astype(str)
    trait = (
        long[["island_id", "accepted_species", "trait_name", "filled_value"]]
        .drop_duplicates(["island_id", "accepted_species", "trait_name"])
        .pivot(index=["island_id", "accepted_species"], columns="trait_name", values="filled_value")
        .reset_index()
    )
    bombus = (
        long[["island_id", "environmental_compatibility_max"]]
        .drop_duplicates("island_id")
        .rename(columns={"environmental_compatibility_max": "bombus_compatibility"})
    )
    return trait.merge(bombus, on="island_id", how="left", validate="many_to_one")


def _classify_color(value: object) -> str | None:
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return None
    if PLAIN_COLOR_RE.search(text):
        return "plain"
    if CONSPICUOUS_COLOR_RE.search(text):
        return "conspicuous"
    return None


def _classify_form(value: object) -> str | None:
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return None
    # Specialized terms take priority for compound labels such as "tubular and zygomorphic".
    if SPECIALIZED_FORM_RE.search(text):
        return "specialized"
    if GENERALIZED_FORM_RE.search(text):
        return "generalized"
    return None


def _count_classified(
    data: pd.DataFrame,
    *,
    column: str,
    classifier,
    positive_label: str,
    prefix: str,
) -> pd.DataFrame:
    x = data[["island_id", "accepted_species", column]].drop_duplicates(
        ["island_id", "accepted_species"]
    ).copy()
    x["contrast_class"] = x[column].map(classifier)
    x = x.loc[x["contrast_class"].notna()].copy()
    x["success"] = x["contrast_class"].eq(positive_label).astype(int)
    out = x.groupby("island_id").agg(
        **{
            f"{prefix}_successes": ("success", "sum"),
            f"{prefix}_trials": ("success", "size"),
        }
    )
    out[f"{prefix}_share"] = out[f"{prefix}_successes"] / out[f"{prefix}_trials"]
    return out.reset_index()


def _count_exact_contrast(
    data: pd.DataFrame,
    *,
    column: str,
    positive: set[str],
    denominator: set[str],
    prefix: str,
) -> pd.DataFrame:
    x = data.loc[
        data[column].astype(str).isin(denominator),
        ["island_id", "accepted_species", column],
    ].drop_duplicates(["island_id", "accepted_species"])
    x = x.copy()
    x["success"] = x[column].astype(str).isin(positive).astype(int)
    out = x.groupby("island_id").agg(
        **{
            f"{prefix}_successes": ("success", "sum"),
            f"{prefix}_trials": ("success", "size"),
        }
    )
    out[f"{prefix}_share"] = out[f"{prefix}_successes"] / out[f"{prefix}_trials"]
    return out.reset_index()


def _guild_counts(data: pd.DataFrame) -> pd.DataFrame:
    x = data.loc[
        data["pollination_functional_guild"].notna(),
        ["island_id", "accepted_species", "pollination_functional_guild"],
    ].drop_duplicates(["island_id", "accepted_species"])
    x = x.copy()
    guild = x["pollination_functional_guild"].astype(str)
    x["showy_alt"] = guild.isin(SHOWY_ALT_GUILDS).astype(int)
    x["any_alt"] = guild.isin(ALL_ALT_GUILDS).astype(int)
    out = x.groupby("island_id").agg(
        showy_alt_successes=("showy_alt", "sum"),
        any_alt_successes=("any_alt", "sum"),
        alternative_guild_trials=("accepted_species", "size"),
    ).reset_index()
    out["showy_alt_share"] = out["showy_alt_successes"] / out["alternative_guild_trials"]
    out["any_alt_share"] = out["any_alt_successes"] / out["alternative_guild_trials"]
    return out


def build(master_csv: Path, island_covariates_csv: Path, output_dir: Path) -> None:
    species = _pivot_species(master_csv)
    cov = pd.read_csv(island_covariates_csv)
    cov["island_id"] = cov["island_id"].astype(str)

    # v1-comparable floral layer: animal-pollinated species only.
    animal = species.loc[
        species["pollen_vector_mode"].astype(str).isin(ANIMAL_VECTOR_VALUES)
    ].copy()

    tables = [
        _count_classified(
            animal,
            column="flower_primary_color",
            classifier=_classify_color,
            positive_label="plain",
            prefix="plain_color",
        ),
        _count_classified(
            animal,
            column="floral_form",
            classifier=_classify_form,
            positive_label="generalized",
            prefix="generalized_form",
        ),
        _count_exact_contrast(
            animal,
            column="self_incompatibility",
            positive=SC_VALUES,
            denominator=SC_VALUES | SI_VALUES,
            prefix="self_compatible",
        ),
        _guild_counts(animal),
    ]

    analysis = cov.copy()
    for table in tables:
        analysis = analysis.merge(table, on="island_id", how="left", validate="one_to_one")

    analysis["log_distance"] = pd.to_numeric(
        analysis["log_distance_to_continent_km"], errors="coerce"
    )
    analysis["log_area"] = pd.to_numeric(analysis["log_island_area_km2"], errors="coerce")
    analysis["bombus_compatibility"] = pd.to_numeric(
        analysis["bombus_channel_state"], errors="coerce"
    )
    analysis["is_northern_midlatitude"] = analysis["analysis_regime"].astype(str).eq(
        "northern_midlatitude"
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    analysis.to_csv(output_dir / "bayesian_m0_m4_island_input.csv", index=False)

    # Preserve the original categorical composition for later multinomial analysis.
    category_rows = []
    for trait_name in (
        "flower_primary_color",
        "floral_form",
        "floral_symmetry",
        "pollination_functional_guild",
    ):
        if trait_name not in animal.columns:
            continue
        x = animal.loc[
            animal[trait_name].notna(),
            ["island_id", "accepted_species", trait_name],
        ].drop_duplicates()
        x = (
            x.groupby(["island_id", trait_name])
            .size()
            .rename("n_species")
            .reset_index()
            .rename(columns={trait_name: "category"})
        )
        x.insert(1, "trait_name", trait_name)
        category_rows.append(x)
    pd.concat(category_rows, ignore_index=True).to_csv(
        output_dir / "bayesian_m0_m4_category_counts.csv", index=False
    )

    manifest = pd.DataFrame(
        {
            "metric": [
                "n_islands",
                "n_species_island_rows",
                "n_animal_species_island_rows",
                "n_northern_midlatitude_islands",
            ],
            "value": [
                int(analysis["island_id"].nunique()),
                int(len(species)),
                int(len(animal)),
                int(analysis["is_northern_midlatitude"].sum()),
            ],
        }
    )
    manifest.to_csv(output_dir / "bayesian_m0_m4_input_manifest.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--island-covariates-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    build(args.master_csv, args.island_covariates_csv, args.output_dir)


if __name__ == "__main__":
    main()
