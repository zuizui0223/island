from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("duckdb")
from island_v2.prepare_v2_bayesian_main_input import prepare_input


def test_main_input_excludes_every_non_direct_fill_tier(tmp_path: Path) -> None:
    master = tmp_path / "master.csv.gz"
    rows = []
    for trait, value in (
        ("pollen_vector_mode", "biotic"),
        ("flower_primary_color", "white"),
        ("floral_form", "open_radial"),
        ("self_incompatibility", "SC"),
        ("pollination_functional_guild", "other_bees"),
    ):
        rows.append(
            {
                "island_id": "island-1",
                "accepted_species": "Directa alba",
                "trait_name": trait,
                "filled_value": value,
                "fill_tier": "species_direct",
            }
        )
    for tier, species in (
        ("genus_inference", "Genusa rubra"),
        ("family_inference", "Familia rubra"),
        ("unresolved_no_evidence", "Unresolveda rubra"),
        ("llm_inference", "Modela rubra"),
    ):
        for trait, value in (
            ("pollen_vector_mode", "biotic"),
            ("flower_primary_color", "red_pink"),
            ("floral_form", "tubular"),
            ("self_incompatibility", "SC"),
        ):
            rows.append(
                {
                    "island_id": "island-1",
                    "accepted_species": species,
                    "trait_name": trait,
                    "filled_value": value,
                    "fill_tier": tier,
                }
            )
    pd.DataFrame(rows).to_csv(master, index=False, compression="gzip")
    covariates = tmp_path / "covariates.csv"
    pd.DataFrame(
        {
            "island_id": ["island-1", "island-2"],
            "bombus_channel_state": [0.25, 0.75],
        }
    ).to_csv(covariates, index=False)

    output = tmp_path / "main.csv"
    prepare_input(master, covariates, output)

    result = pd.read_csv(output).set_index("island_id")
    direct = result.loc["island-1"]
    assert direct["n_species_total"] == 1
    assert direct["n_animal_species"] == 1
    assert direct["color_trials"] == 1
    assert direct["color_plain"] == 1
    assert direct["color_red_pink"] == 0
    assert direct["form_trials"] == 1
    assert direct["form_open_generalized"] == 1
    assert direct["sc_trials"] == 1
    assert direct["sc_successes"] == 1
    assert result.loc["island-2", "bombus_deficit"] == 0.25
