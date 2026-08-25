from pathlib import Path

import pandas as pd
import pytest

from island_v2.purpose_shortest_analysis import (
    build_grounded_trait_composition,
    replace_trait_richness,
)
from island_v2.trait_bombus_analysis import (
    GROUNDED_FILL_TIERS,
    load_trait_rows_for_grounded_analysis,
    require_grounded_trait_rows,
)
from island_v2.v1_v2_bridge_analysis import build_v1_like_animal_counts


def _legacy_sensitivity_master() -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    values = {
        "pollination_functional_guild": "bumblebees",
        "flower_primary_color": "white",
        "floral_form": "open_radial",
    }
    for trait_name, value in values.items():
        rows.append(
            {
                "island_id": "island-1",
                "accepted_species": "Grounded species",
                "trait_name": trait_name,
                "filled_value": value,
                "fill_tier": {
                    "pollination_functional_guild": "species_direct",
                    "flower_primary_color": "genus_inference",
                    "floral_form": "family_inference",
                }[trait_name],
            }
        )

    # These plausible values would change both purpose and bridge counts if the
    # historical sensitivity_all rows were allowed through.
    for species, tier, color, form in (
        ("Global species", "global_fallback", "red_pink", "tubular"),
        ("Unresolved species", "unresolved_no_evidence", "red_pink", "tubular"),
        ("Model species", "llm_inference", "red_pink", "tubular"),
    ):
        for trait_name, value in (
            ("pollination_functional_guild", "bumblebees"),
            ("flower_primary_color", color),
            ("floral_form", form),
        ):
            rows.append(
                {
                    "island_id": "island-1",
                    "accepted_species": species,
                    "trait_name": trait_name,
                    "filled_value": value,
                    "fill_tier": tier,
                }
            )
    return pd.DataFrame(rows)


def test_grounded_filter_is_positive_allow_list_and_fails_closed() -> None:
    master = _legacy_sensitivity_master()
    grounded = require_grounded_trait_rows(master)

    assert set(grounded["accepted_species"]) == {"Grounded species"}
    assert set(grounded["fill_tier"]) == GROUNDED_FILL_TIERS
    assert not grounded["fill_tier"].isin(
        {"global_fallback", "unresolved_no_evidence", "llm_inference"}
    ).any()

    with pytest.raises(Exception, match="lacks fill_tier"):
        require_grounded_trait_rows(master.drop(columns="fill_tier"))


def test_analysis_eligible_false_is_excluded_when_new_flag_exists() -> None:
    master = _legacy_sensitivity_master().loc[
        lambda frame: frame["accepted_species"].eq("Grounded species")
    ].copy()
    rejected = master.copy()
    rejected["accepted_species"] = "Flagged species"
    master["analysis_eligible"] = "true"
    rejected["analysis_eligible"] = "false"

    grounded = require_grounded_trait_rows(pd.concat([master, rejected], ignore_index=True))
    assert set(grounded["accepted_species"]) == {"Grounded species"}


def test_species_master_loader_is_column_bounded_and_requires_fill_tier(tmp_path: Path) -> None:
    path = tmp_path / "master.csv.gz"
    master = _legacy_sensitivity_master()
    master["large_unused_bombus_column"] = "not loaded"
    master.to_csv(path, index=False, compression="gzip")

    loaded = load_trait_rows_for_grounded_analysis(path)
    assert "large_unused_bombus_column" not in loaded.columns
    assert set(loaded.columns) == {
        "island_id",
        "accepted_species",
        "trait_name",
        "filled_value",
        "fill_tier",
    }

    missing = tmp_path / "missing_fill_tier.csv"
    master.drop(columns="fill_tier").to_csv(missing, index=False)
    with pytest.raises(Exception, match="missing columns.*fill_tier"):
        load_trait_rows_for_grounded_analysis(missing)


def test_purpose_rebuilds_composition_and_richness_from_grounded_rows() -> None:
    master = _legacy_sensitivity_master()
    grounded, composition, richness = build_grounded_trait_composition(master)

    assert set(grounded["accepted_species"]) == {"Grounded species"}
    assert set(composition["filled_value"]) == {
        "bumblebees",
        "white",
        "open_radial",
    }
    assert composition["n_species"].eq(1).all()
    assert richness.to_dict("records") == [
        {"island_id": "island-1", "n_trait_species": 1}
    ]

    historical_metrics = pd.DataFrame(
        {"island_id": ["island-1"], "n_trait_species": [4], "metric": [0.5]}
    )
    refreshed = replace_trait_richness(historical_metrics, richness)
    assert refreshed.loc[0, "n_trait_species"] == 1


def test_v1_v2_bridge_counts_never_use_global_or_unresolved_rows() -> None:
    counts = build_v1_like_animal_counts(_legacy_sensitivity_master())
    row = counts.set_index("island_id").loc["island-1"]

    assert row["n_animal_species"] == 1
    assert row["plain_color_trials"] == 1
    assert row["plain_color_successes"] == 1
    assert row["generalized_form_trials"] == 1
    assert row["generalized_form_successes"] == 1


def test_purpose_workflow_passes_species_master_not_aggregated_composition() -> None:
    workflow = Path(".github/workflows/run-purpose-shortest-distance-regime.yml").read_text(
        encoding="utf-8"
    )
    assert 'MASTER="$ROOT/island_species_trait_bombus_master.csv.gz"' in workflow
    assert '--species-master-csv "$MASTER"' in workflow
    assert "--trait-composition-csv" not in workflow
