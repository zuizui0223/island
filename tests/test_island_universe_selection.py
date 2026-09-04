import pandas as pd
import pytest

from island_v2.island_universe_selection import audit_island_universe


def candidate_table() -> pd.DataFrame:
    regimes = [
        "northern_midlatitude",
        "northern_midlatitude",
        "northern_high_latitude",
        "northern_high_latitude",
        "tropical",
        "tropical",
        "southern_extratropical",
        "southern_extratropical",
    ]
    return pd.DataFrame(
        {
            "island_id": [f"i{index}" for index in range(8)],
            "analysis_regime": regimes,
            "area_km2": [5, 6, 7, 8, 9, 10, 11, 12],
            "distance_to_continent_km": [0, 1, 2, 3, 4, 5, 6, 7],
            "climate_pc1": range(8),
            "climate_pc2": range(1, 9),
            "climate_pc3": range(2, 10),
            "climate_pc4": range(3, 11),
        }
    )


def test_audit_reconstructs_nested_attrition_and_selection_weights() -> None:
    candidates = candidate_table()
    effort = pd.DataFrame({"island_id": ["i0", "i1", "i2", "i3", "i4"]})
    occurrences = pd.DataFrame(
        {
            "island_id": ["i0", "i1", "i2", "i3"],
            "species": ["sp0", "sp1", "sp2", "outside"],
        }
    )
    strict = pd.DataFrame({"accepted_species": ["sp0", "sp1", "sp2"]})
    composition = pd.DataFrame({"island_id": ["i0", "i1"]})

    island, stages, regions, quartiles, manifest = audit_island_universe(
        candidates, effort, occurrences, strict, composition
    )

    assert stages["n_islands"].tolist() == [8, 5, 4, 3, 2, 2]
    assert stages["n_lost_from_previous"].tolist() == [0, 3, 1, 1, 1, 0]
    assert int(regions["candidate_islands"].sum()) == 8
    assert int(regions["analysis_islands"].sum()) == 2
    assert len(quartiles) == 8
    assert island["selection_propensity"].between(0, 1).all()
    assert island.loc[island["analysis_included"], "stabilized_ipw"].notna().all()
    assert island.loc[~island["analysis_included"], "stabilized_ipw"].isna().all()
    assert manifest["selection_model"]["outcome_blind_to_trait_values"] is True


def test_audit_rejects_non_nested_composition_islands() -> None:
    candidates = candidate_table()
    effort = pd.DataFrame({"island_id": ["i0"]})
    occurrences = pd.DataFrame({"island_id": ["i0"], "species": ["sp0"]})
    strict = pd.DataFrame({"accepted_species": ["sp0"]})
    composition = pd.DataFrame({"island_id": ["i7"]})

    with pytest.raises(ValueError, match="resolved-trait islands is not nested"):
        audit_island_universe(candidates, effort, occurrences, strict, composition)
