from __future__ import annotations

import pandas as pd
import pytest

from island_v2.chapter1_taxonomic_depth_decomposition import (
    STAGES,
    build_taxonomic_decomposition,
    classify_depth,
    taxonomic_branching_config,
    validate_taxonomy,
)


def test_taxonomy_hash_and_unique_species_fail_closed() -> None:
    taxonomy = pd.DataFrame(
        {
            "accepted_species": ["Alpha one", "Beta two"],
            "genus": ["Alpha", "Beta"],
            "family": ["FamilyA", "FamilyB"],
        }
    )
    observed = validate_taxonomy(taxonomy, "abc", "abc")
    assert observed.equals(taxonomy)
    with pytest.raises(ValueError, match="SHA-256"):
        validate_taxonomy(taxonomy, "abc", "def")
    with pytest.raises(ValueError, match="duplicate"):
        validate_taxonomy(pd.concat([taxonomy, taxonomy.iloc[[0]]]), "abc", "abc")


def _decomposition_inputs() -> tuple[pd.DataFrame, ...]:
    score_rows: list[dict] = []
    taxonomy_rows: list[dict] = []
    flora_rows: list[dict] = []
    position_rows: list[dict] = []
    family_availability: list[dict] = []
    genus_availability: list[dict] = []
    axes = ["generalized_accessible", "selfing_core"]
    for genus_index in range(8):
        family = f"Family{genus_index // 2}"
        genus = f"Genus{genus_index}"
        family_availability.append({"entity_ID": 1, "taxon_group": family})
        genus_availability.append({"entity_ID": 1, "taxon_group": genus})
        for species_index in range(2):
            species = f"{genus} species{species_index}"
            taxonomy_rows.append(
                {"accepted_species": species, "family": family, "genus": genus}
            )
            flora_rows.append(
                {
                    "island_id": "island_1",
                    "accepted_species": species,
                    "origin_status": "native",
                    "floristic_status": "native_nonendemic",
                }
            )
            for axis_index, axis in enumerate(axes):
                score_rows.append(
                    {
                        "accepted_species": species,
                        "syndrome": axis,
                        "syndrome_concordance": genus_index
                        + species_index * 0.2
                        + axis_index * 0.1,
                    }
                )
        for axis_index, axis in enumerate(axes):
            position_rows.extend(
                [
                    {
                        "taxon_group": family,
                        "syndrome": axis,
                        "source_group_position": genus_index + axis_index * 0.1,
                        "n_source_scored_species": 4,
                        "taxonomic_level": "family",
                    },
                    {
                        "taxon_group": genus,
                        "syndrome": axis,
                        "source_group_position": genus_index
                        + 0.05
                        + axis_index * 0.1,
                        "n_source_scored_species": 2,
                        "taxonomic_level": "genus",
                    },
                ]
            )
    return (
        pd.DataFrame(score_rows),
        pd.DataFrame(taxonomy_rows),
        pd.DataFrame(flora_rows),
        pd.DataFrame(position_rows).drop_duplicates(
            ["taxonomic_level", "taxon_group", "syndrome"]
        ),
        pd.DataFrame(family_availability).drop_duplicates(),
        pd.DataFrame(genus_availability).drop_duplicates(),
    )


def test_family_and_genus_nulls_use_identical_observed_species() -> None:
    scores, taxonomy, flora, positions, family_available, genus_available = (
        _decomposition_inputs()
    )
    output = build_taxonomic_decomposition(
        scores,
        taxonomy,
        flora,
        positions,
        {"family": family_available, "genus": genus_available},
        pd.DataFrame(
            {"island_id": ["island_1"], "source_mode": ["geo_k1"], "entity_ID": [1]}
        ),
        pd.DataFrame({"island_id": ["island_1"]}),
        axes=["generalized_accessible", "selfing_core"],
        source_modes=["geo_k1"],
        strata=["all_native", "native_nonendemic"],
        minimum_species=10,
        minimum_families=3,
        minimum_genera=5,
    )
    assert len(output) == 4
    assert output["n_species"].eq(16).all()
    assert output["after_family_residual"].notna().all()
    assert output["after_genus_residual"].notna().all()
    reconstructed = (
        output["family_expected"]
        + output["family_to_genus_increment"]
        + output["after_genus_residual"]
    )
    pd.testing.assert_series_equal(
        reconstructed.reset_index(drop=True),
        output["observed_score"].reset_index(drop=True),
        check_names=False,
    )


def test_missing_genus_availability_removes_species_from_both_nulls() -> None:
    scores, taxonomy, flora, positions, family_available, genus_available = (
        _decomposition_inputs()
    )
    genus_available = genus_available.loc[
        genus_available["taxon_group"].ne("Genus7")
    ]
    output = build_taxonomic_decomposition(
        scores,
        taxonomy,
        flora,
        positions,
        {"family": family_available, "genus": genus_available},
        pd.DataFrame(
            {"island_id": ["island_1"], "source_mode": ["geo_k1"], "entity_ID": [1]}
        ),
        pd.DataFrame({"island_id": ["island_1"]}),
        axes=["generalized_accessible", "selfing_core"],
        source_modes=["geo_k1"],
        strata=["all_native"],
        minimum_species=10,
        minimum_families=3,
        minimum_genera=5,
    )
    assert output["n_species"].eq(14).all()
    assert output["n_genera"].eq(7).all()


def test_depth_classification_requires_all_source_modes() -> None:
    modes = ["geo_k5", "geo_k10", "geo_k20", "geo50_climate10"]
    rows: list[dict] = []
    for stage in STAGES:
        for mode in modes:
            rows.append(
                {
                    "source_mode": mode,
                    "context_layer": "biogeographic_realm",
                    "axis_set": f"taxonomic_stage__{stage}",
                    "stratum": "native_nonendemic",
                    "support_tier": "confirmatory",
                    "context": "Palearctic",
                    "V2_source_mode_family_supported": stage == "observed_score",
                }
            )
    result = classify_depth(pd.DataFrame(rows), modes, ["native_nonendemic"])
    assert result.loc[0, "V2_taxonomic_depth_classification"] == (
        "compatible_with_deep_family_sorting"
    )


def test_taxonomic_branching_has_three_two_axis_vector_stages() -> None:
    base = {
        "alpha": 0.05,
        "context_layers": {
            "analysis_regime": {
                "column": "analysis_regime",
                "contexts": ["northern_midlatitude", "tropical"],
            }
        },
    }
    config = taxonomic_branching_config(
        base, ["generalized_accessible", "selfing_core"]
    )
    assert set(config["axis_sets"]) == {
        f"taxonomic_stage__{stage}" for stage in STAGES
    }
    assert all(len(spec["axes"]) == 2 for spec in config["axis_sets"].values())
