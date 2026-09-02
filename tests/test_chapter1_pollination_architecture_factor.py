from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_pollination_architecture_factor import (
    DEFAULT_TEMPLATES,
    RESIDUAL_NAMES,
    SHARED_FACTOR,
    architecture_branching_config,
    build_genus_decomposition,
    fit_and_project_source_factor,
)
from island_v2.chapter1_pr138_syndrome_analysis import build_island_syndrome_scores


def _species_scores(n_species: int = 90, seed: int = 1424) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    common = rng.normal(size=n_species)
    rows: list[dict] = []
    for index in range(n_species):
        species = f"Genus{index // 3} species{index}"
        values = {
            "large_bee_like": 0.8 * common[index] + rng.normal(0, 0.2),
            "butterfly_like": 0.7 * common[index] + rng.normal(0, 0.3),
            "bird_like": 0.9 * common[index] + rng.normal(0, 0.15),
        }
        for syndrome, score in values.items():
            rows.append(
                {
                    "accepted_species": species,
                    "syndrome": syndrome,
                    "syndrome_concordance": score,
                }
            )
    return pd.DataFrame(rows)


def _gift(n_species: int = 70) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "entity_ID": [index % 5 + 1 for index in range(n_species)],
            "work_species": [
                f"Genus{index // 3} species{index}" for index in range(n_species)
            ],
        }
    )


def test_factor_is_source_only_and_outcome_blind() -> None:
    scores = _species_scores()
    projected, model, audit = fit_and_project_source_factor(scores, _gift())

    perturbed = scores.copy()
    nonsource = perturbed["accepted_species"].str.extract(r"species(\d+)")[0].astype(int) >= 70
    perturbed.loc[nonsource, "syndrome_concordance"] *= -50.0
    _, repeated_model, repeated_audit = fit_and_project_source_factor(
        perturbed, _gift()
    )

    pd.testing.assert_frame_equal(model, repeated_model)
    assert audit["factor_fit_used_island_status"] is False
    assert audit["factor_fit_used_island_distance"] is False
    assert repeated_audit["n_complete_gift_source_species"] == 70
    assert set(projected["syndrome"]) == {
        SHARED_FACTOR,
        *RESIDUAL_NAMES.values(),
    }


def test_source_factor_residuals_are_orthogonal_and_orientation_is_positive() -> None:
    _, model, audit = fit_and_project_source_factor(_species_scores(), _gift())
    assert model["source_factor_correlation"].sum() > 0
    assert model["variance_fraction"].iloc[0] > 0.5
    assert max(
        abs(value)
        for value in audit["source_factor_residual_covariances"].values()
    ) < 1.0e-12


def test_missing_template_is_excluded_not_imputed() -> None:
    scores = _species_scores()
    omitted_species = "Genus29 species89"
    scores = scores.loc[
        ~(
            scores["accepted_species"].eq(omitted_species)
            & scores["syndrome"].eq("bird_like")
        )
    ].copy()
    projected, _, audit = fit_and_project_source_factor(scores, _gift())
    assert omitted_species not in set(projected["accepted_species"])
    assert audit["missing_template_imputation"] is False


def test_architecture_config_is_one_prespecified_four_axis_family() -> None:
    branching = {
        "alpha": 0.05,
        "context_layers": {
            "analysis_regime": {
                "column": "analysis_regime",
                "contexts": ["northern_midlatitude", "tropical"],
            }
        },
    }
    components = [SHARED_FACTOR, *RESIDUAL_NAMES.values()]
    config = architecture_branching_config(branching, components)
    assert config["axis_sets"]["shared_architecture_decomposition"]["axes"] == components
    assert config["axis_sets"]["shared_architecture_decomposition"]["classify"] is False
    for component in components:
        assert config["branch_axes"][component]["components"] == {component: 1.0}


def test_genus_decomposition_emits_beyond_genus_residual() -> None:
    rows: list[dict] = []
    gift_rows: list[dict] = []
    flora_rows: list[dict] = []
    for genus_index in range(8):
        genus = f"Genus{chr(ord('a') + genus_index)}"
        for species_index in range(2):
            species = f"{genus} species{chr(ord('a') + species_index)}"
            for component_index, component in enumerate(
                [SHARED_FACTOR, *RESIDUAL_NAMES.values()]
            ):
                rows.append(
                    {
                        "accepted_species": species,
                        "syndrome": component,
                        "syndrome_concordance": genus_index
                        + species_index * 0.25
                        + component_index * 0.1,
                        "n_informative_traits": len(DEFAULT_TEMPLATES),
                    }
                )
            gift_rows.append({"entity_ID": 1, "work_species": species})
            gift_rows.append({"entity_ID": 2, "work_species": species})
            if species_index == 1:
                flora_rows.append(
                    {
                        "island_id": "island_1",
                        "accepted_species": species,
                        "origin_status": "native",
                        "endemic_status": "nonendemic",
                        "floristic_status": "native_nonendemic",
                    }
                )
    projected = pd.DataFrame(rows)
    flora = pd.DataFrame(flora_rows)
    island_scores = build_island_syndrome_scores(
        flora, projected, ["all_native", "native_nonendemic"]
    )
    assignments = pd.DataFrame(
        {
            "island_id": ["island_1", "island_1"],
            "source_mode": ["geo_k2", "geo_k2"],
            "entity_ID": [1, 2],
        }
    )
    decomposition = build_genus_decomposition(
        projected,
        island_scores,
        flora,
        pd.DataFrame(gift_rows),
        assignments,
        pd.DataFrame({"island_id": ["island_1"]}),
        source_modes=["geo_k2"],
        strata=["all_native", "native_nonendemic"],
        minimum_represented_genera=5,
    )
    assert not decomposition.empty
    assert decomposition["n_represented_genera"].eq(8).all()
    assert decomposition["beyond_genus_residual"].abs().max() > 0
    assert set(decomposition["architecture_component"]) == {
        SHARED_FACTOR,
        *RESIDUAL_NAMES.values(),
    }
