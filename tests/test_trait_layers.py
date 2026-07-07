from pathlib import Path

import pandas as pd

from island_v2.trait_layers import (
    TRACK_EXCLUDED,
    TRACK_MAIN,
    TRACK_SENSITIVITY,
    analysis_track,
    annotate_candidates,
    load_layers,
    trait_to_layer,
)

CONFIG = load_layers(Path("config/trait_layers.yml"))


def test_shipped_config_layers_and_membership():
    mapping = trait_to_layer(CONFIG)
    assert mapping["flower_primary_color"] == "M0_floral_phenotype"
    assert mapping["self_incompatibility"] == "M1_reproductive_assurance_pollen_vector"
    assert mapping["pollination_functional_guild"] == "M2_pollinator_functional_evidence"


def test_m0_admits_species_indirect_to_main():
    # M0 floral phenotype may use species-indirect evidence in the main analysis.
    assert analysis_track("flower_primary_color", "source_backed", "species_indirect", CONFIG) == (
        "M0_floral_phenotype",
        TRACK_MAIN,
    )
    assert analysis_track("floral_symmetry", "source_backed", "species_direct", CONFIG)[1] == TRACK_MAIN


def test_m1_main_requires_species_direct():
    # species-direct M1 -> main; species-indirect M1 -> only sensitivity (below main bar).
    assert analysis_track("self_incompatibility", "source_backed", "species_direct", CONFIG)[1] == TRACK_MAIN
    assert analysis_track("mating_system", "source_backed", "species_indirect", CONFIG)[1] == TRACK_SENSITIVITY


def test_hierarchical_inference_is_always_sensitivity_never_main():
    for trait in ("self_incompatibility", "pollination_functional_guild", "flower_primary_color"):
        layer, track = analysis_track(trait, "hierarchical_inference", "genus_inference", CONFIG)
        assert track == TRACK_SENSITIVITY
        assert layer != "unclassified_trait"


def test_unresolved_and_unknown_trait_are_excluded():
    assert analysis_track("self_incompatibility", "source_backed", "unresolved", CONFIG)[1] == TRACK_EXCLUDED
    assert analysis_track("not_a_trait", "source_backed", "species_direct", CONFIG) == (
        "unclassified_trait",
        TRACK_EXCLUDED,
    )


def test_annotate_candidates_adds_layer_and_track_columns():
    candidates = pd.DataFrame(
        {
            "accepted_species": ["A a", "B b", "C c", "D d"],
            "trait_name": [
                "flower_primary_color",       # M0 species_direct -> main
                "self_incompatibility",        # M1 hierarchical -> sensitivity
                "mating_system",               # M1 species_indirect -> sensitivity
                "pollination_functional_guild",  # M2 species_direct -> main
            ],
            "candidate_kind": ["source_backed", "hierarchical_inference", "source_backed", "source_backed"],
            "evidence_scope": ["species_direct", "genus_inference", "species_indirect", "species_direct"],
        }
    )
    table, summary = annotate_candidates(candidates, CONFIG)
    tracks = dict(zip(table["trait_name"], table["analysis_track"], strict=True))
    assert tracks["flower_primary_color"] == TRACK_MAIN
    assert tracks["self_incompatibility"] == TRACK_SENSITIVITY
    assert tracks["mating_system"] == TRACK_SENSITIVITY
    assert tracks["pollination_functional_guild"] == TRACK_MAIN
    assert summary["n_main_direct"] == 2
    assert summary["n_expanded_sensitivity"] == 2
    assert set(table["trait_layer"]) >= {"M0_floral_phenotype", "M1_reproductive_assurance_pollen_vector"}
