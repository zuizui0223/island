import pandas as pd

from island_v2.chapter1_v13_raw_colour_coupling_audit import (
    build_colour_conditioned_architecture_counts,
)


def test_colour_conditioned_form_counts_use_colour_species_as_denominator():
    axis = pd.DataFrame(
        {
            "accepted_species": ["A a", "A a", "B b", "B b", "C c", "C c", "D d", "D d"],
            "axis": [
                "flower_colour",
                "floral_structural_complexity",
                "flower_colour",
                "floral_structural_complexity",
                "flower_colour",
                "floral_structural_complexity",
                "flower_colour",
                "floral_structural_complexity",
            ],
            "trait_composition": [
                'flower_primary_color=["red_pink"]',
                'floral_form=["tubular"]|tube_depth_class=["deep"]',
                'flower_primary_color=["red_pink"]',
                'floral_form=["open_radial"]|tube_depth_class=["absent_or_open"]',
                'flower_primary_color=["yellow_orange"]',
                'floral_form=["bell_campanulate"]|tube_depth_class=["intermediate"]',
                'flower_primary_color=["yellow_orange"]',
                'floral_form=["open_radial"]|tube_depth_class=["shallow"]',
            ],
            "quality": ["high"] * 8,
        }
    )
    flora = pd.DataFrame(
        {
            "island_id": ["i1"] * 4,
            "accepted_species": ["A a", "B b", "C c", "D d"],
            "origin_status": ["native"] * 4,
            "endemic_status": ["nonendemic"] * 4,
            "floristic_status": ["native_nonendemic"] * 4,
        }
    )

    counts = build_colour_conditioned_architecture_counts(
        axis,
        flora,
        evidence_scope="all",
        strata=["all_observed"],
    ).set_index("combination")

    assert counts.loc["red_pink__butterfly_form_given_colour", "trials"] == 2
    assert counts.loc["red_pink__butterfly_form_given_colour", "successes"] == 1
    assert counts.loc["red_pink__butterfly_deep_tube_given_colour", "trials"] == 2
    assert counts.loc["red_pink__butterfly_deep_tube_given_colour", "successes"] == 1
    assert counts.loc["yellow_orange__bird_form_given_colour", "trials"] == 2
    assert counts.loc["yellow_orange__bird_form_given_colour", "successes"] == 1
    assert counts.loc["yellow_orange__bird_deep_tube_given_colour", "trials"] == 2
    assert counts.loc["yellow_orange__bird_deep_tube_given_colour", "successes"] == 1


def test_colour_conditioned_denominator_excludes_species_without_required_architecture_trait():
    axis = pd.DataFrame(
        {
            "accepted_species": ["A a", "A a", "B b", "B b"],
            "axis": [
                "flower_colour",
                "floral_structural_complexity",
                "flower_colour",
                "floral_structural_complexity",
            ],
            "trait_composition": [
                'flower_primary_color=["blue_purple"]',
                'floral_form=["bilabiate"]',
                'flower_primary_color=["blue_purple"]',
                'tube_depth_class=["deep"]',
            ],
            "quality": ["high"] * 4,
        }
    )
    flora = pd.DataFrame(
        {
            "island_id": ["i1", "i1"],
            "accepted_species": ["A a", "B b"],
            "origin_status": ["native", "native"],
            "endemic_status": ["nonendemic", "nonendemic"],
            "floristic_status": ["native_nonendemic", "native_nonendemic"],
        }
    )

    counts = build_colour_conditioned_architecture_counts(
        axis,
        flora,
        evidence_scope="all",
        strata=["all_observed"],
    ).set_index("combination")

    # Form coupling uses only A a because B b has no floral_form evidence.
    assert counts.loc["blue_purple__large_bee_form_given_colour", "trials"] == 1
    assert counts.loc["blue_purple__large_bee_form_given_colour", "successes"] == 1
    # Tube coupling uses only B b because A a has no tube_depth_class evidence.
    assert counts.loc["blue_purple__large_bee_deep_tube_given_colour", "trials"] == 1
    assert counts.loc["blue_purple__large_bee_deep_tube_given_colour", "successes"] == 1
