import pandas as pd

from island_v2.chapter1_pr138_direct_ledger_pattern import build_direct_ledger_broad_counts


def test_direct_ledger_builder_expands_new_syndrome_axes_and_fails_closed():
    status = pd.DataFrame(
        {
            "island_id": ["i1"] * 8,
            "accepted_species": [
                "sym", "size", "tube", "mate", "auto", "mixed_size", "mixed_mate", "unresolved"
            ],
            "origin_status": ["native"] * 8,
            "endemic_status": ["nonendemic"] * 8,
            "floristic_status": ["native_nonendemic"] * 8,
        }
    )
    ledger = pd.DataFrame(
        {
            "accepted_species": [
                "sym", "size", "tube", "mate", "auto", "mixed_size", "mixed_mate", "unresolved"
            ],
            "trait_name": [
                "floral_symmetry", "flower_size_class", "tube_depth_class", "mating_system",
                "autonomous_selfing_capacity", "flower_size_class", "mating_system", "floral_symmetry"
            ],
            "resolution_status": ["resolved"] * 7 + ["unresolved"],
            "normalized_value": [
                "actinomorphic", "small", "shallow", "predominantly_selfing", "autonomous",
                "small|large", "mixed_mating", "zygomorphic"
            ],
        }
    )
    config = {
        "strata": ["all_native", "native_nonendemic"],
        "broad_outcomes": {
            "actinomorphic_symmetry": {
                "trait_name": "floral_symmetry",
                "positive_states": ["actinomorphic"],
                "negative_states": ["zygomorphic"],
            },
            "small_flower": {
                "trait_name": "flower_size_class",
                "positive_states": ["very_small", "small"],
                "negative_states": ["large", "very_large"],
            },
            "shallow_open_tube": {
                "trait_name": "tube_depth_class",
                "positive_states": ["absent_or_open", "shallow"],
                "negative_states": ["deep"],
            },
            "selfing_mating_system": {
                "trait_name": "mating_system",
                "positive_states": ["predominantly_selfing", "obligate_selfing"],
                "negative_states": ["predominantly_outcrossing"],
            },
            "autonomous_selfing": {
                "trait_name": "autonomous_selfing_capacity",
                "positive_states": ["autonomous", "delayed"],
                "negative_states": ["absent"],
            },
        },
    }
    counts = build_direct_ledger_broad_counts(status, ledger, config)
    native = counts.loc[counts["stratum"].eq("all_native")].set_index("outcome")
    assert set(native.index) == {
        "actinomorphic_symmetry", "small_flower", "shallow_open_tube",
        "selfing_mating_system", "autonomous_selfing",
    }
    assert native["successes"].eq(1).all()
    assert native["trials"].eq(1).all()
    # small|large crosses the predeclared contrast and mixed_mating is not forced to either side.
    assert "mixed_size" not in counts.columns
