from pathlib import Path

import pandas as pd

from island_v2.trait_measurements import annotate_measurements, derive_class, load_rules

RULES = load_rules(Path("config/measurement_classification.yml"))


def test_flower_size_bands():
    assert derive_class("flower_size_class", "1.5", "mm", RULES) == "minute"
    assert derive_class("flower_size_class", "8", "mm", RULES) == "small"
    assert derive_class("flower_size_class", "2", "cm", RULES) == "medium"      # 20 mm
    assert derive_class("flower_size_class", "5", "cm", RULES) == "large"       # 50 mm
    assert derive_class("flower_size_class", "0.2", "m", RULES) == "very_large"  # 200 mm, open top


def test_tube_depth_bands_and_range_midpoint():
    assert derive_class("tube_depth_class", "0.5", "mm", RULES) == "no_tube"
    assert derive_class("tube_depth_class", "3", "mm", RULES) == "shallow"
    # a "10-20 mm" range uses its 15 mm midpoint -> medium band (<= 15)
    assert derive_class("tube_depth_class", "10-20", "mm", RULES) == "medium"


def test_unknown_unit_or_unparseable_is_blank_not_guessed():
    assert derive_class("flower_size_class", "8", "cubits", RULES) == ""
    assert derive_class("flower_size_class", "large-ish", "mm", RULES) == ""
    assert derive_class("not_a_measurement_trait", "8", "mm", RULES) == ""


def test_annotate_appends_class_and_version_without_touching_raw():
    table = pd.DataFrame(
        {
            "accepted_species": ["A a", "B b"],
            "trait_name": ["flower_size_class", "tube_depth_class"],
            "raw_measurement": ["25", "60"],
            "raw_unit": ["mm", "mm"],
            "measurement_structure": ["flower_diameter", "corolla_tube_depth"],
            "measurement_source_text": ["flowers ca. 25 mm", "tube 60 mm long"],
        }
    )
    before_raw = table[["raw_measurement", "raw_unit"]].copy()
    result, summary = annotate_measurements(table, RULES)
    assert result.loc[0, "derived_class"] == "medium"     # 25 mm
    assert result.loc[1, "derived_class"] == "very_deep"  # 60 mm, open top
    assert result.loc[0, "classification_rule_version"] == "size_tube_v1"
    assert summary["n_class_derived"] == 2
    pd.testing.assert_frame_equal(result[["raw_measurement", "raw_unit"]], before_raw)
