import json
from pathlib import Path

from island_v2.chapter1_v8_figure1_inference_map import load_locks


def test_load_locks_accepts_frozen_contract(tmp_path: Path) -> None:
    config = tmp_path / "config"
    config.mkdir()
    (config / "chapter1_v8_figure2_result_lock.json").write_text(
        json.dumps({"contract": "chapter1_v8_figure2_result_lock_v1"})
    )
    (config / "chapter1_v8_figure3_submission_result_lock.json").write_text(
        json.dumps(
            {
                "contract": "chapter1_v8_figure3_submission_result_lock_v1",
                "frozen_results": {
                    "family_attenuation_fraction_range": [0.196, 0.334],
                    "genus_attenuation_fraction_range": [0.788, 0.859],
                    "conditional_family_to_genus_attenuation_range": [0.706, 0.791],
                    "support_ladder": "4/4 -> 4/4 -> 0/4",
                },
            }
        )
    )
    (config / "chapter1_v8_figure4_result_lock.json").write_text(
        json.dumps(
            {
                "contract": "chapter1_v8_figure4_result_lock_v1",
                "frozen_results": {
                    "v6_palearctic_survival": "99/100",
                    "v6_north_tropical_vector_survival": "70/75",
                    "v6_tropical_survival": "35/75",
                    "geometry_promoted": "0/12",
                    "h5c_interaction_estimate": 0.06495254330943623,
                    "h5c_interaction_ci": [-0.09029679363407014, 0.2202018802529426],
                    "h5c_p_value": 0.4122068222871858,
                    "h5d_qualified": "0/8",
                },
            }
        )
    )
    locks = load_locks(tmp_path)
    assert locks["figure3"]["frozen_results"]["support_ladder"] == "4/4 -> 4/4 -> 0/4"
    assert locks["figure4"]["frozen_results"]["v6_palearctic_survival"] == "99/100"
