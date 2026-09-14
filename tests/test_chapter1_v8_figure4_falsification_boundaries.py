import json
from pathlib import Path

import pandas as pd

from island_v2.chapter1_v8_figure4_falsification_boundaries import load_inputs


def test_load_inputs_accepts_frozen_shapes(tmp_path: Path) -> None:
    v6 = tmp_path / "v6"
    geom = tmp_path / "geom"
    h5c = tmp_path / "h5c"
    h5d = tmp_path / "h5d"
    for path in (v6, geom, h5c, h5d):
        path.mkdir()

    rows = []
    for target in [
        "Palearctic_accessibility",
        "tropical_accessibility",
        "north_tropical_vector_difference",
    ]:
        rows.append(
            {
                "target": target,
                "baseline_supported": True,
                "tipping_event": "none_on_frozen_grid",
            }
        )
    pd.DataFrame(rows).to_csv(
        v6 / "species_detection_tipping_surface.csv",
        index=False,
    )

    pd.DataFrame(
        [
            {
                "all_observed_D": 0,
                "all_critical_D": 1,
                "direct_observed_D": 0,
                "direct_critical_D": 1,
                "classification": "monotonic_or_unresolved",
            }
            for _ in range(12)
        ]
    ).to_csv(geom / "observed_geometry_cross_scope.csv", index=False)

    (h5c / "h5c_observed_result.json").write_text(
        json.dumps(
            {
                "classification": "no_pollination_mode_specificity_support",
                "interaction_estimate": 0.06,
                "interaction_ci_low": -0.09,
                "interaction_ci_high": 0.22,
                "interaction_p_value": 0.41,
            }
        )
    )

    pd.DataFrame(
        [
            {
                "classification_accuracy": 0.75,
                "false_distributed_selection_under_smooth": 0.2,
                "qualified": False,
            }
            for _ in range(8)
        ]
    ).to_csv(h5d / "h5d_identifiability_summary.csv", index=False)

    loaded = load_inputs(v6, geom, h5c, h5d)
    assert len(loaded["geometry"]) == 12
    assert len(loaded["h5d"]) == 8
