from __future__ import annotations

import pandas as pd

from island_v2.m0_m4_sem_scoped import select_model_scopes


def test_model_specific_scopes_preserve_m4_structural_comparison() -> None:
    data = pd.DataFrame(
        {
            "island_id": ["core", "weak", "structural", "unresolved"],
            "applicability": [
                "applicable",
                "applicable",
                "structurally_not_applicable",
                "unresolved",
            ],
            "primary_bombus_channel_analysis": [True, True, False, False],
            "bombus_channel_state": [0.8, 0.2, 0.7, 0.1],
        }
    )

    scopes = select_model_scopes(data)
    assert set(scopes["M0"]["island_id"]) == {"core", "weak", "structural"}
    assert set(scopes["M1"]["island_id"]) == {"core", "weak"}
    assert set(scopes["M2"]["island_id"]) == {"core", "weak"}
    assert set(scopes["M3"]["island_id"]) == {"core", "weak"}
    assert set(scopes["M4"]["island_id"]) == {"weak", "structural"}
