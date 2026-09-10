from __future__ import annotations

import pandas as pd

from island_v2.chapter1_database_release import build_inventory


def test_release_inventory_is_fail_closed_for_unreviewed_sources() -> None:
    frame = pd.DataFrame(
        {
            "source_groups": ["open_source|unknown_source", "open_source"],
            "source_lineages": ["a", "b"],
        }
    )
    policy = {
        "default_status": "review_required",
        "rules": [
            {
                "pattern": "^open_source$",
                "status": "redistributable",
                "note": "audited",
            }
        ],
    }
    inventory = build_inventory(frame, policy).set_index("source_token")
    assert inventory.loc["open_source", "redistribution_status"] == "redistributable"
    assert inventory.loc["unknown_source", "redistribution_status"] == "review_required"
