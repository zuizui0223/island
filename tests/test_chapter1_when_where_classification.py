import pandas as pd

from island_v2.chapter1_when_where_classification import classify_when_where_support


def test_classification_distinguishes_supported_null_and_untestable():
    within = pd.DataFrame(
        [
            {"stratum": "all_native", "support_tier": "confirmatory", "context": "north", "where_supported": True},
            {"stratum": "native_nonendemic", "support_tier": "confirmatory", "context": "north", "where_supported": True},
            {"stratum": "all_native", "support_tier": "confirmatory", "context": "tropics", "where_supported": False},
            {"stratum": "all_native", "support_tier": "pilot", "context": "south", "where_supported": True},
            {"stratum": "native_nonendemic", "support_tier": "pilot", "context": "south", "where_supported": True},
        ]
    )
    out = classify_when_where_support(within, ["north", "tropics", "south", "highlat"])
    classes = out.set_index("context")["when_class"].to_dict()
    assert classes["north"] == "confirmatory_persists_in_native_nonendemic"
    assert classes["tropics"] == "confirmatory_tested_null"
    assert classes["south"] == "pilot_signal_confirmatory_not_testable"
    assert classes["highlat"] == "not_testable_current_support"

    south = out.loc[out["context"].eq("south")].iloc[0]
    assert not bool(south["confirmatory_testable_all_native"])
    assert bool(south["pilot_testable_all_native"])
    assert bool(south["pilot_supported_all_native"])
