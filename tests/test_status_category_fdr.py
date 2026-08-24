import importlib

import pandas as pd

mod = importlib.import_module("island_v2.status_category_fdr")


def test_bh_is_applied_within_domain_stratum_regime():
    data = pd.DataFrame(
        {
            "stratum": ["all_native"] * 5,
            "regime": ["tropical"] * 5,
            "domain": ["colour"] * 4 + ["form"],
            "category": ["a", "b", "c", "d", "e"],
            "p_value": [0.001, 0.02, 0.2, 0.8, 0.04],
        }
    )
    out = mod.add_category_fdr(data).set_index("category")
    assert out.loc["a", "q_value_bh"] <= 0.0041
    assert out.loc["e", "q_value_bh"] == 0.04
    assert bool(out.loc["a", "fdr_supported"])
    assert bool(out.loc["e", "fdr_supported"])
