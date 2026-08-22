import importlib

import pandas as pd

mod = importlib.import_module("island_v2.introduced_status_trait_analysis")


def test_confirmed_introduced_excludes_native_and_unresolved():
    status = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i2"],
            "accepted_species": ["A a", "B b", "C c"],
            "origin_status": ["introduced", "native", "unresolved"],
        }
    )
    out = mod.confirmed_introduced_flora(status)
    assert out.to_dict("records") == [
        {"island_id": "i1", "accepted_species": "A a"}
    ]
