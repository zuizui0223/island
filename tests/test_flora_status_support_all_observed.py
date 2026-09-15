import pandas as pd

from island_v2.flora_status_support import STRATA, stratum_mask


def test_all_observed_retains_every_flora_row_regardless_of_status():
    frame = pd.DataFrame(
        {
            "origin_status": ["native", "introduced", "unresolved", "native"],
            "floristic_status": [
                "native_nonendemic",
                "introduced",
                "unresolved",
                "endemic",
            ],
        }
    )
    mask = stratum_mask(frame, "all_observed")
    assert mask.dtype == bool
    assert mask.tolist() == [True, True, True, True]


def test_all_observed_is_explicitly_part_of_supported_strata():
    assert STRATA[0] == "all_observed"
    assert set(STRATA) >= {
        "all_observed",
        "all_native",
        "native_nonendemic",
        "endemic",
        "introduced",
        "unresolved",
    }
