from __future__ import annotations

import pandas as pd
import pytest

from island_v2.within_lineage_population_pilot import exact_island_permutation


def test_arabidopsis_exact_island_permutation() -> None:
    tm = [
        0.934, 0.881, 0.827, 0.957, 0.870, 0.934, 0.056, 0.215, 0.075,
        0.907, 0.974, 0.966, 1.030, 0.878, 0.961, 0.574, 0.134, 0.903,
    ]
    status = ["non_island"] * 18
    for index in (13, 16, 17):
        status[index] = "island"
    table = pd.DataFrame({"tm": tm, "island_status": status})
    result = exact_island_permutation(table)

    assert result["n_exact_labelings"] == 816
    assert result["island_mean_tm"] == pytest.approx(0.6383333333333333)
    assert result["nonisland_mean_tm"] == pytest.approx(0.7440666666666667)
    assert result["island_minus_nonisland_tm"] == pytest.approx(-0.10573333333333335)
    assert result["exact_two_sided_p"] == pytest.approx(0.6703431372549019)
    assert result["exact_one_sided_lower_tm_p"] == pytest.approx(0.2855392156862745)
