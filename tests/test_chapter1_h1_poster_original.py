from pathlib import Path

import pandas as pd

from island_v2.chapter1_h1_poster_original import original_poster_columns


def test_primary_is_exactly_original_including_failure_flags():
    root = Path(__file__).resolve().parents[1]
    got = pd.read_csv(root / "results/h1_poster_original_wcvp_20261004/traitwise_results.csv")
    cols = [
        "context",
        "outcome",
        "geography_slope_log_odds",
        "cluster_robust_se",
        "p_value",
        "optimizer_success",
    ]
    for evidence in ["all", "direct"]:
        original = pd.read_csv(
            root / f"results/geography_20260924/{evidence}/beta_binomial_within_slopes.csv"
        ).query("stratum=='all_observed'")
        actual = got.loc[(got.evidence_scope == evidence) & (got.flora_scope == "broad")]
        pd.testing.assert_frame_equal(
            original[cols].sort_values(cols[:2]).reset_index(drop=True),
            actual[cols].sort_values(cols[:2]).reset_index(drop=True),
            check_exact=True,
        )


def test_display_uses_original_unadjusted_rule():
    x = pd.DataFrame(
        {"geography_slope_log_odds": [0.2], "cluster_robust_se": [0.1], "p_value": [0.049]}
    )
    y = original_poster_columns(x)
    assert y.nominal_supported.iloc[0]
    assert abs(y.ci_low.iloc[0] - 0.004) < 1e-10
