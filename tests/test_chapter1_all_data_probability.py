import numpy as np
import pandas as pd

from island_v2.chapter1_all_data_probability import (
    _fit_single_beta_binomial,
    build_broad_counts,
)


def test_all_observed_retains_unresolved_status_rows():
    flora = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i2", "i2"],
            "accepted_species": ["a", "b", "a", "c"],
            "origin_status": ["native", "unresolved", "native", "introduced"],
            "floristic_status": [
                "native_nonendemic",
                "unresolved",
                "native_nonendemic",
                "introduced",
            ],
        }
    )
    audit = pd.DataFrame(
        {
            "accepted_species": ["a", "b", "c"],
            "trait_name": ["floral_form", "floral_form", "floral_form"],
            "resolved_for_primary": [True, True, True],
            "canonical_signature": ["open_radial", "tubular", "open_radial"],
        }
    )
    config = {
        "strata": ["all_observed", "all_native", "native_nonendemic"],
        "broad_outcomes": {
            "generalized_form": {
                "trait_name": "floral_form",
                "positive_states": ["open_radial"],
                "negative_states": ["tubular"],
            }
        },
    }
    counts = build_broad_counts(flora, audit, config)
    observed = counts.loc[counts["stratum"].eq("all_observed")].set_index("island_id")
    native = counts.loc[counts["stratum"].eq("all_native")].set_index("island_id")
    assert int(observed.loc["i1", "trials"]) == 2
    assert int(observed.loc["i2", "trials"]) == 2
    assert int(native.loc["i1", "trials"]) == 1
    assert int(native.loc["i2", "trials"]) == 1


def test_beta_binomial_recovers_positive_logit_slope():
    x = np.linspace(-2.0, 2.0, 120)
    design = np.column_stack([np.ones(len(x)), x])
    probability = 1.0 / (1.0 + np.exp(-(-0.3 + 1.2 * x)))
    trials = np.full(len(x), 40.0)
    successes = np.round(trials * probability)
    fit = _fit_single_beta_binomial(
        successes,
        trials,
        design,
        ["intercept", "slope"],
        max_iter=500,
    )
    assert fit["success"]
    assert float(fit["theta"][1]) > 0.8
    assert float(fit["kappa"]) > 0.0
