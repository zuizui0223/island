import numpy as np

from island_v2.chapter1_h1_traitwise import cluster_inference, holm_fixed


def test_holm_keeps_missing_tests_in_family():
    got = holm_fixed([0.01, np.nan, 0.04, 0.001])
    np.testing.assert_allclose(got[[0, 2, 3]], [0.03, 0.08, 0.004])
    assert np.isnan(got[1])


def test_cluster_inference_is_two_sided_and_sign_symmetric():
    a = cluster_inference(0.4, np.array([0.1, -0.08, 0.06, -0.04]), 40, 3)
    b = cluster_inference(-0.4, np.array([0.1, -0.08, 0.06, -0.04]), 40, 3)
    assert a["p_two_sided"] == b["p_two_sided"]
    assert np.isclose(a["ci_low"], -b["ci_high"])
    assert a["df"] == 3
    assert a["ci_low"] < 0.4 < a["ci_high"]


def test_cluster_inference_fails_closed_for_single_cluster():
    assert cluster_inference(0.4, np.array([0.1]), 40, 3)["status"] == "not_testable"


def test_each_trait_keeps_its_sign_without_pooling(monkeypatch):
    import pandas as pd

    import island_v2.chapter1_h1_traitwise as module

    cfg = {
        "context_column": "region",
        "contexts": ["one"],
        "model_outcomes": ["a", "b"],
        "response_families": {"domain": ["a", "b"]},
        "alpha": 0.05,
    }
    data = pd.DataFrame({"region": ["one", "one"], "outcome": ["a", "b"]})
    monkeypatch.setattr(module, "build_broad_counts", lambda *args: data)
    monkeypatch.setattr(module, "_prepare", lambda *args: data)
    monkeypatch.setattr(
        module,
        "fit_trait",
        lambda part, outcome, cfg: {
            "status": "ok",
            "estimate": 1 if outcome == "a" else -1,
            "p_two_sided": 0.001,
        },
    )
    got = module.run_scope(None, None, None, cfg, "all_observed")
    assert got.interpretation.tolist() == ["classic_direction", "opposite_direction"]
    assert len(got) == 2
    assert not any("score" in column or "weight" in column for column in got.columns)


def test_failed_optimizer_never_produces_significance(monkeypatch):
    import pandas as pd

    import island_v2.chapter1_h1_traitwise as module

    cfg = {
        "minimum_islands_per_outcome": 2,
        "baseline_covariates": [],
        "geography_column": "distance",
        "cluster_column": "cluster",
        "max_iter": 1,
        "retry_max_iter": 2,
    }
    data = pd.DataFrame(
        {
            "island_id": ["a", "b", "c"],
            "successes": [1, 2, 3],
            "trials": [5, 5, 5],
            "distance": [0, 1, 2],
            "cluster": ["a", "b", "c"],
        }
    )
    monkeypatch.setattr(
        module,
        "_fit_single_beta_binomial",
        lambda *args, **kwargs: {"success": False, "message": "failed"},
    )
    got = module.fit_trait(data, "a", cfg)
    assert got["status"] == "optimizer_failure"
    assert "p_two_sided" not in got


def test_committed_replay_has_all_fixed_families_and_verified_inference():
    from pathlib import Path

    import pandas as pd
    from scipy.stats import t

    root = Path(__file__).resolve().parents[1]
    d = pd.read_csv(root / "results/h1_traitwise_20261004/traitwise_results.csv")
    assert len(d) == 112
    assert not d.duplicated(["evidence_scope", "flora_scope", "context", "outcome"]).any()
    assert d.status.eq("ok").all()
    for _, group in d.groupby(["evidence_scope", "flora_scope"]):
        assert len(group) == 28
        p = 2 * t.sf(abs(group.estimate / group.se), group.df)
        np.testing.assert_allclose(p, group.p_two_sided, rtol=1e-10, atol=1e-14)
        np.testing.assert_allclose(holm_fixed(p), group.p_holm_28, rtol=1e-10, atol=1e-14)
    assert set(d.flora_scope) == {"broad", "wcvp"}


def test_final_nominal_policy_does_not_apply_holm(monkeypatch):
    import pandas as pd

    import island_v2.chapter1_h1_traitwise as module

    cfg = {
        "multiplicity": "none",
        "context_column": "region",
        "contexts": ["one"],
        "model_outcomes": ["a", "b"],
        "response_families": {"domain": ["a", "b"]},
        "alpha": 0.05,
    }
    data = pd.DataFrame({"region": ["one", "one"], "outcome": ["a", "b"]})
    monkeypatch.setattr(module, "build_broad_counts", lambda *args: data)
    monkeypatch.setattr(module, "_prepare", lambda *args: data)
    monkeypatch.setattr(
        module, "fit_trait", lambda *args: {"status": "ok", "estimate": 0.2, "p_two_sided": 0.03}
    )
    got = module.run_scope(None, None, None, cfg, "all_observed")
    assert "p_holm_28" not in got
    assert got.interpretation.eq("classic_direction").all()


def test_final_release_has_nominal_t_p_and_correct_intervals():
    from pathlib import Path

    import pandas as pd
    from scipy.stats import t

    root = Path(__file__).resolve().parents[1]
    d = pd.read_csv(root / "results/h1_final_traitwise_t_20261004/traitwise_results.csv")
    assert len(d) == 112 and d.status.eq("ok").all()
    assert "p_holm_28" not in d
    np.testing.assert_allclose(d.p_two_sided, 2 * t.sf(abs(d.estimate / d.se), d.df), atol=1e-14)
    np.testing.assert_allclose(d.ci_low, d.estimate - t.ppf(0.975, d.df) * d.se, atol=1e-14)
    np.testing.assert_allclose(d.ci_high, d.estimate + t.ppf(0.975, d.df) * d.se, atol=1e-14)
    assert d.interpretation.ne("uncertain_or_not_testable").equals(d.p_two_sided.lt(0.05))
