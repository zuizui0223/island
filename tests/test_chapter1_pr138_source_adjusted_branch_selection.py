from __future__ import annotations

import pandas as pd
import pytest
import typer

import island_v2.chapter1_pr138_source_adjusted_branch_selection as joint


def _scores() -> pd.DataFrame:
    rows = []
    for mode in ["geo_k5", "geo_k10"]:
        for island in ["i1", "i2"]:
            for syndrome in ["generalized_accessible", "selfing_core"]:
                rows.append(
                    {
                        "island_id": island,
                        "stratum": "all_native",
                        "syndrome": syndrome,
                        "syndrome_score": 0.1,
                        "n_species": 3,
                        "source_mode": mode,
                        "source_expectation_eligible": True,
                    }
                )
    rows.append(
        {
            "island_id": "drop_me",
            "stratum": "all_native",
            "syndrome": "selfing_core",
            "syndrome_score": 0.9,
            "n_species": 3,
            "source_mode": "geo_k5",
            "source_expectation_eligible": False,
        }
    )
    return pd.DataFrame(rows)


def test_joint_wrapper_runs_each_source_mode_separately(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[set[str]] = []

    def fake_run(scores, *_args, **_kwargs):
        calls.append(set(scores["source_mode"].astype(str)))
        assert "drop_me" not in set(scores["island_id"].astype(str))
        return tuple(pd.DataFrame({"marker": [len(scores)]}) for _ in range(8))

    monkeypatch.setattr(joint, "run_selection_sensitivity", fake_run)
    outputs = joint.run_joint_sensitivity(
        _scores(),
        pd.DataFrame(),
        pd.DataFrame(),
        {},
        {},
        {},
    )

    assert calls == [{"geo_k10"}, {"geo_k5"}] or calls == [{"geo_k5"}, {"geo_k10"}]
    for frame in outputs:
        assert set(frame["source_mode"]) == {"geo_k5", "geo_k10"}
        assert len(frame) == 2


def test_source_adjusted_scores_fail_closed_on_duplicates() -> None:
    scores = _scores()
    duplicate = pd.concat([scores, scores.iloc[[0]]], ignore_index=True)
    with pytest.raises(typer.BadParameter, match="duplicate source-adjusted"):
        joint._validate_source_adjusted_scores(duplicate)


def test_source_adjusted_scores_require_eligible_source_expectation() -> None:
    scores = _scores()
    scores["source_expectation_eligible"] = False
    with pytest.raises(typer.BadParameter, match="no source-adjusted scores"):
        joint._validate_source_adjusted_scores(scores)
