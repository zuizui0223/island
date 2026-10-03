import pandas as pd
import pytest

from island_v2.chapter1_final_inference_audit import (
    audit_h2,
    audit_h3,
    audit_h4,
    finite_cluster_one_sided_p,
    finite_cluster_two_sided_p,
)


def test_finite_cluster_p_is_more_conservative_than_normal_for_small_g():
    p_small = finite_cluster_two_sided_p(2.0, 1.0, 20)
    p_large = finite_cluster_two_sided_p(2.0, 1.0, 10000)
    assert p_small > p_large


def test_one_sided_direction_is_respected():
    assert finite_cluster_one_sided_p(
        2.0, 1.0, 20, alternative="positive"
    ) < 0.05
    assert finite_cluster_one_sided_p(
        -2.0, 1.0, 20, alternative="positive"
    ) > 0.5
    assert finite_cluster_one_sided_p(
        -2.0, 1.0, 20, alternative="negative"
    ) < 0.05


def test_h2_recomputes_same_eight_test_fdr_family_per_scope():
    rows = []
    for context in ["a", "b", "c", "d"]:
        rows.append(
            {
                "evidence_scope": "all_analysis_eligible",
                "context": context,
                "response": "generalized_accessible",
                "status": "fit",
                "n_clusters": 40,
                "distance_estimate": 0.2,
                "distance_se": 0.05,
            }
        )
        rows.append(
            {
                "evidence_scope": "all_analysis_eligible",
                "context": context,
                "response": "plain_colour",
                "status": "fit",
                "n_clusters": 40,
                "distance_estimate": 0.01,
                "distance_se": 0.05,
            }
        )
        rows.append(
            {
                "evidence_scope": "all_analysis_eligible",
                "context": context,
                "response": "selfing_core",
                "status": "fit",
                "n_clusters": 40,
                "distance_estimate": 0.2,
                "distance_se": 0.05,
            }
        )
    out = audit_h2(pd.DataFrame(rows))
    primary = out.loc[
        out["response"].isin(
            ["generalized_accessible", "plain_colour"]
        )
    ]
    assert len(primary) == 8
    assert (
        primary.loc[
            primary["response"].eq("generalized_accessible"),
            "finite_cluster_H2b_supported",
        ].all()
    )
    assert not (
        primary.loc[
            primary["response"].eq("plain_colour"),
            "finite_cluster_H2b_supported",
        ].any()
    )


def test_h3_uses_publications_as_finite_cluster_df():
    payload = {
        "corrected": {
            "global_gradient": {
                "distance_slope": 0.1,
                "distance_slope_se": 0.04,
                "two_sided_p": 0.01,
                "n_publications": 100,
            },
            "sensitivities": {
                "supplemental_only": {
                    "distance_slope": 0.02,
                    "distance_slope_se": 0.05,
                    "two_sided_p": 0.7,
                    "n_publications": 50,
                }
            },
        }
    }
    out = audit_h3(payload)
    primary = out.loc[out["analysis"].eq("primary")].iloc[0]
    assert primary["finite_publication_t_two_sided_p"] < 0.05
    assert primary[
        "finite_publication_t_one_sided_positive_p"
    ] < 0.025


def test_h4_negative_family_direction_is_preserved():
    frame = pd.DataFrame(
        [
            {
                "family": "reproductive_assurance",
                "analysis": "primary",
                "evaluable": True,
                "estimate": -0.30,
                "se": 0.10,
                "n_publications": 100,
                "two_sided_p": 0.003,
            },
            {
                "family": "accessibility_generalization",
                "analysis": "primary",
                "evaluable": True,
                "estimate": -0.29,
                "se": 0.13,
                "n_publications": 100,
                "two_sided_p": 0.02,
            },
        ]
    )
    out = audit_h4(frame)
    assert (
        out["finite_publication_t_one_sided_negative_p"] < 0.05
    ).all()
