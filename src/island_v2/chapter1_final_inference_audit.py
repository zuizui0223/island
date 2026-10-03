"""Reviewer-facing finite-cluster inference audit for Chapter 1 H2-H4.

This module does not refit or redefine any biological estimand.  It replaces the
large-sample normal reference distribution used for already-computed cluster-robust
standard errors with Student-t reference distributions based on the number of
independent clusters minus one.

H2 uses spatial blocks as clusters.  H3 and H4 use publications.  H2b multiplicity is
recomputed with the same frozen eight-test BH family after replacing normal p-values
with finite-cluster t p-values.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
from scipy.stats import t as student_t

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _bh(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=values.index, dtype=float)
    ok = p.notna()
    if not bool(ok.any()):
        return out
    x = p.loc[ok].to_numpy(float)
    order = np.argsort(x)
    ranked = x[order]
    n = len(ranked)
    adjusted = np.minimum.accumulate(
        (ranked * n / np.arange(1, n + 1))[::-1]
    )[::-1]
    restored = np.empty(n, dtype=float)
    restored[order] = np.clip(adjusted, 0.0, 1.0)
    out.loc[ok] = restored
    return out


def finite_cluster_two_sided_p(
    estimate: float,
    se: float,
    n_clusters: int,
) -> float:
    if (
        not math.isfinite(float(estimate))
        or not math.isfinite(float(se))
        or float(se) <= 0
        or int(n_clusters) <= 1
    ):
        return float("nan")
    statistic = abs(float(estimate) / float(se))
    return float(2.0 * student_t.sf(statistic, df=int(n_clusters) - 1))


def finite_cluster_one_sided_p(
    estimate: float,
    se: float,
    n_clusters: int,
    *,
    alternative: str,
) -> float:
    if (
        not math.isfinite(float(estimate))
        or not math.isfinite(float(se))
        or float(se) <= 0
        or int(n_clusters) <= 1
    ):
        return float("nan")
    statistic = float(estimate) / float(se)
    if alternative == "positive":
        return float(student_t.sf(statistic, df=int(n_clusters) - 1))
    if alternative == "negative":
        return float(student_t.cdf(statistic, df=int(n_clusters) - 1))
    raise ValueError(f"unsupported alternative: {alternative}")


def audit_h2(frame: pd.DataFrame) -> pd.DataFrame:
    required = {
        "evidence_scope",
        "context",
        "response",
        "status",
        "n_clusters",
        "distance_estimate",
        "distance_se",
    }
    if missing := required - set(frame.columns):
        raise typer.BadParameter(
            f"H2 table missing columns: {sorted(missing)}"
        )
    out = frame.copy()
    out["distance_p_finite_cluster_t"] = np.nan
    fit = out["status"].astype(str).eq("fit")
    for idx in out.index[fit]:
        out.loc[idx, "distance_p_finite_cluster_t"] = (
            finite_cluster_two_sided_p(
                float(out.loc[idx, "distance_estimate"]),
                float(out.loc[idx, "distance_se"]),
                int(out.loc[idx, "n_clusters"]),
            )
        )
    out["primary_H2b_q_finite_cluster_t"] = np.nan
    primary = (
        fit
        & out["response"].astype(str).isin(
            ["plain_colour", "generalized_accessible"]
        )
    )
    for scope in out.loc[primary, "evidence_scope"].dropna().unique():
        mask = primary & out["evidence_scope"].astype(str).eq(str(scope))
        out.loc[mask, "primary_H2b_q_finite_cluster_t"] = _bh(
            out.loc[mask, "distance_p_finite_cluster_t"]
        )
    out["finite_cluster_H2b_supported"] = (
        out["primary_H2b_q_finite_cluster_t"].le(0.05).fillna(False)
    )
    return out


def audit_h3(payload: dict[str, Any]) -> pd.DataFrame:
    corrected = payload["corrected"]
    analyses = {
        "primary": corrected["global_gradient"],
        **corrected["sensitivities"],
    }
    rows: list[dict[str, Any]] = []
    for analysis, item in analyses.items():
        estimate = float(item["distance_slope"])
        se = float(item["distance_slope_se"])
        clusters = int(item["n_publications"])
        rows.append(
            {
                "analysis": analysis,
                "estimate": estimate,
                "se": se,
                "n_publications": clusters,
                "asymptotic_two_sided_p": float(item["two_sided_p"]),
                "finite_publication_t_two_sided_p":
                    finite_cluster_two_sided_p(
                        estimate, se, clusters
                    ),
                "finite_publication_t_one_sided_positive_p":
                    finite_cluster_one_sided_p(
                        estimate,
                        se,
                        clusters,
                        alternative="positive",
                    ),
            }
        )
    return pd.DataFrame(rows)


def audit_h4(frame: pd.DataFrame) -> pd.DataFrame:
    required = {
        "family",
        "analysis",
        "evaluable",
        "estimate",
        "se",
        "n_publications",
        "two_sided_p",
    }
    if missing := required - set(frame.columns):
        raise typer.BadParameter(
            f"H4 table missing columns: {sorted(missing)}"
        )
    out = frame.copy()
    out["finite_publication_t_two_sided_p"] = np.nan
    out["finite_publication_t_one_sided_negative_p"] = np.nan
    evaluable = (
        out["evaluable"]
        .fillna(False)
        .astype(str)
        .str.lower()
        .isin({"true", "1", "yes"})
    )
    for idx in out.index[evaluable]:
        estimate = float(out.loc[idx, "estimate"])
        se = float(out.loc[idx, "se"])
        clusters = int(out.loc[idx, "n_publications"])
        out.loc[idx, "finite_publication_t_two_sided_p"] = (
            finite_cluster_two_sided_p(
                estimate,
                se,
                clusters,
            )
        )
        out.loc[
            idx,
            "finite_publication_t_one_sided_negative_p",
        ] = finite_cluster_one_sided_p(
            estimate,
            se,
            clusters,
            alternative="negative",
        )
    return out


def final_decision_summary(
    h2_all: pd.DataFrame,
    h2_direct: pd.DataFrame,
    h3: pd.DataFrame,
    h4: pd.DataFrame,
) -> dict[str, Any]:
    h2 = pd.concat([h2_all, h2_direct], ignore_index=True)
    h2_primary = h2.loc[
        h2["response"].isin(
            ["plain_colour", "generalized_accessible"]
        )
    ].copy()
    h3_primary = h3.loc[h3["analysis"].eq("primary")].iloc[0]
    h4_primary = h4.loc[h4["analysis"].eq("primary")].copy()
    return {
        "H2": {
            "n_primary_H2b_tests": int(len(h2_primary)),
            "n_finite_cluster_FDR_supported": int(
                h2_primary["finite_cluster_H2b_supported"].sum()
            ),
            "supported_cells": sorted(
                (
                    h2_primary.loc[
                        h2_primary["finite_cluster_H2b_supported"],
                        ["evidence_scope", "context", "response"],
                    ]
                    .astype(str)
                    .agg("|".join, axis=1)
                    .tolist()
                )
            ),
            "claim_boundary":
                "conditional decomposition, not causal mediation",
        },
        "H3": {
            "estimate": float(h3_primary["estimate"]),
            "finite_publication_t_two_sided_p": float(
                h3_primary["finite_publication_t_two_sided_p"]
            ),
            "supported": bool(
                h3_primary["finite_publication_t_two_sided_p"] <= 0.05
                and h3_primary["estimate"] > 0
            ),
            "claim_boundary":
                "independent ecological-pressure correlate, not historical mediation",
        },
        "H4": {
            "families": {
                str(row["family"]): {
                    "estimate": float(row["estimate"]),
                    "finite_publication_t_two_sided_p": float(
                        row["finite_publication_t_two_sided_p"]
                    ),
                    "supported": bool(
                        row["estimate"] < 0
                        and row["finite_publication_t_two_sided_p"] <= 0.05
                    ),
                }
                for _, row in h4_primary.iterrows()
            },
            "claim_boundary":
                "post-hoc exact-species functional triangulation only",
        },
    }


@app.command("run")
def run(
    h2_all_csv: Path = typer.Option(..., exists=True),
    h2_direct_csv: Path = typer.Option(..., exists=True),
    h3_json: Path = typer.Option(..., exists=True),
    h4_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    h2_all = audit_h2(pd.read_csv(h2_all_csv))
    h2_direct = audit_h2(pd.read_csv(h2_direct_csv))
    h3 = audit_h3(
        json.loads(h3_json.read_text(encoding="utf-8"))
    )
    h4 = audit_h4(pd.read_csv(h4_csv))
    summary = final_decision_summary(
        h2_all,
        h2_direct,
        h3,
        h4,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    h2_all.to_csv(
        output_dir / "h2_all_finite_cluster_audit.csv",
        index=False,
    )
    h2_direct.to_csv(
        output_dir / "h2_direct_finite_cluster_audit.csv",
        index=False,
    )
    h3.to_csv(
        output_dir / "h3_finite_publication_audit.csv",
        index=False,
    )
    h4.to_csv(
        output_dir / "h4_finite_publication_audit.csv",
        index=False,
    )
    (output_dir / "h2_h4_final_decision_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
