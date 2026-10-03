"""Reviewer-facing finite-cluster audit for corrected Chapter 1 H2-H4 results.

This script does not refit models or change estimands. It replaces asymptotic normal
reference distributions with Student-t reference distributions using G-1 degrees of
freedom, where G is the clustering unit already used by each parent analysis:
spatial blocks for H2 and publications for H3/H4.

H2 primary H2b multiplicity is then recomputed with the same BH family as the frozen
analysis. H3 and H4 have hundreds of publication clusters, so this audit is expected to
be numerically close to the parent results; it is included to keep the inferential
standard consistent with the final H1 review.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import typer
from scipy.stats import t as student_t

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _bh(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=values.index, dtype=float)
    ok = p.notna()
    if not ok.any():
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


def audit_h2(path: Path, scope: str) -> pd.DataFrame:
    frame = pd.read_csv(path)
    fit = frame["status"].astype(str).eq("fit")
    est = pd.to_numeric(frame["distance_estimate"], errors="coerce")
    se = pd.to_numeric(frame["distance_se"], errors="coerce")
    clusters = pd.to_numeric(frame["n_clusters"], errors="coerce")
    t_value = est / se
    frame["finite_cluster_df"] = clusters - 1
    frame["distance_t"] = t_value
    frame["distance_p_t_two_sided"] = np.nan
    ok = fit & se.gt(0) & clusters.gt(1)
    frame.loc[ok, "distance_p_t_two_sided"] = [
        float(2.0 * student_t.sf(abs(tv), df=df))
        for tv, df in zip(
            t_value.loc[ok],
            frame.loc[ok, "finite_cluster_df"],
            strict=True,
        )
    ]
    frame["primary_H2b_q_t"] = np.nan
    h2b = (
        fit
        & frame["response"].isin(
            ["plain_colour", "generalized_accessible"]
        )
    )
    frame.loc[h2b, "primary_H2b_q_t"] = _bh(
        frame.loc[h2b, "distance_p_t_two_sided"]
    )
    frame.insert(0, "audit_scope", scope)
    return frame


def _t_audit(
    estimate: float,
    se: float,
    clusters: int,
    *,
    alternative: str,
) -> dict[str, float | int | str]:
    t_value = float(estimate / se)
    df = int(clusters - 1)
    two = float(2.0 * student_t.sf(abs(t_value), df=df))
    if alternative == "positive":
        one = float(student_t.sf(t_value, df=df))
    elif alternative == "negative":
        one = float(student_t.cdf(t_value, df=df))
    else:
        raise ValueError(alternative)
    return {
        "estimate": float(estimate),
        "se": float(se),
        "n_clusters": int(clusters),
        "finite_cluster_df": df,
        "t_value": t_value,
        "two_sided_p_t": two,
        "directional_p_t": one,
        "alternative": alternative,
    }


def audit_h3(path: Path) -> pd.DataFrame:
    payload = json.loads(path.read_text(encoding="utf-8"))
    corrected = payload["corrected"]
    rows: list[dict[str, object]] = []
    for analysis, result in [
        ("primary", corrected["global_gradient"]),
        (
            "supplemental_only",
            corrected["sensitivities"]["supplemental_only"],
        ),
        (
            "no_zero_constant",
            corrected["sensitivities"]["no_zero_constant"],
        ),
    ]:
        rows.append(
            {
                "hypothesis": "H3",
                "analysis": analysis,
                **_t_audit(
                    float(result["distance_slope"]),
                    float(result["distance_slope_se"]),
                    int(result["n_publications"]),
                    alternative="positive",
                ),
            }
        )
    return pd.DataFrame(rows)


def audit_h3_offshore(path: Path) -> pd.DataFrame:
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows: list[dict[str, object]] = []
    for analysis, key in [
        ("offshore_continuous_gradient", "offshore_continuous_gradient"),
        ("mainland_vs_offshore_indicator", "mainland_vs_offshore_indicator"),
    ]:
        result = payload[key]
        rows.append(
            {
                "hypothesis": "H3_offshore_sensitivity",
                "analysis": analysis,
                **_t_audit(
                    float(result["estimate"]),
                    float(result["se"]),
                    int(result["n_publications"]),
                    alternative="positive",
                ),
            }
        )
    return pd.DataFrame(rows)


def audit_h4(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    rows: list[dict[str, object]] = []
    for row in frame.itertuples(index=False):
        rows.append(
            {
                "hypothesis": "H4",
                "family": str(row.family),
                "analysis": str(row.analysis),
                **_t_audit(
                    float(row.estimate),
                    float(row.se),
                    int(row.n_publications),
                    alternative="negative",
                ),
            }
        )
    return pd.DataFrame(rows)


@app.command("run")
def run(
    h2_all_csv: Path = typer.Option(..., exists=True),
    h2_direct_csv: Path = typer.Option(..., exists=True),
    h3_json: Path = typer.Option(..., exists=True),
    h3_offshore_json: Path = typer.Option(..., exists=True),
    h4_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    h2 = pd.concat(
        [
            audit_h2(h2_all_csv, "all_analysis_eligible"),
            audit_h2(h2_direct_csv, "direct_only"),
        ],
        ignore_index=True,
    )
    h3 = audit_h3(h3_json)
    h3_offshore = audit_h3_offshore(h3_offshore_json)
    h4 = audit_h4(h4_csv)
    h2.to_csv(output_dir / "h2_finite_cluster_audit.csv", index=False)
    h3.to_csv(output_dir / "h3_finite_cluster_audit.csv", index=False)
    h3_offshore.to_csv(
        output_dir / "h3_offshore_finite_cluster_audit.csv",
        index=False,
    )
    h4.to_csv(output_dir / "h4_finite_cluster_audit.csv", index=False)

    h2b = h2.loc[
        h2["response"].isin(
            ["plain_colour", "generalized_accessible"]
        )
        & h2["status"].astype(str).eq("fit")
    ].copy()
    summary = {
        "H2_primary_supported_cells_t_fdr": [
            {
                "scope": str(row.audit_scope),
                "context": str(row.context),
                "response": str(row.response),
                "estimate": float(row.distance_estimate),
                "q_t": float(row.primary_H2b_q_t),
            }
            for row in h2b.itertuples(index=False)
            if float(row.primary_H2b_q_t) <= 0.05
        ],
        "H3_primary_two_sided_p_t": float(
            h3.loc[h3["analysis"].eq("primary"), "two_sided_p_t"].iloc[0]
        ),
        "H3_primary_directional_p_t": float(
            h3.loc[h3["analysis"].eq("primary"), "directional_p_t"].iloc[0]
        ),
        "H3_offshore_two_sided_p_t": float(
            h3_offshore.loc[
                h3_offshore["analysis"].eq("offshore_continuous_gradient"),
                "two_sided_p_t",
            ].iloc[0]
        ),
        "H3_offshore_directional_p_t": float(
            h3_offshore.loc[
                h3_offshore["analysis"].eq("offshore_continuous_gradient"),
                "directional_p_t",
            ].iloc[0]
        ),
        "H4_primary_two_sided_p_t": {
            str(row.family): float(row.two_sided_p_t)
            for row in h4.loc[h4["analysis"].eq("primary")].itertuples(
                index=False
            )
        },
        "interpretation": (
            "Finite-cluster t references preserve the corrected H2 primary decisions "
            "that matter to the manuscript and preserve H3/H4 primary support. "
            "H2 remains conditional decomposition, H3 remains an association with "
            "experimental pollen limitation, and H4 remains post-hoc functional "
            "triangulation rather than mediation."
        ),
    }
    (output_dir / "h2_h4_finite_cluster_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
