"""Freeze Chapter 1 regional trait vectors before pollination-syndrome interpretation.

Formal biogeographic contingency is defined by the cluster-robust joint Wald test
of all isolation x context interactions, with FDR correction across atomic states
within each stratum x trait domain. A significant slope in one context and a
non-significant slope in another is never, by itself, labelled a context difference.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _bh(p: pd.Series) -> pd.Series:
    values = pd.to_numeric(p, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    ok = values.notna()
    if not ok.any():
        return out
    x = values.loc[ok].to_numpy(dtype=float)
    order = np.argsort(x)
    ranked = x[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    restored = np.empty(n, dtype=float)
    restored[order] = adjusted
    out.loc[ok] = restored
    return out


def _direction(estimate: float, q: float, p: float) -> str:
    if not np.isfinite(estimate):
        return "unresolved"
    sign = "positive" if estimate > 0 else "negative" if estimate < 0 else "zero"
    if np.isfinite(q) and q <= 0.05:
        return f"fdr_{sign}"
    if np.isfinite(p) and p <= 0.05:
        return f"nominal_{sign}"
    return f"uncertain_{sign}"


def freeze_trait_vectors(
    slopes: pd.DataFrame,
    coefficients: pd.DataFrame,
    lineage: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    required_slopes = {
        "stratum",
        "trait_name",
        "trait_state",
        "context",
        "isolation_slope_log_odds_per_sd",
        "p_value",
        "n_context_islands",
        "interaction_joint_wald_chisq",
        "interaction_joint_df",
        "interaction_joint_p",
    }
    missing = required_slopes - set(slopes.columns)
    if missing:
        raise typer.BadParameter(f"slope table missing columns: {sorted(missing)}")

    vector = slopes.copy()
    vector["q_within_stratum_trait_context"] = vector.groupby(
        ["stratum", "trait_name", "context"], group_keys=False
    )["p_value"].apply(_bh)
    vector["direction_class"] = [
        _direction(est, q, p)
        for est, q, p in zip(
            pd.to_numeric(vector["isolation_slope_log_odds_per_sd"], errors="coerce"),
            pd.to_numeric(vector["q_within_stratum_trait_context"], errors="coerce"),
            pd.to_numeric(vector["p_value"], errors="coerce"),
            strict=True,
        )
    ]
    vector["fdr_supported"] = vector["q_within_stratum_trait_context"].le(0.05)

    # One joint interaction test per atomic trait category. FDR is applied across
    # states within a stratum x trait domain; this is the formal when/where test.
    joint = (
        vector[
            [
                "stratum",
                "trait_name",
                "trait_state",
                "interaction_joint_wald_chisq",
                "interaction_joint_df",
                "interaction_joint_p",
            ]
        ]
        .drop_duplicates(["stratum", "trait_name", "trait_state"])
        .reset_index(drop=True)
    )
    joint["q_joint_within_stratum_trait"] = joint.groupby(
        ["stratum", "trait_name"], group_keys=False
    )["interaction_joint_p"].apply(_bh)
    joint["contingency_fdr_supported"] = joint["q_joint_within_stratum_trait"].le(0.05)
    vector = vector.merge(
        joint[
            [
                "stratum",
                "trait_name",
                "trait_state",
                "q_joint_within_stratum_trait",
                "contingency_fdr_supported",
            ]
        ],
        on=["stratum", "trait_name", "trait_state"],
        how="left",
        validate="many_to_one",
    )

    # Individual interaction contrasts are secondary decomposition of a joint test.
    interaction = coefficients.loc[
        coefficients.get("model", pd.Series(index=coefficients.index, dtype=object)).eq("M2")
        & coefficients.get("predictor", pd.Series(index=coefficients.index, dtype=object))
        .astype(str)
        .str.startswith("z_isolation:context[")
    ].copy()
    if not interaction.empty:
        interaction["contrast_context"] = interaction["predictor"].str.extract(
            r"z_isolation:context\[(.*)\]", expand=False
        )
        interaction["q_contrast_within_stratum_trait"] = interaction.groupby(
            ["stratum", "trait_name"], group_keys=False
        )["p_value"].apply(_bh)
        interaction["contrast_fdr_supported"] = interaction[
            "q_contrast_within_stratum_trait"
        ].le(0.05)
        interaction = interaction.merge(
            joint,
            on=["stratum", "trait_name", "trait_state"],
            how="left",
            validate="many_to_one",
            suffixes=("", "_joint"),
        )

    # M3 is a broad lineage guardrail, not an atomic-category replacement.
    lineage_out = lineage.copy()
    if not lineage_out.empty and "p_value" in lineage_out.columns:
        lineage_out["q_within_stratum_context"] = lineage_out.groupby(
            ["stratum", "context"], group_keys=False
        )["p_value"].apply(_bh)
        lineage_out["lineage_direction_class"] = [
            _direction(est, q, p)
            for est, q, p in zip(
                pd.to_numeric(lineage_out["estimate"], errors="coerce"),
                pd.to_numeric(lineage_out["q_within_stratum_context"], errors="coerce"),
                pd.to_numeric(lineage_out["p_value"], errors="coerce"),
                strict=True,
            )
        ]

    pattern_rows: list[dict[str, object]] = []
    for keys, group in vector.groupby(["stratum", "trait_name", "trait_state"], sort=True):
        supported = group.loc[group["fdr_supported"]].copy()
        nominal = group.loc[pd.to_numeric(group["p_value"], errors="coerce").le(0.05)]
        j = joint.loc[
            joint["stratum"].eq(keys[0])
            & joint["trait_name"].eq(keys[1])
            & joint["trait_state"].eq(keys[2])
        ].iloc[0]
        if bool(j["contingency_fdr_supported"]):
            cls = "fdr_supported_biogeographic_contingency"
        elif pd.to_numeric(pd.Series([j["interaction_joint_p"]]), errors="coerce").iloc[0] <= 0.05:
            cls = "nominal_biogeographic_contingency"
        elif len(supported) >= 1:
            cls = "regional_slope_signal_without_supported_heterogeneity"
        elif len(nominal) >= 1:
            cls = "nominal_slope_signal_without_supported_heterogeneity"
        else:
            cls = "no_supported_isolation_or_contingency_signal"
        pattern_rows.append(
            {
                "stratum": keys[0],
                "trait_name": keys[1],
                "trait_state": keys[2],
                "pattern_class": cls,
                "interaction_joint_wald_chisq": float(j["interaction_joint_wald_chisq"]),
                "interaction_joint_df": int(j["interaction_joint_df"]),
                "interaction_joint_p": float(j["interaction_joint_p"]),
                "q_joint_within_stratum_trait": float(j["q_joint_within_stratum_trait"]),
                "contingency_fdr_supported": bool(j["contingency_fdr_supported"]),
                "n_contexts_tested": int(group["context"].nunique()),
                "n_fdr_supported_contexts": int(len(supported)),
                "n_nominal_supported_contexts": int(len(nominal)),
                "contexts": "|".join(sorted(group["context"].astype(str).unique())),
            }
        )
    patterns = pd.DataFrame(pattern_rows)
    return vector, interaction, lineage_out, patterns


@app.command("run")
def run(
    slopes_csv: Path = typer.Option(..., exists=True),
    coefficients_csv: Path = typer.Option(..., exists=True),
    lineage_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    vector, interaction, lineage, patterns = freeze_trait_vectors(
        pd.read_csv(slopes_csv),
        pd.read_csv(coefficients_csv),
        pd.read_csv(lineage_csv),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    vector.to_csv(output_dir / "regional_trait_vector_frozen.csv", index=False)
    interaction.to_csv(output_dir / "regional_context_interactions_frozen.csv", index=False)
    lineage.to_csv(output_dir / "lineage_guardrail_frozen.csv", index=False)
    patterns.to_csv(output_dir / "regional_pattern_classification_frozen.csv", index=False)
    manifest = {
        "contract": "chapter1_regional_trait_vector_freeze_v2_joint_wald",
        "atomic_trait_categories": True,
        "pollination_syndrome_labels_included": False,
        "primary_contingency_test": "cluster-robust joint Wald; BH across atomic states within stratum x trait",
        "simple_slope_fdr": "BH within stratum x trait x context; descriptive within-context isolation signal",
        "individual_interaction_fdr": "secondary contrast decomposition only",
        "lineage_guardrail": "M3 genus-composition-preserving broad outcomes retained separately",
        "n_trait_vector_rows": int(len(vector)),
        "n_interaction_rows": int(len(interaction)),
        "n_pattern_rows": int(len(patterns)),
        "n_fdr_supported_contingency_categories": (
            int(patterns["contingency_fdr_supported"].sum()) if not patterns.empty else 0
        ),
    }
    (output_dir / "regional_trait_vector_freeze_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
