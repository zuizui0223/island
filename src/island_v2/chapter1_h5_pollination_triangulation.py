"""Integrate plant-side syndrome consistency with independent H5 pollination tests."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _prediction_match(estimate: float, prediction: str) -> bool:
    if prediction == "negative":
        return estimate < 0
    if prediction in {"positive", "nonnegative"}:
        return estimate >= 0
    return True


def build_template_consistency(
    all_slopes: pd.DataFrame,
    direct_slopes: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    spec = config["plant_side_consistency"]
    rows: list[dict[str, Any]] = []
    scopes = {
        "all_analysis_eligible": all_slopes,
        "direct_only": direct_slopes,
    }
    for scope, frame in scopes.items():
        target = frame.loc[
            frame["context_layer"].astype(str).eq(str(spec["context_layer"]))
            & frame["axis_set"].astype(str).eq("sampled_guild_concordance")
            & frame["support_tier"].astype(str).eq("confirmatory")
        ].copy()
        for context, predictions in spec["predictions"].items():
            for syndrome, prediction in predictions.items():
                if prediction == "no_primary_prediction":
                    continue
                for stratum in [str(x) for x in spec["strata"]]:
                    hit = target.loc[
                        target["context"].astype(str).eq(str(context))
                        & target["stratum"].astype(str).eq(stratum)
                        & target["syndrome"].astype(str).eq(str(syndrome))
                    ]
                    if len(hit) != 1:
                        rows.append(
                            {
                                "evidence_scope": scope,
                                "context": context,
                                "stratum": stratum,
                                "syndrome": syndrome,
                                "prediction": prediction,
                                "status": "not_testable",
                            }
                        )
                        continue
                    record = hit.iloc[0]
                    estimate = float(record["distance_slope"])
                    q = float(record["q_axis_family"])
                    rows.append(
                        {
                            "evidence_scope": scope,
                            "context": context,
                            "stratum": stratum,
                            "syndrome": syndrome,
                            "prediction": prediction,
                            "status": "fit",
                            "distance_slope": estimate,
                            "q_value": q,
                            "sign_consistent": _prediction_match(estimate, str(prediction)),
                            "fdr_supported": bool(q <= 0.05),
                            "prediction_supported": bool(
                                _prediction_match(estimate, str(prediction)) and q <= 0.05
                            ),
                            "n_islands": int(record["n_islands"]),
                        }
                    )
    return pd.DataFrame(rows)


def architecture_specificity(
    all_factor: pd.DataFrame,
    direct_factor: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    ceiling = float(
        config["plant_side_consistency"]["identity_specificity_guardrail"][
            "architecture_factor_variance_fraction_ceiling"
        ]
    )
    rows = []
    for scope, frame in {
        "all_analysis_eligible": all_factor,
        "direct_only": direct_factor,
    }.items():
        fraction = float(pd.to_numeric(frame["variance_fraction"], errors="coerce").dropna().max())
        rows.append(
            {
                "evidence_scope": scope,
                "common_architecture_variance_fraction": fraction,
                "specificity_ceiling": ceiling,
                "named_pollinator_identity_promotable": bool(fraction <= ceiling),
            }
        )
    return pd.DataFrame(rows)


def integrated_decision(
    consistency: pd.DataFrame,
    specificity: pd.DataFrame,
    globi_classification: pd.DataFrame,
    config: dict[str, Any],
) -> dict[str, Any]:
    primary_scope = str(config["plant_side_consistency"]["primary_evidence_scope"])
    primary = consistency.loc[
        consistency["evidence_scope"].eq(primary_scope) & consistency["status"].eq("fit")
    ].copy()
    predicted_pairs = int(primary[["context", "syndrome"]].drop_duplicates().shape[0])
    robust_pairs = 0
    for _, group in primary.groupby(["context", "syndrome"]):
        if len(group) and group["prediction_supported"].all():
            robust_pairs += 1
    plant_side_consistent = bool(predicted_pairs > 0 and robust_pairs == predicted_pairs)
    identity_specific = bool(specificity["named_pollinator_identity_promotable"].all())
    globi_promoted = bool(
        not globi_classification.empty and globi_classification["promoted"].astype(bool).any()
    )
    existing = config["independent_existing_checks"]
    independent_promoted = bool(
        globi_promoted
        or str(existing["N1_channel_heterogeneity"]["status"]) == "promoted"
        or str(existing["H5c_biotic_vs_wind"]["status"]) == "promoted"
    )
    if independent_promoted:
        mechanism_status = "independent_pollination_association_supported"
    else:
        mechanism_status = "pollination_plausible_but_not_identified"
    return {
        "plant_side_prediction_pairs": predicted_pairs,
        "plant_side_robust_pairs": robust_pairs,
        "plant_side_pollination_architecture_consistent": plant_side_consistent,
        "named_pollinator_identity_specificity_passed": identity_specific,
        "new_globi_area_bridge_promoted": globi_promoted,
        "independent_pollination_support_present": independent_promoted,
        "mechanism_status": mechanism_status,
        "claim": (
            "Plant-side floral responses can be consistent with pollination-associated architecture, "
            "but named syndromes do not identify realized pollinators and no single global pollinator "
            "mechanism is promoted without independent support."
        ),
    }


@app.command("run")
def run(
    all_branch_slopes_csv: Path = typer.Option(..., exists=True),
    direct_branch_slopes_csv: Path = typer.Option(..., exists=True),
    all_factor_csv: Path = typer.Option(..., exists=True),
    direct_factor_csv: Path = typer.Option(..., exists=True),
    globi_classification_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if config.get("contract") != "chapter1_h5_pollination_triangulation_v2":
        raise typer.BadParameter("unexpected H5 triangulation contract")
    consistency = build_template_consistency(
        pd.read_csv(all_branch_slopes_csv),
        pd.read_csv(direct_branch_slopes_csv),
        config,
    )
    specificity = architecture_specificity(
        pd.read_csv(all_factor_csv), pd.read_csv(direct_factor_csv), config
    )
    globi = pd.read_csv(globi_classification_csv)
    decision = integrated_decision(consistency, specificity, globi, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    consistency.to_csv(output_dir / "h5_template_consistency.csv", index=False)
    specificity.to_csv(output_dir / "h5_architecture_specificity.csv", index=False)
    (output_dir / "h5_pollination_triangulation_decision.json").write_text(
        json.dumps(decision, indent=2) + "\n", encoding="utf-8"
    )
    lines = [
        "# Chapter 1 H5 pollination triangulation",
        "",
        f"- plant-side predeclared pairs robust in direct evidence: {decision['plant_side_robust_pairs']}/{decision['plant_side_prediction_pairs']}",
        f"- named pollinator identity specificity passed: {decision['named_pollinator_identity_specificity_passed']}",
        f"- new GloBI area bridge promoted: {decision['new_globi_area_bridge_promoted']}",
        f"- mechanism status: **{decision['mechanism_status']}**",
        "",
        decision["claim"],
    ]
    (output_dir / "H5_RESULT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo(json.dumps(decision, indent=2))


if __name__ == "__main__":
    app()
