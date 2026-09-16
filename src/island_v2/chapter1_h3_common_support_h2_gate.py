"""Re-fit the original all-data beta-binomial H2 model on H3 common taxonomic support."""
from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.chapter1_all_data_probability import run_probability_analysis

app = typer.Typer(add_completion=False, no_args_is_help=True)


def counts_from_common_scores(scores: pd.DataFrame) -> pd.DataFrame:
    required = {"island_id", "outcome", "observed_score", "n_species"}
    if missing := required - set(scores.columns):
        raise typer.BadParameter(f"H3 island scores missing columns: {sorted(missing)}")
    out = scores[["island_id", "outcome", "observed_score", "n_species"]].copy()
    out["observed_score"] = pd.to_numeric(out["observed_score"], errors="coerce")
    out["n_species"] = pd.to_numeric(out["n_species"], errors="coerce")
    out = out.dropna(subset=["observed_score", "n_species"])
    out = out.loc[out["n_species"].gt(0)].copy()
    out["trials"] = out["n_species"]
    out["successes"] = out["observed_score"] * out["trials"]
    out["stratum"] = "all_observed"
    return out[["island_id", "successes", "trials", "outcome", "stratum"]]


def gate_config(probability_config: dict[str, Any], depth_config: dict[str, Any]) -> dict[str, Any]:
    config = deepcopy(probability_config)
    contexts = [str(x) for x in depth_config["contexts"]]
    config["contexts"] = contexts
    config["between_contexts"] = [contexts]
    config["strata"] = ["all_observed"]
    config["model_outcomes"] = [str(x) for x in depth_config["atomic_outcomes"]]
    return config


@app.command("run")
def run(
    island_scores_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    probability_config_path: Path = typer.Option(..., exists=True),
    depth_config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    scores = pd.read_csv(island_scores_csv)
    covariates = pd.read_csv(covariates_csv)
    probability_config = yaml.safe_load(probability_config_path.read_text(encoding="utf-8"))
    depth_config = yaml.safe_load(depth_config_path.read_text(encoding="utf-8"))
    counts = counts_from_common_scores(scores)
    config = gate_config(probability_config, depth_config)
    within_slopes, between_slopes, within, between = run_probability_analysis(
        counts, covariates, config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output_dir / "h3_common_support_counts.csv.gz", index=False, compression="gzip")
    within_slopes.to_csv(output_dir / "h3_common_support_within_slopes.csv", index=False)
    between_slopes.to_csv(output_dir / "h3_common_support_between_slopes.csv", index=False)
    within.to_csv(output_dir / "h3_common_support_within_omnibus.csv", index=False)
    between.to_csv(output_dir / "h3_common_support_between_omnibus.csv", index=False)
    row = between.iloc[0].to_dict() if len(between) else {}
    manifest = {
        "contract": "chapter1_h3_common_support_h2_gate_v1",
        "evidence_scope": evidence_scope,
        "model_family": "same_beta_binomial_as_all_data_H2",
        "flora_scope": "all_observed",
        "species_scope": "H3_family_and_genus_LOO_common_support",
        "status": str(row.get("status", "not_testable")),
        "p_value": row.get("p_value"),
        "n_unique_islands": row.get("n_unique_islands"),
        "n_clusters": row.get("n_clusters"),
        "gate_passed": bool(row.get("status") == "fit" and float(row.get("p_value", 1.0)) < 0.05),
    }
    (output_dir / "h3_common_support_h2_gate.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
