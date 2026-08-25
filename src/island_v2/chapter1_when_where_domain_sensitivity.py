"""Domain-level sensitivity for the Chapter 1 when/where response-vector tests."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import typer
import yaml

from island_v2.chapter1_when_where_omnibus import run_when_where_omnibus

app = typer.Typer(add_completion=False, no_args_is_help=True)


def run_domain_sensitivity(
    composition: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if "trait_name" not in composition.columns:
        raise typer.BadParameter("composition missing trait_name")
    within_parts: list[pd.DataFrame] = []
    between_parts: list[pd.DataFrame] = []
    persistence_parts: list[pd.DataFrame] = []
    for domain in sorted(composition["trait_name"].dropna().astype(str).unique()):
        subset = composition.loc[composition["trait_name"].astype(str).eq(domain)].copy()
        within, between, persistence = run_when_where_omnibus(subset, covariates, config)
        for table, collector in (
            (within, within_parts),
            (between, between_parts),
            (persistence, persistence_parts),
        ):
            if not table.empty:
                table.insert(0, "trait_domain", domain)
                collector.append(table)
    return (
        pd.concat(within_parts, ignore_index=True) if within_parts else pd.DataFrame(),
        pd.concat(between_parts, ignore_index=True) if between_parts else pd.DataFrame(),
        pd.concat(persistence_parts, ignore_index=True) if persistence_parts else pd.DataFrame(),
    )


@app.command("run")
def run(
    composition_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    within, between, persistence = run_domain_sensitivity(
        pd.read_csv(composition_csv),
        pd.read_csv(covariates_csv),
        config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    within.to_csv(output_dir / "when_where_domain_within.csv", index=False)
    between.to_csv(output_dir / "when_where_domain_between.csv", index=False)
    persistence.to_csv(output_dir / "when_where_domain_persistence.csv", index=False)


if __name__ == "__main__":
    app()
