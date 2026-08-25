"""Lower-dimensional domain sensitivity for Chapter 1 when/where inference."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import typer
import yaml

from island_v2.chapter1_when_where_omnibus import _bh, run_when_where_omnibus

app = typer.Typer(add_completion=False, no_args_is_help=True)


def run_domain_omnibus(
    composition: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if "trait_name" not in composition.columns:
        raise typer.BadParameter("composition missing trait_name")
    within_parts: list[pd.DataFrame] = []
    between_parts: list[pd.DataFrame] = []
    for domain in sorted(composition["trait_name"].dropna().astype(str).unique()):
        subset = composition.loc[composition["trait_name"].astype(str).eq(domain)].copy()
        within, between, _ = run_when_where_omnibus(subset, covariates, config)
        if not within.empty:
            within.insert(0, "trait_domain", domain)
            within_parts.append(within)
        if not between.empty:
            between.insert(0, "trait_domain", domain)
            between_parts.append(between)

    within_out = pd.concat(within_parts, ignore_index=True) if within_parts else pd.DataFrame()
    between_out = pd.concat(between_parts, ignore_index=True) if between_parts else pd.DataFrame()
    if not within_out.empty:
        within_out["q_across_domains_contexts"] = within_out.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["p_value"].apply(_bh)
        within_out["domain_where_supported"] = within_out["q_across_domains_contexts"].le(0.05)
    if not between_out.empty:
        between_out["q_across_domains_pairs"] = between_out.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["p_value"].apply(_bh)
        between_out["domain_between_supported"] = between_out["q_across_domains_pairs"].le(0.05)
    return within_out, between_out


@app.command("run")
def run(
    composition_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    within, between = run_domain_omnibus(
        pd.read_csv(composition_csv), pd.read_csv(covariates_csv), config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    within.to_csv(output_dir / "when_where_domain_within.csv", index=False)
    between.to_csv(output_dir / "when_where_domain_between.csv", index=False)


if __name__ == "__main__":
    app()
