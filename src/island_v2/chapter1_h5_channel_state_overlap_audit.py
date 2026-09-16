"""Reporting-only identifiability audit for canonical exact-island channel states."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


def load_config(path: Path) -> dict[str, Any]:
    cfg = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(cfg, dict) or cfg.get("contract") != "chapter1_h5_channel_state_overlap_audit_v1":
        raise typer.BadParameter("unexpected channel-state overlap contract")
    return cfg


def audit_channel(
    table: pd.DataFrame, covariates: pd.DataFrame, channel: str, cfg: dict[str, Any]
) -> pd.DataFrame:
    required = {"island_id", "channel_id", "observation_state"}
    if missing := required - set(table.columns):
        raise typer.BadParameter(f"{channel} table missing columns: {sorted(missing)}")
    context = str(cfg["context_column"])
    if not {"island_id", context}.issubset(covariates.columns):
        raise typer.BadParameter("covariates missing island/context columns")
    mapping = cfg["state_projection"]
    work = table.loc[table["channel_id"].astype(str).eq(channel)].copy()
    work["state"] = work["observation_state"].astype(str).map(mapping)
    work = work.merge(
        covariates[["island_id", context]].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    rows: list[dict[str, Any]] = []
    minimum = int(cfg["primary_overlap_reference"])
    for value in sorted(work[context].dropna().astype(str).unique()):
        part = work.loc[work[context].astype(str).eq(value)]
        counts = part["state"].value_counts()
        retained = int(counts.get("retained", 0))
        disrupted = int(counts.get("disrupted", 0))
        rows.append(
            {
                "channel": channel,
                "context": value,
                "retained": retained,
                "disrupted": disrupted,
                "strict_evaluable": retained + disrupted,
                "insufficient_effort": int(part["observation_state"].eq("insufficient_effort").sum()),
                "unresolved": int(part["observation_state"].eq("unresolved").sum()),
                "passes_reference_overlap": retained >= minimum and disrupted >= minimum,
            }
        )
    return pd.DataFrame(rows)


def summarize(matrix: pd.DataFrame, cfg: dict[str, Any]) -> dict[str, Any]:
    primary = [str(x) for x in cfg["primary_contexts"]]
    channels = [str(x) for x in cfg["channels"]]
    flags: dict[str, bool] = {}
    for channel in channels:
        part = matrix.loc[
            matrix["channel"].eq(channel) & matrix["context"].isin(primary)
        ]
        flags[channel] = bool(
            set(part["context"]) == set(primary)
            and part["passes_reference_overlap"].all()
        )
    return {
        "contract": cfg["contract"],
        "primary_contexts": primary,
        "primary_overlap_reference": int(cfg["primary_overlap_reference"]),
        "channel_primary_overlap_pass": flags,
        "n_channels_passing_both_primary_contexts": int(sum(flags.values())),
        "n_channels_total": len(channels),
        "mechanism_promoted": False,
        "interpretation": "Strict channel-state overlap is an identifiability audit only.",
    }


@app.command("run")
def run(
    bombus_csv: Path = typer.Option(..., exists=True),
    non_bombus_bees_csv: Path = typer.Option(..., exists=True),
    lepidoptera_csv: Path = typer.Option(..., exists=True),
    flower_visiting_birds_csv: Path = typer.Option(..., exists=True),
    diptera_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    cfg = load_config(config_path)
    cov = pd.read_csv(covariates_csv, dtype={"island_id": str})
    paths = {
        "bombus": bombus_csv,
        "non_bombus_bees": non_bombus_bees_csv,
        "lepidoptera": lepidoptera_csv,
        "flower_visiting_birds": flower_visiting_birds_csv,
        "diptera": diptera_csv,
    }
    parts = [
        audit_channel(pd.read_csv(paths[ch], dtype={"island_id": str}), cov, ch, cfg)
        for ch in cfg["channels"]
    ]
    matrix = pd.concat(parts, ignore_index=True)
    decision = summarize(matrix, cfg)
    output_dir.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(output_dir / "channel_state_context_overlap.csv", index=False)
    (output_dir / "channel_state_overlap_decision.json").write_text(
        json.dumps(decision, indent=2) + "\n", encoding="utf-8"
    )
    primary = matrix.loc[matrix["context"].isin(cfg["primary_contexts"])].copy()
    lines = [
        "# Five-channel strict state-overlap audit",
        "",
        f"Reference overlap = at least {cfg['primary_overlap_reference']} retained and disrupted islands per context.",
        "",
        "| channel | context | retained | disrupted | passes |",
        "|---|---|---:|---:|---|",
    ]
    for row in primary.to_dict("records"):
        lines.append(
            f"| {row['channel']} | {row['context']} | {row['retained']} | {row['disrupted']} | {bool(row['passes_reference_overlap'])} |"
        )
    lines.extend([
        "",
        f"**Channels passing both primary contexts: {decision['n_channels_passing_both_primary_contexts']}/{decision['n_channels_total']}.**",
        "",
        "> Failure means the current strict occurrence states cannot identify the channel-disruption mechanism for the primary contrast; it does not mean the channel is biologically irrelevant.",
    ])
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo(json.dumps(decision, indent=2))


if __name__ == "__main__":
    app()
