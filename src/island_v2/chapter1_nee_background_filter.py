"""Filter broad GBIF acquisition proxies to the frozen N1 observation backgrounds.

Acquisition may deliberately use a broader GBIF taxon than the analytical background
when that reduces the number of global downloads.  This module is the boundary that
prevents the acquisition proxy from silently redefining observation effort.  In
particular, non-Bombus bee acquisition may use Hymenoptera, but only the seven bee
families frozen in ``chapter1_nee_channel_observation_policy.yml`` count as background.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


def load_policy(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("channel observation policy must be a mapping")
    if payload.get("contract") != "chapter1_nee_channel_observation_policy_v1":
        raise typer.BadParameter("unexpected channel observation policy contract")
    return payload


def filter_frozen_background(
    records: pd.DataFrame,
    channel_id: str,
    policy: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Return only records belonging to the predeclared analytical background group."""
    channels = policy.get("channels", {})
    if channel_id not in channels:
        raise ValueError(f"unregistered N1 channel: {channel_id}")

    before = int(len(records))
    channel = channels[channel_id]
    if channel_id == "non_bombus_bees":
        if "family" not in records.columns:
            raise ValueError("non-Bombus bee background filtering requires family")
        allowed = {str(value).strip() for value in channel.get("background_taxa", []) if str(value).strip()}
        if not allowed:
            raise ValueError("frozen non-Bombus bee family set is empty")
        family = records["family"].fillna("").astype(str).str.strip()
        result = records.loc[family.isin(allowed)].copy()
        rule = "family_in_frozen_seven_bee_families"
        allowed_values = sorted(allowed)
    elif channel_id == "bombus":
        if "family" not in records.columns:
            raise ValueError("Bombus background filtering requires family")
        family = records["family"].fillna("").astype(str).str.strip()
        result = records.loc[family.eq("Apidae")].copy()
        rule = "family_exact_Apidae"
        allowed_values = ["Apidae"]
    else:
        # Lepidoptera, Aves and Diptera are acquired under the exact broad taxon
        # declared in the policy. Their order/class fields need not be inferred from
        # focal function and are not required by the downstream effort classifier.
        result = records.copy()
        rule = "exact_acquisition_taxon_is_frozen_background"
        allowed_values = [str(channel.get("background_group", ""))]

    result = result.reset_index(drop=True)
    receipt = {
        "contract": "chapter1_nee_background_filter_v1",
        "channel_id": channel_id,
        "background_group": channel.get("background_group"),
        "filter_rule": rule,
        "allowed_values": allowed_values,
        "n_input_exact_island_records": before,
        "n_frozen_background_records": int(len(result)),
        "n_records_excluded_by_background_definition": int(before - len(result)),
        "uses_functional_target_catalog": False,
        "uses_focal_plant_traits": False,
    }
    return result, receipt


@app.command("filter")
def filter_command(
    records_csv: Path = typer.Option(..., exists=True),
    channel_id: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
    policy_path: Path = typer.Option(Path("config/chapter1_nee_channel_observation_policy.yml")),
) -> None:
    policy = load_policy(policy_path)
    records = pd.read_csv(records_csv, dtype=str).fillna("")
    filtered, receipt = filter_frozen_background(records, channel_id, policy)
    output_dir.mkdir(parents=True, exist_ok=True)
    filtered.to_csv(
        output_dir / f"{channel_id}_frozen_background_occurrences.csv.gz",
        index=False,
        compression="gzip",
    )
    (output_dir / f"{channel_id}_background_filter_receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(receipt))


if __name__ == "__main__":
    app()
