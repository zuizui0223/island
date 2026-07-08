"""Compress discovered literature leads into a bounded review packet.

Stage-1 discovery yields many leads per species. Trait relevance says a lead is
about flowers/pollination; it does not say the lead is about *this* species, so
the raw dump is dominated by off-target and off-species hits. This module applies
the **taxon relevance gate** (drop leads whose genus never appears in the
title/abstract) and then keeps only the strongest few leads per
(species, trait-group), spreading evenly across them, so the packet lands in a
human-reviewable size band (default 150-300). Only then is it worth a human
evaluating ``irrelevant_literature_rate``.

This decides nothing biological: it filters and ranks unreviewed leads. No trait
value, native/establishment status, Bombus applicability, or analysis inclusion
is assigned.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_COLUMNS = {
    "query_taxon",
    "query_trait_group",
    "taxon_relevance_score",
    "title_relevance_score",
    "abstract_relevance_score",
}

# Prohibited finalized-value columns must never ride along in a review packet.
PROHIBITED_COLUMNS = {
    "flower_primary_color", "floral_symmetry", "floral_form", "trait", "trait_value",
    "trait_state", "native_status", "establishment_means", "is_native",
    "bombus_applicability", "applicability", "analysis_included", "accepted",
    "review_status_accepted",
}


@app.callback()
def main() -> None:
    """Compress discovered leads into a bounded, taxon-gated review packet."""


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise typer.BadParameter("lead review-packet config must be a mapping")
    return config


def compress_review_packet(
    leads: pd.DataFrame, config: dict[str, Any]
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Gate on taxon relevance, then round-robin the strongest leads per cell.

    Fail-closed: if the taxon-gated set already fits within target_max, all of it is
    kept; otherwise leads are taken best-first, one per (species, trait-group) cell
    per pass, so no single species or group monopolises the packet.
    """
    missing = REQUIRED_COLUMNS.difference(leads.columns)
    if missing:
        raise typer.BadParameter(f"leads table missing columns: {', '.join(sorted(missing))}")
    prohibited = PROHIBITED_COLUMNS.intersection(leads.columns)
    if prohibited:
        raise typer.BadParameter(
            f"leads table carries prohibited finalized-value columns: {sorted(prohibited)}"
        )

    target_max = int(config.get("target_max", 300))
    target_min = int(config.get("target_min", 150))
    min_taxon = int(config.get("min_taxon_relevance_score", 1))

    table = leads.copy()
    for col in ("taxon_relevance_score", "title_relevance_score", "abstract_relevance_score"):
        table[col] = pd.to_numeric(table[col], errors="coerce").fillna(0).astype(int)
    if "priority_rank" in table.columns:
        table["priority_rank"] = pd.to_numeric(table["priority_rank"], errors="coerce").fillna(0).astype(int)
    else:
        table["priority_rank"] = 0

    n_input = int(len(table))
    gated = table.loc[table["taxon_relevance_score"] >= min_taxon].copy()
    n_after_gate = int(len(gated))

    # Best-first within each (species, trait-group) cell.
    gated = gated.sort_values(
        ["query_taxon", "query_trait_group", "taxon_relevance_score",
         "title_relevance_score", "abstract_relevance_score", "priority_rank"],
        ascending=[True, True, False, False, False, True],
    )
    buckets: dict[tuple[str, str], list[int]] = {}
    for key, group in gated.groupby(["query_taxon", "query_trait_group"], sort=True):
        buckets[key] = group.index.tolist()

    order = sorted(buckets, key=lambda k: (-len(buckets[k]), k))
    selected: list[int] = []
    if n_after_gate <= target_max:
        selected = gated.index.tolist()  # already tractable; keep the whole gated set
    else:
        exhausted = False
        while len(selected) < target_max and not exhausted:
            exhausted = True
            for key in order:
                if buckets[key]:
                    selected.append(buckets[key].pop(0))
                    exhausted = False
                    if len(selected) >= target_max:
                        break

    packet = gated.loc[selected].reset_index(drop=True)
    summary = {
        "note": (
            "Bounded, taxon-gated review packet of unreviewed leads. Leads whose genus "
            "never appears (taxon_relevance_score < min) are gated out as not credibly "
            "about the target species. No biological value is decided; this only filters "
            "and ranks so a human can evaluate irrelevant_literature_rate on a tractable set."
        ),
        "min_taxon_relevance_score": min_taxon,
        "target_min": target_min,
        "target_max": target_max,
        "n_input_leads": n_input,
        "n_after_taxon_gate": n_after_gate,
        "n_gated_out_zero_taxon": n_input - n_after_gate,
        "n_selected": int(len(packet)),
        "n_species_in_packet": int(packet["query_taxon"].nunique()) if len(packet) else 0,
        "within_target_band": bool(target_min <= len(packet) <= target_max) or n_after_gate < target_min,
        "n_by_trait_group": {
            str(k): int(v) for k, v in packet["query_trait_group"].value_counts().sort_index().items()
        } if len(packet) else {},
    }
    return packet, summary


@app.command("compress")
def compress(
    leads_csv: Path = typer.Option(..., exists=True, help="Discovered leads CSV (from source discovery)."),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/lead_review_packet.yml")),
) -> None:
    """Write a bounded, taxon-gated review packet plus a compression summary."""
    config = load_config(config_path)
    leads = pd.read_csv(leads_csv, dtype=str).fillna("")
    packet, summary = compress_review_packet(leads, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    packet.to_csv(output_dir / "lead_review_packet.csv", index=False)
    (output_dir / "lead_review_packet_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"Compressed {summary['n_input_leads']} leads -> {summary['n_selected']} in packet "
        f"({summary['n_gated_out_zero_taxon']} gated out on taxon relevance; "
        f"{summary['n_species_in_packet']} species)."
    )


if __name__ == "__main__":
    app()
