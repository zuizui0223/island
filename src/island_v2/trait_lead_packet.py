"""Compress discovered literature leads into tiered review outputs.

Stage-1A discovery yields many leads per species. Trait relevance says a lead is
about flowers/pollination; it does not say the lead is about *this* species. So
leads are split by taxon-relevance **tier** and routed, not merged:

- **Tier S** (species named with confidence) -> the bounded review packet, compressed
  round-robin per (species, trait-group) so no species/group monopolises it.
- **Tier A** (genus-only / flora-monograph rescue) -> a separate quarantined CSV, kept
  for optional follow-up but not the main review target.
- **Tier B** (species never named; background literature) -> retained in its own CSV,
  never a review target.

Only the Tier-S packet is what a human evaluates for ``irrelevant_literature_rate``.
This decides nothing biological: it filters and ranks unreviewed leads. No trait
value, native/establishment status, Bombus applicability, or analysis inclusion is
assigned.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.trait_source_discovery import (
    TIER_BACKGROUND,
    TIER_GENUS_FLORA,
    TIER_SPECIES,
    taxon_relevance_tier,
)

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

_SORT_COLUMNS = [
    "query_taxon", "query_trait_group", "taxon_relevance_score",
    "title_relevance_score", "abstract_relevance_score", "priority_rank",
]
_SORT_ASCENDING = [True, True, False, False, False, True]


@app.callback()
def main() -> None:
    """Split discovered leads into a Tier-S packet, a Tier-A rescue, and Tier-B background."""


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise typer.BadParameter("lead review-packet config must be a mapping")
    return config


def _prepare(leads: pd.DataFrame) -> pd.DataFrame:
    missing = REQUIRED_COLUMNS.difference(leads.columns)
    if missing:
        raise typer.BadParameter(f"leads table missing columns: {', '.join(sorted(missing))}")
    prohibited = PROHIBITED_COLUMNS.intersection(leads.columns)
    if prohibited:
        raise typer.BadParameter(
            f"leads table carries prohibited finalized-value columns: {sorted(prohibited)}"
        )
    table = leads.copy()
    for col in ("taxon_relevance_score", "title_relevance_score", "abstract_relevance_score"):
        table[col] = pd.to_numeric(table[col], errors="coerce").fillna(0).astype(int)
    if "priority_rank" in table.columns:
        table["priority_rank"] = pd.to_numeric(table["priority_rank"], errors="coerce").fillna(0).astype(int)
    else:
        table["priority_rank"] = 0
    # Derive the tier from the score when the column is absent (older lead tables).
    if "taxon_relevance_tier" not in table.columns:
        table["taxon_relevance_tier"] = table["taxon_relevance_score"].map(taxon_relevance_tier)
    return table.sort_values(_SORT_COLUMNS, ascending=_SORT_ASCENDING).reset_index(drop=True)


def _round_robin_cap(tier_s: pd.DataFrame, target_max: int) -> pd.DataFrame:
    """Take Tier-S leads best-first, one per (species, trait-group) per pass, up to cap."""
    if len(tier_s) <= target_max:
        return tier_s.reset_index(drop=True)
    buckets: dict[tuple[str, str], list[int]] = {}
    for key, group in tier_s.groupby(["query_taxon", "query_trait_group"], sort=True):
        buckets[key] = group.index.tolist()
    order = sorted(buckets, key=lambda k: (-len(buckets[k]), k))
    selected: list[int] = []
    exhausted = False
    while len(selected) < target_max and not exhausted:
        exhausted = True
        for key in order:
            if buckets[key]:
                selected.append(buckets[key].pop(0))
                exhausted = False
                if len(selected) >= target_max:
                    break
    return tier_s.loc[selected].reset_index(drop=True)


def compress_review_packet(
    leads: pd.DataFrame, config: dict[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Split leads by taxon-relevance tier into (packet, rescue, background, summary).

    Tier S is compressed round-robin into the bounded packet; Tier A (genus/flora) is
    isolated for optional follow-up; Tier B (background) is retained but never reviewed.
    """
    target_max = int(config.get("target_max", 300))
    target_min = int(config.get("target_min", 150))

    table = _prepare(leads)
    tier = table["taxon_relevance_tier"].astype(str)
    tier_s = table.loc[tier == TIER_SPECIES]
    rescue = table.loc[tier == TIER_GENUS_FLORA].reset_index(drop=True)
    background = table.loc[tier == TIER_BACKGROUND].reset_index(drop=True)

    packet = _round_robin_cap(tier_s, target_max)

    summary = {
        "note": (
            "Tier-routed, unreviewed leads. Tier S (species named) is the review packet; "
            "Tier A (genus-only / flora rescue) is quarantined for optional follow-up; "
            "Tier B (species never named) is background literature, retained but not "
            "reviewed. No biological value is decided; only the Tier-S packet feeds a "
            "human irrelevant_literature_rate evaluation."
        ),
        "target_min": target_min,
        "target_max": target_max,
        "n_input_leads": int(len(table)),
        "n_tier_S_species_confident": int(len(tier_s)),
        "n_tier_A_genus_flora": int(len(rescue)),
        "n_tier_B_background": int(len(background)),
        "n_selected_packet": int(len(packet)),
        "n_species_in_packet": int(packet["query_taxon"].nunique()) if len(packet) else 0,
        "packet_within_target_band": bool(target_min <= len(packet) <= target_max)
        or int(len(tier_s)) < target_min,
        "packet_by_trait_group": {
            str(k): int(v) for k, v in packet["query_trait_group"].value_counts().sort_index().items()
        } if len(packet) else {},
    }
    return packet, rescue, background, summary


@app.command("compress")
def compress(
    leads_csv: Path = typer.Option(..., exists=True, help="Discovered leads CSV (from source discovery)."),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/lead_review_packet.yml")),
) -> None:
    """Write the Tier-S review packet, the Tier-A rescue, the Tier-B background, and a summary."""
    config = load_config(config_path)
    leads = pd.read_csv(leads_csv, dtype=str).fillna("")
    packet, rescue, background, summary = compress_review_packet(leads, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    packet.to_csv(output_dir / "lead_review_packet.csv", index=False)
    rescue.to_csv(output_dir / "lead_rescue_genus_flora.csv", index=False)
    background.to_csv(output_dir / "lead_background_literature.csv", index=False)
    (output_dir / "lead_review_packet_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"{summary['n_input_leads']} leads -> Tier S {summary['n_tier_S_species_confident']} "
        f"(packet {summary['n_selected_packet']}, {summary['n_species_in_packet']} species), "
        f"Tier A rescue {summary['n_tier_A_genus_flora']}, Tier B background {summary['n_tier_B_background']}."
    )


if __name__ == "__main__":
    app()
