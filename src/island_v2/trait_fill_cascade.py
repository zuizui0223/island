"""Fill-first taxonomic cascade to drive 9-column trait unknown toward zero.

The operational priority is trait yield, not measurement. Direct source-backed
evidence covers a small fraction of the master species, so this cascade fills
every remaining species-trait by descending a source-blind fallback ladder:

    species_direct -> synonym_direct -> genus_inference -> family_inference
    -> global_fallback

Every fill records ``fill_tier``, ``evidence_scope`` and ``confidence`` so a
low-resolution fill is always separable for sensitivity analysis and never
masquerades as accepted direct evidence. Fills are written to staging, never to
curated. A cascade fill is candidate coverage, never a biological absence, and
one trait is never derived from another (self-compatibility never fills
autonomous selfing). Inference tiers also retain the supporting value
distribution so a modal fill can be checked against a distribution-aware draw.
"""

from __future__ import annotations

import glob as globlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.search_enabled_llm_campaign import _parse_csv_row as parse_llm_result

app = typer.Typer(add_completion=False, no_args_is_help=True)

EMPTY_VALUES = {"", "unknown", "na", "n/a", "none", "unresolved"}

LLM_TRAIT_MAP = {
    "flower_color": "flower_primary_color",
    "flower_shape": "floral_form",
    "pollination_guild": "pollination_functional_guild",
    "mating_system": "mating_system",
    "self_incompatibility": "self_incompatibility",
}

OUTPUT_COLUMNS = [
    "accepted_species",
    "genus",
    "family",
    "trait_name",
    "filled_value",
    "fill_tier",
    "evidence_scope",
    "confidence",
    "support_n",
    "value_distribution",
]


@app.callback()
def main() -> None:
    """Fill every master species-trait by taxonomic cascade, tagged by resolution."""


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _is_present(value: object) -> bool:
    return _text(value).lower() not in EMPTY_VALUES


def _write_gzip(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})


def load_config(path: Path) -> dict[str, Any]:
    """Load and minimally validate the versioned cascade configuration."""
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    required = {"master_taxa_csv", "species_column", "genus_column", "family_column",
                "target_traits", "evidence_sources", "tier_labels"}
    if not isinstance(config, dict) or not required.issubset(config):
        raise typer.BadParameter(
            "cascade config must contain master_taxa_csv, species/genus/family columns, "
            "target_traits, evidence_sources, and tier_labels"
        )
    return config


# --- evidence loading: return (accepted_species, trait_name, value, weight) ---

def _evidence_candidate_index(source: dict[str, Any]) -> pd.DataFrame:
    path = Path(source["path"])
    if not path.exists():
        return pd.DataFrame(columns=["accepted_species", "trait_name", "value", "weight"])
    frame = pd.read_csv(path, dtype=str).fillna("")
    weight = pd.to_numeric(frame.get("n_wave_observations", 1), errors="coerce").fillna(1.0)
    return pd.DataFrame(
        {
            "accepted_species": frame["accepted_species"].map(_text),
            "trait_name": frame["trait_name"].map(_text),
            "value": frame["candidate_value"].map(_text),
            "weight": weight.to_numpy(),
        }
    )


def _evidence_candidate_long(source: dict[str, Any]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in sorted(globlib.glob(source["glob"], recursive=True)):
        frame = pd.read_csv(path, dtype=str).fillna("")
        if frame.empty:
            continue
        species_col = next((c for c in ("accepted_species", "scientific_name", "species") if c in frame), None)
        value_col = next((c for c in ("candidate_value", "trait_value", "value") if c in frame), None)
        if species_col is None or value_col is None or "trait_name" not in frame:
            continue
        frames.append(
            pd.DataFrame(
                {
                    "accepted_species": frame[species_col].map(_text),
                    "trait_name": frame["trait_name"].map(_text),
                    "value": frame[value_col].map(_text),
                    "weight": 1.0,
                }
            )
        )
    if not frames:
        return pd.DataFrame(columns=["accepted_species", "trait_name", "value", "weight"])
    return pd.concat(frames, ignore_index=True)


def _evidence_curated(source: dict[str, Any]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in sorted(globlib.glob(source["glob"], recursive=True)):
        frame = pd.read_csv(path, dtype=str).fillna("")
        if frame.empty:
            continue
        if "review_status" in frame.columns:
            frame = frame.loc[frame["review_status"].map(_text).str.lower().eq("accepted")]
        value_col = next((c for c in ("trait_value", "candidate_value", "value") if c in frame), None)
        if value_col is None or "trait_name" not in frame or "accepted_species" not in frame:
            continue
        frames.append(
            pd.DataFrame(
                {
                    "accepted_species": frame["accepted_species"].map(_text),
                    "trait_name": frame["trait_name"].map(_text),
                    "value": frame[value_col].map(_text),
                    # curated evidence outweighs machine candidates in the modal vote.
                    "weight": 1000.0,
                }
            )
        )
    if not frames:
        return pd.DataFrame(columns=["accepted_species", "trait_name", "value", "weight"])
    return pd.concat(frames, ignore_index=True)


def _evidence_search_enabled_llm(source: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for path in sorted(globlib.glob(source["glob"], recursive=True)):
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                parsed = parse_llm_result(_text(json.loads(line).get("result")))
            except (ValueError, KeyError, json.JSONDecodeError):
                continue
            species = _text(parsed.get("species", ""))
            for llm_column, trait_name in LLM_TRAIT_MAP.items():
                value = _text(parsed.get(llm_column, ""))
                if species and _is_present(value):
                    rows.append(
                        {"accepted_species": species, "trait_name": trait_name, "value": value, "weight": 1.0}
                    )
    return pd.DataFrame(rows, columns=["accepted_species", "trait_name", "value", "weight"])


EVIDENCE_ADAPTERS = {
    "candidate_index": _evidence_candidate_index,
    "candidate_long": _evidence_candidate_long,
    "curated_evidence": _evidence_curated,
    "search_enabled_llm": _evidence_search_enabled_llm,
}


def load_direct_evidence(config: dict[str, Any], traits: set[str]) -> pd.DataFrame:
    """Union all evidence channels into (species, trait, value, weight) rows."""
    frames: list[pd.DataFrame] = []
    for key, source in config["evidence_sources"].items():
        adapter = source.get("adapter")
        if adapter not in EVIDENCE_ADAPTERS:
            raise typer.BadParameter(f"evidence source {key!r} has unknown adapter {adapter!r}")
        frame = EVIDENCE_ADAPTERS[adapter](source)
        if not frame.empty:
            frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=["accepted_species", "trait_name", "value", "weight"])
    evidence = pd.concat(frames, ignore_index=True)
    evidence = evidence.loc[
        evidence["trait_name"].isin(traits)
        & evidence["accepted_species"].map(_is_present)
        & evidence["value"].map(_is_present)
    ]
    return evidence.reset_index(drop=True)


def _modal(frame: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    """Return the max-weight value per group, ties broken alphabetically."""
    if frame.empty:
        return pd.DataFrame(columns=[*group_cols, "value", "weight", "value_distribution", "support_n"])
    totals = frame.groupby([*group_cols, "value"], as_index=False)["weight"].sum()
    totals = totals.sort_values([*group_cols, "weight", "value"], ascending=[*([True] * len(group_cols)), False, True])
    winners = totals.drop_duplicates(group_cols, keep="first").rename(columns={"value": "value"})
    dist = (
        totals.groupby(group_cols)
        .apply(lambda g: json.dumps({str(v): round(float(w), 3) for v, w in zip(g["value"], g["weight"], strict=True)}), include_groups=False)
        .rename("value_distribution")
        .reset_index()
    )
    support = frame.groupby(group_cols)["value"].size().rename("support_n").reset_index()
    return winners.merge(dist, on=group_cols).merge(support, on=group_cols)


def build_fills(master: pd.DataFrame, evidence: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Fill every master species-trait by descending the taxonomic cascade."""
    traits = list(config["target_traits"])
    labels = config["tier_labels"]
    min_genus = int(config.get("min_genus_support", 1))
    min_family = int(config.get("min_family_support", 3))

    # Species-level modal direct value (species_direct tier source).
    species_direct = _modal(evidence, ["accepted_species", "trait_name"])
    species_value = {
        (r.accepted_species, r.trait_name): r for r in species_direct.itertuples(index=False)
    }

    # Attach each species-direct value to its genus/family for inference tiers.
    taxo = master[["accepted_species", "genus", "family"]].drop_duplicates("accepted_species")
    direct_taxo = species_direct.merge(taxo, on="accepted_species", how="left")

    genus_modal = _modal(
        direct_taxo.rename(columns={"value": "value"})[["genus", "trait_name", "value", "weight"]],
        ["genus", "trait_name"],
    )
    genus_lookup = {(r.genus, r.trait_name): r for r in genus_modal.itertuples(index=False)}
    family_modal = _modal(
        direct_taxo[["family", "trait_name", "value", "weight"]], ["family", "trait_name"]
    )
    family_lookup = {(r.family, r.trait_name): r for r in family_modal.itertuples(index=False)}
    global_modal = _modal(evidence, ["trait_name"])
    global_lookup = {r.trait_name: r for r in global_modal.itertuples(index=False)}

    rows: list[dict[str, Any]] = []
    for sp in master.itertuples(index=False):
        species = _text(sp.accepted_species)
        genus = _text(sp.genus)
        family = _text(sp.family)
        for trait in traits:
            direct = species_value.get((species, trait))
            if direct is not None:
                tier, value, support, dist = "species_direct", direct.value, int(direct.support_n), direct.value_distribution
            else:
                gen = genus_lookup.get((genus, trait))
                fam = family_lookup.get((family, trait))
                glob = global_lookup.get(trait)
                if gen is not None and int(gen.support_n) >= min_genus:
                    tier, value, support, dist = "genus_inference", gen.value, int(gen.support_n), gen.value_distribution
                elif fam is not None and int(fam.support_n) >= min_family:
                    tier, value, support, dist = "family_inference", fam.value, int(fam.support_n), fam.value_distribution
                elif glob is not None:
                    tier, value, support, dist = "global_fallback", glob.value, int(glob.support_n), glob.value_distribution
                else:
                    continue  # trait has zero evidence anywhere; leave genuinely unknown
            label = labels[tier]
            rows.append(
                {
                    "accepted_species": species,
                    "genus": genus,
                    "family": family,
                    "trait_name": trait,
                    "filled_value": value,
                    "fill_tier": tier,
                    "evidence_scope": label["evidence_scope"],
                    "confidence": label["confidence"],
                    "support_n": support,
                    "value_distribution": dist if tier != "species_direct" else "",
                }
            )
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def build_coverage_summary(fills: pd.DataFrame, config: dict[str, Any], n_master: int) -> dict[str, Any]:
    """Per-trait fill coverage and tier composition against the master."""
    traits = list(config["target_traits"])
    by_trait: dict[str, Any] = {}
    for trait in traits:
        sub = fills.loc[fills["trait_name"].eq(trait)]
        n_filled = int(sub["accepted_species"].nunique())
        tiers = {str(k): int(v) for k, v in sub["fill_tier"].value_counts().items()}
        direct = int(sub.loc[sub["fill_tier"].eq("species_direct"), "accepted_species"].nunique())
        by_trait[trait] = {
            "n_filled": n_filled,
            "fill_rate": n_filled / n_master if n_master else 0.0,
            "n_species_direct": direct,
            "species_direct_rate": direct / n_master if n_master else 0.0,
            "n_unknown_remaining": n_master - n_filled,
            "by_tier": tiers,
        }
    overall_tiers = {str(k): int(v) for k, v in fills["fill_tier"].value_counts().items()}
    return {
        "version": "1.0",
        "n_master_species": n_master,
        "n_target_traits": len(traits),
        "by_trait": by_trait,
        "fills_by_tier": overall_tiers,
        "interpretation": (
            "Fill-first cascade coverage. filled_value at genus/family/global tiers is "
            "low-resolution inference tagged by fill_tier and confidence, separable for "
            "sensitivity analysis, never accepted evidence or biological absence."
        ),
    }


def build_benchmark(fills: pd.DataFrame, master: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Wide per-field fill for a deterministic species sample (the 100-sp benchmark)."""
    size = int(config.get("benchmark_sample_size", 100))
    species = sorted(master["accepted_species"].map(_text).loc[lambda s: s.ne("")].unique())[:size]
    sub = fills.loc[fills["accepted_species"].isin(species)].copy()
    sub["cell"] = sub["filled_value"] + " [" + sub["fill_tier"] + "]"
    wide = sub.pivot_table(
        index="accepted_species", columns="trait_name", values="cell", aggfunc="first"
    ).reindex(species)
    return wide.reset_index()


@app.command("run")
def run(
    config_path: Path = typer.Option(Path("config/trait_fill_cascade.yml"), exists=True),
    output_dir: Path = typer.Option(..., help="Directory for cascade fill outputs."),
) -> None:
    """Write the full fill table, per-trait coverage summary, and benchmark sample."""
    config = load_config(config_path)
    master = pd.read_csv(config["master_taxa_csv"], dtype=str).fillna("").rename(
        columns={
            config["species_column"]: "accepted_species",
            config["genus_column"]: "genus",
            config["family_column"]: "family",
        }
    )
    master = master.loc[master["accepted_species"].map(_is_present)].drop_duplicates("accepted_species")
    n_master = int(len(master))

    traits = set(config["target_traits"])
    evidence = load_direct_evidence(config, traits)
    fills = build_fills(master, evidence, config)
    summary = build_coverage_summary(fills, config, n_master)
    benchmark = build_benchmark(fills, master, config)

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_gzip(fills, output_dir / "trait_fills.csv.gz")
    benchmark.to_csv(output_dir / "benchmark_sample.csv", index=False)
    (output_dir / "fill_coverage_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    colour = summary["by_trait"].get("flower_primary_color", {})
    typer.echo(
        f"Filled {len(fills)} species-trait cells over {n_master} species x "
        f"{summary['n_target_traits']} traits. Colour direct "
        f"{colour.get('n_species_direct', 0)} -> filled {colour.get('n_filled', 0)} "
        f"({colour.get('fill_rate', 0):.1%}); tiers {summary['fills_by_tier']}."
    )


if __name__ == "__main__":
    app()
