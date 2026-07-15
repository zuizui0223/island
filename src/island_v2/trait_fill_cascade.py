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
import hashlib
import json
import os
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.angiosperm_scope import classify_scope, load_config as load_scope_config
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
                "target_traits", "evidence_sources", "tier_labels", "angiosperm_scope_config"}
    if not isinstance(config, dict) or not required.issubset(config):
        raise typer.BadParameter(
            "cascade config must contain master_taxa_csv, species/genus/family columns, "
            "target_traits, evidence_sources, tier_labels, and angiosperm_scope_config"
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
        if "candidate_kind" in frame.columns:
            frame = frame.loc[frame["candidate_kind"].map(_text).eq("source_backed")].copy()
        if "evidence_scope" in frame.columns:
            frame = frame.loc[
                frame["evidence_scope"].map(_text).isin({"species_direct", "species_indirect"})
            ].copy()
        if frame.empty:
            continue
        species_col = next((c for c in ("accepted_species", "scientific_name", "species") if c in frame), None)
        value_col = next(
            (
                c
                for c in ("candidate_value", "standardized_value", "trait_value", "value")
                if c in frame
            ),
            None,
        )
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


def build_model(master: pd.DataFrame, evidence: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    """Build the shared cascade model once from all evidence, for per-shard reuse.

    Genus/family/global distributions are computed over the full eligible master
    and evidence, so applying the model to any species subset yields identical
    fills to a single whole-universe pass.
    """
    # Species-level modal direct value (species_direct tier source).
    species_direct = _modal(evidence, ["accepted_species", "trait_name"])
    species_value = {
        (r.accepted_species, r.trait_name): r for r in species_direct.itertuples(index=False)
    }

    # Attach each species-direct value to its genus/family for inference tiers.
    taxo = master[["accepted_species", "genus", "family"]].drop_duplicates("accepted_species")
    direct_taxo = species_direct.merge(taxo, on="accepted_species", how="left")

    genus_modal = _modal(
        direct_taxo[["genus", "trait_name", "value", "weight"]], ["genus", "trait_name"]
    )
    genus_lookup = {(r.genus, r.trait_name): r for r in genus_modal.itertuples(index=False)}
    family_modal = _modal(
        direct_taxo[["family", "trait_name", "value", "weight"]], ["family", "trait_name"]
    )
    family_lookup = {(r.family, r.trait_name): r for r in family_modal.itertuples(index=False)}
    global_modal = _modal(evidence, ["trait_name"])
    global_lookup = {r.trait_name: r for r in global_modal.itertuples(index=False)}
    return {
        "species_value": species_value,
        "genus_lookup": genus_lookup,
        "family_lookup": family_lookup,
        "global_lookup": global_lookup,
    }


def fill_species_frame(
    species_frame: pd.DataFrame, model: dict[str, Any], config: dict[str, Any]
) -> pd.DataFrame:
    """Fill one species subset by descending the cascade with a prebuilt model."""
    traits = list(config["target_traits"])
    labels = config["tier_labels"]
    min_genus = int(config.get("min_genus_support", 1))
    min_family = int(config.get("min_family_support", 3))
    species_value = model["species_value"]
    genus_lookup = model["genus_lookup"]
    family_lookup = model["family_lookup"]
    global_lookup = model["global_lookup"]

    rows: list[dict[str, Any]] = []
    for sp in species_frame.itertuples(index=False):
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


def build_fills(master: pd.DataFrame, evidence: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Fill every master species-trait by descending the taxonomic cascade."""
    model = build_model(master, evidence, config)
    return fill_species_frame(master, model, config)


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


# --- sharding, checkpoint, resume ------------------------------------------

CONTRACT_VERSION = "trait_fill_cascade_shards_v1"


def _sha_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def stable_shard(species: str, shard_count: int) -> int:
    """Return a stable shard for a species, independent of input row order."""
    if shard_count < 1:
        raise ValueError("shard_count must be positive")
    return int(_sha_text(_text(species))[:16], 16) % shard_count


def _now() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _atomic_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def _atomic_write_gzip(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    _write_gzip(frame, tmp)
    os.replace(tmp, path)


def model_fingerprint(evidence: pd.DataFrame, config: dict[str, Any]) -> str:
    """Fingerprint the shared model: any evidence or cascade-parameter change flips it."""
    ordered = (
        evidence[["accepted_species", "trait_name", "value", "weight"]]
        .astype(str)
        .sort_values(["accepted_species", "trait_name", "value", "weight"])
        if not evidence.empty
        else evidence
    )
    payload = "\n".join(
        "\t".join(row) for row in ordered.itertuples(index=False, name=None)
    ) if not evidence.empty else ""
    params = json.dumps(
        {
            "traits": list(config["target_traits"]),
            "min_genus_support": int(config.get("min_genus_support", 1)),
            "min_family_support": int(config.get("min_family_support", 3)),
            "tier_labels": config["tier_labels"],
            "contract": CONTRACT_VERSION,
        },
        sort_keys=True,
    )
    return _sha_text(params + "\n" + payload)


def shard_species_fingerprint(species: list[str]) -> str:
    """Fingerprint the exact species membership of a shard."""
    return _sha_text("\n".join(sorted(species)))


def shard_summary_from_fills(shard_fills: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    """Compact per-shard trait counts, disjoint across shards so global sums are exact."""
    traits = list(config["target_traits"])
    by_trait: dict[str, Any] = {}
    for trait in traits:
        sub = shard_fills.loc[shard_fills["trait_name"].eq(trait)]
        by_trait[trait] = {
            "n_filled": int(sub["accepted_species"].nunique()),
            "n_species_direct": int(
                sub.loc[sub["fill_tier"].eq("species_direct"), "accepted_species"].nunique()
            ),
            "by_tier": {str(k): int(v) for k, v in sub["fill_tier"].value_counts().items()},
        }
    return {
        "n_species": int(shard_fills["accepted_species"].nunique()) if not shard_fills.empty else 0,
        "by_trait": by_trait,
        "fills_by_tier": {str(k): int(v) for k, v in shard_fills["fill_tier"].value_counts().items()},
    }


def aggregate_shard_summaries(
    summaries: list[dict[str, Any]], config: dict[str, Any], n_eligible: int
) -> dict[str, Any]:
    """Sum disjoint per-shard summaries into the global coverage summary."""
    traits = list(config["target_traits"])
    by_trait: dict[str, Any] = {}
    fills_by_tier: dict[str, int] = {}
    for trait in traits:
        n_filled = sum(s["by_trait"].get(trait, {}).get("n_filled", 0) for s in summaries)
        direct = sum(s["by_trait"].get(trait, {}).get("n_species_direct", 0) for s in summaries)
        tiers: dict[str, int] = {}
        for s in summaries:
            for tier, count in s["by_trait"].get(trait, {}).get("by_tier", {}).items():
                tiers[tier] = tiers.get(tier, 0) + count
        by_trait[trait] = {
            "n_filled": n_filled,
            "fill_rate": n_filled / n_eligible if n_eligible else 0.0,
            "n_species_direct": direct,
            "species_direct_rate": direct / n_eligible if n_eligible else 0.0,
            "n_unknown_remaining": n_eligible - n_filled,
            "by_tier": tiers,
        }
    for s in summaries:
        for tier, count in s.get("fills_by_tier", {}).items():
            fills_by_tier[tier] = fills_by_tier.get(tier, 0) + count
    return {
        "version": "1.0",
        "n_master_species": n_eligible,
        "n_target_traits": len(traits),
        "by_trait": by_trait,
        "fills_by_tier": fills_by_tier,
        "interpretation": (
            "Fill-first cascade coverage. filled_value at genus/family/global tiers is "
            "low-resolution inference tagged by fill_tier and confidence, separable for "
            "sensitivity analysis, never accepted evidence or biological absence."
        ),
    }


def load_scoped_master_and_evidence(
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Load master, apply the angiosperm gate, and load eligible-only evidence."""
    master = pd.read_csv(config["master_taxa_csv"], dtype=str).fillna("").rename(
        columns={
            config["species_column"]: "accepted_species",
            config["genus_column"]: "genus",
            config["family_column"]: "family",
        }
    )
    master = master.loc[master["accepted_species"].map(_is_present)].drop_duplicates("accepted_species")

    scope_config = load_scope_config(Path(config["angiosperm_scope_config"]))
    scope = classify_scope(master, scope_config)
    eligible_species = set(scope.loc[scope["angiosperm_analysis_eligible"], "accepted_species"])
    eligible_master = (
        master.loc[master["accepted_species"].isin(eligible_species)]
        .sort_values("accepted_species")
        .reset_index(drop=True)
    )

    traits = set(config["target_traits"])
    evidence = load_direct_evidence(config, traits)
    evidence = evidence.loc[evidence["accepted_species"].isin(eligible_species)].reset_index(drop=True)
    return master, eligible_master, evidence, {"scope": scope}


@app.command("run")
def run(
    config_path: Path = typer.Option(Path("config/trait_fill_cascade.yml"), exists=True),
    output_dir: Path = typer.Option(..., help="Directory for cascade fill outputs."),
    shard_count: int = typer.Option(128, min=1, help="Number of deterministic species shards."),
    max_shards: int = typer.Option(
        0, min=0, help="Max stale shards to process this run (0 = all remaining)."
    ),
    no_resume: bool = typer.Option(False, help="Recompute every shard, ignoring checkpoints."),
) -> None:
    """Fill eligible species over resumable shards, then finalize when all are current.

    Each shard is checkpointed with the model and species fingerprints. A shard is
    recomputed only when its inputs changed or it is missing, so repeat runs do not
    recompute the whole universe. --max-shards bounds work per invocation for
    scheduled lanes; rerun to resume and finalize.
    """
    config = load_config(config_path)
    master, eligible_master, evidence, extra = load_scoped_master_and_evidence(config)
    scope = extra["scope"]
    n_master = int(len(master))
    n_eligible = int(len(eligible_master))

    model = build_model(eligible_master, evidence, config)
    fingerprint = model_fingerprint(evidence, config)

    assignments = eligible_master["accepted_species"].map(lambda s: stable_shard(s, shard_count))
    shards_dir = output_dir / "shards"
    shards_dir.mkdir(parents=True, exist_ok=True)

    processed, skipped, stale = 0, 0, 0
    for shard_index in range(shard_count):
        shard_species = eligible_master.loc[assignments.eq(shard_index)].reset_index(drop=True)
        species_fp = shard_species_fingerprint(list(shard_species["accepted_species"]))
        shard_dir = shards_dir / f"shard_{shard_index:05d}"
        checkpoint_path = shard_dir / "shard_checkpoint.json"

        current = None
        if not no_resume and checkpoint_path.exists():
            try:
                current = json.loads(checkpoint_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                current = None
        up_to_date = (
            current is not None
            and current.get("model_fingerprint") == fingerprint
            and current.get("species_fingerprint") == species_fp
            and (shard_dir / "fills.csv.gz").exists()
        )
        if up_to_date:
            skipped += 1
            continue

        stale += 1
        if max_shards and processed >= max_shards:
            continue

        shard_fills = fill_species_frame(shard_species, model, config)
        summary = shard_summary_from_fills(shard_fills, config)
        _atomic_write_gzip(shard_fills, shard_dir / "fills.csv.gz")
        _atomic_json(
            {
                "contract_version": CONTRACT_VERSION,
                "shard_index": shard_index,
                "shard_count": shard_count,
                "model_fingerprint": fingerprint,
                "species_fingerprint": species_fp,
                "summary": summary,
                "updated_at": _now(),
            },
            checkpoint_path,
        )
        processed += 1

    remaining = stale - processed
    complete = remaining == 0
    _finalize(output_dir, shards_dir, shard_count, config, master, eligible_master, scope, fingerprint,
              n_master, n_eligible, complete)

    typer.echo(
        f"Shards {shard_count}: processed {processed}, skipped up-to-date {skipped}, "
        f"remaining stale {remaining}. Eligible {n_eligible}/{n_master}. "
        f"{'Finalized (all shards current).' if complete else 'Partial; rerun to resume.'}"
    )


def _finalize(
    output_dir: Path,
    shards_dir: Path,
    shard_count: int,
    config: dict[str, Any],
    master: pd.DataFrame,
    eligible_master: pd.DataFrame,
    scope: pd.DataFrame,
    fingerprint: str,
    n_master: int,
    n_eligible: int,
    complete: bool,
) -> None:
    """Aggregate current shard summaries into the coverage summary and benchmark."""
    summaries: list[dict[str, Any]] = []
    shards_done = 0
    for shard_index in range(shard_count):
        checkpoint_path = shards_dir / f"shard_{shard_index:05d}" / "shard_checkpoint.json"
        if not checkpoint_path.exists():
            continue
        payload = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        if payload.get("model_fingerprint") != fingerprint:
            continue
        summaries.append(payload["summary"])
        shards_done += 1

    coverage = aggregate_shard_summaries(summaries, config, n_eligible)
    coverage["n_master_species_all"] = n_master
    coverage["n_angiosperm_eligible"] = n_eligible
    coverage["n_out_of_scope"] = n_master - n_eligible
    coverage["out_of_scope_by_group"] = {
        str(k): int(v)
        for k, v in scope.loc[~scope["angiosperm_analysis_eligible"], "taxonomic_group"].value_counts().items()
    }
    coverage["denominator"] = "angiosperm_eligible"
    coverage["shards_complete"] = complete
    coverage["n_shards_current"] = shards_done
    coverage["shard_count"] = shard_count

    output_dir.mkdir(parents=True, exist_ok=True)
    scope.to_csv(output_dir / "angiosperm_scope_by_species.csv", index=False)
    (output_dir / "fill_coverage_summary.json").write_text(
        json.dumps(coverage, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    _atomic_json(
        {
            "contract_version": CONTRACT_VERSION,
            "shard_count": shard_count,
            "model_fingerprint": fingerprint,
            "n_master_species_all": n_master,
            "n_angiosperm_eligible": n_eligible,
            "n_shards_current": shards_done,
            "shards_complete": complete,
            "updated_at": _now(),
        },
        output_dir / "cascade_manifest.json",
    )

    # Benchmark reads only the shards holding the deterministic sample species.
    size = int(config.get("benchmark_sample_size", 100))
    sample = sorted(eligible_master["accepted_species"].map(_text).loc[lambda s: s.ne("")].unique())[:size]
    sample_shards = {stable_shard(s, shard_count) for s in sample}
    frames = []
    for shard_index in sorted(sample_shards):
        path = shards_dir / f"shard_{shard_index:05d}" / "fills.csv.gz"
        if path.exists():
            frames.append(pd.read_csv(path, dtype=str).fillna(""))
    if frames:
        benchmark = build_benchmark(pd.concat(frames, ignore_index=True), eligible_master, config)
    else:
        benchmark = pd.DataFrame(columns=["accepted_species"])
    benchmark.to_csv(output_dir / "benchmark_sample.csv", index=False)


@app.command("status")
def status(
    output_dir: Path = typer.Option(..., help="Cascade output directory to inspect."),
    config_path: Path = typer.Option(Path("config/trait_fill_cascade.yml"), exists=True),
    shard_count: int = typer.Option(128, min=1, help="Shard count used for the run."),
) -> None:
    """Report how many shards are current for the present inputs, without filling."""
    config = load_config(config_path)
    _, eligible_master, evidence, _ = load_scoped_master_and_evidence(config)
    fingerprint = model_fingerprint(evidence, config)
    assignments = eligible_master["accepted_species"].map(lambda s: stable_shard(s, shard_count))
    shards_dir = output_dir / "shards"

    current, stale, missing = 0, 0, 0
    for shard_index in range(shard_count):
        shard_species = list(eligible_master.loc[assignments.eq(shard_index), "accepted_species"])
        species_fp = shard_species_fingerprint(shard_species)
        checkpoint_path = shards_dir / f"shard_{shard_index:05d}" / "shard_checkpoint.json"
        fills_path = shards_dir / f"shard_{shard_index:05d}" / "fills.csv.gz"
        if not checkpoint_path.exists() or not fills_path.exists():
            missing += 1
            continue
        payload = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        if payload.get("model_fingerprint") == fingerprint and payload.get("species_fingerprint") == species_fp:
            current += 1
        else:
            stale += 1
    typer.echo(
        f"Shards {shard_count}: current {current}, stale {stale}, missing {missing}. "
        f"{'All current.' if current == shard_count else 'Rerun run to resume.'}"
    )


if __name__ == "__main__":
    app()
