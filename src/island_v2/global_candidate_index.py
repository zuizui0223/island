"""Build cumulative, machine-only indexes from global trait-campaign waves.

The outputs are discovery and review aids. They never accept a trait value, infer a
biological absence, or enter confirmatory models. Every indexed row retains the
source excerpt and campaign provenance that produced it.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_CANDIDATE_COLUMNS = {
    "accepted_species",
    "trait_name",
    "candidate_value",
    "source_url",
    "source_excerpt",
    "source_type",
    "evidence_scope",
    "campaign_task",
    "campaign_phase",
    "target_for_task",
}

INDEX_KEY = [
    "accepted_species",
    "trait_name",
    "candidate_value",
    "source_url",
    "source_excerpt",
    "source_type",
    "evidence_scope",
    "matched_term",
    "campaign_task",
    "campaign_phase",
]

ALTERNATIVE_GUILD_VALUES = {
    "other_bees",
    "flies",
    "butterflies",
    "moths",
    "birds",
    "bats",
    "beetles_wasps_ants",
    "mixed_or_generalist",
}
NON_BEE_GUILD_VALUES = {
    "flies",
    "butterflies",
    "moths",
    "birds",
    "bats",
    "beetles_wasps_ants",
}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _bool_series(series: pd.Series) -> pd.Series:
    return series.astype("string").fillna("").str.lower().isin({"true", "1", "yes", "y"})


def _write_gzip(frame: pd.DataFrame, path: Path) -> None:
    """Write deterministic gzip bytes so no-op rebuilds do not create commits."""
    frame.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})


def load_wave_candidates(campaign_dir: Path) -> tuple[pd.DataFrame, list[str]]:
    """Read all completed wave candidate files and attach their wave IDs."""
    paths = sorted((campaign_dir / "waves").glob("wave_*/machine_candidates.csv.gz"))
    frames: list[pd.DataFrame] = []
    wave_ids: list[str] = []
    for path in paths:
        frame = pd.read_csv(path, dtype=str).fillna("")
        missing = REQUIRED_CANDIDATE_COLUMNS.difference(frame.columns)
        if missing:
            raise typer.BadParameter(f"candidate wave {path} missing columns: {sorted(missing)}")
        if "matched_term" not in frame.columns:
            frame["matched_term"] = ""
        frame["wave_id"] = path.parent.name
        frames.append(frame)
        wave_ids.append(path.parent.name)

    columns = sorted(REQUIRED_CANDIDATE_COLUMNS | {"matched_term", "wave_id"})
    if not frames:
        return pd.DataFrame(columns=columns), wave_ids
    return pd.concat(frames, ignore_index=True, sort=False).fillna(""), wave_ids


def build_candidate_index(raw: pd.DataFrame) -> pd.DataFrame:
    """Deduplicate identical source-backed candidates while retaining wave history."""
    if raw.empty:
        return pd.DataFrame(
            columns=[
                *INDEX_KEY,
                "target_for_task",
                "first_wave_id",
                "last_wave_id",
                "n_wave_observations",
            ]
        )

    frame = raw.copy()
    for column in INDEX_KEY:
        if column not in frame.columns:
            frame[column] = ""
        frame[column] = frame[column].map(_text)
    frame["target_for_task"] = _bool_series(frame["target_for_task"])
    frame["wave_id"] = frame["wave_id"].map(_text)
    frame = frame.sort_values(["wave_id", *INDEX_KEY]).reset_index(drop=True)

    rows: list[dict[str, Any]] = []
    for key, group in frame.groupby(INDEX_KEY, dropna=False, sort=True):
        record = dict(zip(INDEX_KEY, key, strict=True))
        waves = sorted({value for value in group["wave_id"].astype(str) if value})
        record.update(
            {
                "target_for_task": bool(group["target_for_task"].astype(bool).any()),
                "first_wave_id": waves[0] if waves else "",
                "last_wave_id": waves[-1] if waves else "",
                "n_wave_observations": int(len(group)),
            }
        )
        rows.append(record)
    return pd.DataFrame(rows).sort_values(
        ["accepted_species", "campaign_phase", "trait_name", "candidate_value", "source_url"]
    ).reset_index(drop=True)


def load_campaign_ledger(campaign_dir: Path) -> pd.DataFrame:
    path = campaign_dir / "campaign_ledger.csv.gz"
    if not path.exists():
        raise typer.BadParameter(f"campaign ledger not found: {path}")
    ledger = pd.read_csv(path, dtype=str).fillna("")
    required = {
        "accepted_species",
        "family",
        "n_islands",
        "n_records",
        "machine_biotic_candidate",
        "reproductive_wikimedia_status",
        "reproductive_openalex_status",
        "floral_access_wikimedia_status",
        "alternative_guild_openalex_status",
    }
    missing = required.difference(ledger.columns)
    if missing:
        raise typer.BadParameter(f"campaign ledger missing columns: {sorted(missing)}")
    ledger["machine_biotic_candidate"] = _bool_series(ledger["machine_biotic_candidate"])
    return ledger


def _phase_count(index: pd.DataFrame, phase: str) -> pd.Series:
    if index.empty:
        return pd.Series(dtype=int)
    return index.loc[index["campaign_phase"].eq(phase)].groupby("accepted_species").size()


def build_counterfactual_screening(index: pd.DataFrame, ledger: pd.DataFrame) -> pd.DataFrame:
    """Build a machine-only queue for alternative-pollinator evidence review."""
    base = ledger.loc[ledger["machine_biotic_candidate"].astype(bool)].copy()
    base = base[
        [
            "accepted_species",
            "family",
            "n_islands",
            "n_records",
            "reproductive_wikimedia_status",
            "reproductive_openalex_status",
            "floral_access_wikimedia_status",
            "alternative_guild_openalex_status",
        ]
    ].drop_duplicates("accepted_species")

    output_columns = [
        *base.columns,
        "n_reproductive_candidates",
        "n_floral_access_candidates",
        "n_alternative_guild_candidates",
        "guild_candidate_values",
        "machine_has_bumblebee_candidate",
        "machine_has_other_bee_candidate",
        "machine_has_nonbee_candidate",
        "machine_has_mixed_generalist_candidate",
        "machine_has_alternative_guild_candidate",
        "machine_screening_label",
        "evidence_use_tier",
    ]
    if base.empty:
        return pd.DataFrame(columns=output_columns)

    reproductive = _phase_count(index, "reproductive_and_pollen_vector")
    floral = _phase_count(index, "floral_access")
    alternative = _phase_count(index, "alternative_pollinator_function")

    guild = index.loc[index["trait_name"].eq("pollination_functional_guild")].copy()
    guild_sets = (
        guild.groupby("accepted_species")["candidate_value"]
        .agg(lambda values: sorted({_text(value) for value in values if _text(value)}))
        .to_dict()
        if not guild.empty
        else {}
    )

    rows: list[dict[str, Any]] = []
    for record in base.to_dict("records"):
        species = _text(record["accepted_species"])
        values = set(guild_sets.get(species, []))
        has_bumblebee = "bumblebees" in values
        has_other_bee = "other_bees" in values
        has_nonbee = bool(values.intersection(NON_BEE_GUILD_VALUES))
        has_mixed = "mixed_or_generalist" in values
        has_alternative = bool(values.intersection(ALTERNATIVE_GUILD_VALUES))
        if has_alternative:
            label = "alternative_or_generalist_candidate"
        elif has_bumblebee:
            label = "bumblebee_candidate_without_alternative_yet"
        else:
            label = "biotic_vector_candidate_without_guild_yet"
        rows.append(
            {
                **record,
                "n_reproductive_candidates": int(reproductive.get(species, 0)),
                "n_floral_access_candidates": int(floral.get(species, 0)),
                "n_alternative_guild_candidates": int(alternative.get(species, 0)),
                "guild_candidate_values": "|".join(sorted(values)),
                "machine_has_bumblebee_candidate": has_bumblebee,
                "machine_has_other_bee_candidate": has_other_bee,
                "machine_has_nonbee_candidate": has_nonbee,
                "machine_has_mixed_generalist_candidate": has_mixed,
                "machine_has_alternative_guild_candidate": has_alternative,
                "machine_screening_label": label,
                "evidence_use_tier": "machine_unreviewed_screening_only",
            }
        )
    return pd.DataFrame(rows, columns=output_columns).sort_values(
        ["machine_screening_label", "accepted_species"]
    ).reset_index(drop=True)


def build_status(
    index: pd.DataFrame,
    screening: pd.DataFrame,
    wave_ids: list[str],
) -> dict[str, Any]:
    def counts(column: str) -> dict[str, int]:
        if index.empty or column not in index.columns:
            return {}
        values = index[column].value_counts().sort_index()
        return {str(key): int(value) for key, value in values.items()}

    labels = (
        screening["machine_screening_label"].value_counts().sort_index()
        if not screening.empty
        else pd.Series(dtype=int)
    )
    return {
        "version": "1.0",
        "n_waves_indexed": int(len(set(wave_ids))),
        "latest_wave_id": sorted(set(wave_ids))[-1] if wave_ids else "",
        "n_unique_candidate_rows": int(len(index)),
        "n_species_with_candidates": int(index["accepted_species"].nunique()) if not index.empty else 0,
        "n_machine_biotic_screening_species": int(len(screening)),
        "n_by_campaign_task": counts("campaign_task"),
        "n_by_campaign_phase": counts("campaign_phase"),
        "n_by_trait_name": counts("trait_name"),
        "screening_label_counts": {str(key): int(value) for key, value in labels.items()},
        "interpretation": (
            "Cumulative machine-only source-backed candidates and screening labels. "
            "They are not accepted trait values, biological absences, or confirmatory-model inputs."
        ),
    }


@app.command("build")
def build(
    campaign_dir: Path = typer.Option(..., exists=True),
    output_dir: Path | None = typer.Option(None),
) -> None:
    """Write cumulative candidate and counterfactual-screening indexes."""
    destination = output_dir or campaign_dir
    destination.mkdir(parents=True, exist_ok=True)
    raw, wave_ids = load_wave_candidates(campaign_dir)
    index = build_candidate_index(raw)
    ledger = load_campaign_ledger(campaign_dir)
    screening = build_counterfactual_screening(index, ledger)
    status = build_status(index, screening, wave_ids)

    _write_gzip(index, destination / "machine_candidate_index.csv.gz")
    _write_gzip(screening, destination / "counterfactual_screening.csv.gz")
    (destination / "candidate_index_status.json").write_text(
        json.dumps(status, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(
        f"Indexed {status['n_unique_candidate_rows']} machine candidate row(s) from "
        f"{status['n_waves_indexed']} wave(s); {len(screening)} biotic screening species."
    )


if __name__ == "__main__":
    app()
