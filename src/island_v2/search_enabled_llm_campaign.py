"""Build and ingest resumable search-enabled LLM trait jobs at global scale.

This module does not call a specific model provider. It creates deterministic jobs for a
search-enabled LLM, validates returned 9-column CSV rows, and maintains a checkpoint ledger.
The same job manifest can be processed by any provider that can browse/search the web.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

PROMPT_VERSION = "search_enabled_v1_style_v1"
OUTPUT_COLUMNS = [
    "species",
    "flower_color",
    "flower_shape",
    "pollination_guild",
    "pollination_notes",
    "mating_system",
    "self_incompatibility",
    "evidence_type",
    "confidence",
]
POLLINATION = {"bees", "bumblebees", "flies", "butterflies", "moths", "birds", "wind", "self", "mixed", "unknown"}
MATING = {"obligate_outcrossing", "mixed_mating", "mainly_selfing", "obligate_selfing", "unknown"}
SELF_INCOMPATIBILITY = {"SI", "SC", "likely_SI", "likely_SC", "unknown"}
EVIDENCE_TYPES = {"field_study", "review", "flora", "horticulture", "inference"}
CONFIDENCE = {"high", "medium", "low"}
TERMINAL = {"completed", "invalid"}

PROMPT_TEMPLATE = """You are an expert in plant reproductive biology and pollination ecology.

Search the web broadly for the focal species below. Information source type is not restricted.
Use direct species-level information when available. When evidence is incomplete, reasonable
inference from genus/family biology is allowed, but mark evidence_type=inference and use
likely_SI/likely_SC when self-incompatibility is inferred rather than directly reported.

Focal species: {species}
Family: {family}
Genus: {genus}

Extract:
- typical flower color(s)
- flower shape / symmetry
- main pollination guild: bees, bumblebees, flies, butterflies, moths, birds, wind, self, mixed, unknown
- mating system: obligate_outcrossing, mixed_mating, mainly_selfing, obligate_selfing, unknown
- self-incompatibility: SI, SC, likely_SI, likely_SC, unknown
- short free-text notes on pollination biology
- evidence type: field_study, review, flora, horticulture, inference
- confidence: high, medium, low

Return exactly one CSV row and no prose, with columns in this order:
species,flower_color,flower_shape,pollination_guild,pollination_notes,mating_system,self_incompatibility,evidence_type,confidence

Rules:
- Fields containing commas must be enclosed in double quotes.
- Do not equate SC with autonomous selfing.
- If direct species-level evidence is absent, inferred values are allowed but must remain visibly lower-confidence.
- If information is insufficient, use unknown and confidence=low.
"""


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _stable_shard(species: str, shard_count: int) -> int:
    return int(_sha(species)[:16], 16) % shard_count


def _load_master(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str).fillna("")
    if "accepted_species" not in frame.columns:
        raise typer.BadParameter("master CSV must contain accepted_species")
    if "genus" not in frame.columns:
        frame["genus"] = frame["accepted_species"].map(lambda x: _text(x).split()[0] if _text(x) else "")
    if "family" not in frame.columns:
        frame["family"] = ""
    for column in ("accepted_species", "genus", "family"):
        frame[column] = frame[column].map(_text)
    frame = frame.loc[frame["accepted_species"].ne("")].copy()
    frame = frame.sort_values("accepted_species").drop_duplicates("accepted_species", keep="first")
    return frame[["accepted_species", "genus", "family"]].reset_index(drop=True)


def _load_ledger(path: Path, master: pd.DataFrame) -> pd.DataFrame:
    base = master.rename(columns={"accepted_species": "species"}).copy()
    if path.exists():
        old = pd.read_csv(path, dtype=str).fillna("")
        keep = [column for column in old.columns if column not in {"genus", "family"}]
        base = base.merge(old[keep], on="species", how="left", validate="one_to_one")
    defaults = {
        "status": "pending",
        "attempts": "0",
        "job_id": "",
        "prompt_version": PROMPT_VERSION,
        "prompt_sha256": "",
        "result_sha256": "",
        "last_error": "",
    }
    for column, value in defaults.items():
        if column not in base.columns:
            base[column] = value
        base[column] = base[column].fillna("").replace("", value)
    base["attempts"] = pd.to_numeric(base["attempts"], errors="coerce").fillna(0).astype(int)
    return base.sort_values("species").reset_index(drop=True)


def _job_record(row: pd.Series, shard_index: int, shard_count: int) -> dict[str, str]:
    species = _text(row["species"])
    genus = _text(row["genus"])
    family = _text(row["family"])
    prompt = PROMPT_TEMPLATE.format(species=species, genus=genus, family=family)
    job_id = f"trait_{_sha(species)[:20]}"
    return {
        "job_id": job_id,
        "species": species,
        "genus": genus,
        "family": family,
        "shard_index": str(shard_index),
        "shard_count": str(shard_count),
        "prompt_version": PROMPT_VERSION,
        "prompt_sha256": _sha(prompt),
        "prompt": prompt,
    }


def _parse_csv_row(text: str) -> dict[str, str]:
    rows = list(csv.reader(io.StringIO(text.strip())))
    if len(rows) != 1 or len(rows[0]) != len(OUTPUT_COLUMNS):
        raise ValueError("result must contain exactly one 9-column CSV row")
    return dict(zip(OUTPUT_COLUMNS, [_text(value) for value in rows[0]], strict=True))


def _validate_row(row: dict[str, str], expected_species: str) -> None:
    if row["species"] != expected_species:
        raise ValueError(f"species mismatch: expected {expected_species!r}, got {row['species']!r}")
    if row["pollination_guild"] not in POLLINATION:
        raise ValueError("invalid pollination_guild")
    if row["mating_system"] not in MATING:
        raise ValueError("invalid mating_system")
    if row["self_incompatibility"] not in SELF_INCOMPATIBILITY:
        raise ValueError("invalid self_incompatibility")
    if row["evidence_type"] not in EVIDENCE_TYPES:
        raise ValueError("invalid evidence_type")
    if row["confidence"] not in CONFIDENCE:
        raise ValueError("invalid confidence")


@app.command("prepare")
def prepare(
    master_csv: Path = typer.Option(..., exists=True),
    campaign_dir: Path = typer.Option(...),
    shard_index: int = typer.Option(0, min=0),
    shard_count: int = typer.Option(128, min=1, max=4096),
    batch_size: int = typer.Option(1000, min=1, max=10000),
    retry_invalid: bool = typer.Option(False),
) -> None:
    if shard_index >= shard_count:
        raise typer.BadParameter("shard_index must be smaller than shard_count")
    campaign_dir.mkdir(parents=True, exist_ok=True)
    ledger_path = campaign_dir / "campaign_ledger.csv"
    master = _load_master(master_csv)
    ledger = _load_ledger(ledger_path, master)
    ledger["shard"] = ledger["species"].map(lambda value: _stable_shard(value, shard_count))
    allowed_status = {"pending", "retry"}
    if retry_invalid:
        allowed_status.add("invalid")
    eligible = ledger.loc[
        ledger["shard"].eq(shard_index) & ledger["status"].isin(allowed_status)
    ].head(batch_size).copy()

    jobs: list[dict[str, str]] = []
    for index, row in eligible.iterrows():
        job = _job_record(row, shard_index, shard_count)
        jobs.append(job)
        ledger.loc[index, "status"] = "queued"
        ledger.loc[index, "attempts"] = int(ledger.loc[index, "attempts"]) + 1
        ledger.loc[index, "job_id"] = job["job_id"]
        ledger.loc[index, "prompt_version"] = PROMPT_VERSION
        ledger.loc[index, "prompt_sha256"] = job["prompt_sha256"]
        ledger.loc[index, "last_error"] = ""

    job_path = campaign_dir / f"jobs_shard_{shard_index:04d}.jsonl"
    job_path.write_text(
        "\n".join(json.dumps(job, ensure_ascii=False, sort_keys=True) for job in jobs)
        + ("\n" if jobs else ""),
        encoding="utf-8",
    )
    ledger.drop(columns=["shard"], errors="ignore").to_csv(ledger_path, index=False)
    status = {
        "n_global_species": len(master),
        "shard_index": shard_index,
        "shard_count": shard_count,
        "n_species_in_shard": int(ledger["shard"].eq(shard_index).sum()),
        "n_jobs_prepared": len(jobs),
        "n_completed_total": int(ledger["status"].eq("completed").sum()),
        "n_pending_total": int(ledger["status"].isin({"pending", "retry"}).sum()),
        "job_manifest": str(job_path),
        "prompt_version": PROMPT_VERSION,
    }
    (campaign_dir / "campaign_status.json").write_text(
        json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(status, ensure_ascii=False))


@app.command("ingest")
def ingest(
    campaign_dir: Path = typer.Option(..., exists=True),
    results_jsonl: Path = typer.Option(..., exists=True),
    model_provider: str = typer.Option(...),
    model_id: str = typer.Option(...),
) -> None:
    ledger_path = campaign_dir / "campaign_ledger.csv"
    if not ledger_path.exists():
        raise typer.BadParameter("campaign ledger does not exist; run prepare first")
    ledger = pd.read_csv(ledger_path, dtype=str).fillna("")
    by_job = {row["job_id"]: index for index, row in ledger.iterrows() if _text(row.get("job_id"))}
    accepted: list[dict[str, str]] = []
    rejected: list[dict[str, str]] = []

    for line_number, raw in enumerate(results_jsonl.read_text(encoding="utf-8").splitlines(), start=1):
        if not raw.strip():
            continue
        try:
            payload = json.loads(raw)
            job_id = _text(payload.get("job_id"))
            result_text = _text(payload.get("result"))
            if job_id not in by_job:
                raise ValueError("unknown job_id")
            index = by_job[job_id]
            expected_species = _text(ledger.loc[index, "species"])
            row = _parse_csv_row(result_text)
            _validate_row(row, expected_species)
            result_sha = _sha(result_text)
            accepted.append(
                {
                    **row,
                    "job_id": job_id,
                    "model_provider": _text(model_provider),
                    "model_id": _text(model_id),
                    "prompt_version": _text(ledger.loc[index, "prompt_version"]),
                    "prompt_sha256": _text(ledger.loc[index, "prompt_sha256"]),
                    "result_sha256": result_sha,
                }
            )
            ledger.loc[index, "status"] = "completed"
            ledger.loc[index, "result_sha256"] = result_sha
            ledger.loc[index, "last_error"] = ""
        except Exception as exc:
            rejected.append({"line_number": str(line_number), "error": str(exc), "raw": raw})
            try:
                payload = json.loads(raw)
                job_id = _text(payload.get("job_id"))
                if job_id in by_job:
                    index = by_job[job_id]
                    ledger.loc[index, "status"] = "invalid"
                    ledger.loc[index, "last_error"] = str(exc)
            except Exception:
                pass

    results_path = campaign_dir / "trait_results.csv"
    new = pd.DataFrame(accepted)
    if results_path.exists() and not new.empty:
        old = pd.read_csv(results_path, dtype=str).fillna("")
        new = pd.concat([old, new], ignore_index=True)
        new = new.drop_duplicates("species", keep="last")
    if not new.empty:
        new.to_csv(results_path, index=False)
    pd.DataFrame(rejected, columns=["line_number", "error", "raw"]).to_csv(
        campaign_dir / "rejected_results.csv", index=False
    )
    ledger.to_csv(ledger_path, index=False)

    coverage = {
        column: int(new[column].ne("unknown").sum()) if not new.empty else 0
        for column in (
            "flower_color",
            "flower_shape",
            "pollination_guild",
            "mating_system",
            "self_incompatibility",
        )
    }
    report = {
        "n_accepted": len(accepted),
        "n_rejected": len(rejected),
        "n_completed_total": int(ledger["status"].eq("completed").sum()),
        "coverage": coverage,
        "model_provider": model_provider,
        "model_id": model_id,
    }
    (campaign_dir / "ingest_status.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
