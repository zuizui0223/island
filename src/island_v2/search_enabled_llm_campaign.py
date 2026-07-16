"""Build and ingest resumable search-enabled LLM trait jobs at global scale.

This module does not call a specific model provider. It creates deterministic jobs for a
search-enabled LLM, validates returned 9-column CSV rows, and maintains a checkpoint ledger.
The same job manifest can be processed by any provider that can browse/search the web.
Returned payloads must echo the job's ``prompt_sha256``. Only rows that match the
checkpointed prompt and pass schema validation enter the validated-results file.
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

PROMPT_VERSION = "search_enabled_v1_style_v2"
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
POLLINATION = {
    "bees",
    "bumblebees",
    "flies",
    "butterflies",
    "moths",
    "birds",
    "wind",
    "self",
    "mixed",
    "unknown",
}
MATING = {"obligate_outcrossing", "mixed_mating", "mainly_selfing", "obligate_selfing", "unknown"}
SELF_INCOMPATIBILITY = {"SI", "SC", "likely_SI", "likely_SC", "unknown"}
EVIDENCE_TYPES = {"field_study", "review", "flora", "horticulture", "inference"}
CONFIDENCE = {"high", "medium", "low"}
TERMINAL = {"completed", "invalid"}
VALIDATED_RESULTS_FILENAME = "validated_trait_results.csv"
VALIDATED_RESULT_COLUMNS = [
    *OUTPUT_COLUMNS,
    "job_id",
    "job_sha256",
    "target_fields",
    "model_provider",
    "model_id",
    "prompt_version",
    "prompt_sha256",
    "result_sha256",
]
JOB_HASH_COLUMNS = [
    "species",
    "genus",
    "family",
    "target_fields",
    "prompt_version",
    "prompt_sha256",
    "prompt_policy_sha256",
]
JOB_HISTORY_FILENAME = "job_history.csv"
JOB_HISTORY_COLUMNS = [
    "job_id",
    "job_sha256",
    "species",
    "target_fields",
    "prompt_version",
    "prompt_sha256",
    "result_sha256",
    "status",
]
REJECTED_RESULT_COLUMNS = [
    "results_jsonl_sha256",
    "line_number",
    "job_id",
    "error",
    "raw",
]
TRAIT_FIELDS = [
    "flower_color",
    "flower_shape",
    "pollination_guild",
    "mating_system",
    "self_incompatibility",
]

PROMPT_TEMPLATE = """You are an expert in plant reproductive biology and pollination ecology.

Search the web broadly for the focal species below. Information source type is not restricted.
Use direct species-level information when available. When evidence is incomplete, reasonable
inference from genus/family biology is allowed, but mark evidence_type=inference and use
likely_SI/likely_SC when self-incompatibility is inferred rather than directly reported.

Focal species: {species}
Family: {family}
Genus: {genus}
Fields still unresolved and eligible for this job: {target_fields}

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
- Set every trait field not listed in "Fields still unresolved" to unknown.
- If information is insufficient, use unknown and confidence=low.
"""


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_result_sha(row: object) -> str:
    """Hash the normalized nine-column result contract in its declared order."""
    getter = getattr(row, "get", None)
    if getter is None:
        raise TypeError("row must provide a get method")
    material = json.dumps(
        {column: _text(getter(column, "")) for column in OUTPUT_COLUMNS},
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=False,
    )
    return _sha(material)


def canonical_job_sha(row: object) -> str:
    """Hash the provider-independent job identity used by jobs and ledgers."""
    getter = getattr(row, "get", None)
    if getter is None:
        raise TypeError("row must provide a get method")
    material = json.dumps(
        {column: _text(getter(column, "")) for column in JOB_HASH_COLUMNS},
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=False,
    )
    return _sha(material)


def _prompt_policy_sha() -> str:
    material = {
        "prompt_version": PROMPT_VERSION,
        "prompt_template": PROMPT_TEMPLATE,
        "output_columns": OUTPUT_COLUMNS,
        "trait_fields": TRAIT_FIELDS,
        "pollination": sorted(POLLINATION),
        "mating": sorted(MATING),
        "self_incompatibility": sorted(SELF_INCOMPATIBILITY),
        "evidence_types": sorted(EVIDENCE_TYPES),
        "confidence": sorted(CONFIDENCE),
    }
    return _sha(
        json.dumps(
            material,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
    )


def _frame_fingerprint(frame: pd.DataFrame, columns: list[str]) -> str:
    records = [
        {column: _text(record.get(column)) for column in columns}
        for record in frame[columns].to_dict("records")
    ]
    return _sha(
        json.dumps(
            records,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=False,
        )
    )


def _stable_shard(species: str, shard_count: int) -> int:
    return int(_sha(species)[:16], 16) % shard_count


def _load_master(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str).fillna("")
    if "accepted_species" not in frame.columns:
        raise typer.BadParameter("master CSV must contain accepted_species")
    if "genus" not in frame.columns:
        frame["genus"] = frame["accepted_species"].map(
            lambda x: _text(x).split()[0] if _text(x) else ""
        )
    if "family" not in frame.columns:
        frame["family"] = ""
    for column in ("accepted_species", "genus", "family"):
        frame[column] = frame[column].map(_text)
    frame = frame.loc[frame["accepted_species"].ne("")].copy()
    frame = frame.sort_values("accepted_species").drop_duplicates("accepted_species", keep="first")
    return frame[["accepted_species", "genus", "family"]].reset_index(drop=True)


def _scope_master(
    master: pd.DataFrame,
    applicability_csv: Path | None,
    current_traits_csv: Path | None,
) -> tuple[pd.DataFrame, dict[str, int], set[str]]:
    master_species = set(master["accepted_species"])
    eligible_species = set(master_species)
    if applicability_csv is not None:
        applicability = pd.read_csv(applicability_csv, dtype=str).fillna("")
        required = {"accepted_species", "angiosperm_analysis_eligible"}
        missing = required.difference(applicability.columns)
        if missing:
            raise typer.BadParameter(f"applicability CSV is missing columns: {sorted(missing)}")
        if applicability["accepted_species"].duplicated().any():
            raise typer.BadParameter("applicability CSV contains duplicate species")
        if set(applicability["accepted_species"].map(_text)) != master_species:
            raise typer.BadParameter("applicability CSV species do not exactly match the master")
        eligible = applicability["angiosperm_analysis_eligible"].map(
            lambda value: _text(value).casefold() in {"true", "1", "yes"}
        )
        eligible_species = set(applicability.loc[eligible, "accepted_species"].map(_text))

    target_fields_by_species = {species: list(TRAIT_FIELDS) for species in eligible_species}
    if current_traits_csv is not None:
        current = pd.read_csv(current_traits_csv, dtype=str).fillna("")
        required = {"species", *TRAIT_FIELDS}
        missing = required.difference(current.columns)
        if missing:
            raise typer.BadParameter(f"current traits CSV is missing columns: {sorted(missing)}")
        if current["species"].duplicated().any():
            raise typer.BadParameter("current traits CSV contains duplicate species")
        if set(current["species"].map(_text)) != master_species:
            raise typer.BadParameter("current traits CSV species do not exactly match the master")
        target_fields_by_species = {}
        for record in current.to_dict("records"):
            species = _text(record.get("species"))
            if species not in eligible_species:
                continue
            gaps = [
                field
                for field in TRAIT_FIELDS
                if _text(record.get(field)).casefold() in {"", "unknown"}
            ]
            if gaps:
                target_fields_by_species[species] = gaps

    scoped = master.loc[master["accepted_species"].isin(target_fields_by_species)].copy()
    scoped["target_fields"] = scoped["accepted_species"].map(
        lambda species: "|".join(target_fields_by_species[_text(species)])
    )
    scoped = scoped.sort_values("accepted_species").reset_index(drop=True)
    return (
        scoped,
        {
            "n_master_species": len(master),
            "n_trait_applicable_species": len(eligible_species),
            "n_trait_not_applicable_species": len(master) - len(eligible_species),
            "n_species_with_unresolved_fields": len(scoped),
            "n_species_all_fields_known": len(eligible_species) - len(scoped),
        },
        eligible_species,
    )


def _job_record(row: object, shard_index: int, shard_count: int) -> dict[str, str]:
    getter = getattr(row, "get", None)
    if getter is None:
        raise TypeError("row must provide a get method")
    species = _text(getter("species"))
    genus = _text(getter("genus"))
    family = _text(getter("family"))
    target_fields = _text(getter("target_fields"))
    prompt = PROMPT_TEMPLATE.format(
        species=species,
        genus=genus,
        family=family,
        target_fields=", ".join(target_fields.split("|")),
    )
    job = {
        "species": species,
        "genus": genus,
        "family": family,
        "target_fields": target_fields,
        "shard_index": str(shard_index),
        "shard_count": str(shard_count),
        "prompt_version": PROMPT_VERSION,
        "prompt_sha256": _sha(prompt),
        "prompt_policy_sha256": _prompt_policy_sha(),
        "prompt": prompt,
    }
    job_sha = canonical_job_sha(job)
    return {"job_id": f"trait_{job_sha[:20]}", "job_sha256": job_sha, **job}


def _load_ledger(
    path: Path,
    master: pd.DataFrame,
    fingerprints: dict[str, str],
    completed_job_ids: set[str],
) -> pd.DataFrame:
    old_by_species: dict[str, dict[str, str]] = {}
    if path.exists():
        old = pd.read_csv(path, dtype=str).fillna("")
        if "species" not in old.columns or old["species"].duplicated().any():
            raise typer.BadParameter("campaign ledger must contain one row per species")
        old_by_species = {_text(record.get("species")): record for record in old.to_dict("records")}

    records: list[dict[str, object]] = []
    for source in master.rename(columns={"accepted_species": "species"}).to_dict("records"):
        expected = _job_record(source, 0, 1)
        old_record = old_by_species.get(_text(source.get("species")), {})
        same_job = all(
            _text(old_record.get(column)) == _text(expected.get(column))
            for column in (
                "job_id",
                "job_sha256",
                "target_fields",
                "prompt_version",
                "prompt_sha256",
            )
        )
        old_status = _text(old_record.get("status"))
        status = (
            old_status
            if same_job and old_status in {*TERMINAL, "pending", "queued", "retry"}
            else "pending"
        )
        if status == "completed" and expected["job_id"] not in completed_job_ids:
            status = "pending"
        old_attempts = pd.to_numeric(_text(old_record.get("attempts")), errors="coerce")
        attempts = int(old_attempts) if same_job and pd.notna(old_attempts) else 0
        records.append(
            {
                "species": expected["species"],
                "genus": expected["genus"],
                "family": expected["family"],
                "target_fields": expected["target_fields"],
                "status": status,
                "attempts": attempts,
                "job_id": expected["job_id"],
                "job_sha256": expected["job_sha256"],
                "prompt_version": expected["prompt_version"],
                "prompt_sha256": expected["prompt_sha256"],
                "result_sha256": (
                    _text(old_record.get("result_sha256")) if status == "completed" else ""
                ),
                "last_error": _text(old_record.get("last_error")) if same_job else "",
                **fingerprints,
            }
        )
    columns = [
        "species",
        "genus",
        "family",
        "target_fields",
        "status",
        "attempts",
        "job_id",
        "job_sha256",
        "prompt_version",
        "prompt_sha256",
        "result_sha256",
        "last_error",
        *fingerprints,
    ]
    return pd.DataFrame(records, columns=columns).sort_values("species").reset_index(drop=True)


def _load_job_history(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=JOB_HISTORY_COLUMNS)
    frame = pd.read_csv(path, dtype=str).fillna("")
    for column in JOB_HISTORY_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""
    frame = frame[JOB_HISTORY_COLUMNS]
    if frame["job_id"].duplicated().any():
        raise typer.BadParameter("job history contains duplicate job_id values")
    return frame


def _load_validated_results(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=VALIDATED_RESULT_COLUMNS)
    frame = pd.read_csv(path, dtype=str).fillna("")
    for column in VALIDATED_RESULT_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""
    return frame[VALIDATED_RESULT_COLUMNS]


def _upsert_history(history: pd.DataFrame, record: dict[str, str]) -> pd.DataFrame:
    matches = history.index[history["job_id"].eq(record["job_id"])].tolist()
    if not matches:
        return pd.concat(
            [history, pd.DataFrame([record], columns=JOB_HISTORY_COLUMNS)],
            ignore_index=True,
        )
    index = matches[0]
    immutable = (
        "job_sha256",
        "species",
        "target_fields",
        "prompt_version",
        "prompt_sha256",
    )
    for column in immutable:
        existing = _text(history.loc[index, column])
        incoming = _text(record.get(column))
        if existing and existing != incoming:
            raise typer.BadParameter(
                f"job history identity mismatch for {record['job_id']}: {column}"
            )
        history.loc[index, column] = incoming
    history.loc[index, "result_sha256"] = _text(record.get("result_sha256"))
    history.loc[index, "status"] = _text(record.get("status"))
    return history


def _prune_stale_validated_results(
    *,
    campaign_dir: Path,
    full_master: pd.DataFrame,
    eligible_species: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    results_path = campaign_dir / VALIDATED_RESULTS_FILENAME
    history_path = campaign_dir / JOB_HISTORY_FILENAME
    validated = _load_validated_results(results_path)
    history = _load_job_history(history_path)
    if validated.empty:
        completed_without_result = history["status"].eq("completed")
        history.loc[completed_without_result, "status"] = "stale"
        return validated, history, 0

    duplicate_job_ids = set(
        validated.loc[validated["job_id"].duplicated(keep=False), "job_id"].map(_text)
    )
    history_by_job = {_text(record.get("job_id")): record for record in history.to_dict("records")}
    master_by_species = {
        _text(record.get("accepted_species")): record for record in full_master.to_dict("records")
    }
    keep: list[bool] = []
    stale_job_ids: set[str] = set()
    active_job_ids: set[str] = set()
    for record in validated.to_dict("records"):
        job_id = _text(record.get("job_id"))
        species = _text(record.get("species"))
        source = master_by_species.get(species)
        historical = history_by_job.get(job_id)
        target_fields = _text(record.get("target_fields"))
        fields = [field for field in target_fields.split("|") if field]
        expected = None
        if source is not None:
            expected = _job_record(
                {
                    "species": species,
                    "genus": source.get("genus"),
                    "family": source.get("family"),
                    "target_fields": target_fields,
                },
                0,
                1,
            )
        valid = bool(
            job_id
            and job_id not in duplicate_job_ids
            and species in eligible_species
            and historical is not None
            and _text(historical.get("status")) == "completed"
            and expected is not None
            and fields
            and set(fields).issubset(TRAIT_FIELDS)
            and _text(record.get("job_sha256")) == _text(expected.get("job_sha256"))
            and job_id == _text(expected.get("job_id"))
            and _text(record.get("prompt_version")) == PROMPT_VERSION
            and _text(record.get("prompt_sha256")) == _text(expected.get("prompt_sha256"))
            and _text(record.get("result_sha256")) == canonical_result_sha(record)
            and all(
                _text(record.get(column)) == _text(historical.get(column))
                for column in (
                    "job_id",
                    "job_sha256",
                    "species",
                    "target_fields",
                    "prompt_version",
                    "prompt_sha256",
                    "result_sha256",
                )
            )
            and all(
                field in fields or _text(record.get(field)) == "unknown" for field in TRAIT_FIELDS
            )
        )
        keep.append(valid)
        if valid:
            active_job_ids.add(job_id)
        elif job_id:
            stale_job_ids.add(job_id)

    history.loc[history["job_id"].isin(stale_job_ids), "status"] = "stale"
    completed_without_result = history["status"].eq("completed") & ~history["job_id"].isin(
        active_job_ids
    )
    history.loc[completed_without_result, "status"] = "stale"
    pruned = validated.loc[keep, VALIDATED_RESULT_COLUMNS].reset_index(drop=True)
    return pruned, history, len(validated) - len(pruned)


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
    applicability_csv: Path | None = typer.Option(None, exists=True),
    current_traits_csv: Path | None = typer.Option(None, exists=True),
    shard_index: int = typer.Option(0, min=0),
    shard_count: int = typer.Option(128, min=1, max=4096),
    batch_size: int = typer.Option(1000, min=1, max=10000),
    retry_invalid: bool = typer.Option(False),
    resume_queued: bool = typer.Option(
        False,
        help="Re-emit queued jobs after an interrupted provider submission.",
    ),
) -> None:
    if shard_index >= shard_count:
        raise typer.BadParameter("shard_index must be smaller than shard_count")
    campaign_dir.mkdir(parents=True, exist_ok=True)
    ledger_path = campaign_dir / "campaign_ledger.csv"
    full_master = _load_master(master_csv)
    master, scope_report, eligible_species = _scope_master(
        full_master,
        applicability_csv=applicability_csv,
        current_traits_csv=current_traits_csv,
    )
    validated, history, n_stale_results_removed = _prune_stale_validated_results(
        campaign_dir=campaign_dir,
        full_master=full_master,
        eligible_species=eligible_species,
    )
    fingerprints = {
        "prompt_policy_sha256": _prompt_policy_sha(),
        "master_csv_sha256": _file_sha(master_csv),
        "master_fingerprint_sha256": _frame_fingerprint(
            full_master, ["accepted_species", "genus", "family"]
        ),
        "applicability_csv_sha256": (
            _file_sha(applicability_csv) if applicability_csv is not None else ""
        ),
        "current_traits_csv_sha256": (
            _file_sha(current_traits_csv) if current_traits_csv is not None else ""
        ),
        "scope_fingerprint_sha256": _frame_fingerprint(
            master, ["accepted_species", "genus", "family", "target_fields"]
        ),
    }
    completed_job_ids = set(history.loc[history["status"].eq("completed"), "job_id"].map(_text))
    ledger = _load_ledger(
        ledger_path,
        master,
        fingerprints,
        completed_job_ids,
    )
    ledger["shard"] = ledger["species"].map(lambda value: _stable_shard(value, shard_count))
    allowed_status = {"pending", "retry"}
    if retry_invalid:
        allowed_status.add("invalid")
    if resume_queued:
        allowed_status.add("queued")
    eligible = ledger.loc[
        ledger["shard"].eq(shard_index) & ledger["status"].isin(allowed_status)
    ].copy()
    if resume_queued:
        eligible["resume_priority"] = eligible["status"].ne("queued").astype(int)
        eligible = eligible.sort_values(["resume_priority", "species"])
    eligible = eligible.head(batch_size)

    jobs: list[dict[str, str]] = []
    for index, row in eligible.iterrows():
        job = _job_record(row, shard_index, shard_count)
        jobs.append(job)
        ledger.loc[index, "status"] = "queued"
        ledger.loc[index, "attempts"] = int(ledger.loc[index, "attempts"]) + 1
        ledger.loc[index, "last_error"] = ""
        history = _upsert_history(
            history,
            {
                "job_id": job["job_id"],
                "job_sha256": job["job_sha256"],
                "species": job["species"],
                "target_fields": job["target_fields"],
                "prompt_version": job["prompt_version"],
                "prompt_sha256": job["prompt_sha256"],
                "result_sha256": "",
                "status": "queued",
            },
        )

    job_path = campaign_dir / f"jobs_shard_{shard_index:04d}.jsonl"
    job_path.write_text(
        "\n".join(json.dumps(job, ensure_ascii=False, sort_keys=True) for job in jobs)
        + ("\n" if jobs else ""),
        encoding="utf-8",
    )
    ledger.drop(columns=["shard"], errors="ignore").to_csv(ledger_path, index=False)
    validated.to_csv(campaign_dir / VALIDATED_RESULTS_FILENAME, index=False)
    history.to_csv(campaign_dir / JOB_HISTORY_FILENAME, index=False)
    input_manifest = {
        "contract_version": "2.0",
        **fingerprints,
        "n_master_species": len(full_master),
        "n_trait_applicable_species": len(eligible_species),
        "n_campaign_species": len(master),
        "master_csv": str(master_csv),
        "applicability_csv": str(applicability_csv) if applicability_csv else "",
        "current_traits_csv": str(current_traits_csv) if current_traits_csv else "",
    }
    (campaign_dir / "campaign_input_manifest.json").write_text(
        json.dumps(input_manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    status = {
        **scope_report,
        "n_global_species": len(full_master),
        "n_campaign_species": len(master),
        "shard_index": shard_index,
        "shard_count": shard_count,
        "n_species_in_shard": int(ledger["shard"].eq(shard_index).sum()),
        "n_jobs_prepared": len(jobs),
        "n_completed_total": int(ledger["status"].eq("completed").sum()),
        "n_pending_total": int(ledger["status"].isin({"pending", "retry"}).sum()),
        "n_queued_total": int(ledger["status"].eq("queued").sum()),
        "n_validated_job_results": len(validated),
        "n_stale_validated_results_removed": n_stale_results_removed,
        "resume_queued": resume_queued,
        "job_manifest": str(job_path),
        "prompt_version": PROMPT_VERSION,
        **fingerprints,
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
    required_ledger_columns = {
        "species",
        "genus",
        "family",
        "target_fields",
        "status",
        "job_id",
        "job_sha256",
        "prompt_version",
        "prompt_sha256",
        "result_sha256",
        "last_error",
    }
    missing_ledger_columns = required_ledger_columns.difference(ledger.columns)
    if missing_ledger_columns:
        raise typer.BadParameter(
            "campaign ledger predates hash-locked prompts; run prepare to migrate it: "
            f"{sorted(missing_ledger_columns)}"
        )
    if ledger["species"].duplicated().any() or ledger["job_id"].duplicated().any():
        raise typer.BadParameter("campaign ledger species and job_id values must be unique")
    for record in ledger.to_dict("records"):
        expected = _job_record(record, 0, 1)
        for column in ("job_id", "job_sha256", "prompt_version", "prompt_sha256"):
            if _text(record.get(column)) != _text(expected.get(column)):
                raise typer.BadParameter(
                    f"campaign ledger has stale {column} for {_text(record.get('species'))}"
                )

    by_job = {row["job_id"]: index for index, row in ledger.iterrows() if _text(row.get("job_id"))}
    history_path = campaign_dir / JOB_HISTORY_FILENAME
    history = _load_job_history(history_path)
    history_by_job = {_text(record.get("job_id")): record for record in history.to_dict("records")}
    for record in ledger.loc[ledger["status"].eq("queued")].to_dict("records"):
        historical = history_by_job.get(_text(record.get("job_id")))
        if historical is None or any(
            _text(historical.get(column)) != _text(record.get(column))
            for column in (
                "job_id",
                "job_sha256",
                "species",
                "target_fields",
                "prompt_version",
                "prompt_sha256",
            )
        ):
            raise typer.BadParameter(
                f"queued job is not bound to job history: {_text(record.get('job_id'))}"
            )

    accepted_by_job: dict[str, dict[str, str]] = {}
    rejected: list[dict[str, str]] = []
    errors_by_job: dict[str, list[str]] = {}
    results_sha = _file_sha(results_jsonl)

    for line_number, raw in enumerate(
        results_jsonl.read_text(encoding="utf-8").splitlines(), start=1
    ):
        if not raw.strip():
            continue
        job_id = ""
        try:
            payload = json.loads(raw)
            if not isinstance(payload, dict):
                raise ValueError("result payload must be a JSON object")
            job_id = _text(payload.get("job_id"))
            result_text = _text(payload.get("result"))
            if job_id not in by_job:
                raise ValueError("unknown job_id")
            index = by_job[job_id]
            job_status = _text(ledger.loc[index, "status"])
            if job_status == "completed":
                raise ValueError("job already completed")
            if job_status != "queued":
                raise ValueError(f"job status must be queued, got {job_status!r}")
            expected_prompt_sha = _text(ledger.loc[index, "prompt_sha256"])
            payload_prompt_sha = _text(payload.get("prompt_sha256"))
            if not payload_prompt_sha:
                raise ValueError("missing prompt_sha256")
            if not expected_prompt_sha:
                raise ValueError("campaign ledger is missing prompt_sha256")
            if payload_prompt_sha != expected_prompt_sha:
                raise ValueError("prompt_sha256 mismatch")
            payload_job_sha = _text(payload.get("job_sha256"))
            expected_job_sha = _text(ledger.loc[index, "job_sha256"])
            if payload_job_sha and payload_job_sha != expected_job_sha:
                raise ValueError("job_sha256 mismatch")
            expected_species = _text(ledger.loc[index, "species"])
            row = _parse_csv_row(result_text)
            _validate_row(row, expected_species)
            target_fields = set(_text(ledger.loc[index, "target_fields"]).split("|"))
            for field in TRAIT_FIELDS:
                if field not in target_fields and row[field] != "unknown":
                    raise ValueError(f"non-target field {field} must be unknown")
            if job_id in accepted_by_job:
                raise ValueError("duplicate valid result for job_id in the same JSONL")
            accepted_record = {
                **row,
                "job_id": job_id,
                "job_sha256": expected_job_sha,
                "target_fields": _text(ledger.loc[index, "target_fields"]),
                "model_provider": _text(model_provider),
                "model_id": _text(model_id),
                "prompt_version": _text(ledger.loc[index, "prompt_version"]),
                "prompt_sha256": payload_prompt_sha,
                "result_sha256": "",
            }
            accepted_record["result_sha256"] = canonical_result_sha(accepted_record)
            accepted_by_job[job_id] = accepted_record
        except Exception as exc:
            error = str(exc)
            rejected.append(
                {
                    "results_jsonl_sha256": results_sha,
                    "line_number": str(line_number),
                    "job_id": job_id,
                    "error": error,
                    "raw": raw,
                }
            )
            if job_id in by_job:
                errors_by_job.setdefault(job_id, []).append(error)

    for job_id, accepted_record in accepted_by_job.items():
        index = by_job[job_id]
        result_sha = accepted_record["result_sha256"]
        ledger.loc[index, "status"] = "completed"
        ledger.loc[index, "result_sha256"] = result_sha
        ledger.loc[index, "last_error"] = ""
        history = _upsert_history(
            history,
            {
                "job_id": job_id,
                "job_sha256": accepted_record["job_sha256"],
                "species": accepted_record["species"],
                "target_fields": accepted_record["target_fields"],
                "prompt_version": accepted_record["prompt_version"],
                "prompt_sha256": accepted_record["prompt_sha256"],
                "result_sha256": result_sha,
                "status": "completed",
            },
        )
    for job_id, errors in errors_by_job.items():
        if job_id in accepted_by_job:
            continue
        index = by_job[job_id]
        if _text(ledger.loc[index, "status"]) != "queued":
            continue
        ledger.loc[index, "status"] = "invalid"
        ledger.loc[index, "result_sha256"] = ""
        ledger.loc[index, "last_error"] = errors[-1]
        history = _upsert_history(
            history,
            {
                "job_id": job_id,
                "job_sha256": _text(ledger.loc[index, "job_sha256"]),
                "species": _text(ledger.loc[index, "species"]),
                "target_fields": _text(ledger.loc[index, "target_fields"]),
                "prompt_version": _text(ledger.loc[index, "prompt_version"]),
                "prompt_sha256": _text(ledger.loc[index, "prompt_sha256"]),
                "result_sha256": "",
                "status": "invalid",
            },
        )

    results_path = campaign_dir / VALIDATED_RESULTS_FILENAME
    old = _load_validated_results(results_path)
    new = pd.DataFrame(accepted_by_job.values(), columns=VALIDATED_RESULT_COLUMNS)
    validated = pd.concat([old, new], ignore_index=True)
    if not validated.empty:
        validated = validated.drop_duplicates("job_id", keep="last")
    validated.to_csv(results_path, index=False)
    rejected_path = campaign_dir / "rejected_results.csv"
    if rejected_path.exists():
        old_rejected = pd.read_csv(rejected_path, dtype=str).fillna("")
        for column in REJECTED_RESULT_COLUMNS:
            if column not in old_rejected.columns:
                old_rejected[column] = ""
        old_rejected = old_rejected[REJECTED_RESULT_COLUMNS]
    else:
        old_rejected = pd.DataFrame(columns=REJECTED_RESULT_COLUMNS)
    rejected_frame = pd.concat(
        [old_rejected, pd.DataFrame(rejected, columns=REJECTED_RESULT_COLUMNS)],
        ignore_index=True,
    )
    rejected_frame.to_csv(rejected_path, index=False)
    history.to_csv(history_path, index=False)
    ledger.to_csv(ledger_path, index=False)

    coverage = {
        column: (
            int(
                validated.loc[
                    validated[column].ne("unknown") & validated[column].ne(""), "species"
                ].nunique()
            )
            if not validated.empty
            else 0
        )
        for column in (
            "flower_color",
            "flower_shape",
            "pollination_guild",
            "mating_system",
            "self_incompatibility",
        )
    }
    report = {
        "n_accepted": len(accepted_by_job),
        "n_rejected": len(rejected),
        "n_validated_total": len(validated),
        "n_validated_species": int(validated["species"].nunique()) if len(validated) else 0,
        "n_completed_total": int(ledger["status"].eq("completed").sum()),
        "coverage": coverage,
        "model_provider": model_provider,
        "model_id": model_id,
        "validated_results_file": str(results_path),
    }
    (campaign_dir / "ingest_status.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
