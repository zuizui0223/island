"""Acquire the Meyer et al. CC0 SC/SI synthesis for the island species master.

The Dryad output calls its sc/si column ``mating_system``, while the dataset
metadata defines those states as largely self-compatible versus largely
self-incompatible.  This importer therefore maps the source only to the v2
``self_incompatibility`` trait.  It validates the exact published file hash,
removes audited DIO-only false-SI rows, retains exact accepted-name matches, and
optionally reconciles unmatched names through strict GBIF exact synonyms.

Species-direct candidates and low-confidence genus inferences are deliberately
written to separate files.  A genus inference requires at least five direct
congeners with a unanimous state and is never represented as source-backed
species evidence.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Iterable, Mapping

import pandas as pd
import requests
import typer
import yaml

from island_v2.angiosperm_scope import classify_scope, load_config as load_scope_config

app = typer.Typer(add_completion=False, no_args_is_help=True)

CANDIDATE_COLUMNS = [
    "evidence_id",
    "source_name",
    "source_record_id",
    "accepted_species",
    "genus",
    "family",
    "name_match_method",
    "trait_name",
    "raw_value",
    "standardized_value",
    "candidate_kind",
    "evidence_scope",
    "source_type",
    "source_reliability",
    "source_citation",
    "source_url",
    "evidence_excerpt",
    "supporting_taxa",
    "inference_rule",
    "confidence",
    "inference_note",
    "needs_human_review",
    "review_status",
    "analysis_role",
]


@app.callback()
def main() -> None:
    """Prepare audited species-direct and genus-inferred SC/SI candidates."""


def _text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _key(value: Any) -> str:
    return _text(value).casefold()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise typer.BadParameter(f"Expected a mapping in {path}")
    return value


def _write_csv_gzip(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})


def _read_jsonl_gzip(path: Path) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    records: dict[str, dict[str, Any]] = {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            name = _text(record.get("source_scientific_name"))
            if name:
                records[name] = record
    return records


def _write_jsonl_gzip(records: Iterable[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(records, key=lambda row: _key(row.get("source_scientific_name")))
    if path.suffix == ".gz":
        with path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
                for record in ordered:
                    line = json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n"
                    zipped.write(line.encode("utf-8"))
    else:
        path.write_text(
            "".join(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n" for row in ordered),
            encoding="utf-8",
        )


def _fetch_gbif_record(
    name: str,
    family: str,
    endpoint: str,
    timeout_seconds: float = 30.0,
) -> dict[str, Any]:
    try:
        response = requests.get(
            endpoint,
            params={"name": name, "rank": "SPECIES"},
            headers={"User-Agent": "island-floral-v2/0.1 (academic trait reconciliation)"},
            timeout=timeout_seconds,
        )
        response.raise_for_status()
        return {
            "source_scientific_name": name,
            "source_family": family,
            "http_status": response.status_code,
            "payload": response.json(),
            "error": "",
        }
    except (requests.RequestException, ValueError) as error:
        return {
            "source_scientific_name": name,
            "source_family": family,
            "http_status": "",
            "payload": {},
            "error": f"{type(error).__name__}: {error}",
        }


def fetch_gbif_records(
    source_rows: pd.DataFrame,
    *,
    endpoint: str,
    cache_path: Path,
    workers: int = 8,
) -> dict[str, dict[str, Any]]:
    """Fetch missing GBIF name matches and persist every response for audit."""
    cached = _read_jsonl_gzip(cache_path)
    requested = {
        _text(row.genus_species): _text(row.Family)
        for row in source_rows[["genus_species", "Family"]].drop_duplicates().itertuples(index=False)
    }
    missing = [
        (name, family)
        for name, family in requested.items()
        if name not in cached or _text(cached[name].get("error"))
    ]
    if missing:
        with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
            jobs = {
                pool.submit(_fetch_gbif_record, name, family, endpoint): name
                for name, family in missing
            }
            for future in as_completed(jobs):
                record = future.result()
                cached[_text(record["source_scientific_name"])] = record
        _write_jsonl_gzip(cached.values(), cache_path)
    return {name: cached[name] for name in requested if name in cached}


def strict_gbif_resolution(
    *,
    source_name: str,
    source_family: str,
    payload: Mapping[str, Any],
    master_by_name: Mapping[str, Mapping[str, str]],
    reconciliation: Mapping[str, Any],
) -> tuple[str, str]:
    """Return (accepted master species, audit reason) for a safe GBIF synonym.

    Only exact species-rank Plantae synonyms are transferable.  The source,
    GBIF and master families must agree whenever they are present, and the GBIF
    accepted species must occur uniquely in the current master.
    """
    if _text(payload.get("matchType")) != _text(reconciliation["allowed_match_type"]):
        return "", "reject_non_exact_match"
    if _text(payload.get("rank")) != _text(reconciliation["allowed_rank"]):
        return "", "reject_non_species_rank"
    if _text(payload.get("kingdom")) != _text(reconciliation["allowed_kingdom"]):
        return "", "reject_non_plantae"
    if _text(payload.get("status")) not in set(reconciliation.get("allowed_status") or []):
        return "", "reject_not_exact_synonym"
    confidence = pd.to_numeric(payload.get("confidence"), errors="coerce")
    if pd.isna(confidence) or float(confidence) < float(reconciliation["minimum_confidence"]):
        return "", "reject_low_confidence"
    if not _text(payload.get("acceptedUsageKey")):
        return "", "reject_synonym_without_accepted_usage_key"

    accepted_name = _text(payload.get("species"))
    master_row = master_by_name.get(_key(accepted_name))
    if not accepted_name or master_row is None:
        return "", "reject_accepted_species_absent_from_master"
    master_family = _text(master_row.get("family"))
    gbif_family = _text(payload.get("family"))
    families = {_key(value) for value in (source_family, gbif_family, master_family) if _text(value)}
    if len(families) > 1:
        return "", "reject_family_conflict"
    if _key(accepted_name) == _key(source_name):
        return "", "reject_no_synonym_change"
    return _text(master_row["accepted_species"]), "accepted_exact_gbif_synonym"


def _source_record_id(row: Mapping[str, Any], row_number: int) -> str:
    for column in ("Unnamed: 0", "", "index"):
        value = _text(row.get(column))
        if value:
            return value
    return str(row_number)


def _direct_candidate(
    row: Mapping[str, Any],
    *,
    row_number: int,
    accepted_species: str,
    genus: str,
    family: str,
    match_method: str,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    source = config["source"]
    raw_value = _text(row["mating_system"]).casefold()
    value = _text(config["value_map"][raw_value])
    record_id = _source_record_id(row, row_number)
    source_name = _text(row["genus_species"])
    reconciliation_note = (
        ""
        if match_method == "accepted_name_exact"
        else f" GBIF exact synonym reconciliation maps {source_name} to {accepted_species}."
    )
    state = "largely self-compatible" if value == "SC" else "largely self-incompatible"
    return {
        "evidence_id": f"meyer2026:{record_id}:{accepted_species}",
        "source_name": source["name"],
        "source_record_id": f"Output_Data_MS.csv:{record_id}",
        "accepted_species": accepted_species,
        "genus": genus,
        "family": family,
        "name_match_method": match_method,
        "trait_name": "self_incompatibility",
        "raw_value": raw_value,
        "standardized_value": value,
        "candidate_kind": "source_backed",
        "evidence_scope": "species_direct",
        "source_type": "curated_trait_database",
        "source_reliability": "B_curated_database_or_institution",
        "source_citation": source["citation"],
        "source_url": source["dataset_url"],
        "evidence_excerpt": (
            f"The CC0 literature synthesis classifies {source_name} as {state} ({value})."
            f"{reconciliation_note}"
        ),
        "supporting_taxa": "[]",
        "inference_rule": "",
        "confidence": "medium",
        "inference_note": (
            "Species-level literature synthesis; mapped to self_incompatibility, not mating_system. "
            "DIO-only records are excluded before matching."
        ),
        "needs_human_review": True,
        "review_status": "pending",
        "analysis_role": "evidence_rich_subset",
    }


def build_genus_inferences(
    master: pd.DataFrame,
    direct_candidates: pd.DataFrame,
    config: Mapping[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build a separately labelled unanimous-genus SC/SI completion layer."""
    policy = config["genus_inference"]
    if not policy.get("enabled", False) or direct_candidates.empty:
        return pd.DataFrame(columns=CANDIDATE_COLUMNS), {
            "n_qualifying_genera": 0,
            "n_inferred_species": 0,
            "leave_one_out_n": 0,
            "leave_one_out_accuracy": None,
        }

    minimum = int(policy["minimum_direct_species"])
    purity_threshold = float(policy["minimum_modal_purity"])
    direct = direct_candidates[
        ["accepted_species", "genus", "family", "standardized_value"]
    ].drop_duplicates("accepted_species")

    stats_rows: list[dict[str, Any]] = []
    for genus, group in direct.groupby("genus", sort=True):
        counts = group["standardized_value"].value_counts()
        support = int(counts.sum())
        winning_value = _text(counts.index[0])
        purity = float(counts.iloc[0] / support)
        if support >= minimum and purity >= purity_threshold:
            stats_rows.append(
                {
                    "genus": genus,
                    "value": winning_value,
                    "support_n": support,
                    "purity": purity,
                    "supporting_taxa": sorted(group["accepted_species"].tolist()),
                }
            )
    qualifying = {row["genus"]: row for row in stats_rows}

    # Cross-validation is evaluated with the same requirements after withholding
    # one direct species. It therefore measures transfer to an unseen congener.
    loo_total = 0
    loo_correct = 0
    for genus, group in direct.groupby("genus", sort=False):
        if len(group) <= minimum:
            continue
        for held_out in group.itertuples(index=False):
            remaining = group.loc[group["accepted_species"].ne(held_out.accepted_species)]
            counts = remaining["standardized_value"].value_counts()
            support = int(counts.sum())
            purity = float(counts.iloc[0] / support) if support else 0.0
            if support < minimum or purity < purity_threshold:
                continue
            loo_total += 1
            loo_correct += _text(counts.index[0]) == _text(held_out.standardized_value)

    direct_species = set(direct["accepted_species"])
    rows: list[dict[str, Any]] = []
    source = config["source"]
    for target in master.itertuples(index=False):
        species = _text(target.accepted_species)
        genus = _text(target.genus)
        if species in direct_species or genus not in qualifying:
            continue
        stat = qualifying[genus]
        value = stat["value"]
        support = int(stat["support_n"])
        supporting_taxa = stat["supporting_taxa"]
        rows.append(
            {
                "evidence_id": f"meyer2026:genus:{genus}:{species}",
                "source_name": source["name"],
                "source_record_id": f"genus:{genus}:{value}:{species}",
                "accepted_species": species,
                "genus": genus,
                "family": _text(target.family),
                "name_match_method": "accepted_name_target_for_genus_inference",
                "trait_name": "self_incompatibility",
                "raw_value": value,
                "standardized_value": value,
                "candidate_kind": "hierarchical_inference",
                "evidence_scope": "genus_inference",
                "source_type": "curated_trait_database",
                "source_reliability": "B_curated_database_or_institution",
                "source_citation": source["citation"],
                "source_url": source["dataset_url"],
                "evidence_excerpt": (
                    f"All {support} directly classified {genus} species in the source-backed "
                    f"island overlap are {value}; this target species itself was not reported."
                ),
                "supporting_taxa": json.dumps(supporting_taxa, ensure_ascii=False),
                "inference_rule": (
                    f"Unanimous genus transfer with support_n={support}, "
                    f"minimum_direct_species={minimum}, minimum_modal_purity={purity_threshold:.3f}."
                ),
                "confidence": _text(policy["confidence"]),
                "inference_note": (
                    "Low-confidence genus inference for completion/sensitivity only; never "
                    "species-direct evidence and superseded by any species-level record."
                ),
                "needs_human_review": True,
                "review_status": "pending",
                "analysis_role": _text(policy["analysis_role"]),
            }
        )
    report = {
        "minimum_direct_species": minimum,
        "minimum_modal_purity": purity_threshold,
        "n_qualifying_genera": len(qualifying),
        "n_inferred_species": len(rows),
        "leave_one_out_n": loo_total,
        "leave_one_out_correct": loo_correct,
        "leave_one_out_accuracy": loo_correct / loo_total if loo_total else None,
    }
    return pd.DataFrame(rows, columns=CANDIDATE_COLUMNS), report


def prepare_meyer_2026(
    *,
    source_csv: Path,
    master_csv: Path,
    output_dir: Path,
    config_path: Path = Path("config/meyer_2026_reproductive_traits.yml"),
    scope_config_path: Path = Path("config/angiosperm_scope.yml"),
    resolve_gbif: bool = True,
    gbif_workers: int = 8,
) -> dict[str, Any]:
    """Validate, reconcile and write the complete auditable source overlap."""
    config = _load_yaml(config_path)
    observed_hash = _sha256(source_csv)
    expected_hash = _text(config["source"]["expected_output_sha256"])
    if observed_hash != expected_hash:
        raise typer.BadParameter(
            f"Source SHA-256 mismatch: expected {expected_hash}, observed {observed_hash}"
        )

    source = pd.read_csv(source_csv, dtype=str).fillna("")
    required_source = {"genus_species", "Family", "mating_system"}
    missing_source = required_source.difference(source.columns)
    if missing_source:
        raise typer.BadParameter(f"Meyer source lacks columns: {sorted(missing_source)}")
    expected_rows = int(config["source"]["expected_output_rows"])
    if len(source) != expected_rows:
        raise typer.BadParameter(f"Expected {expected_rows} source rows, found {len(source)}")

    master = pd.read_csv(master_csv, dtype=str).fillna("")
    required_master = {"accepted_species", "genus", "family"}
    missing_master = required_master.difference(master.columns)
    if missing_master:
        raise typer.BadParameter(f"Master lacks columns: {sorted(missing_master)}")
    master["accepted_species"] = master["accepted_species"].map(_text)
    master = master.loc[master["accepted_species"].ne("")].drop_duplicates("accepted_species")
    scope = classify_scope(master, load_scope_config(scope_config_path))
    eligible_names = set(
        scope.loc[scope["angiosperm_analysis_eligible"], "accepted_species"].map(_text)
    )
    eligible_master = master.loc[master["accepted_species"].isin(eligible_names)].copy()
    master_by_name = {
        _key(record["accepted_species"]): record
        for record in eligible_master[["accepted_species", "genus", "family"]].to_dict("records")
    }

    value_map = {_key(key): _text(value) for key, value in config["value_map"].items()}
    excluded_names = {_text(name) for name in config["excluded_dio_only_names"]}
    usable = source.loc[source["mating_system"].map(_key).isin(value_map)].copy()
    exclusions = usable.loc[usable["genus_species"].map(_text).isin(excluded_names)].copy()
    usable = usable.loc[~usable["genus_species"].map(_text).isin(excluded_names)].copy()

    resolution_rows: list[dict[str, Any]] = []
    resolved: list[tuple[int, dict[str, Any], str, str]] = []
    unmatched_indices: list[int] = []
    for row_number, (index, row) in enumerate(usable.iterrows(), start=2):
        record = row.to_dict()
        source_name = _text(record["genus_species"])
        master_row = master_by_name.get(_key(source_name))
        if master_row is not None:
            accepted = _text(master_row["accepted_species"])
            resolved.append((row_number, record, accepted, "accepted_name_exact"))
            resolution_rows.append(
                {
                    "source_scientific_name": source_name,
                    "source_family": _text(record["Family"]),
                    "raw_value": _text(record["mating_system"]),
                    "accepted_species": accepted,
                    "name_match_method": "accepted_name_exact",
                    "match_status": "matched",
                    "reason": "accepted_name_exact",
                    "gbif_match_type": "",
                    "gbif_status": "",
                    "gbif_confidence": "",
                    "gbif_usage_key": "",
                    "gbif_accepted_usage_key": "",
                }
            )
        else:
            unmatched_indices.append(index)

    gbif_cache_path = output_dir / "gbif_species_match_cache.jsonl.gz"
    gbif_records: dict[str, dict[str, Any]] = {}
    if resolve_gbif and unmatched_indices:
        gbif_records = fetch_gbif_records(
            usable.loc[unmatched_indices],
            endpoint=_text(config["gbif_reconciliation"]["endpoint"]),
            cache_path=gbif_cache_path,
            workers=gbif_workers,
        )

    source_row_number = {index: number for number, index in enumerate(usable.index, start=2)}
    for index in unmatched_indices:
        record = usable.loc[index].to_dict()
        source_name = _text(record["genus_species"])
        cached = gbif_records.get(source_name, {})
        payload = cached.get("payload") if isinstance(cached, dict) else {}
        payload = payload if isinstance(payload, dict) else {}
        if resolve_gbif and not _text(cached.get("error")):
            accepted, reason = strict_gbif_resolution(
                source_name=source_name,
                source_family=_text(record["Family"]),
                payload=payload,
                master_by_name=master_by_name,
                reconciliation=config["gbif_reconciliation"],
            )
        elif resolve_gbif:
            accepted, reason = "", "gbif_request_error"
        else:
            accepted, reason = "", "gbif_resolution_disabled"
        if accepted:
            resolved.append(
                (
                    source_row_number[index],
                    record,
                    accepted,
                    "gbif_exact_synonym_to_accepted",
                )
            )
        resolution_rows.append(
            {
                "source_scientific_name": source_name,
                "source_family": _text(record["Family"]),
                "raw_value": _text(record["mating_system"]),
                "accepted_species": accepted,
                "name_match_method": "gbif_exact_synonym_to_accepted" if accepted else "unmatched",
                "match_status": "matched" if accepted else "unmatched",
                "reason": reason,
                "gbif_match_type": _text(payload.get("matchType")),
                "gbif_status": _text(payload.get("status")),
                "gbif_confidence": _text(payload.get("confidence")),
                "gbif_usage_key": _text(payload.get("usageKey")),
                "gbif_accepted_usage_key": _text(payload.get("acceptedUsageKey")),
            }
        )

    # A many-to-one synonym collapse is permitted only when every source row
    # agrees. Conflicting mappings are quarantined rather than voted.
    resolved_frame = pd.DataFrame(
        [
            {
                "row_number": row_number,
                "record": record,
                "accepted_species": accepted,
                "match_method": method,
                "standardized_value": value_map[_key(record["mating_system"])],
            }
            for row_number, record, accepted, method in resolved
        ]
    )
    conflict_species: set[str] = set()
    if not resolved_frame.empty:
        conflict_counts = resolved_frame.groupby("accepted_species")["standardized_value"].nunique()
        conflict_species = set(conflict_counts.loc[conflict_counts.gt(1)].index)
        resolved_frame = resolved_frame.loc[
            ~resolved_frame["accepted_species"].isin(conflict_species)
        ].copy()
        resolved_frame["match_priority"] = resolved_frame["match_method"].map(
            {"accepted_name_exact": 0, "gbif_exact_synonym_to_accepted": 1}
        )
        resolved_frame = (
            resolved_frame.sort_values(["accepted_species", "match_priority", "row_number"])
            .drop_duplicates("accepted_species", keep="first")
            .reset_index(drop=True)
        )

    direct_rows: list[dict[str, Any]] = []
    for row in resolved_frame.to_dict("records"):
        master_row = master_by_name[_key(row["accepted_species"])]
        direct_rows.append(
            _direct_candidate(
                row["record"],
                row_number=int(row["row_number"]),
                accepted_species=_text(row["accepted_species"]),
                genus=_text(master_row["genus"]),
                family=_text(master_row["family"]),
                match_method=_text(row["match_method"]),
                config=config,
            )
        )
    direct = pd.DataFrame(direct_rows, columns=CANDIDATE_COLUMNS).sort_values(
        ["accepted_species", "standardized_value"]
    )
    inference, inference_report = build_genus_inferences(eligible_master, direct, config)

    output_dir.mkdir(parents=True, exist_ok=True)
    direct_path = output_dir / "trait_candidates.csv.gz"
    inference_path = output_dir / "genus_inference_candidates.csv.gz"
    match_path = output_dir / "name_match_audit.csv.gz"
    unmatched_path = output_dir / "unmatched_taxa.csv"
    exclusion_path = output_dir / "source_exclusions.csv"
    _write_csv_gzip(direct, direct_path)
    _write_csv_gzip(inference, inference_path)
    match_audit = pd.DataFrame(resolution_rows).sort_values("source_scientific_name")
    _write_csv_gzip(match_audit, match_path)
    match_audit.loc[match_audit["match_status"].eq("unmatched")].to_csv(
        unmatched_path, index=False
    )
    exclusions.to_csv(exclusion_path, index=False)

    match_methods = {
        str(key): int(value) for key, value in direct["name_match_method"].value_counts().items()
    }
    value_counts = {
        str(key): int(value) for key, value in direct["standardized_value"].value_counts().items()
    }
    report = {
        "version": "1.0",
        "source_name": config["source"]["name"],
        "source_path": str(source_csv),
        "source_sha256": observed_hash,
        "source_license": config["source"]["license"],
        "source_url": config["source"]["dataset_url"],
        "n_source_rows": int(len(source)),
        "n_source_sc_si_rows": int(source["mating_system"].map(_key).isin(value_map).sum()),
        "n_dio_only_rows_excluded": int(len(exclusions)),
        "n_usable_source_rows": int(len(usable)),
        "n_direct_candidate_species": int(direct["accepted_species"].nunique()),
        "direct_candidates_by_value": value_counts,
        "direct_candidates_by_name_match_method": match_methods,
        "n_conflicting_master_mappings_quarantined": len(conflict_species),
        "conflicting_master_species": sorted(conflict_species),
        "n_unmatched_usable_source_rows": int(match_audit["match_status"].eq("unmatched").sum()),
        "n_angiosperm_master_species": int(len(eligible_master)),
        "genus_inference": inference_report,
        "outputs": {
            "trait_candidates": str(direct_path),
            "genus_inference_candidates": str(inference_path),
            "name_match_audit": str(match_path),
            "unmatched_taxa": str(unmatched_path),
            "source_exclusions": str(exclusion_path),
            "gbif_match_cache": str(gbif_cache_path) if gbif_cache_path.exists() else "",
        },
        "policy": (
            "SC/SI are species-direct only after exact accepted-name or strict GBIF exact-synonym "
            "reconciliation. Genus inference is a separate low-confidence layer requiring the "
            "configured support and purity; fuzzy and family-conflicting matches remain unmatched."
        ),
    }
    manifest_path = output_dir / "acquisition_manifest.json"
    manifest_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    report["manifest"] = str(manifest_path)
    return report


@app.command("prepare")
def prepare(
    source_csv: Path = typer.Option(..., exists=True),
    master_csv: Path = typer.Option(
        Path("data/v2/staging/gbif/collected/island_taxa.csv"), exists=True
    ),
    output_dir: Path = typer.Option(
        Path("data/v2/staging/traits/bulk/meyer_2026_reproductive_traits")
    ),
    config_path: Path = typer.Option(
        Path("config/meyer_2026_reproductive_traits.yml"), exists=True
    ),
    scope_config_path: Path = typer.Option(Path("config/angiosperm_scope.yml"), exists=True),
    resolve_gbif: bool = typer.Option(True, "--resolve-gbif/--no-resolve-gbif"),
    gbif_workers: int = typer.Option(8, min=1, max=32),
) -> None:
    """Acquire exact and strict-synonym SC/SI evidence plus genus completion."""
    report = prepare_meyer_2026(
        source_csv=source_csv,
        master_csv=master_csv,
        output_dir=output_dir,
        config_path=config_path,
        scope_config_path=scope_config_path,
        resolve_gbif=resolve_gbif,
        gbif_workers=gbif_workers,
    )
    typer.echo(
        f"Prepared {report['n_direct_candidate_species']:,} direct SC/SI species and "
        f"{report['genus_inference']['n_inferred_species']:,} labelled genus inferences; "
        f"{report['n_unmatched_usable_source_rows']:,} usable source rows remain unmatched."
    )


if __name__ == "__main__":
    app()
