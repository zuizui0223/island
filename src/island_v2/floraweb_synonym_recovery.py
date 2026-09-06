"""Recover FloraWeb/BIOLFLOR evidence hidden behind strict synonyms.

The first stage uses the frozen WFO 2026-06 Darwin Core archive only as a
local prefilter.  It never accepts a mapping: it merely identifies FloraWeb
names whose WFO synonym concept can reach a species in the fixed island
universe.  A later stage must confirm every candidate independently with the
live WFO exact matcher and the GBIF exact species matcher before evidence is
promoted.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import shutil
import subprocess
import zipfile
from collections import defaultdict
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.floraweb_trait_source import (
    INDEPENDENT_EVIDENCE_COLUMNS,
    STRICT_EVIDENCE_COLUMNS,
    _coalesce_states,
    build_source_package,
)

WFO_RELEASE = "2026-06"
WFO_MATCH_URL = "https://list.worldfloraonline.org/matching_rest.php"
GBIF_MATCH_URL = "https://api.gbif.org/v1/species/match"
SOURCE_GROUP = "floraweb_biolflor_synonym_recovery"
SOURCE_RUN_ID = "floraweb-biolflor-two-backbone-synonym-wave39-20260828"
SOURCE_ARTIFACT = "floraweb-biolflor-two-backbone-synonym-wave39"


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _tsv_rows(archive: zipfile.ZipFile, member: str) -> Iterable[dict[str, str]]:
    with (
        archive.open(member) as raw,
        io.TextIOWrapper(raw, encoding="utf-8", newline="") as handle,
    ):
        yield from csv.DictReader(handle, delimiter="\t")


def build_wfo_local_prefilter(
    *,
    wfo_zip: Path,
    floraweb_name_audit: pd.DataFrame,
    master: pd.DataFrame,
    universe_species: set[str],
) -> pd.DataFrame:
    """Return candidates for independent WFO+GBIF confirmation.

    The accepted concept is reconstructed through ``synonym.taxonID`` and
    ``taxon.nameID``.  This is important: WFO name identifiers cannot safely
    be treated as accepted concept identifiers, and doing so can produce
    biologically impossible cross-family mappings.
    """

    required_audit = {"source_name", "accepted"}
    required_master = {"accepted_species", "family"}
    if missing := required_audit.difference(floraweb_name_audit.columns):
        raise ValueError(f"FloraWeb name audit missing columns: {sorted(missing)}")
    if missing := required_master.difference(master.columns):
        raise ValueError(f"master missing columns: {sorted(missing)}")

    audit = floraweb_name_audit.fillna("").copy()
    accepted_flag = audit["accepted"].astype(str).str.casefold().eq("true")
    source_names = {
        _text(value)
        for value in audit.loc[~accepted_flag, "source_name"]
        if len(_text(value).split()) == 2
    }
    master_work = master.fillna("").copy()
    master_work = master_work.loc[master_work["accepted_species"].isin(universe_species)].copy()
    if master_work["accepted_species"].nunique() != len(universe_species):
        raise ValueError("master does not cover the fixed universe exactly")
    master_family = dict(zip(master_work["accepted_species"], master_work["family"], strict=True))
    master_names = set(master_family)

    source_to_name_ids: dict[str, set[str]] = defaultdict(set)
    name_id_to_source: dict[str, str] = {}
    source_name_id_rank: dict[str, str] = {}

    with zipfile.ZipFile(wfo_zip) as archive:
        for row in _tsv_rows(archive, "name.tsv"):
            scientific_name = _text(row.get("scientificName"))
            if scientific_name not in source_names:
                continue
            name_id = _text(row.get("ID"))
            rank = _text(row.get("rank")).casefold()
            if not name_id:
                continue
            source_to_name_ids[scientific_name].add(name_id)
            name_id_to_source[name_id] = scientific_name
            source_name_id_rank[name_id] = rank

        source_name_ids = set(name_id_to_source)
        name_id_to_taxon_ids: dict[str, set[str]] = defaultdict(set)
        synonym_target_taxon_ids: set[str] = set()
        for row in _tsv_rows(archive, "synonym.tsv"):
            name_id = _text(row.get("nameID"))
            if name_id not in source_name_ids:
                continue
            taxon_id = _text(row.get("taxonID"))
            if not taxon_id:
                continue
            name_id_to_taxon_ids[name_id].add(taxon_id)
            synonym_target_taxon_ids.add(taxon_id)

        taxon_to_accepted_name_id: dict[str, str] = {}
        for row in _tsv_rows(archive, "taxon.tsv"):
            taxon_id = _text(row.get("ID"))
            name_id = _text(row.get("nameID"))
            if taxon_id in synonym_target_taxon_ids or name_id in source_name_ids:
                taxon_to_accepted_name_id[taxon_id] = name_id
            if name_id in source_name_ids and taxon_id:
                name_id_to_taxon_ids[name_id].add(taxon_id)

        accepted_name_ids = {
            taxon_to_accepted_name_id[taxon_id]
            for taxon_ids in name_id_to_taxon_ids.values()
            for taxon_id in taxon_ids
            if taxon_id in taxon_to_accepted_name_id
        }
        accepted_name_by_id: dict[str, tuple[str, str]] = {}
        for row in _tsv_rows(archive, "name.tsv"):
            name_id = _text(row.get("ID"))
            if name_id not in accepted_name_ids:
                continue
            accepted_name_by_id[name_id] = (
                _text(row.get("scientificName")),
                _text(row.get("rank")).casefold(),
            )

    rows: list[dict[str, object]] = []
    for source_name in sorted(source_names):
        source_ids = source_to_name_ids.get(source_name, set())
        routes: set[tuple[str, str, str, str]] = set()
        for source_id in source_ids:
            if source_name_id_rank.get(source_id) != "species":
                continue
            for taxon_id in name_id_to_taxon_ids.get(source_id, set()):
                accepted_name_id = taxon_to_accepted_name_id.get(taxon_id, "")
                accepted_name, accepted_rank = accepted_name_by_id.get(accepted_name_id, ("", ""))
                if accepted_rank != "species" or accepted_name not in master_names:
                    continue
                routes.add((source_id, taxon_id, accepted_name_id, accepted_name))
        local_candidates = sorted({route[3] for route in routes})
        if not local_candidates:
            continue
        rows.append(
            {
                "source_name": source_name,
                "local_candidate_accepted_species": "|".join(local_candidates),
                "local_candidate_count": len(local_candidates),
                "local_mapping_ambiguous": len(local_candidates) != 1,
                "master_family": "|".join(
                    sorted({master_family[name] for name in local_candidates})
                ),
                "wfo_source_name_ids": "|".join(sorted(source_ids)),
                "wfo_route_count": len(routes),
                "wfo_routes_json": json.dumps(
                    [
                        {
                            "source_name_id": source_id,
                            "accepted_taxon_id": taxon_id,
                            "accepted_name_id": accepted_name_id,
                            "accepted_species": accepted_name,
                        }
                        for source_id, taxon_id, accepted_name_id, accepted_name in sorted(routes)
                    ],
                    ensure_ascii=False,
                    sort_keys=True,
                ),
                "prefilter_only": True,
                "acceptance_status": "pending_live_two_backbone_confirmation",
            }
        )
    return pd.DataFrame(rows)


def _curl_json(
    *,
    curl_path: str,
    endpoint: str,
    params: dict[str, object],
    timeout_seconds: int,
) -> tuple[int, dict[str, Any], str]:
    command = [
        curl_path,
        "--location",
        "--silent",
        "--show-error",
        "--retry",
        "2",
        "--retry-delay",
        "1",
        "--max-time",
        str(timeout_seconds),
        "--get",
        endpoint,
    ]
    for key, value in params.items():
        command.extend(["--data-urlencode", f"{key}={value}"])
    command.extend(["--write-out", "\n%{http_code}"])
    completed = subprocess.run(
        command,
        capture_output=True,
        check=False,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    body, separator, status_text = completed.stdout.rpartition("\n")
    status = int(status_text) if separator and status_text.isdigit() else 0
    error = " ".join(completed.stderr.strip().split())
    if completed.returncode != 0:
        return status, {}, error or f"curl_exit_{completed.returncode}"
    try:
        payload = json.loads(body)
    except json.JSONDecodeError as exc:
        return status, {}, f"json_decode_error:{exc.msg}"
    if not isinstance(payload, dict):
        return status, {}, "json_root_not_object"
    return status, payload, error


def _fetch_one(row: dict[str, str], *, curl_path: str, timeout_seconds: int) -> dict[str, Any]:
    source_name = row["source_name"]
    wfo_params: dict[str, object] = {
        "input_string": source_name,
        "fuzzy_names": 0,
        "fuzzy_authors": 0,
        "check_homonyms": "true",
        "check_rank": "true",
        "accept_single_candidate": "false",
    }
    gbif_params: dict[str, object] = {
        "name": source_name,
        "rank": "SPECIES",
        "kingdom": "Plantae",
        "strict": "true",
        "verbose": "true",
    }
    wfo_status, wfo, wfo_error = _curl_json(
        curl_path=curl_path,
        endpoint=WFO_MATCH_URL,
        params=wfo_params,
        timeout_seconds=timeout_seconds,
    )
    gbif_status, gbif, gbif_error = _curl_json(
        curl_path=curl_path,
        endpoint=GBIF_MATCH_URL,
        params=gbif_params,
        timeout_seconds=timeout_seconds,
    )
    return {
        **row,
        "retrieved_at_utc": datetime.now(UTC).isoformat(),
        "wfo_endpoint": WFO_MATCH_URL,
        "wfo_request_params": wfo_params,
        "wfo_status": wfo_status,
        "wfo": wfo,
        "wfo_error": wfo_error,
        "gbif_endpoint": GBIF_MATCH_URL,
        "gbif_request_params": gbif_params,
        "gbif_status": gbif_status,
        "gbif": gbif,
        "gbif_error": gbif_error,
    }


def fetch_two_backbone_responses(
    *,
    prefilter: pd.DataFrame,
    checkpoint_jsonl: Path,
    workers: int,
    timeout_seconds: int,
    curl_path: str,
    max_records: int | None = None,
) -> dict[str, int]:
    required = {"source_name", "local_candidate_accepted_species", "master_family"}
    if missing := required.difference(prefilter.columns):
        raise ValueError(f"prefilter missing columns: {sorted(missing)}")
    if shutil.which(curl_path) is None:
        raise ValueError(f"curl executable not found: {curl_path}")
    checkpoint_jsonl.parent.mkdir(parents=True, exist_ok=True)
    completed_names: set[str] = set()
    if checkpoint_jsonl.exists():
        with checkpoint_jsonl.open(encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    completed_names.add(_text(json.loads(line).get("source_name")))
    rows = [
        {key: _text(value) for key, value in row.items()}
        for row in prefilter.fillna("").to_dict(orient="records")
        if _text(row.get("source_name")) not in completed_names
    ]
    rows.sort(key=lambda row: row["source_name"])
    if max_records is not None:
        rows = rows[:max_records]
    successes = 0
    failures = 0
    with (
        checkpoint_jsonl.open("a", encoding="utf-8", newline="\n") as output,
        ThreadPoolExecutor(max_workers=workers) as executor,
    ):
        futures = {
            executor.submit(
                _fetch_one,
                row,
                curl_path=curl_path,
                timeout_seconds=timeout_seconds,
            ): row["source_name"]
            for row in rows
        }
        for future in as_completed(futures):
            result = future.result()
            output.write(json.dumps(result, ensure_ascii=False, sort_keys=True) + "\n")
            output.flush()
            if result["wfo_status"] == 200 and result["gbif_status"] == 200:
                successes += 1
            else:
                failures += 1
    return {
        "already_completed": len(completed_names),
        "attempted_this_run": len(rows),
        "both_http_200": successes,
        "http_or_decode_failures": failures,
        "checkpoint_records": len(completed_names) + len(rows),
    }


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    if path.suffix == ".gz":
        with gzip.open(path, mode="rt", encoding="utf-8") as handle:
            return [json.loads(line) for line in handle if line.strip()]
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _wfo_accepted_species(record: dict[str, Any]) -> str:
    match = (record.get("wfo") or {}).get("match") or {}
    placement = _text(match.get("placement"))
    if not placement.startswith("Code/Plantae/"):
        return ""
    accepted_placement = placement.split("$", 1)[0]
    parts = [part for part in accepted_placement.split("/") if part]
    return " ".join(parts[-2:]) if len(parts) >= 2 else ""


def audit_two_backbone_responses(
    records: Iterable[dict[str, Any]],
    master: pd.DataFrame,
    universe_species: set[str],
) -> pd.DataFrame:
    work = master.fillna("").copy()
    work = work.loc[work["accepted_species"].isin(universe_species)]
    master_family = dict(zip(work["accepted_species"], work["family"], strict=True))
    if len(master_family) != len(universe_species):
        raise ValueError("master does not cover the fixed universe exactly")
    rows: list[dict[str, object]] = []
    for record in records:
        source_name = _text(record.get("source_name"))
        wfo = record.get("wfo") or {}
        gbif = record.get("gbif") or {}
        match = wfo.get("match") or {}
        accepted = _wfo_accepted_species(record)
        reason = "accepted_strict_two_backbone"
        if record.get("wfo_status") != 200 or not match:
            reason = "wfo_no_unambiguous_exact_match"
        elif wfo.get("parsedName", {}).get("rank") != "species":
            reason = "wfo_not_species_rank"
        elif _text(wfo.get("parsedName", {}).get("canonical_form")) != source_name:
            reason = "wfo_canonical_input_mismatch"
        elif wfo.get("params", {}).get("fuzzyNameParts") != 0:
            reason = "wfo_fuzzy_match_rejected"
        elif not wfo.get("params", {}).get("checkHomonyms"):
            reason = "wfo_homonym_check_not_enabled"
        elif not wfo.get("params", {}).get("checkRank"):
            reason = "wfo_rank_check_not_enabled"
        elif record.get("gbif_status") != 200:
            reason = "gbif_no_response"
        elif gbif.get("matchType") != "EXACT":
            reason = "gbif_not_exact"
        elif gbif.get("rank") != "SPECIES" or gbif.get("kingdom") != "Plantae":
            reason = "gbif_not_plant_species"
        elif accepted != _text(gbif.get("species")):
            reason = "backbones_disagree"
        elif accepted not in master_family:
            reason = "agreed_name_not_in_fixed_master"
        elif accepted not in _text(record.get("local_candidate_accepted_species")).split("|"):
            reason = "live_mapping_not_in_frozen_wfo_prefilter"
        else:
            family = master_family[accepted]
            placement_parts = set(_text(match.get("placement")).split("/"))
            if family not in placement_parts or _text(gbif.get("family")) != family:
                reason = "family_conflict"
        rows.append(
            {
                "source_name": source_name,
                "accepted_species": accepted if reason == "accepted_strict_two_backbone" else "",
                "mapping_reason": reason,
                "master_family": master_family.get(accepted, ""),
                "wfo_match_id": _text(match.get("wfo_id")),
                "wfo_accepted_species": accepted,
                "wfo_family_in_placement": bool(
                    master_family.get(accepted, "")
                    and master_family[accepted] in set(_text(match.get("placement")).split("/"))
                ),
                "wfo_classification_version": _text(
                    wfo.get("params", {}).get("classificationVersion")
                ),
                "gbif_usage_key": _text(gbif.get("usageKey")),
                "gbif_accepted_usage_key": _text(
                    gbif.get("acceptedUsageKey") or gbif.get("usageKey")
                ),
                "gbif_accepted_species": _text(gbif.get("species")),
                "gbif_family": _text(gbif.get("family")),
                "retrieved_at_utc": _text(record.get("retrieved_at_utc")),
            }
        )
    result = pd.DataFrame(rows)
    if len(result):
        result.sort_values("source_name", inplace=True)
        result.reset_index(drop=True, inplace=True)
    return result


def _stable_id(*parts: object) -> str:
    payload = "\x1f".join(_text(part) for part in parts)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


def _name_usage_id(source_record_id: object) -> str:
    parts = _text(source_record_id).split(":")
    return parts[2] if len(parts) >= 3 and parts[:2] == ["floraweb", "name-use-id"] else ""


def _collapse_recovered_evidence(
    frame: pd.DataFrame,
    *,
    original_names_by_usage: dict[str, set[str]],
    mapping_by_name: dict[str, dict[str, str]],
    exact_provider_keys: set[tuple[str, str]],
    independent: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    kept: list[dict[str, object]] = []
    provenance: list[dict[str, object]] = []
    conflicts: list[dict[str, object]] = []
    for (species, trait), group in frame.groupby(["accepted_species", "trait_name"], sort=True):
        if (species, trait) in exact_provider_keys:
            provenance.append(
                {
                    "accepted_species": species,
                    "trait_name": trait,
                    "decision": "excluded_existing_exact_floraweb_species_trait",
                    "supporting_rows": len(group),
                }
            )
            continue
        values = sorted({_text(value) for value in group["normalized_value"] if _text(value)})
        if not values:
            continue
        if len(values) == 1:
            normalized_value = values[0]
        elif trait == "flower_primary_color":
            normalized_value = "|".join(
                sorted({token for value in values for token in value.split("|") if token})
            )
        elif independent:
            normalized_value = _coalesce_states(
                trait,
                {token for value in values for token in value.split("|") if token},
            )
        else:
            conflicts.append(
                {
                    "accepted_species": species,
                    "trait_name": trait,
                    "normalized_values": json.dumps(values, ensure_ascii=False),
                    "decision": "excluded_source_taxon_concept_conflict",
                    "supporting_rows": len(group),
                }
            )
            continue

        group = group.sort_values(["source_record_id", "source_url"])
        representative = group.iloc[0].to_dict()
        usage_ids = sorted({_name_usage_id(value) for value in group["source_record_id"]})
        usage_ids = [value for value in usage_ids if value]
        searched_names = sorted(
            {
                name
                for usage_id in usage_ids
                for name in original_names_by_usage.get(usage_id, set())
            }
        )
        resolution_lineages = sorted(
            {
                (
                    f"wfo:{mapping_by_name[name]['wfo_classification_version']}:"
                    f"{mapping_by_name[name]['wfo_match_id']};gbif:usage:"
                    f"{mapping_by_name[name]['gbif_accepted_usage_key']}"
                )
                for name in searched_names
                if name in mapping_by_name
            }
        )
        source_prefix = "rothmaler" if trait == "flower_primary_color" else "biolflor"
        representative.update(
            {
                "normalized_value": normalized_value,
                "evidence_scope": "synonym_direct",
                "name_match_method": "synonym_exact_two_backbone",
                "source_group": SOURCE_GROUP,
                "source_record_id": "floraweb-synonym:" + _stable_id(species, trait),
                "source_lineage": (
                    f"{source_prefix}:floraweb:accepted-species:{_stable_id(species)}:trait:{trait}"
                ),
                "lineage_method": "underlying_source_accepted_species_trait",
                "source_run_id": SOURCE_RUN_ID,
                "source_artifact": SOURCE_ARTIFACT,
                "acceptance_contract": ("floraweb_species_trait_strict_wfo_gbif_synonym_v1"),
                "inference_rule": "",
            }
        )
        kept.append(representative)
        provenance.append(
            {
                "accepted_species": species,
                "trait_name": trait,
                "normalized_value": normalized_value,
                "decision": "accepted_synonym_direct",
                "searched_names": "|".join(searched_names),
                "name_usage_ids": "|".join(usage_ids),
                "supporting_rows": len(group),
                "source_urls": "|".join(sorted(set(group["source_url"]))),
                "source_record_ids": "|".join(sorted(set(group["source_record_id"]))),
                "source_excerpts": " || ".join(
                    sorted({_text(value) for value in group["source_excerpt"] if _text(value)})
                ),
                "name_resolution_lineages": "|".join(resolution_lineages),
            }
        )

    columns = INDEPENDENT_EVIDENCE_COLUMNS if independent else STRICT_EVIDENCE_COLUMNS
    evidence = pd.DataFrame(kept, columns=columns)
    if len(evidence):
        evidence.sort_values(["accepted_species", "trait_name"], inplace=True)
        evidence.reset_index(drop=True, inplace=True)
    return evidence, pd.DataFrame(provenance), pd.DataFrame(conflicts)


def build_recovered_source_package(
    *,
    raw: pd.DataFrame,
    mapping_audit: pd.DataFrame,
    master: pd.DataFrame,
    universe_species: set[str],
    exact_strict: pd.DataFrame,
    exact_independent: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    accepted = mapping_audit.loc[
        mapping_audit["mapping_reason"].eq("accepted_strict_two_backbone")
    ].copy()
    if accepted["source_name"].duplicated().any():
        raise ValueError("mapping audit contains duplicate accepted source names")
    mapping = dict(zip(accepted["source_name"], accepted["accepted_species"], strict=True))
    fixed_master = master.loc[master["accepted_species"].isin(universe_species)].copy()
    if fixed_master["accepted_species"].nunique() != len(universe_species):
        raise ValueError("master does not cover the fixed universe exactly")

    recovered_raw = raw.loc[raw["canonical_name"].isin(mapping)].copy()
    recovered_raw["original_source_name"] = recovered_raw["canonical_name"]
    recovered_raw["canonical_name"] = recovered_raw["canonical_name"].map(mapping)
    strict, independent, extraction, _ = build_source_package(recovered_raw, fixed_master)
    original_names_by_usage = {
        str(usage_id): set(group["original_source_name"])
        for usage_id, group in recovered_raw.groupby("name_usage_id")
    }
    mapping_by_name = {
        str(row.source_name): {column: _text(getattr(row, column)) for column in accepted.columns}
        for row in accepted.itertuples(index=False)
    }
    exact_strict_keys = set(
        exact_strict[["accepted_species", "trait_name"]].itertuples(index=False, name=None)
    )
    exact_independent_keys = set(
        exact_independent[["accepted_species", "trait_name"]].itertuples(index=False, name=None)
    )
    strict_out, strict_provenance, strict_conflicts = _collapse_recovered_evidence(
        strict,
        original_names_by_usage=original_names_by_usage,
        mapping_by_name=mapping_by_name,
        exact_provider_keys=exact_strict_keys,
        independent=False,
    )
    independent_out, independent_provenance, independent_conflicts = _collapse_recovered_evidence(
        independent,
        original_names_by_usage=original_names_by_usage,
        mapping_by_name=mapping_by_name,
        exact_provider_keys=exact_independent_keys,
        independent=True,
    )
    provenance = pd.concat(
        [
            strict_provenance.assign(strict_three_axis_included=True),
            independent_provenance.assign(strict_three_axis_included=False),
        ],
        ignore_index=True,
    )
    conflicts = pd.concat(
        [
            strict_conflicts.assign(strict_three_axis_included=True),
            independent_conflicts.assign(strict_three_axis_included=False),
        ],
        ignore_index=True,
    )
    return {
        "strict": strict_out,
        "independent": independent_out,
        "provenance": provenance,
        "conflicts": conflicts,
        "extraction": extraction,
        "recovered_raw": recovered_raw,
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})


def _write_sorted_jsonl_gzip(records: list[dict[str, Any]], path: Path) -> None:
    with (
        path.open("wb") as raw,
        gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as compressed,
        io.TextIOWrapper(compressed, encoding="utf-8", newline="\n") as text,
    ):
        for record in sorted(records, key=lambda value: _text(value.get("source_name"))):
            text.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")


def _prefilter(args: argparse.Namespace) -> None:
    audit = pd.read_csv(args.floraweb_name_audit, dtype=str).fillna("")
    master = pd.read_csv(args.master_csv, dtype=str).fillna("")
    universe = set(
        pd.read_csv(
            args.universe_coverage,
            usecols=["accepted_species"],
            dtype=str,
        )["accepted_species"]
    )
    result = build_wfo_local_prefilter(
        wfo_zip=args.wfo_zip,
        floraweb_name_audit=audit,
        master=master,
        universe_species=universe,
    )
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_csv, index=False)
    print(
        json.dumps(
            {
                "contract": "floraweb_wfo_2026_06_local_prefilter_v1",
                "fixed_universe_species": len(universe),
                "pending_live_confirmation": len(result),
                "ambiguous_local_candidates": int(
                    result.get("local_mapping_ambiguous", pd.Series(dtype=bool)).sum()
                ),
            },
            indent=2,
            sort_keys=True,
        )
    )


def _fetch(args: argparse.Namespace) -> None:
    prefilter = pd.read_csv(args.prefilter_csv, dtype=str).fillna("")
    summary = fetch_two_backbone_responses(
        prefilter=prefilter,
        checkpoint_jsonl=args.checkpoint_jsonl,
        workers=args.workers,
        timeout_seconds=args.timeout_seconds,
        curl_path=args.curl_path,
        max_records=args.max_records,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


def _audit(args: argparse.Namespace) -> None:
    records = read_jsonl(args.responses_jsonl)
    master = pd.read_csv(args.master_csv, dtype=str).fillna("")
    universe = set(
        pd.read_csv(
            args.universe_coverage,
            usecols=["accepted_species"],
            dtype=str,
        )["accepted_species"]
    )
    audit = audit_two_backbone_responses(records, master, universe)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(args.output_csv, index=False)
    print(
        json.dumps(
            {
                "contract": "floraweb_wfo_gbif_two_backbone_mapping_audit_v1",
                "records": len(audit),
                "accepted": int(audit["mapping_reason"].eq("accepted_strict_two_backbone").sum()),
                "by_reason": {
                    str(key): int(value)
                    for key, value in audit["mapping_reason"].value_counts().items()
                },
            },
            indent=2,
            sort_keys=True,
        )
    )


def _build(args: argparse.Namespace) -> None:
    raw = pd.read_csv(args.raw_csv, dtype=str).fillna("")
    prefilter = pd.read_csv(args.prefilter_csv, dtype=str).fillna("")
    mapping_audit = pd.read_csv(args.mapping_audit_csv, dtype=str).fillna("")
    master = pd.read_csv(args.master_csv, dtype=str).fillna("")
    exact_strict = pd.read_csv(args.exact_strict, dtype=str).fillna("")
    exact_independent = pd.read_csv(args.exact_independent, dtype=str).fillna("")
    universe = set(
        pd.read_csv(
            args.universe_coverage,
            usecols=["accepted_species"],
            dtype=str,
        )["accepted_species"]
    )
    response_records = read_jsonl(args.responses_jsonl)
    if len(prefilter) != len(response_records) or len(mapping_audit) != len(prefilter):
        raise ValueError("prefilter, response, and mapping-audit row counts differ")
    package = build_recovered_source_package(
        raw=raw,
        mapping_audit=mapping_audit,
        master=master,
        universe_species=universe,
        exact_strict=exact_strict,
        exact_independent=exact_independent,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frames = {
        "floraweb_synonym_reviewed_direct_evidence.csv.gz": package["strict"],
        "floraweb_synonym_independent_trait_evidence.csv.gz": package["independent"],
        "floraweb_synonym_provenance_audit.csv.gz": package["provenance"],
        "floraweb_synonym_conflict_audit.csv.gz": package["conflicts"],
        "floraweb_synonym_extraction_audit.csv.gz": package["extraction"],
        "floraweb_two_backbone_mapping_audit.csv.gz": mapping_audit,
        "floraweb_wfo_local_prefilter.csv.gz": prefilter,
    }
    for name, frame in frames.items():
        _write_gzip_csv(frame, args.output_dir / name)
    response_name = "floraweb_two_backbone_responses.jsonl.gz"
    _write_sorted_jsonl_gzip(response_records, args.output_dir / response_name)
    strict = package["strict"]
    independent = package["independent"]
    summary: dict[str, object] = {
        "contract": "floraweb_biolflor_two_backbone_synonym_source_package_v1",
        "source_run_id": SOURCE_RUN_ID,
        "fixed_universe_species": len(universe),
        "fixed_species_axis_denominator": len(universe) * 3,
        "wfo_release": WFO_RELEASE,
        "web_requests": {
            "wfo": len(response_records),
            "gbif": len(response_records),
            "query_cost_usd": 0,
        },
        "mapping": {
            "local_prefilter_candidates": len(prefilter),
            "strict_two_backbone_accepted_names": int(
                mapping_audit["mapping_reason"].eq("accepted_strict_two_backbone").sum()
            ),
            "strict_two_backbone_accepted_species": int(
                mapping_audit.loc[
                    mapping_audit["mapping_reason"].eq("accepted_strict_two_backbone"),
                    "accepted_species",
                ].nunique()
            ),
            "by_reason": {
                str(key): int(value)
                for key, value in mapping_audit["mapping_reason"].value_counts().items()
            },
        },
        "strict_direct": {
            "rows": len(strict),
            "species": int(strict["accepted_species"].nunique()),
            "species_trait": int(
                strict[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
            ),
            "by_trait": {
                str(key): int(value)
                for key, value in strict["trait_name"].value_counts().sort_index().items()
            },
            "by_quality": {
                str(key): int(value)
                for key, value in strict["quality"].value_counts().sort_index().items()
            },
        },
        "independent_traits": {
            "rows": len(independent),
            "species_trait": int(
                independent[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
            ),
            "by_trait": {
                str(key): int(value)
                for key, value in independent["trait_name"].value_counts().sort_index().items()
            },
            "strict_three_axis_included": False,
        },
        "checks": {
            "wfo_fuzzy_matches_accepted": 0,
            "gbif_non_exact_matches_accepted": 0,
            "backbone_disagreements_accepted": 0,
            "family_conflicts_accepted": 0,
            "existing_exact_floraweb_species_trait_recounted": 0,
            "family_inference_used": False,
            "global_fallback_used": False,
            "pollen_vector_in_strict_axis": False,
            "reward_type_in_strict_axis": False,
        },
    }
    artifact_hashes = {name: _sha256(args.output_dir / name) for name in [*frames, response_name]}
    summary["artifact_hashes"] = artifact_hashes
    summary_path = args.output_dir / "floraweb_synonym_source_summary.json"
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    prefilter = subparsers.add_parser("prefilter")
    prefilter.add_argument("--wfo-zip", required=True, type=Path)
    prefilter.add_argument("--floraweb-name-audit", required=True, type=Path)
    prefilter.add_argument("--master-csv", required=True, type=Path)
    prefilter.add_argument("--universe-coverage", required=True, type=Path)
    prefilter.add_argument("--output-csv", required=True, type=Path)
    prefilter.set_defaults(handler=_prefilter)
    fetch = subparsers.add_parser("fetch")
    fetch.add_argument("--prefilter-csv", required=True, type=Path)
    fetch.add_argument("--checkpoint-jsonl", required=True, type=Path)
    fetch.add_argument("--workers", type=int, default=4)
    fetch.add_argument("--timeout-seconds", type=int, default=45)
    fetch.add_argument("--curl-path", default="curl.exe" if shutil.which("curl.exe") else "curl")
    fetch.add_argument("--max-records", type=int)
    fetch.set_defaults(handler=_fetch)
    audit = subparsers.add_parser("audit")
    audit.add_argument("--responses-jsonl", required=True, type=Path)
    audit.add_argument("--master-csv", required=True, type=Path)
    audit.add_argument("--universe-coverage", required=True, type=Path)
    audit.add_argument("--output-csv", required=True, type=Path)
    audit.set_defaults(handler=_audit)
    build = subparsers.add_parser("build")
    build.add_argument("--raw-csv", required=True, type=Path)
    build.add_argument("--prefilter-csv", required=True, type=Path)
    build.add_argument("--responses-jsonl", required=True, type=Path)
    build.add_argument("--mapping-audit-csv", required=True, type=Path)
    build.add_argument("--master-csv", required=True, type=Path)
    build.add_argument("--universe-coverage", required=True, type=Path)
    build.add_argument("--exact-strict", required=True, type=Path)
    build.add_argument("--exact-independent", required=True, type=Path)
    build.add_argument("--output-dir", required=True, type=Path)
    build.set_defaults(handler=_build)
    args = parser.parse_args()
    args.handler(args)


if __name__ == "__main__":
    main()
