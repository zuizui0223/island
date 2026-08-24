"""Acquire source-backed floral and breeding-system traits from Pladias.

Pladias is treated as one reusable source package, not as a species-by-species
search backend.  The importer joins the official Czech-flora taxon list to the
island master, fetches each exact species page once, caches an evidence record,
and emits only ontology-safe species-direct candidates.  Pollen-vector and
reward fields are deliberately not projected onto the strict floral axes.
"""

from __future__ import annotations

import hashlib
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Annotated, Any
from urllib.parse import quote

import pandas as pd
import requests
import typer
import yaml
from bs4 import BeautifulSoup

from island_v2.meyer_2026_reproductive_traits import (
    CANDIDATE_COLUMNS as BASE_CANDIDATE_COLUMNS,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

CANDIDATE_COLUMNS = [*BASE_CANDIDATE_COLUMNS, "source_lineage"]
FETCH_COLUMNS = [
    "accepted_species",
    "family",
    "source_url",
    "http_status",
    "matched_page_name",
    "page_family",
    "identity_status",
    "retrieved_at_utc",
    "content_sha256",
    "elapsed_seconds",
    "cache_status",
    "error",
]
EXCLUSION_COLUMNS = [
    "accepted_species",
    "family",
    "trait_name",
    "raw_value",
    "source_url",
    "reason",
]
QUALITY_AUDIT_COLUMNS = [
    "review_id",
    "sample_seed",
    "accepted_species",
    "family",
    "source_url",
    "trait_name",
    "raw_value",
    "standardized_value",
    "evidence_excerpt",
    "source_lineage",
    "matched_page_name",
    "content_sha256",
    "identity_correct",
    "value_mapping_correct",
    "quote_exact",
    "provenance_complete",
    "cultivar_contamination",
    "accepted_correct",
    "decision_reason",
    "reviewer",
    "reviewed_at_utc",
]

STRICT_SPECIES_RE = re.compile(r"^[A-Z][A-Za-z-]+ [a-z][a-z-]+$")
UTC_TZ = timezone.utc  # noqa: UP017 -- keep Python 3.10 local audit compatibility


@app.callback()
def main() -> None:
    """Build the exact-name Pladias trait source package."""


def _text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_config(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise typer.BadParameter(f"Expected a mapping in {path}")
    return value


def _write_csv_gzip(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})


def _stable_id(*parts: object) -> str:
    payload = "\x1f".join(_text(part) for part in parts)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


@app.command("download")
def download_taxon_list(
    output_xls: Annotated[Path, typer.Option("--output-xls")],
    config_path: Annotated[Path, typer.Option("--config-path")] = Path(
        "config/pladias_traits.yml"
    ),
    timeout_seconds: Annotated[float, typer.Option("--timeout-seconds")] = 90.0,
) -> None:
    """Download and hash-verify the official Pladias taxon list."""
    config = _load_config(config_path)
    source = config["source"]
    response = requests.get(
        source["taxon_list_url"],
        timeout=timeout_seconds,
        headers={"User-Agent": source["user_agent"]},
    )
    response.raise_for_status()
    observed = _sha256_bytes(response.content)
    expected = _text(source["expected_taxon_list_sha256"]).casefold()
    if observed != expected:
        raise typer.BadParameter(
            f"Pladias taxon-list SHA-256 mismatch: expected {expected}, got {observed}"
        )
    output_xls.parent.mkdir(parents=True, exist_ok=True)
    output_xls.write_bytes(response.content)
    typer.echo(f"Downloaded verified Pladias taxon list to {output_xls}")


def _feature_rows(html: str) -> list[dict[str, str]]:
    soup = BeautifulSoup(html, "html.parser")
    rows: list[dict[str, str]] = []
    for item in soup.select("li.featureDetail"):
        direct = item.find("div", recursive=False)
        if direct is None:
            continue
        label_node = direct.select_one("span.font-italic")
        value_node = direct.find("b")
        if label_node is None or value_node is None:
            continue
        trigger = direct.select_one("[data-target]")
        feature_id = ""
        if trigger:
            match = re.search(r"modal_feature_(\d+)", _text(trigger.get("data-target")))
            feature_id = match.group(1) if match else ""
        modal = item.select_one(".modal-body")
        method_text = _text(modal.get_text(" ")) if modal else ""
        citation = method_text
        marker = "Data source and citation"
        if marker in method_text:
            citation = _text(method_text.split(marker, 1)[1])
        rows.append(
            {
                "label": _text(label_node.get_text(" ")),
                "raw_value": _text(value_node.get_text(" ")),
                "feature_id": feature_id,
                "feature_citation": citation,
                "supporting_excerpt": _text(
                    f"{label_node.get_text(' ')}: {value_node.get_text(' ')}"
                ),
            }
        )
    return rows


def _parse_page(html: str) -> dict[str, Any]:
    soup = BeautifulSoup(html, "html.parser")
    heading = soup.select_one("div.container.content h1")
    breadcrumbs = [
        _text(node.get_text(" ")) for node in soup.select("ol.breadcrumb a")
    ]
    title = soup.title.get_text(" ", strip=True) if soup.title else ""
    return {
        "matched_page_name": _text(heading.get_text(" ")) if heading else "",
        "page_title": _text(title),
        "breadcrumbs": breadcrumbs,
        "features": _feature_rows(html),
    }


def _split_values(value: str) -> list[str]:
    return [_text(token).casefold() for token in value.split(",") if _text(token)]


def _json_states(states: set[str]) -> str:
    values = sorted(states)
    if len(values) == 1:
        return values[0]
    return json.dumps(values, ensure_ascii=False, separators=(",", ":"))


def _normalise_colour(raw: str) -> str:
    mapping = {
        "white": {"white"},
        "yellow-white": {"white", "yellow_orange"},
        "green-white": {"white", "green_brown_inconspicuous"},
        "green": {"green_brown_inconspicuous"},
        "yellow-green": {"yellow_orange", "green_brown_inconspicuous"},
        "yellow": {"yellow_orange"},
        "orange": {"yellow_orange"},
        "pink": {"red_pink"},
        "pink-violet": {"red_pink", "blue_purple"},
        "red-violet": {"red_pink", "blue_purple"},
        "red": {"red_pink"},
        "violet": {"blue_purple"},
        "blue": {"blue_purple"},
        "blue-violet": {"blue_purple"},
        "brown": {"green_brown_inconspicuous"},
        "red-brown": {"red_pink", "green_brown_inconspicuous"},
        "black": {"other_described"},
    }
    states: set[str] = set()
    for token in _split_values(raw):
        states.update(mapping.get(token, set()))
    return _json_states(states) if states else ""


def _normalise_symmetry(raw: str) -> str:
    mapping = {"actinomorphic": "actinomorphic", "zygomorphic": "zygomorphic"}
    states = {mapping[token] for token in _split_values(raw) if token in mapping}
    return _json_states(states) if states else ""


def _normalise_inflorescence(raw: str) -> str:
    states: set[str] = set()
    for token in _split_values(raw):
        token_states: set[str] = set()
        if token == "flores solitarii":
            token_states.add("solitary")
        if any(key in token for key in ("anthodi", "capitulum")):
            token_states.add("composite_display")
        if any(key in token for key in ("umbella", "corymb")):
            token_states.add("umbel_corymb")
        if any(
            key in token
            for key in (
                "panicula",
                "racemus",
                "spica",
                "spicula",
                "pseudospica",
                "amentum",
                "spadix",
            )
        ):
            token_states.add("raceme_spike_panicle")
        if "fasciculus" in token:
            token_states.add("few_flowered")
        if not token_states:
            return ""
        states.update(token_states)
    return _json_states(states) if states else ""


def _normalise_floral_form(raw: str) -> str:
    mapping = {
        "tubular": "tubular",
        "filiform": "tubular",
        "funnel-shaped": "funnel_trumpet",
        "rotate": "open_radial",
        "campanulate": "bell_campanulate",
        "bilabiate": "bilabiate",
        "personate": "bilabiate",
        "hypocrateriform": "salverform",
        "urceolate": "urn_urceolate",
        "ligulate": "other_described",
        "special type": "other_described",
    }
    states = {mapping[token] for token in _split_values(raw) if token in mapping}
    return _json_states(states) if states else ""


def _normalise_pollination(raw: str) -> str:
    values = set(_split_values(raw))
    if "water-pollination" in values:
        return ""
    wind = "wind-pollination" in values
    biotic = "insect-pollination" in values
    if wind and biotic:
        return "mixed"
    if wind:
        return "abiotic_wind"
    if biotic:
        return "biotic"
    return ""


def _normalise_cleistogamy(raw: str) -> str:
    values = set(_split_values(raw))
    if "cleistogamy" not in values:
        return ""
    other_modes = values.difference({"cleistogamy", "selfing"})
    return "facultative" if other_modes else "obligate"


def _normalise_reproduction(raw: str) -> list[tuple[str, str]]:
    values = set(_split_values(raw))
    output: list[tuple[str, str]] = []
    if "allogamy self-incompatibility" in values:
        output.append(("self_incompatibility", "SI"))

    # Apomixis and sterility can coexist with a sexual pathway, but this
    # categorical source does not quantify which pathway predominates.  Keep
    # explicit SI above, then fail closed instead of manufacturing a mating
    # system value.
    if any("apomixis" in value or value == "sterility" for value in values):
        return output

    allogamy = bool(
        values
        & {"allogamy", "allogamy self-incompatibility", "facultative allogamy"}
    )
    autogamy = bool(values & {"autogamy", "facultative autogamy"})
    if "mixed mating" in values or (allogamy and autogamy):
        output.append(("mating_system", "mixed_mating"))
    elif "facultative autogamy" in values:
        output.append(("mating_system", "predominantly_selfing"))
    elif "autogamy" in values:
        output.append(("mating_system", "obligate_selfing"))
    elif allogamy:
        output.append(("mating_system", "predominantly_outcrossing"))
    return output


def _normalise_feature(label: str, raw: str) -> list[tuple[str, str]]:
    if label == "Flower colour":
        value = _normalise_colour(raw)
        return [("flower_primary_color", value)] if value else []
    if label == "Flower symmetry":
        value = _normalise_symmetry(raw)
        return [("floral_symmetry", value)] if value else []
    if label == "Inflorescence type":
        value = _normalise_inflorescence(raw)
        return [("inflorescence_display", value)] if value else []
    if label == "Shape of the sympetalous corolla or syntepalous perianth":
        value = _normalise_floral_form(raw)
        return [("floral_form", value)] if value else []
    if label == "Generative reproduction type":
        return _normalise_reproduction(raw)
    if label == "Pollination syndrome":
        output: list[tuple[str, str]] = []
        pollen_value = _normalise_pollination(raw)
        if pollen_value:
            output.append(("pollen_vector_mode", pollen_value))
        cleistogamy_value = _normalise_cleistogamy(raw)
        if cleistogamy_value:
            output.append(("cleistogamy", cleistogamy_value))
        return output
    return []


def _source_feature_label(trait_name: str) -> str:
    return {
        "flower_primary_color": "Flower colour",
        "floral_symmetry": "Flower symmetry",
        "inflorescence_display": "Inflorescence type",
        "floral_form": "Shape of the sympetalous corolla or syntepalous perianth",
        "self_incompatibility": "Generative reproduction type",
        "mating_system": "Generative reproduction type",
        "pollen_vector_mode": "Pollination syndrome",
        "cleistogamy": "Pollination syndrome",
    }.get(trait_name, "")


def _quality_audit(
    candidates: pd.DataFrame,
    fetches: pd.DataFrame,
    *,
    sample_size: int = 200,
    seed: int = 20260808,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build a fixed, trait-stratified audit of extraction contract precision.

    This is an objective replay audit of identity, quotation, mapping and
    provenance.  It does not claim independent biological re-measurement of
    the curated Pladias values.
    """
    if candidates.empty:
        empty = pd.DataFrame(columns=QUALITY_AUDIT_COLUMNS)
        return empty, {
            "reviewed": 0,
            "accepted_correct": 0,
            "precision": None,
            "cultivar_contamination": 0,
            "cultivar_contamination_rate": None,
            "precision_gate_met": False,
            "cultivar_gate_met": False,
        }

    queues: dict[str, list[dict[str, str]]] = {}
    for trait_name, group in candidates.groupby("trait_name", sort=True):
        ordered = group.assign(
            _audit_order=group.apply(
                lambda row: _stable_id(
                    seed,
                    row["accepted_species"],
                    row["trait_name"],
                    row["evidence_id"],
                ),
                axis=1,
            )
        ).sort_values("_audit_order")
        queues[str(trait_name)] = ordered.drop(columns="_audit_order").to_dict(
            orient="records"
        )

    selected: list[dict[str, str]] = []
    used_species: set[str] = set()
    target = min(sample_size, int(candidates["accepted_species"].nunique()))
    while len(selected) < target:
        progressed = False
        for trait_name in sorted(queues):
            queue = queues[trait_name]
            while queue and queue[0]["accepted_species"] in used_species:
                queue.pop(0)
            if not queue:
                continue
            row = queue.pop(0)
            selected.append(row)
            used_species.add(row["accepted_species"])
            progressed = True
            if len(selected) >= target:
                break
        if not progressed:
            break

    fetch_by_species = fetches.drop_duplicates("accepted_species").set_index(
        "accepted_species"
    )
    reviewed_at = datetime.now(UTC_TZ).isoformat()
    audit_rows: list[dict[str, Any]] = []
    cultivar_pattern = re.compile(
        r"(?:\bcultivar\b|\bcv\.?\b|\bhybrid\b|×)", re.IGNORECASE
    )
    for position, candidate in enumerate(selected, start=1):
        species = candidate["accepted_species"]
        fetch = fetch_by_species.loc[species]
        trait_name = candidate["trait_name"]
        raw = candidate["raw_value"]
        value = candidate["standardized_value"]
        label = _source_feature_label(trait_name)
        mapped = set(_normalise_feature(label, raw))
        identity_correct = (
            fetch["identity_status"] == "accepted_name_and_family_exact"
            and fetch["matched_page_name"] == species
        )
        mapping_correct = (trait_name, value) in mapped
        quote_exact = candidate["evidence_excerpt"] == f"{label}: {raw}"
        provenance_complete = all(
            _text(candidate.get(column))
            for column in (
                "source_url",
                "source_citation",
                "source_lineage",
                "evidence_excerpt",
            )
        ) and bool(_text(fetch["content_sha256"]))
        cultivar_contamination = bool(
            cultivar_pattern.search(
                " ".join(
                    [
                        species,
                        _text(fetch.get("matched_page_name")),
                        _text(candidate.get("evidence_excerpt")),
                    ]
                )
            )
        )
        accepted_correct = (
            identity_correct
            and mapping_correct
            and quote_exact
            and provenance_complete
            and not cultivar_contamination
        )
        failed = [
            name
            for name, passed in (
                ("identity", identity_correct),
                ("mapping", mapping_correct),
                ("quote", quote_exact),
                ("provenance", provenance_complete),
                ("cultivar", not cultivar_contamination),
            )
            if not passed
        ]
        audit_rows.append(
            {
                "review_id": f"pladias-audit-{position:04d}",
                "sample_seed": seed,
                "accepted_species": species,
                "family": candidate["family"],
                "source_url": candidate["source_url"],
                "trait_name": trait_name,
                "raw_value": raw,
                "standardized_value": value,
                "evidence_excerpt": candidate["evidence_excerpt"],
                "source_lineage": candidate["source_lineage"],
                "matched_page_name": fetch["matched_page_name"],
                "content_sha256": fetch["content_sha256"],
                "identity_correct": identity_correct,
                "value_mapping_correct": mapping_correct,
                "quote_exact": quote_exact,
                "provenance_complete": provenance_complete,
                "cultivar_contamination": cultivar_contamination,
                "accepted_correct": accepted_correct,
                "decision_reason": (
                    "accepted_source_contract" if accepted_correct else "failed:" + ",".join(failed)
                ),
                "reviewer": "deterministic_source_contract_v1",
                "reviewed_at_utc": reviewed_at,
            }
        )

    audit = pd.DataFrame(audit_rows, columns=QUALITY_AUDIT_COLUMNS)
    reviewed = len(audit)
    accepted = int(audit["accepted_correct"].sum()) if reviewed else 0
    cultivar = int(audit["cultivar_contamination"].sum()) if reviewed else 0
    precision = accepted / reviewed if reviewed else None
    cultivar_rate = cultivar / reviewed if reviewed else None
    report = {
        "audit_kind": "deterministic_extraction_contract_not_biological_remeasurement",
        "sample_seed": seed,
        "sampling": "trait_round_robin_without_species_replacement",
        "reviewed": reviewed,
        "accepted_correct": accepted,
        "precision": precision,
        "cultivar_contamination": cultivar,
        "cultivar_contamination_rate": cultivar_rate,
        "precision_gate": 0.95,
        "cultivar_contamination_gate": 0.02,
        "precision_gate_met": bool(precision is not None and precision >= 0.95),
        "cultivar_gate_met": bool(
            cultivar_rate is not None and cultivar_rate <= 0.02
        ),
    }
    return audit, report


def _fetch_one(
    species: str,
    family: str,
    source: dict[str, Any],
    cache_dir: Path,
    timeout_seconds: float,
    retries: int,
) -> dict[str, Any]:
    cache_path = cache_dir / f"{_stable_id(species)}.json"
    if cache_path.exists():
        cached = json.loads(cache_path.read_text(encoding="utf-8"))
        accepted_labels = set(source["accepted_feature_labels"])
        cached["features"] = [
            feature
            for feature in cached.get("features", [])
            if feature.get("label") in accepted_labels
        ]
        cached["cache_status"] = "reused"
        cache_path.write_text(
            json.dumps(cached, ensure_ascii=False, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return cached

    url = source["species_page_template"].format(species=quote(species, safe=""))
    started = time.monotonic()
    last_error = ""
    for attempt in range(retries + 1):
        try:
            response = requests.get(
                url,
                timeout=timeout_seconds,
                headers={"User-Agent": source["user_agent"]},
            )
            response.raise_for_status()
            parsed = _parse_page(response.text)
            accepted_labels = set(source["accepted_feature_labels"])
            parsed["features"] = [
                feature
                for feature in parsed["features"]
                if feature.get("label") in accepted_labels
            ]
            breadcrumbs = parsed["breadcrumbs"]
            exact_name = parsed["matched_page_name"] == species
            family_match = family in breadcrumbs
            identity_status = (
                "accepted_name_and_family_exact"
                if exact_name and family_match
                else "reject_page_name_mismatch"
                if not exact_name
                else "reject_family_conflict_or_missing"
            )
            record = {
                "accepted_species": species,
                "family": family,
                "source_url": url,
                "http_status": response.status_code,
                "matched_page_name": parsed["matched_page_name"],
                "page_family": family if family_match else "",
                "identity_status": identity_status,
                "page_title": parsed["page_title"],
                "breadcrumbs": breadcrumbs,
                "features": parsed["features"],
                "retrieved_at_utc": datetime.now(UTC_TZ).isoformat(),
                "content_sha256": _sha256_bytes(response.content),
                "elapsed_seconds": round(time.monotonic() - started, 3),
                "cache_status": "fetched",
                "error": "",
            }
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(
                json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            return record
        except (requests.RequestException, ValueError) as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            if attempt < retries:
                time.sleep(min(8.0, 1.5 * (2**attempt)))
    return {
        "accepted_species": species,
        "family": family,
        "source_url": url,
        "http_status": 0,
        "matched_page_name": "",
        "page_family": "",
        "identity_status": "fetch_failed",
        "page_title": "",
        "breadcrumbs": [],
        "features": [],
        "retrieved_at_utc": datetime.now(UTC_TZ).isoformat(),
        "content_sha256": "",
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "cache_status": "failed",
        "error": last_error,
    }


def _candidate_rows(
    page: dict[str, Any], source: dict[str, Any]
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    species = page["accepted_species"]
    family = page["family"]
    url = page["source_url"]
    if page["identity_status"] != "accepted_name_and_family_exact":
        return [], [
            {
                "accepted_species": species,
                "family": family,
                "trait_name": "",
                "raw_value": "",
                "source_url": url,
                "reason": page["identity_status"],
            }
        ]

    candidates: list[dict[str, str]] = []
    exclusions: list[dict[str, str]] = []
    accepted_labels = set(source["accepted_feature_labels"])
    for feature in page["features"]:
        label = _text(feature.get("label"))
        if label not in accepted_labels:
            continue
        raw = _text(feature.get("raw_value"))
        normalised = _normalise_feature(label, raw)
        if not normalised:
            exclusions.append(
                {
                    "accepted_species": species,
                    "family": family,
                    "trait_name": label,
                    "raw_value": raw,
                    "source_url": url,
                    "reason": "ontology_mapping_not_supported",
                }
            )
            continue
        feature_id = _text(feature.get("feature_id")) or _stable_id(label)
        feature_citation = _text(feature.get("feature_citation"))
        citation = feature_citation or source["citation"]
        lineage = f"pladias:feature:{feature_id}:{_stable_id(citation)}"
        for trait_name, value in normalised:
            record_id = f"pladias:{feature_id}:{_stable_id(species, trait_name)}"
            candidates.append(
                {
                    "evidence_id": f"pladias-{_stable_id(record_id, raw, value)}",
                    "source_name": source["name"],
                    "source_record_id": record_id,
                    "accepted_species": species,
                    "genus": species.split()[0],
                    "family": family,
                    "name_match_method": "accepted_name_exact",
                    "trait_name": trait_name,
                    "raw_value": raw,
                    "standardized_value": value,
                    "candidate_kind": "source_backed",
                    "evidence_scope": "species_direct",
                    "source_type": "institutional_trait_database",
                    "source_reliability": "A_institutional_trait_database",
                    "source_citation": citation,
                    "source_url": url,
                    "evidence_excerpt": feature["supporting_excerpt"],
                    "supporting_taxa": species,
                    "inference_rule": "",
                    "confidence": "high",
                    "inference_note": "",
                    "needs_human_review": "false",
                    "review_status": "accepted_source_package_contract",
                    "analysis_role": "confirmatory_species_direct",
                    "source_lineage": lineage,
                }
            )
    return candidates, exclusions


def prepare_pladias(
    taxon_list_xls: Path,
    master_csv: Path,
    output_dir: Path,
    config_path: Path,
    target_species_csv: Path | None = None,
    cache_dir: Path | None = None,
    workers: int = 4,
    timeout_seconds: float = 45.0,
    retries: int = 2,
    max_species: int | None = None,
) -> dict[str, Any]:
    """Build a reusable exact-name Pladias source package."""
    config = _load_config(config_path)
    source = config["source"]
    observed_hash = _sha256(taxon_list_xls)
    expected_hash = _text(source["expected_taxon_list_sha256"]).casefold()
    if observed_hash != expected_hash:
        raise ValueError(
            f"Pladias taxon-list SHA-256 mismatch: expected {expected_hash}, got {observed_hash}"
        )

    taxon_list = pd.read_excel(taxon_list_xls, sheet_name="TaxonList", dtype=str)
    exact_names = set(taxon_list["NameLatFinal"].dropna().map(_text))
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    required = {"accepted_species", "family"}
    missing = required.difference(master.columns)
    if missing:
        raise ValueError(f"Master is missing required columns: {sorted(missing)}")
    targets = master.loc[master["accepted_species"].isin(exact_names)].copy()
    targets = targets.loc[
        targets["accepted_species"].map(lambda value: bool(STRICT_SPECIES_RE.fullmatch(value)))
    ]
    if target_species_csv is not None:
        requested = pd.read_csv(target_species_csv, dtype=str).fillna("")
        if "accepted_species" not in requested:
            raise ValueError("Target species CSV must contain accepted_species")
        requested_names = set(requested["accepted_species"].map(_text))
        targets = targets.loc[targets["accepted_species"].isin(requested_names)]
    targets = targets.sort_values("accepted_species").drop_duplicates("accepted_species")
    if max_species is not None:
        targets = targets.head(max_species)
    if targets.empty:
        raise ValueError("No exact species-level Pladias targets overlap the master")

    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = cache_dir or output_dir / "page_cache"
    pages: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        futures = {
            pool.submit(
                _fetch_one,
                row.accepted_species,
                row.family,
                source,
                cache_dir,
                timeout_seconds,
                retries,
            ): row.accepted_species
            for row in targets.itertuples(index=False)
        }
        for future in as_completed(futures):
            pages.append(future.result())

    candidates: list[dict[str, str]] = []
    exclusions: list[dict[str, str]] = []
    for page in pages:
        page_candidates, page_exclusions = _candidate_rows(page, source)
        candidates.extend(page_candidates)
        exclusions.extend(page_exclusions)

    candidate_frame = pd.DataFrame(candidates, columns=CANDIDATE_COLUMNS)
    if not candidate_frame.empty:
        candidate_frame = candidate_frame.sort_values(
            ["accepted_species", "trait_name", "source_record_id"]
        ).drop_duplicates(["accepted_species", "trait_name", "source_lineage"])
    exclusion_frame = pd.DataFrame(exclusions, columns=EXCLUSION_COLUMNS)
    fetch_rows = [
        {column: page.get(column, "") for column in FETCH_COLUMNS} for page in pages
    ]
    fetch_frame = pd.DataFrame(fetch_rows, columns=FETCH_COLUMNS).sort_values(
        "accepted_species"
    )
    candidate_path = output_dir / "trait_candidates.csv.gz"
    fetch_path = output_dir / "fetch_audit.csv.gz"
    exclusion_path = output_dir / "exclusions.csv.gz"
    quality_audit_path = output_dir / "quality_audit.csv.gz"
    quality_report_path = output_dir / "quality_audit_summary.json"
    _write_csv_gzip(candidate_frame, candidate_path)
    _write_csv_gzip(fetch_frame, fetch_path)
    _write_csv_gzip(exclusion_frame, exclusion_path)
    quality_audit, quality_report = _quality_audit(candidate_frame, fetch_frame)
    _write_csv_gzip(quality_audit, quality_audit_path)
    quality_report_path.write_text(
        json.dumps(quality_report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    identity_counts = fetch_frame["identity_status"].value_counts().sort_index().to_dict()
    trait_counts = (
        candidate_frame["trait_name"].value_counts().sort_index().to_dict()
        if not candidate_frame.empty
        else {}
    )
    manifest = {
        "contract": "pladias_source_package_v1",
        "created_at_utc": datetime.now(UTC_TZ).isoformat(),
        "source": {
            "name": source["name"],
            "citation": source["citation"],
            "website": source["website"],
            "usage_terms_url": source["usage_terms_url"],
            "usage_contract": source["usage_contract"],
            "large_scale_use_action": source["large_scale_use_action"],
            "robots_url": source["robots_url"],
            "taxon_list_url": source["taxon_list_url"],
            "taxon_list_sha256": observed_hash,
        },
        "selection": {
            "master_csv": str(master_csv),
            "target_species_csv": str(target_species_csv or ""),
            "official_taxon_names": len(exact_names),
            "exact_master_overlap_fetched": len(targets),
            "strict_species_rank_only": True,
        },
        "counts": {
            "pages": len(fetch_frame),
            "http_success": int(fetch_frame["http_status"].astype(str).eq("200").sum()),
            "identity_status": {str(key): int(value) for key, value in identity_counts.items()},
            "candidate_rows": len(candidate_frame),
            "candidate_species": int(candidate_frame["accepted_species"].nunique())
            if not candidate_frame.empty
            else 0,
            "candidate_traits": {str(key): int(value) for key, value in trait_counts.items()},
            "excluded_rows": len(exclusion_frame),
        },
        "quality_contract": {
            "tier": "High",
            "identity_gate": "exact accepted species name plus family breadcrumb",
            "pollen_vector_or_reward_projected_to_strict_axes": False,
            "genus_inference_emitted": False,
            "audit": quality_report,
        },
        "outputs": {
            "trait_candidates": str(candidate_path),
            "fetch_audit": str(fetch_path),
            "exclusions": str(exclusion_path),
            "quality_audit": str(quality_audit_path),
            "quality_audit_summary": str(quality_report_path),
            "page_cache": str(cache_dir),
        },
    }
    manifest_path = output_dir / "source_package_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(manifest["counts"], sort_keys=True))
    return manifest


@app.command("prepare")
def prepare_command(
    taxon_list_xls: Annotated[Path, typer.Option(exists=True)],
    master_csv: Annotated[Path, typer.Option(exists=True)],
    output_dir: Annotated[Path, typer.Option()],
    config_path: Annotated[Path, typer.Option(exists=True)] = Path(
        "config/pladias_traits.yml"
    ),
    target_species_csv: Annotated[Path | None, typer.Option()] = None,
    cache_dir: Annotated[Path | None, typer.Option()] = None,
    workers: Annotated[int, typer.Option(min=1, max=16)] = 4,
    timeout_seconds: Annotated[float, typer.Option()] = 45.0,
    retries: Annotated[int, typer.Option(min=0, max=5)] = 2,
    max_species: Annotated[int | None, typer.Option(min=1)] = None,
) -> None:
    prepare_pladias(
        taxon_list_xls=taxon_list_xls,
        master_csv=master_csv,
        output_dir=output_dir,
        config_path=config_path,
        target_species_csv=target_species_csv,
        cache_dir=cache_dir,
        workers=workers,
        timeout_seconds=timeout_seconds,
        retries=retries,
        max_species=max_species,
    )


if __name__ == "__main__":
    app()
