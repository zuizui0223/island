"""Acquire six visible floral traits from the complete PlantNET NSW Flora index.

PlantNET is an official National Herbarium of New South Wales flora.  The lane
downloads its public species index once, fetches only exact fixed-universe
species profiles, caches every completed page, and extracts evidence only from
the profile's Description section.  Search-result text is never evidence.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import json
import re
import time
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urljoin, urlparse

import pandas as pd
import requests
from bs4 import BeautifulSoup

from analysis.acquire_flora_of_australia_traits_20260809 import (
    AXIS,
    CULTIVAR,
    extract_description,
)
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS
from island_v2.trait_measurements import load_rules

BASE = "https://plantnet.rbgsyd.nsw.gov.au"
INDEX_URL = f"{BASE}/PlantNet/NSWflora/fullindex.html"
PROFILE_PATH = "/cgi-bin/NSWfl.pl?page=nswfl&lvl=sp&name={name}"
USER_AGENT = (
    "island-trait-evidence-research/1.0 "
    "(academic provenance-preserving acquisition; "
    "https://github.com/zuizui0223/island)"
)
SOURCE_RUN_ID = "plantnet-nsw-source-scale-20260809"
SOURCE_ARTIFACT = "committed-plantnet-nsw-source-package-20260809"
SOURCE_FILE = "plantnet_nsw_evidence.csv.gz"
ACCEPTANCE_CONTRACT = "plantnet_nsw_official_description_v1"
AUDIT_SEED = "plantnet-nsw-audit-20260809"
AUDIT_SIZE = 200


def _text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _read_gzip(path: Path) -> str:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        return handle.read()


def _write_gzip(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write(value)


def _request_text(url: str, cache_path: Path, attempts: int = 5) -> tuple[str, str]:
    if cache_path.exists():
        return "cached", _read_gzip(cache_path)
    error = ""
    for attempt in range(attempts):
        try:
            response = requests.get(
                url,
                timeout=120,
                headers={"User-Agent": USER_AGENT},
            )
            response.raise_for_status()
            if len(response.content) < 1_000:
                raise RuntimeError(f"unexpected short response: {len(response.content)}")
            _write_gzip(cache_path, response.text)
            return "fetched", response.text
        except (requests.RequestException, RuntimeError) as exc:
            error = repr(exc)
            if attempt + 1 < attempts:
                time.sleep(2**attempt)
    raise RuntimeError(f"failed to fetch {url}: {error}")


def parse_index(html: str) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    soup = BeautifulSoup(html, "html.parser")
    for anchor in soup.find_all("a", href=True):
        query = parse_qs(urlparse(anchor["href"]).query)
        if query.get("lvl") != ["sp"] or not query.get("name"):
            continue
        name_param = query["name"][0]
        species = _text(name_param.replace("~", " "))
        if len(species.split()) != 2 or "×" in species:
            continue
        rows.append(
            {
                "plantnet_species": species,
                "name_parameter": name_param,
                "profile_url": urljoin(BASE, PROFILE_PATH.format(name=name_param)),
            }
        )
    index = pd.DataFrame(rows).drop_duplicates("plantnet_species")
    if len(index) < 7_500:
        raise RuntimeError(f"unexpected PlantNET species index size: {len(index)}")
    return index.sort_values("plantnet_species").reset_index(drop=True)


def acquire_index(cache_dir: Path) -> tuple[pd.DataFrame, str]:
    status, html = _request_text(INDEX_URL, cache_dir / "plantnet_fullindex.html.gz")
    return parse_index(html), status


def parse_profile(html: str) -> dict[str, str]:
    soup = BeautifulSoup(html, "html.parser")
    heading = soup.select_one("td.head2 i")
    identity = _text(heading.get_text(" ", strip=True) if heading else "")
    family = ""
    for cell in soup.select('td[align="right"]'):
        text = _text(cell.get_text(" ", strip=True))
        match = re.fullmatch(r"Family\s+([A-Za-z-]+)", text, re.IGNORECASE)
        if match:
            family = match.group(1)
            break

    description_label = soup.find(
        "b",
        string=lambda value: bool(value and _text(value).casefold() == "description:"),
    )
    description = ""
    main_text = ""
    if description_label:
        container = description_label.find_parent("td") or soup.body
        main_text = _text(container.get_text(" ", strip=True) if container else "")
        if "Description:" in main_text:
            description = main_text.split("Description:", maxsplit=1)[1]
            for delimiter in (
                "Distribution and occurrence:",
                "Distribution:",
                "Text by ",
            ):
                description = description.split(delimiter, maxsplit=1)[0]
            description = _text(description)

    source_text = main_text or _text(soup.get_text(" ", strip=True))
    attribution = ""
    match = re.search(
        r"\bText by\s+(.+?)(?=\s+Taxon concept:|\s+N\.S\.W\. occurrences|$)",
        source_text,
        re.IGNORECASE,
    )
    if match:
        attribution = f"Text by {_text(match.group(1))}"
    return {
        "profile_scientific_name": identity,
        "profile_family": family,
        "description": description,
        "source_attribution": attribution,
    }


def _fetch_profile(row: dict[str, str], cache_dir: Path) -> dict[str, str]:
    cache_path = cache_dir / "profiles" / f"{_sha(row['plantnet_species'])[:24]}.html.gz"
    status, html = _request_text(row["profile_url"], cache_path)
    return {**row, "fetch_status": status, "content_sha256": _sha(html), **parse_profile(html)}


def acquire_profiles(
    selected: pd.DataFrame,
    cache_dir: Path,
    *,
    workers: int,
) -> pd.DataFrame:
    records = selected.to_dict("records")
    fetched: list[dict[str, str]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(_fetch_profile, row, cache_dir) for row in records]
        for completed, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            fetched.append(future.result())
            if completed % 250 == 0 or completed == len(futures):
                print(f"profiles {completed}/{len(futures)}", flush=True)
    return pd.DataFrame(fetched).sort_values("plantnet_species").reset_index(drop=True)


def _source_lineage(citation: str, quote: str, species: str) -> tuple[str, str]:
    normalized = re.sub(
        r"[^a-z0-9]+",
        " ",
        f"{citation}|{quote}".casefold(),
    ).strip()
    if normalized:
        return (
            f"plantnet-text:{_sha(normalized)[:32]}",
            "citation_plus_normalized_excerpt_fingerprint",
        )
    return f"plantnet-profile:{_sha(species)[:32]}", "official_profile_name_fallback"


def build_evidence(
    inventory: pd.DataFrame,
    master: pd.DataFrame,
    baseline_pairs: set[tuple[str, str]],
    rules: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    family_map = dict(zip(master["accepted_species"], master["family"], strict=True))
    rows: list[dict[str, str]] = []
    treatment_rows: list[dict[str, Any]] = []
    for item in inventory.itertuples(index=False):
        species = _text(item.plantnet_species)
        expected_family = _text(family_map.get(species))
        family = _text(item.profile_family)
        identity_ok = bool(
            species in family_map
            and _text(item.profile_scientific_name) == species
            and family
            and family == expected_family
        )
        description = _text(item.description)
        cultivar_excluded = bool(CULTIVAR.search(f"{species} {description}"))
        extracted = (
            extract_description(description, rules)
            if identity_ok and description and not cultivar_excluded
            else {}
        )
        novel = {
            trait: result
            for trait, result in extracted.items()
            if (species, trait) not in baseline_pairs
        }
        source_citation = "PlantNET NSW FloraOnline"
        if _text(item.source_attribution):
            source_citation += f"; {_text(item.source_attribution)}"
        for trait, result in novel.items():
            value = _text(result["value"])
            quote = _text(result["quote"])
            lineage, lineage_method = _source_lineage(source_citation, quote, species)
            record_id = f"plantnet-nsw-{_sha(f'{species}|{trait}|{value}|{quote}')[:24]}"
            rows.append(
                {
                    "accepted_species": species,
                    "axis": AXIS[trait],
                    "trait_name": trait,
                    "normalized_value": value,
                    "quality": "high",
                    "source_group": "latest_public_web",
                    "source_provider": "plantnet_nsw_flora_source_scale",
                    "source_url": _text(item.profile_url),
                    "source_record_id": record_id,
                    "source_citation": source_citation,
                    "source_excerpt": quote,
                    "evidence_scope": "species_direct",
                    "name_match_method": "accepted_name_exact_plus_family_exact",
                    "source_lineage": lineage,
                    "lineage_method": lineage_method,
                    "source_run_id": SOURCE_RUN_ID,
                    "source_artifact": SOURCE_ARTIFACT,
                    "source_file": SOURCE_FILE,
                    "acceptance_contract": ACCEPTANCE_CONTRACT,
                }
            )
        treatment_rows.append(
            {
                "accepted_species": species,
                "profile_scientific_name": _text(item.profile_scientific_name),
                "profile_family": family,
                "master_family": expected_family,
                "identity_gate_passed": identity_ok,
                "cultivar_excluded": cultivar_excluded,
                "description_chars": len(description),
                "all_extracted_traits": "|".join(sorted(extracted)),
                "novel_extracted_traits": "|".join(sorted(novel)),
                "source_attribution": _text(item.source_attribution),
                "source_url": _text(item.profile_url),
                "content_sha256": _text(item.content_sha256),
                "fetch_status": _text(item.fetch_status),
            }
        )
    evidence = pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)
    if not evidence.empty:
        evidence = evidence.drop_duplicates("source_record_id").sort_values(
            ["trait_name", "accepted_species", "source_record_id"],
            kind="stable",
        )
    return evidence.reset_index(drop=True), pd.DataFrame(treatment_rows)


def deterministic_audit_queue(
    evidence: pd.DataFrame,
    size: int = AUDIT_SIZE,
) -> pd.DataFrame:
    if len(evidence) < size:
        raise ValueError(f"only {len(evidence)} candidates available for {size}-row audit")
    pool = evidence.copy()
    pool["_audit_hash"] = pool["source_record_id"].map(
        lambda value: _sha(f"{AUDIT_SEED}|{value}")
    )
    selected = [
        group.sort_values("_audit_hash").head(min(20, len(group)))
        for _, group in pool.groupby("trait_name", sort=True)
    ]
    audit = pd.concat(selected).drop_duplicates("source_record_id")
    remaining = pool.loc[~pool.index.isin(audit.index)].sort_values("_audit_hash")
    if len(audit) < size:
        audit = pd.concat([audit, remaining.head(size - len(audit))])
    audit = audit.sort_values(["trait_name", "_audit_hash"], kind="stable").head(size)
    return pd.DataFrame(
        {
            "candidate_id": audit["source_record_id"],
            "accepted_species": audit["accepted_species"],
            "trait_name": audit["trait_name"],
            "normalized_value": audit["normalized_value"],
            "provider": audit["source_provider"],
            "source_url": audit["source_url"],
            "source_lineage": audit["source_lineage"],
            "source_citation": audit["source_citation"],
            "exact_supporting_quote": audit["source_excerpt"],
            "name_match_method": audit["name_match_method"],
            "cultivar_status": "wild_or_naturalised_flora_treatment",
            "accepted_correct": "",
            "reviewer": "",
            "reviewed_at_utc": "",
            "audit_reason": "",
        }
    ).reset_index(drop=True)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(
    *,
    master_csv: Path,
    baseline_evidence_csvs: list[Path],
    universe_coverage_csv: Path,
    cache_dir: Path,
    output_dir: Path,
    measurement_config: Path,
    workers: int,
) -> dict[str, Any]:
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    universe = pd.read_csv(
        universe_coverage_csv,
        dtype=str,
        usecols=["accepted_species"],
    ).fillna("")
    universe_names = set(universe["accepted_species"])
    if len(universe_names) != 106_295:
        raise ValueError("coverage universe must contain exactly 106,295 species")
    master = master.loc[master["accepted_species"].isin(universe_names)].copy()
    if len(master) != 106_295 or master["accepted_species"].nunique() != 106_295:
        raise ValueError("master must uniquely cover the fixed universe")
    baseline_pairs: set[tuple[str, str]] = set()
    for path in baseline_evidence_csvs:
        baseline = pd.read_csv(path, dtype=str).fillna("")
        baseline_pairs.update(zip(baseline["accepted_species"], baseline["trait_name"]))

    index, index_status = acquire_index(cache_dir)
    selected = index.loc[index["plantnet_species"].isin(universe_names)].copy()
    inventory = acquire_profiles(selected, cache_dir, workers=workers)
    rules = load_rules(measurement_config)
    evidence, treatments = build_evidence(
        inventory,
        master,
        baseline_pairs,
        rules,
    )
    audit = deterministic_audit_queue(evidence)
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(output_dir / SOURCE_FILE, index=False, compression="gzip")
    treatments.to_csv(
        output_dir / "plantnet_nsw_treatment_inventory.csv.gz",
        index=False,
        compression="gzip",
    )
    index.to_csv(
        output_dir / "plantnet_nsw_species_index.csv.gz",
        index=False,
        compression="gzip",
    )
    audit.to_csv(output_dir / "plantnet_nsw_audit_queue_200.csv", index=False)
    summary = {
        "contract": "plantnet_nsw_source_scale_acquisition_v1",
        "source_run_id": SOURCE_RUN_ID,
        "source_landing_url": INDEX_URL,
        "collection_species_profiles": len(index),
        "exact_master_profiles": len(selected),
        "identity_gate_passed": int(treatments["identity_gate_passed"].sum()),
        "profiles_with_description": int(treatments["description_chars"].gt(0).sum()),
        "index_fetch_status": index_status,
        "profile_fetch_status": treatments["fetch_status"].value_counts().to_dict(),
        "search_queries": 0,
        "api_cost_usd": 0.0,
        "baseline_species_trait_pairs": len(baseline_pairs),
        "novel_candidate_rows": len(evidence),
        "novel_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "novel_species_axis": int(
            evidence[["accepted_species", "axis"]].drop_duplicates().shape[0]
        ),
        "unique_species_with_novel_evidence": int(evidence["accepted_species"].nunique()),
        "by_trait": (
            evidence.groupby("trait_name")
            .agg(rows=("accepted_species", "size"), species=("accepted_species", "nunique"))
            .reset_index()
            .to_dict("records")
        ),
        "audit_queue_rows": len(audit),
        "audit_seed": AUDIT_SEED,
        "measurement_classification_rule_version": rules.get(
            "classification_rule_version"
        ),
        "guardrails": {
            "species_identity": "accepted_name_exact_plus_family_exact",
            "rank": "species_index_only",
            "description_field": "Description_section_only",
            "cultivar_hybrid_excluded": True,
            "snippet_as_evidence": False,
            "family_inference": False,
            "global_fallback": False,
            "axis_level_extraction": False,
        },
    }
    (output_dir / "plantnet_nsw_source_scale_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    manifest = []
    for path in sorted(output_dir.iterdir()):
        if path.is_file() and path.name != "file_manifest.sha256":
            manifest.append(f"{_sha256_file(path)}  {path.name}")
    (output_dir / "file_manifest.sha256").write_text(
        "\n".join(manifest) + "\n",
        encoding="ascii",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--master-csv",
        type=Path,
        default=Path("data/v2/staging/gbif/collected/island_taxa.csv"),
    )
    parser.add_argument(
        "--baseline-evidence-csv",
        type=Path,
        action="append",
        required=True,
    )
    parser.add_argument("--universe-coverage-csv", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--measurement-config",
        type=Path,
        default=Path("config/measurement_classification.yml"),
    )
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if not 1 <= args.workers <= 16:
        raise ValueError("workers must be between 1 and 16")
    print(
        json.dumps(
            run(
                master_csv=args.master_csv,
                baseline_evidence_csvs=args.baseline_evidence_csv,
                universe_coverage_csv=args.universe_coverage_csv,
                cache_dir=args.cache_dir,
                output_dir=args.output_dir,
                measurement_config=args.measurement_config,
                workers=args.workers,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
